#ifndef __itkCylinderMatchingImageFilter_cxx
#define __itkCylinderMatchingImageFilter_cxx

#include "itkCylinderMatchingImageFilter.h"

#include "itkOtsuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

namespace itk
{
//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
CylinderMatchingImageFilter< TInputImage, TOutputImage >
::CylinderMatchingImageFilter() :
  m_SmoothingSigma(1.0), m_Volume(0.0)
{
  m_Center.Fill(0);
  m_Direction.Fill(0);
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
void CylinderMatchingImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  //----------------------------------------------------------------------------
  // processing steps: smoothing, find and apply Otsu threshold, label objects, sort
  // them by size, use calculate centroid and eigenvector of largest component
  // as cylinder parameters

  typename InputImageType::ConstPointer image = this->GetInput();
  
  typedef itk::SmoothingRecursiveGaussianImageFilter<InputImageType, InputImageType> SmoothingImageFilterType;
  typename SmoothingImageFilterType::Pointer smoothingFilter = SmoothingImageFilterType::New();
  smoothingFilter->SetSigma(m_SmoothingSigma);
  smoothingFilter->SetInput(image);

  typedef itk::OtsuThresholdImageFilter<InputImageType, LabelImageType > ThresholdImageFilterType;
  typename ThresholdImageFilterType::Pointer otsuFilter = ThresholdImageFilterType::New();
  otsuFilter->SetInput(smoothingFilter->GetOutput());
  otsuFilter->SetInsideValue(0);
  otsuFilter->SetOutsideValue(1);
  otsuFilter->Update();
  //double otsuThreshold = otsuFilter->GetThreshold(); // just for debugging

  typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType> LabelingImageFilterType;
  typename LabelingImageFilterType::Pointer labeingFilter = LabelingImageFilterType::New();
  labeingFilter->SetInput(otsuFilter->GetOutput());

  typedef itk::RelabelComponentImageFilter<LabelImageType, OutputImageType> RelabelingImageFilterType;
  typename RelabelingImageFilterType::Pointer relablingFilter = RelabelingImageFilterType::New();
  relablingFilter->SetInput(labeingFilter->GetOutput());
  relablingFilter->Update();

  typename OutputImageType::Pointer segmentation = relablingFilter->GetOutput();

  PointType tmpPoint;
  typedef std::vector< PointType > PointList;
  PointList segPoints;
  ImageRegionConstIteratorWithIndex<OutputImageType> it(segmentation,
    segmentation->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    if (it.Value()==1)
    {
      segmentation->TransformIndexToPhysicalPoint(it.GetIndex(), tmpPoint);
      segPoints.push_back(tmpPoint);
    }
  }
  VectorType centerVector; centerVector.Fill(0);
  for (typename PointList::const_iterator it=segPoints.begin(); it!=segPoints.end(); ++it)
    centerVector += it->GetVectorFromOrigin();
  m_Center = (centerVector/double(segPoints.size()));

  vnl_matrix< double > normalizedSecondOrderCentralMoments(3, 3, 0.0);
  for (typename PointList::const_iterator it=segPoints.begin(); it!=segPoints.end(); ++it)
  {
    vnl_vector<double> v = (*it - m_Center).GetVnlVector();
    normalizedSecondOrderCentralMoments += outer_product(v,v);
  }
  vnl_symmetric_eigensystem< double > eig(normalizedSecondOrderCentralMoments);

  m_Direction[0] = eig.get_eigenvector(2)[0];
  m_Direction[1] = eig.get_eigenvector(2)[1];
  m_Direction[2] = eig.get_eigenvector(2)[2];

  // get volume (for plausibility check)
  typename InputImageType::SpacingType spacing = image->GetSpacing();
  double voxelVolume = spacing[0]*spacing[1]*spacing[2];
  m_Volume = double(segPoints.size())*voxelVolume;

  this->GetOutput()->Graft(relablingFilter->GetOutput());
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
void CylinderMatchingImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Center: " << this->m_Center << std::endl;
  os << indent << "Direction: " << this->m_Direction << std::endl;
  os << indent << "Volume: " << this->m_Volume << std::endl;
}

}// end namespace

#endif
