#ifndef __itkCylinderMatchingImageFilter_cxx
#define __itkCylinderMatchingImageFilter_cxx

#include "itkCylinderMatchingImageFilter.h"

#include "itkOtsuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include <math.h>

namespace itk
{
//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
CylinderMatchingImageFilter< TInputImage, TOutputImage >
::CylinderMatchingImageFilter()
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

  using SmoothingImageFilterType = itk::SmoothingRecursiveGaussianImageFilter<InputImageType, InputImageType>;
  auto smoothingFilter = SmoothingImageFilterType::New();
  smoothingFilter->SetSigma(m_SmoothingSigma);
  smoothingFilter->SetInput(image);

  using ThresholdImageFilterType = itk::OtsuThresholdImageFilter<InputImageType, LabelImageType > ;
  auto otsuFilter = ThresholdImageFilterType::New();
  otsuFilter->SetInput(smoothingFilter->GetOutput());
  otsuFilter->SetInsideValue(0);
  otsuFilter->SetOutsideValue(1);
  otsuFilter->SetReturnBinMidpoint(false);
  otsuFilter->Update();
  //double otsuThreshold = otsuFilter->GetThreshold(); // just for debugging

  using LabelingImageFilterType = itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType>;
  auto labeingFilter = LabelingImageFilterType::New();
  labeingFilter->SetInput(otsuFilter->GetOutput());

  using RelabelingImageFilterType = itk::RelabelComponentImageFilter<LabelImageType, OutputImageType>;
  auto relablingFilter = RelabelingImageFilterType::New();
  relablingFilter->SetInput(labeingFilter->GetOutput());
  relablingFilter->Update();

  typename OutputImageType::Pointer segmentation = relablingFilter->GetOutput();
  this->GetOutput()->Graft(relablingFilter->GetOutput());

/*
 // non-robust version
  PointType tmpPoint;
  using PointList = std::vector< PointType >;
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
  this->FitLine(segPoints, m_Center, m_Direction);
*/

  double sliceSpacing = std::max(segmentation->GetSpacing()[0], segmentation->GetSpacing()[1]);
  PointList initialAxisPoints = this->GetCylinderAxisPoints(segmentation,2.0);
  //this->FitLine(initialAxisPoints, m_Center, m_Direction);
  this->RobustFitLine(initialAxisPoints, m_Center, m_Direction, 2.0*sliceSpacing);

}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
typename CylinderMatchingImageFilter< TInputImage, TOutputImage >::PointList
CylinderMatchingImageFilter< TInputImage, TOutputImage >
::GetCylinderAxisPoints(typename OutputImageType::Pointer segmentation, double threshold) const
{
  size_t nSlices = segmentation->GetLargestPossibleRegion().GetSize()[2];
  std::vector<double> sliceNumSegmentationVoxels(nSlices, 0);
  PointType tmpPoint; tmpPoint.Fill(0);
  VectorType tmpVector; tmpVector.Fill(0);
  std::vector<VectorType> sliceCenterOfGravity(nSlices, tmpVector);
  ImageRegionConstIteratorWithIndex<OutputImageType> it(segmentation,
    segmentation->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    if (it.Value()==1)
    {
      size_t sliceId = it.GetIndex()[2];
      segmentation->TransformIndexToPhysicalPoint(it.GetIndex(), tmpPoint);
      sliceNumSegmentationVoxels[sliceId] += 1;
      sliceCenterOfGravity[sliceId] += tmpPoint.GetVectorFromOrigin();
    }
  }
  std::vector<double> areaSizes;
  for (size_t i=0; i<nSlices; ++i)
  {
    double n = sliceNumSegmentationVoxels[i];
    if (n>0) sliceCenterOfGravity[i] /= n;
    if (n>10) areaSizes.push_back(n);
  }
  for (size_t i=0; i<areaSizes.size(); ++i) areaSizes[i] = sqrt(areaSizes[i]/3.14);
  std::sort(areaSizes.begin(), areaSizes.end());
  double median = areaSizes[size_t(round((areaSizes.size()-1.0)*0.5))];
  using PointList = std::vector< PointType >;
  PointList axisPoints;
  for (size_t i=0; i<nSlices; ++i)
  {
    double a = sliceNumSegmentationVoxels[i];
    double r = sqrt(a/3.14);
    if (fabs(r-median)<threshold)
      axisPoints.push_back(sliceCenterOfGravity[i]);
  }
  return axisPoints;
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
std::pair<double, double> CylinderMatchingImageFilter< TInputImage, TOutputImage >
::GetThresholds(std::vector<double> values, double qLow, double qHigh, double sigmaMul) const
{
  std::sort(values.begin(), values.end());
  double qLowValue = values[size_t(round((values.size()-1.0)*qLow))];
  double qHighValue = values[size_t(round((values.size()-1.0)*qHigh))];
  double mean = 0.0;
  double tmpCount = 0.0;
  for (size_t i=0; i<values.size(); ++i)
  {
    double a = values[i];
    if (a>=qLowValue && a<=qHighValue)
    {
      mean+=a;
      tmpCount+=1;
    }
  }
  mean/=tmpCount;
  tmpCount = 0.0;
  double std = 0.0;
  for (size_t i=0; i<values.size(); ++i)
  {
    double a = values[i];
    if (a>=qLowValue && a<=qHighValue)
    {
      std+=pow(a-mean,2.0);
      tmpCount+=1;
    }
  }
  std = sqrt(std/tmpCount);
  double tLow = mean-sigmaMul*std;
  double tHigh = mean+sigmaMul*std;
  return std::pair<double, double>(tLow, tHigh);
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
void CylinderMatchingImageFilter< TInputImage, TOutputImage >
::FitLine(const PointList& axisPoints, PointType& center, VectorType& direction) const
{
  VectorType centerVector; centerVector.Fill(0);
  for (typename PointList::const_iterator it=axisPoints.begin(); it!=axisPoints.end(); ++it)
    centerVector += it->GetVectorFromOrigin();
  center = (centerVector/double(axisPoints.size()));

  vnl_matrix< double > normalizedSecondOrderCentralMoments(3, 3, 0.0);
  for (typename PointList::const_iterator it=axisPoints.begin(); it!=axisPoints.end(); ++it)
  {
    vnl_vector<double> v = (*it - center).GetVnlVector();
    normalizedSecondOrderCentralMoments += outer_product(v,v);
  }
  vnl_symmetric_eigensystem< double > eig(normalizedSecondOrderCentralMoments);

  direction[0] = eig.get_eigenvector(2)[0];
  direction[1] = eig.get_eigenvector(2)[1];
  direction[2] = eig.get_eigenvector(2)[2];
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
typename CylinderMatchingImageFilter< TInputImage, TOutputImage >::PointList
CylinderMatchingImageFilter< TInputImage, TOutputImage >
::RobustFitLine(const PointList& initialAxisPoints, PointType& center, VectorType& direction, double threshold) const
{
  PointList currentAxisPoints = initialAxisPoints;
  bool hasOutlier = true;
  while (hasOutlier)
  {
    this->FitLine(currentAxisPoints, center, direction);
    size_t n = currentAxisPoints.size();
    std::vector<double> distances(n,0.0);
    for (size_t i=0; i<n; ++i)
      distances[i] = this->DistanceFromAxis(currentAxisPoints[i], center, direction);
    std::vector<double>::const_iterator maxIt = std::max_element(distances.begin(), distances.end());
    if (*maxIt>threshold)
      currentAxisPoints.erase(currentAxisPoints.begin()+(maxIt-distances.begin())); // remove farthest outlier
    else
      hasOutlier = false;
  }
  return currentAxisPoints;
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
double CylinderMatchingImageFilter< TInputImage, TOutputImage >
::DistanceFromAxis(const PointType& p, PointType& center, VectorType& direction) const
{
  double dCenter = (p-center).GetNorm();
  double dAlongAxis = (p-center)*direction;
  return std::sqrt(dCenter*dCenter-dAlongAxis*dAlongAxis);
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
void CylinderMatchingImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Center: " << this->m_Center << std::endl;
  os << indent << "Direction: " << this->m_Direction << std::endl;
}

}// end namespace

#endif
