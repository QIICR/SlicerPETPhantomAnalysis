#ifndef __itkCylinderUniformityMeasurementImageFilter_cxx
#define __itkCylinderUniformityMeasurementImageFilter_cxx

#include "itkCylinderUniformityMeasurementImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <numeric>

namespace itk
{
//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
CylinderUniformityMeasurementImageFilter< TInputImage, TOutputImage >
::CylinderUniformityMeasurementImageFilter()
{
  m_Center.Fill(0);
  m_Direction.Fill(0);
  m_SliceMeasurements = MeasurementVectorType::New();
  m_SliceOffsets = MeasurementVectorType::New();
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
void CylinderUniformityMeasurementImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  typename InputImageType::ConstPointer image = this->GetInput();
  size_t nSlices = image->GetLargestPossibleRegion().GetSize()[2];

  // initialize output measurements
  m_SliceMeasurements->Initialize();
  m_SliceMeasurements->Reserve(nSlices);
  m_SliceOffsets->Initialize();
  m_SliceOffsets->Reserve(nSlices);

  // initialize output image
  typename OutputImageType::Pointer measurementRegion = this->GetOutput();
  measurementRegion->CopyInformation(image);
  measurementRegion->SetRegions(image->GetLargestPossibleRegion());
  measurementRegion->Allocate();
  measurementRegion->FillBuffer(0);

  // temporary variables
  std::vector<double> cylinderVoxels;
  std::vector<size_t> sliceVoxelsInRegion(nSlices, 0.0);

  // collect measurements and create output measurement region
  itk::ImageRegionConstIteratorWithIndex<InputImageType> imageIt(image, image->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<OutputImageType> measurementRegionIt(measurementRegion, measurementRegion->GetLargestPossibleRegion());
  imageIt.GoToBegin();
  measurementRegionIt.GoToBegin();
  PointType p;
  while (!imageIt.IsAtEnd())
  {
    image->TransformIndexToPhysicalPoint(imageIt.GetIndex(), p);
    double distanceFromAxis = this->DistanceFromAxis(p);
    double distanceAlongAxis = this->DistanceAlongAxis(p);
    double isInside = this->IsInside(p);
    if (distanceFromAxis<m_Radius)
    {
      size_t z = imageIt.GetIndex()[2];
      m_SliceMeasurements->SetElement(z,
        m_SliceMeasurements->GetElement(z)+imageIt.Value());
      sliceVoxelsInRegion[z]+=1;
    }
    if (isInside)
    {
      cylinderVoxels.push_back(imageIt.Value());
    }

    OutputPixelType regionLabel = m_LabelOutside;
    if (isInside)
      regionLabel = m_LabelInside;
    else if (distanceFromAxis<m_Radius)
      regionLabel = m_LabelInsideRadiusLimit;
    else if (fabs(distanceAlongAxis)<m_Height/2.0)
      regionLabel = m_LabelInsideHeightLimit;
    measurementRegionIt.Set(regionLabel);

    ++imageIt;
    ++measurementRegionIt;
  }

  // normalize slice based measurements based on number of measured voxels
  for (size_t z=0; z<nSlices; z++)
    if (sliceVoxelsInRegion[z]>0)
      m_SliceMeasurements->SetElement(z,
        m_SliceMeasurements->GetElement(z)/double(sliceVoxelsInRegion[z]));

  // calculate mean and std
  m_CylinderMean = std::accumulate(cylinderVoxels.begin(),
    cylinderVoxels.end(), 0.0)/double(cylinderVoxels.size());
  m_CylinderStd = 0.0;
  for (const auto& it : cylinderVoxels)
    m_CylinderStd += pow(it-m_CylinderMean,2.0);
  m_CylinderStd = sqrt(m_CylinderStd/double(cylinderVoxels.size()));

  // calculate slice offsets
  IndexType idx; idx.Fill(0);
  for (size_t z=0; z<nSlices; ++z)
  {
    idx[2] = z;
    image->TransformIndexToPhysicalPoint(idx, p);
    m_SliceOffsets->SetElement(z, p[2]-m_Center[2]);
  }

  // calculate maximum relative difference
  m_MaxRelativeDifference = 0.0;
  for (size_t z=0; z<nSlices; z++)
  {
    double oz = m_SliceOffsets->GetElement(z);
    if (fabs(oz)<(m_Direction[2]*m_Height/2.0))
    {
      double relativeDifference =
        (m_SliceMeasurements->GetElement(z)-m_CylinderMean)/m_CylinderMean;
      if (fabs(relativeDifference)>fabs(m_MaxRelativeDifference))
        m_MaxRelativeDifference = relativeDifference;
    }
  }
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
bool CylinderUniformityMeasurementImageFilter< TInputImage, TOutputImage >
::IsInside(PointType p) const
{
  return this->DistanceFromAxis(p)<m_Radius &&
    fabs(this->DistanceAlongAxis(p))<m_Height/2.0;
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
double CylinderUniformityMeasurementImageFilter< TInputImage, TOutputImage >
::DistanceFromAxis(PointType p) const
{
  double dCenter = (p-m_Center).GetNorm();
  double dAlongAxis = DistanceAlongAxis(p);
  return std::sqrt(dCenter*dCenter-dAlongAxis*dAlongAxis);
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
double CylinderUniformityMeasurementImageFilter< TInputImage, TOutputImage >
::DistanceAlongAxis(PointType p) const
{
  return (p-m_Center)*m_Direction;
}

//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage >
void CylinderUniformityMeasurementImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Center: " << this->m_Center << std::endl;
  os << indent << "Direction: " << this->m_Direction << std::endl;
  os << indent << "Radius: " << this->m_Radius << std::endl;
  os << indent << "Height: " << this->m_Height << std::endl;
}

}// end namespace

#endif
