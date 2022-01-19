#ifndef __itkCylinderUniformityMeasurementImageFilter_h
#define __itkCylinderUniformityMeasurementImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorContainer.h"
#include "itkVector.h"

namespace itk
{

/**\class CylinderUniformityMeasurementImageFilter
 * \brief short description
 * More detailed desciption
 */
template< class TInputImage, class TOutputImage=TInputImage >
class CylinderUniformityMeasurementImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef CylinderUniformityMeasurementImageFilter                         Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage >  Superclass;
  typedef SmartPointer< Self >                            Pointer;

  /** Useful class typedefs*/
  typedef TInputImage                                      InputImageType;
  typedef TOutputImage                                     OutputImageType;
  typedef typename TInputImage::PixelType                  PixelType;
  typedef typename TOutputImage::PixelType                 OutputPixelType;
  typedef typename TInputImage::PointType                  PointType;
  typedef typename TInputImage::IndexType                  IndexType;
  typedef Vector<double, TInputImage::ImageDimension > VectorType;
  typedef VectorContainer<size_t, double>                  MeasurementVectorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CylinderUniformityMeasurementImageFilter, ImageToImageFilter);

  /** Dimension of the underlying image. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Set/Get Macros */
  itkSetMacro(Center, PointType);
  itkGetMacro(Center, PointType);
  itkSetMacro(Direction, VectorType);
  itkGetMacro(Direction, VectorType);
  itkSetMacro(Radius, double);
  itkGetMacro(Radius, double);
  itkSetMacro(Height, double);
  itkGetMacro(Height, double);
  itkGetMacro(CylinderMean, double);
  itkGetMacro(CylinderStd, double);
  itkGetMacro(MaxRelativeDifference, double);
  itkGetObjectMacro(SliceMeasurements, MeasurementVectorType);
  itkGetObjectMacro(SliceOffsets, MeasurementVectorType);

  itkSetMacro(LabelInside, OutputPixelType);
  itkGetMacro(LabelInside, OutputPixelType);
  itkSetMacro(LabelInsideRadiusLimit, OutputPixelType);
  itkGetMacro(LabelInsideRadiusLimit, OutputPixelType);
  itkSetMacro(LabelInsideHeightLimit, OutputPixelType);
  itkGetMacro(LabelInsideHeightLimit, OutputPixelType);
  itkSetMacro(LabelOutside, OutputPixelType);
  itkGetMacro(LabelOutside, OutputPixelType);

  bool IsInside(PointType p) const;
  double DistanceFromAxis(PointType p) const;
  double DistanceAlongAxis(PointType p) const;

protected:
  CylinderUniformityMeasurementImageFilter();
  ~CylinderUniformityMeasurementImageFilter() override = default;

  void PrintSelf(std::ostream& os, Indent indent) const override;

  /** Does the real work. */
  void GenerateData() override;

private:
  CylinderUniformityMeasurementImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  PointType m_Center;
  VectorType m_Direction;
  double m_Radius;
  double m_Height;
  MeasurementVectorType::Pointer m_SliceMeasurements;
  MeasurementVectorType::Pointer m_SliceOffsets;
  double m_CylinderMean;
  double m_CylinderStd;
  double m_MaxRelativeDifference;
  OutputPixelType m_LabelInside;
  OutputPixelType m_LabelInsideRadiusLimit;
  OutputPixelType m_LabelInsideHeightLimit;
  OutputPixelType m_LabelOutside;
};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCylinderUniformityMeasurementImageFilter.cxx"
#endif

#endif // __itkCylinderUniformityMeasurementImageFilter_h
