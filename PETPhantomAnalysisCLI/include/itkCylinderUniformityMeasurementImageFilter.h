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
  /** Standard class type aliases. */
  using Self = CylinderUniformityMeasurementImageFilter;
  using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer = SmartPointer< Self >;

  /** Useful ctype aliases. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using PixelType = typename TInputImage::PixelType;
  using OutputPixelType = typename TOutputImage::PixelType;
  using PointType = typename TInputImage::PointType;
  using IndexType = typename TInputImage::IndexType;
  using VectorType = Vector<double, TInputImage::ImageDimension > ;
  using MeasurementVectorType = VectorContainer<size_t, double>;

  ITK_DISALLOW_COPY_AND_ASSIGN(CylinderUniformityMeasurementImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CylinderUniformityMeasurementImageFilter, ImageToImageFilter);

  /** Dimension of the underlying image. */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

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
  PointType m_Center;
  VectorType m_Direction;
  double m_Radius{ 0.0 };
  double m_Height{ 0.0 };
  MeasurementVectorType::Pointer m_SliceMeasurements;
  MeasurementVectorType::Pointer m_SliceOffsets;
  double m_CylinderMean{ 0.0 };
  double m_CylinderStd{ 0.0 };
  double m_MaxRelativeDifference{ 0.0 };
  OutputPixelType m_LabelInside{ 1 };
  OutputPixelType m_LabelInsideRadiusLimit{ 2 };
  OutputPixelType m_LabelInsideHeightLimit{ 3 };
  OutputPixelType m_LabelOutside{ 0 };
};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCylinderUniformityMeasurementImageFilter.cxx"
#endif

#endif // __itkCylinderUniformityMeasurementImageFilter_h
