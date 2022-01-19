#ifndef __itkCylinderMatchingImageFilter_h
#define __itkCylinderMatchingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorContainer.h"
#include "itkVector.h"

namespace itk
{

/**\class CylinderMatchingImageFilter
 * \brief short description
 * More detailed desciption
 */
template< class TInputImage, class TOutputImage=TInputImage >
class CylinderMatchingImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class type aliases. */
  using Self = CylinderMatchingImageFilter;
  using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer = SmartPointer< Self >;

  /** Useful type aliases. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using LabelImageType = Image<short, TInputImage::ImageDimension>;
  using PixelType = typename TInputImage::PixelType;
  using PointType = typename TInputImage::PointType;
  using IndexType = typename TInputImage::IndexType;
  using VectorType = Vector<double, TInputImage::ImageDimension>;

  ITK_DISALLOW_COPY_AND_ASSIGN(CylinderMatchingImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CylinderMatchingImageFilter, ImageToImageFilter);

  /** Dimension of the underlying image. */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  /** Set/Get Macros */
  itkSetMacro(SmoothingSigma, double);
  itkGetMacro(SmoothingSigma, double);
  itkGetMacro(Center, PointType);
  itkGetMacro(Direction, VectorType);

protected:
  CylinderMatchingImageFilter();
  ~CylinderMatchingImageFilter() override = default;

  void PrintSelf(std::ostream& os, Indent indent) const override;

  /** Does the real work. */
  void GenerateData() override;

private:
  double m_SmoothingSigma{ 1.0 };
  PointType m_Center;
  VectorType m_Direction;

  using PointList = std::vector< PointType >;
  PointList GetCylinderAxisPoints(typename OutputImageType::Pointer segmentation, double threshold=1.0) const;
  void FitLine(const PointList& axisPoints, PointType& center, VectorType& direction) const;
  PointList RobustFitLine(const PointList& axisPoints, PointType& center, VectorType& direction, double threshold) const;
  double DistanceFromAxis(const PointType& p, PointType& center, VectorType& direction) const;
  std::pair<double, double> GetThresholds(std::vector<double> values, double qLow=0.25, double qHigh=0.75, double sigmaMul=2.0) const;
};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCylinderMatchingImageFilter.cxx"
#endif

#endif // __itkCylinderMatchingImageFilter_h
