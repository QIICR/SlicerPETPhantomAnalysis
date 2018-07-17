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
  /** Standard class typedefs. */
  typedef CylinderMatchingImageFilter                     Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;

  /** Useful class typedefs*/
  typedef TInputImage                                      InputImageType;
  typedef TOutputImage                                     OutputImageType;
  typedef Image<short, TInputImage::ImageDimension>        LabelImageType;
  typedef typename TInputImage::PixelType                  PixelType;
  typedef typename TInputImage::PointType                  PointType;
  typedef typename TInputImage::IndexType                  IndexType;
  typedef Vector<double, TInputImage::ImageDimension>      VectorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CylinderMatchingImageFilter, ImageToImageFilter);

  /** Dimension of the underlying image. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Set/Get Macros */
  itkSetMacro(SmoothingSigma, double);
  itkGetMacro(SmoothingSigma, double);
  itkGetMacro(Center, PointType);
  itkGetMacro(Direction, VectorType);

protected:
  CylinderMatchingImageFilter();
  ~CylinderMatchingImageFilter(){}

  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  /** Does the real work. */
  virtual void GenerateData();

private:
  CylinderMatchingImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  double m_SmoothingSigma;
  PointType m_Center;
  VectorType m_Direction;
  double m_Volume;

  typedef std::vector< PointType > PointList;
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
