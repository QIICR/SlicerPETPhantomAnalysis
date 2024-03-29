/*==============================================================================
 
 Program: PETPhantomAnlysisCLI
 
 (c) Copyright University of Iowa All Rights Reserved.
 
 See COPYRIGHT.txt
 or http://www.slicer.org/copyright/copyright.txt for details.
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 
 ==============================================================================*/

#ifndef __itkCylinderUniformityMeasurementImageFilter_h
#define __itkCylinderUniformityMeasurementImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorContainer.h"
#include "itkVector.h"

namespace itk
{

/**\class CylinderUniformityMeasurementImageFilter
 * \brief Provide unformity measurements for a cylinder in an image
 * \author	Chrsitian Bauer
 * Obtains unformity measurements for a cylindrical object in an image, given
 * the meausrement cylinder's center, main-axis direction, radius, and height.
 * The measurements include measurements (and offsets) in individual slices,
 * mean, standard deviation and maximum relative difference.
 * In addition, the filter outputs a label image with the used measurement
 * regions.
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

  /** Set the cylinder center. */
  itkSetMacro(Center, PointType);

  /** Get the cylinder center. */
  itkGetMacro(Center, PointType);

  /** Set the cylinder orientation. */
  itkSetMacro(Direction, VectorType);

  /** Get the cylinder orientation. */
  itkGetMacro(Direction, VectorType);

  /** Set the cylinder radius. */
  itkSetMacro(Radius, double);

  /** Get the cylinder radius. */
  itkGetMacro(Radius, double);

  /** Set the cylinder height. */
  itkSetMacro(Height, double);

  /** Get the cylinder height'. */
  itkGetMacro(Height, double);

  /** Get mean of slice measurements. */
  itkGetMacro(CylinderMean, double);

  /** Get standard deviation of slice measurements. */
  itkGetMacro(CylinderStd, double);

  /** Get the maximum relative difference among all slices measurements. */
  itkGetMacro(MaxRelativeDifference, double);

  /** Get measurements for individual slices. */
  itkGetObjectMacro(SliceMeasurements, MeasurementVectorType);

  /** Get offsets for the invidivual slices. */
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
