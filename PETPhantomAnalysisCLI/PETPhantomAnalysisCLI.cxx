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
 
#include "itkImageFileWriter.h"

#include "itkPluginUtilities.h"

#include "PETPhantomAnalysisCLICLP.h"

#include "itkCylinderMatchingImageFilter.h"
#include "itkCylinderUniformityMeasurementImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include <iostream>

// rapidjson
#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/ostreamwrapper.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module. Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template<typename TValue>
itk::FixedArray<TValue, 3> lps2ras(const itk::FixedArray<TValue, 3>& lps)
{
  itk::FixedArray<TValue, 3> ras;
  ras[0] = -lps[0];
  ras[1] = -lps[1];
  ras[2] = lps[2];
  return ras;
}

template<typename TValue, unsigned int VLength>
void jsonAddMember(rapidjson::Value& element, rapidjson::Document::AllocatorType& allocator,
  const std::string& name, const itk::FixedArray<TValue, VLength>& values)
{
  using namespace rapidjson;
  Value r_values(kArrayType);
  for (size_t i=0; i<VLength; ++i)
    r_values.PushBack(values[i], allocator);
  Value elementName(name.c_str(), allocator);
  element.AddMember(elementName.Move(), r_values.Move(), allocator);
}

template<class T>
void jsonAddMember(rapidjson::Value& element, rapidjson::Document::AllocatorType& allocator,
  const std::string& name, const std::vector<T>& values)
{
  using namespace rapidjson;
  Value r_values(kArrayType);
  //for (std::vector<double>::const_iterator it=values.begin(); it!=values.end(); ++it)
  //  r_values.PushBack(*it, allocator);
  for (const auto& it : values)
    r_values.PushBack(it, allocator);
  Value elementName(name.c_str(), allocator);
  element.AddMember(elementName.Move(), r_values.Move(), allocator);
}

} // end of anonymous namespace


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  try
  {
    using PixelType = double;
    using LabelPixelType = short;
    using ImageType = itk::Image<PixelType, 3>;
    using LabelImageType = itk::Image<LabelPixelType, 3>;
    
    // read PET scan
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName( inputVolume.c_str() );
    reader->Update();
    ImageType::Pointer petScan = reader->GetOutput();
        
    // normalize input data, if necessary
    if (normalizationFactor!=1.0)
    {
      using NormalizerFilterType = itk::ShiftScaleImageFilter<ImageType, ImageType>;
      auto normalizer = NormalizerFilterType::New();
      normalizer->SetShift(0.0);
      normalizer->SetScale(normalizationFactor);
      normalizer->SetInput(petScan);
      normalizer->Update();
      petScan = normalizer->GetOutput();
    }
    
    // find cylinder center and orientation
    using CylinderMatchingType = itk::CylinderMatchingImageFilter<ImageType, LabelImageType>;
    auto cylinderMatching = CylinderMatchingType::New();
    cylinderMatching->SetInput(petScan);
    cylinderMatching->SetSmoothingSigma(12.0);
    cylinderMatching->Update();
    CylinderMatchingType::PointType cylinderCenter = cylinderMatching->GetCenter();
    CylinderMatchingType::VectorType cylinderDirection = cylinderMatching->GetDirection();
    if (cylinderDirection[2]<0) cylinderDirection *= -1.0;
    
    // measure calibration and uniformity
    using UniformityMeasurementFilterType = itk::CylinderUniformityMeasurementImageFilter<ImageType, LabelImageType>;
    auto uniformityMeasurements = UniformityMeasurementFilterType::New();
    uniformityMeasurements->SetInput(petScan);
    uniformityMeasurements->SetCenter(cylinderCenter);
    uniformityMeasurements->SetDirection(cylinderDirection);
    uniformityMeasurements->SetRadius(73.0); // todo: or make configurable?
    uniformityMeasurements->SetHeight(160.0); // todo: or make configurable?
    uniformityMeasurements->SetLabelInsideRadiusLimit(0);
    uniformityMeasurements->SetLabelInsideHeightLimit(0);
    uniformityMeasurements->Update();
    LabelImageType::Pointer measurementRegion = uniformityMeasurements->GetOutput();
        
    // write measurements
    std::ofstream writeFile;
    writeFile.open( returnParameterFile.c_str() );
    writeFile << "Mean_s = " << uniformityMeasurements->GetCylinderMean() << std::endl;
    writeFile << "Std_s = " << uniformityMeasurements->GetCylinderStd() << std::endl;
    writeFile << "MaxRelDiff_s = " << uniformityMeasurements->GetMaxRelativeDifference() << std::endl;
    writeFile.close();
    
    // write measurement region volume, if requested
    if (outputVolume.size()>0)
    {
      using WriterType = itk::ImageFileWriter<LabelImageType>;
      auto writer = WriterType::New();
      writer->SetFileName( outputVolume.c_str() );
      writer->SetInput( measurementRegion );
      writer->SetUseCompression(1);
      writer->Update();
    }
        
    // write JSON report, if requested
    if (measurementsData.size()>0)
    {
      // prepare JSON document
      using namespace rapidjson;
      Document document;
      document.SetObject();
      Document::AllocatorType& allocator = document.GetAllocator();
      
      // add information
      document.AddMember("NormalizationFactor", normalizationFactor, allocator);
      jsonAddMember(document, allocator, "CylinderCenter", lps2ras(uniformityMeasurements->GetCenter()));
      jsonAddMember(document, allocator, "CylinderDirection", lps2ras(uniformityMeasurements->GetDirection()));
      document.AddMember("CylinderRadius", uniformityMeasurements->GetRadius(), allocator);
      document.AddMember("CylinderHeight", uniformityMeasurements->GetHeight(), allocator);
      document.AddMember("CylinderMean", uniformityMeasurements->GetCylinderMean(), allocator);
      document.AddMember("CylinderStd", uniformityMeasurements->GetCylinderStd(), allocator);
      document.AddMember("MaxRelDiff", uniformityMeasurements->GetMaxRelativeDifference(), allocator);
      jsonAddMember(document, allocator, "SliceMeasurements",
        uniformityMeasurements->GetSliceMeasurements()->CastToSTLConstContainer());
      jsonAddMember(document, allocator, "SliceOffsets",
        uniformityMeasurements->GetSliceOffsets()->CastToSTLConstContainer());
      
      // write JSON document to file
      std::ofstream outputFile(measurementsData.c_str());
      OStreamWrapper osw(outputFile);
      PrettyWriter<OStreamWrapper> writer(osw);
      document.Accept(writer);
    }
        
  }
  catch(const itk::ExceptionObject & excep )
  {
    std::cerr << argv[0] << ": exception caught!" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
