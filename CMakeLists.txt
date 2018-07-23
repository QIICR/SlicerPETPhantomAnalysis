cmake_minimum_required(VERSION 2.8.9)

project(PETCPhantom)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/PETPhantomAnalysis")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Christian Bauer (University of Iowa), Ethan Ulrich (University of Iowa), Andrey Fedorov (SPL), Reinhard R. Beichel (University of Iowa), John Buatti (University of Iowa)")
set(EXTENSION_DESCRIPTION "The PET Phantom Analysis Extension allows measurement of calibration and uniformity in a cylinder phantom PET scan.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/QIICR/SlicerPETPhantomAnalysis/master/PETCphantom.png")
set(EXTENSION_SCREENSHOTURLS "https://www.slicer.org/w/images/thumb/0/01/PETCylinderPhantomAnalysisModule.png/800px-PETCylinderPhantomAnalysisModule.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#find_package(DCMQI REQUIRED)
if(CMAKE_CONFIGURATION_TYPES)
  set(DEPENDENCIES_ADDITIONAL_MODULE_PATHS "")
  foreach(config ${CMAKE_CONFIGURATION_TYPES})
    list(APPEND DEPENDENCIES_ADDITIONAL_MODULE_PATHS "${CMAKE_BINARY_DIR}/${Slicer_CLIMODULES_LIB_DIR}/${config}")
  endforeach()
else()
  set(DEPENDENCIES_ADDITIONAL_MODULE_PATHS
    ${CMAKE_BINARY_DIR}/${Slicer_CLIMODULES_LIB_DIR}
    )
endif()


#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PETPhantomAnalysisCLI)
add_subdirectory(PETPhantomAnalysis)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
