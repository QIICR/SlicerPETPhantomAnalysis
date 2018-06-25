cmake_minimum_required(VERSION 2.8.9)

project(PETCPhantom)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/PETCPhantom")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Christian Bauer (University of Iowa), Ethan Ulrich (University of Iowa), Andrey Fedorov (SPL), Reinhard R. Beichel (University of Iowa), John Buatti (University of Iowa)")
set(EXTENSION_DESCRIPTION "The PETCPhantom Extension allows measurement of calibration and uniformity in a cylinder phantom PET scan.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/QIICR/PETPhantomAnalysis/master/PETCPhantom.png")
set(EXTENSION_SCREENSHOTURLS "")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PETPhantomAnalysisCLI)
add_subdirectory(PETPhantomAnalysis)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
