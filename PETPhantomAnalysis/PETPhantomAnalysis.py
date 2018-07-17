import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import getpass
from datetime import datetime
import json
import shutil
import tempfile
import dicom

#
# PETPhantomAnalysis
#

class PETPhantomAnalysis(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "PET Cylinder Phantom Analysis"
    self.parent.categories = ["Quantification"]
    self.parent.dependencies = []
    self.parent.contributors = ["Christian Bauer (University of Iowa)"]
    self.parent.helpText = """
    Measurement of calibration and uniformity in a cylinder phantom PET scan. \
    <a href="https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/PETCPhantomAnalysis">Documentation.</a>
    """
    self.parent.acknowledgementText = """
    This work was partially funded by NIH grants U01-CA140206 and U24-CA180918.
""" # replace with organization, grant and thanks.

#
# PETPhantomAnalysisWidget
#

class PETPhantomAnalysisWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self,parent=None):
    ScriptedLoadableModuleWidget.__init__(self, parent)
    self.slicerTempDir = slicer.util.tempDirectory()

  def setup(self):
    self.measurementsLogic = PETPhantomAnalysisLogic()

    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "PET Cylinder Phantom Analysis"
    self.layout.addWidget(parametersCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    #self.inputSelector.addAttribute("vtkMRMLScalarVolumeNode", "DICOM.instanceUIDs", None) # add this to restrict to dicom datasets
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Input SUVbw normalized DICOM PET volume")
    self.inputSelector.connect("currentNodeChanged(bool)",self.refreshUIElements)
    self.inputSelector.connect("currentNodeChanged(bool)",self.inputVolumeSelected)
    parametersFormLayout.addRow("Input Volume", self.inputSelector)
    
    self.normalizationFactorLineEdit = qt.QLineEdit("1.0")
    self.normalizationFactorValidator = qt.QDoubleValidator()
    self.normalizationFactorValidator.bottom = 0.0
    self.normalizationFactorLineEdit.setValidator(self.normalizationFactorValidator)
    self.normalizationFactorLineEdit.setToolTip( "Normalization factor for input dataset")
    self.normalizationFactorLineEdit.connect("textChanged(QString)",self.refreshUIElements)
    parametersFormLayout.addRow("Normalization Factor",self.normalizationFactorLineEdit)

    self.segmentationSelector = slicer.qMRMLNodeComboBox()
    self.segmentationSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.segmentationSelector.selectNodeUponCreation = True
    self.segmentationSelector.addEnabled = True
    self.segmentationSelector.removeEnabled = True
    self.segmentationSelector.noneEnabled = False
    self.segmentationSelector.showHidden = False
    self.segmentationSelector.showChildNodeTypes = False
    self.segmentationSelector.setMRMLScene( slicer.mrmlScene )
    self.segmentationSelector.setToolTip( "Output cylinder reference region volume")
    self.segmentationSelector.connect("currentNodeChanged(bool)",self.refreshUIElements)
    parametersFormLayout.addRow("Output Volume", self.segmentationSelector)

    self.segmentButton = qt.QPushButton("Analyze cylinder phantom")
    self.segmentButton.toolTip = "Analyze cylinder phantom and show results"
    self.segmentButton.connect('clicked(bool)', self.onApplyButton)
    parametersFormLayout.addRow("Analysis",self.segmentButton)

    # Normalization form
    normalizationCollapsibleButton = ctk.ctkCollapsibleButton()
    normalizationCollapsibleButton.text = "Normalization Metadata"
    self.layout.addWidget(normalizationCollapsibleButton)
    normalizationFormLayout = qt.QFormLayout(normalizationCollapsibleButton)

    self.inputTypeLabel = qt.QLabel('unspecified')
    normalizationFormLayout.addRow("Input Volume Type",self.inputTypeLabel)

    self.halfLifeLineEdit = qt.QLineEdit(float('inf'))
    self.halfLifeLineEditValidator = qt.QDoubleValidator()
    self.halfLifeLineEditValidator.bottom = 0.0
    self.halfLifeLineEdit.setValidator(self.halfLifeLineEditValidator)
    self.halfLifeLineEdit.setToolTip( "Radionuclide half life time")
    self.halfLifeLineEdit.connect("textChanged(QString)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow("Radionuclide Half Life (s)",self.halfLifeLineEdit)
    
    self.activityLineEdit = qt.QLineEdit("1.0")
    self.activityLineEditValidator = qt.QDoubleValidator()
    self.activityLineEditValidator.bottom = 0.0
    self.activityLineEdit.setValidator(self.activityLineEditValidator)
    self.activityLineEdit.setToolTip( "Normalization factor for input dataset")
    self.activityLineEdit.connect("textChanged(QString)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow(u"Syringe Activity (\u00B5Ci)",self.activityLineEdit)
    
    self.activityTimeEdit = qt.QTimeEdit()
    self.activityTimeEdit.connect("timeChanged(QTime)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow("Syringe Assay Time",self.activityTimeEdit)
    
    self.residualActivityLineEdit = qt.QLineEdit("0.0")
    self.residualActivityLineEditValidator = qt.QDoubleValidator()
    self.residualActivityLineEditValidator.bottom = 0.0
    self.residualActivityLineEdit.setValidator(self.residualActivityLineEditValidator)
    self.residualActivityLineEdit.setToolTip( "Normalization factor for input dataset")
    self.residualActivityLineEdit.connect("textChanged(QString)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow(u"Syringe Residual Activity (\u00B5Ci)",self.residualActivityLineEdit)
    
    self.residualactivityTimeEdit = qt.QTimeEdit()
    self.residualactivityTimeEdit.connect("timeChanged(QTime)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow("Syringe Residual Activity Time",self.residualactivityTimeEdit)
    
    self.scanTimeEdit = qt.QTimeEdit()
    self.scanTimeEdit.connect("timeChanged(QTime)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow("Scan Time",self.scanTimeEdit)
    
    self.weightLineEdit = qt.QLineEdit("1000.0")
    self.weightLineEditValidator = qt.QDoubleValidator()
    self.weightLineEditValidator.bottom = 0.0
    self.weightLineEdit.setValidator(self.weightLineEditValidator)
    self.weightLineEdit.setToolTip( "Normalization factor for input dataset")
    self.weightLineEdit.connect("textChanged(QString)",self.updateNormalizationFactor)
    normalizationFormLayout.addRow("Fill Weight (g)",self.weightLineEdit)
    
    self.inputVolumeSUVNormalizationLabel = qt.QLabel("1.0")
    self.inputVolumeSUVNormalizationLabel.setToolTip( "Normalization factor applied to input volume")
    normalizationFormLayout.addRow("Input Volume Normalization Factor",self.inputVolumeSUVNormalizationLabel)

    # Measurement results
    measurementsCollapsibleButton = ctk.ctkCollapsibleButton()
    measurementsCollapsibleButton.text = "Measurements"
    self.layout.addWidget(measurementsCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(measurementsCollapsibleButton)

    self.meanValueLineEdit = qt.QLineEdit("-")
    self.meanValueLineEdit.setReadOnly(True)
    self.meanValueLineEdit.setToolTip( "Mean value in measurement region")
    parametersFormLayout.addRow("Mean",self.meanValueLineEdit)
    self.stdValueLineEdit = qt.QLineEdit("-")
    self.stdValueLineEdit.setReadOnly(True)
    self.stdValueLineEdit.setToolTip( "Standard deviation in measurement region")
    parametersFormLayout.addRow("Standard Deviation",self.stdValueLineEdit)
    self.maxRelDiffValueLineEdit = qt.QLineEdit("-")
    self.maxRelDiffValueLineEdit.setReadOnly(True)
    self.maxRelDiffValueLineEdit.setToolTip( "Maximum relative difference  in axial direction in reference region")
    parametersFormLayout.addRow("Maximum Relative Difference",self.maxRelDiffValueLineEdit)

    # visualize matched cylinder
    self.cylinderModelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode")
    self.cylinderModelNode.CreateDefaultDisplayNodes()
    self.cylinderModelNode.SetName("MatchedPhantomCylinderModelNode")
    self.cylinderModelNode.GetDisplayNode().HideFromEditorsOn()
    self.cylinderModelNode.GetDisplayNode().SetColor(0,1,0)
    self.cylinderModelNode.GetDisplayNode().SliceIntersectionVisibilityOn()
    self.cylinderModelNode.HideFromEditorsOn()

    # Add vertical spacer
    self.layout.addStretch(1)

    self.inputVolumeSelected()
    self.refreshUIElements()

  def enter(self):
    if self.cylinderModelNode: self.cylinderModelNode.GetDisplayNode().VisibilityOn()
    if self.cylinderModelNode: self.cylinderModelNode.GetDisplayNode().SetSliceIntersectionVisibility(True)
    self.inputVolumeSelected()

  def exit(self):
    if self.cylinderModelNode: self.cylinderModelNode.GetDisplayNode().VisibilityOff()
    if self.cylinderModelNode: self.cylinderModelNode.GetDisplayNode().SetSliceIntersectionVisibility(False)
    pass

  def inputVolumeSelected(self):
    self.inputTypeLabel.text = 'unspecified'
    self.halfLifeLineEdit.text = float('inf')
    self.activityLineEdit.text = 1 
    self.activityTimeEdit.time = qt.QTime(0,0)
    self.residualActivityLineEdit.text = 0
    self.residualactivityTimeEdit.time = qt.QTime(0,0)
    self.scanTimeEdit.time = qt.QTime(0,0)
    self.weightLineEdit.text = 1000.0
    self.inputVolumeSUVNormalizationLabel.text = 1.0

    inputVolume = self.inputSelector.currentNode()
    if not inputVolume: return
    dicomInstanceUIDs = inputVolume.GetAttribute('DICOM.instanceUIDs')
    rwvmInstanceUID = inputVolume.GetAttribute('DICOM.RWV.instanceUID')
    if dicomInstanceUIDs:
      self.inputTypeLabel.text = "DICOM volume"
      sourceFileName = slicer.dicomDatabase.fileForInstance(dicomInstanceUIDs.split(" ")[0])
      d = dicom.read_file(sourceFileName) if sourceFileName else None
      if d and d.Modality=='PT':
        self.inputTypeLabel.text = "PET DICOM volume"
        self.halfLifeLineEdit.text = d.RadiopharmaceuticalInformationSequence[0].RadionuclideHalfLife
        self.activityLineEdit.text = self._doseInMicroCi(d.RadiopharmaceuticalInformationSequence[0].RadionuclideTotalDose, d.Units)
        self.activityTimeEdit.time = self._getTime(d.RadiopharmaceuticalInformationSequence[0].RadiopharmaceuticalStartTime)
        self.residualactivityTimeEdit.time = self.activityTimeEdit.time
        self.scanTimeEdit.time = self._getTime(d.SeriesTime)
        self.weightLineEdit.text = d.PatientWeight*1000.0
        if inputVolume.GetVoxelValueUnits() and inputVolume.GetVoxelValueUnits().GetCodeValue()=='{SUVbw}g/ml': # calculate SUV normalization factor based on PET DICOM
          self.inputTypeLabel.text = "SUVbw normalized PET DICOM volume"
          halfLife = d.RadiopharmaceuticalInformationSequence[0].RadionuclideHalfLife
          injectedDose = self._doseInMicroCi(d.RadiopharmaceuticalInformationSequence[0].RadionuclideTotalDose, d.Units)
          weight = d.PatientWeight
          decayTime = 0.001*(\
            self._getTime(d.SeriesTime).msecsSinceStartOfDay()-\
            self._getTime(d.RadiopharmaceuticalInformationSequence[0].RadiopharmaceuticalStartTime).msecsSinceStartOfDay())
          decayedDose = injectedDose * pow(2.0, -decayTime/halfLife)
          Ci2BqFactor = 37000000000.0
          decayedDose = decayedDose * Ci2BqFactor * 1e-9 # convert to kBq/mL
          SUVNormalizationFactor = weight / decayedDose
          self.inputVolumeSUVNormalizationLabel.text = SUVNormalizationFactor
        elif rwvmInstanceUID: # load SUV normalization factor from RWVM file
          self.inputTypeLabel.text = "SUV normalized PET DICOM volume"
          sourceFileName = slicer.dicomDatabase.fileForInstance(rwvmInstanceUID)
          d = dicom.read_file(sourceFileName)
          self.inputVolumeSUVNormalizationLabel.text = d.ReferencedImageRealWorldValueMappingSequence[0].RealWorldValueMappingSequence[0].RealWorldValueSlope
    self.updateNormalizationFactor()

  def updateNormalizationFactor(self):
    if self.inputTypeLabel.text == 'unspecified':
      return
    try:
      halfLife = float(self.halfLifeLineEdit.text)
      seriesTime = float(self.scanTimeEdit.time.msecsSinceStartOfDay())/1000.0
      activity = float(self.activityLineEdit.text)
      assayTime = float(self.activityTimeEdit.time.msecsSinceStartOfDay())/1000.0
      residualActivity = float(self.residualActivityLineEdit.text)
      residualActivityTime =  float(self.residualactivityTimeEdit.time.msecsSinceStartOfDay())/1000.0
      weight = float(self.weightLineEdit.text)/1000.0
      t1 = residualActivityTime-assayTime
      injectedDose = activity * pow(2.0, -t1/halfLife) - residualActivity
      t2 = seriesTime - residualActivityTime
      decayedDose = injectedDose * pow(2.0, -t2/halfLife)
      Ci2BqFactor = 37000000000.0
      decayedDose = decayedDose * Ci2BqFactor * 1e-9 # convert to kBq/mL
      SUVNormalizationFactor = weight / decayedDose
      normalizationFactor = SUVNormalizationFactor/float(self.inputVolumeSUVNormalizationLabel.text)
      self.normalizationFactorLineEdit.text = normalizationFactor
    except:
      pass

  def _getTime(self, time):
    dot = time.find('.') if time.find('.')>=0 else len(time)
    hours = int(time[:dot-4])
    minutes = int(time[dot-4:dot-2])
    seconds = int(round(float(time[dot-2:dot])))
    mseconds = int(round(float(time[dot:-1])*1000.0)) if dot<len(time) else 0
    return qt.QTime(hours, minutes, seconds, mseconds)

  def _doseInMicroCi(self, dose, units='BQML'):    
    Ci2BqFactor = 37000000000.0
    magnitudeFactors = {"M":1000.0, "k":1.0, "m":1e-6, "u":1e-9}
    if units[-2:].lower!='ci':
      dose = dose/Ci2BqFactor
    if len(units)>2 and units[0] in magnitudeFactors:
      dose = dose*magnitudeFactors[units[0]]
    return dose*1e6

  def refreshUIElements(self):
    self.meanValueLineEdit.text = "-"
    self.stdValueLineEdit.text = "-"
    self.maxRelDiffValueLineEdit.text = "-"
    self.segmentButton.enabled = self.inputSelector.currentNode()!=None

  def onApplyButton(self):
    inputVolume = self.inputSelector.currentNode()
    outputVolume = self.segmentationSelector.currentNode()
    jsonFile = tempfile.NamedTemporaryFile()
    normalizationFactor = 1.0
    try: 
      normalizationFactor = float(self.normalizationFactorLineEdit.text)
    except:
     pass

    cliParams = {'inputVolume': inputVolume.GetID(), \
                 'normalizationFactor': normalizationFactor, \
                 'measurementsData': jsonFile.name,
                 }
    if outputVolume: cliParams['outputVolume'] =  outputVolume.GetID()

    pd = qt.QProgressDialog('Running PET Cylinder Phantom Analysis...', 'Cancel', 0, 100, slicer.util.mainWindow())
    pd.setModal(True)
    pd.setMinimumDuration(0)
    pd.show()
    pd.setValue(30)
    cliNode = None
    cliNode = slicer.cli.run(slicer.modules.petphantomanalysiscli, cliNode, cliParams, wait_for_completion=False)
    while cliNode.IsBusy():
      slicer.app.processEvents()
      if pd.wasCanceled:
        cliNode.Cancel()
    pd.setValue(100)
    if pd.wasCanceled:
      return

    if cliNode.GetStatusString() != 'Completed':
      qt.QMessageBox().warning(None,"Warning","Analysis of PET cylinder phantom was not successful.")
      return

    for i in range(cliNode.GetNumberOfParametersInGroup(1)):
      name = cliNode.GetParameterName(1,i)
      value = cliNode.GetParameterDefault(1,i)
      if name=='Mean_s': self.meanValueLineEdit.setText(value)
      if name=='Std_s': self.stdValueLineEdit.setText(value)
      if name=='MaxRelDiff_s': self.maxRelDiffValueLineEdit.setText(value)

    if self.meanValueLineEdit.text=='--':
      qt.QMessageBox().warning(None,"Warning","Analysis of PET cylinder phantom was not successful.")
      return

    try:
      measurements = json.load(open(jsonFile.name))
    except:
      qt.QMessageBox().warning(None,"Warning","Analysis of PET cylinder phantom was not successful.")
      return

    # plot normalized slice measurements
    mean = measurements['CylinderMean']
    sliceOffsets = measurements['SliceOffsets']
    sliceMeasurements = measurements['SliceMeasurements']
    normalizedSliceMeasurements = [v/mean for v in sliceMeasurements]
    self._createPlot(normalizedSliceMeasurements, sliceOffsets)
    
    # visualize matched cylinder
    self.cylinderModelNode.SetAndObserveMesh(self._createCylinderMesh(
      measurements['CylinderCenter'], measurements['CylinderDirection']))   

  def _createPlot(self, normalizedSliceMeasurements, sliceOffsets):    
    tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
    table = tableNode.GetTable()

    arrX = vtk.vtkFloatArray()
    arrX.SetName("offset")
    table.AddColumn(arrX)

    arrY = vtk.vtkFloatArray()
    arrY.SetName("intensity")
    table.AddColumn(arrY)
    
    numPoints = len(normalizedSliceMeasurements)
    table.SetNumberOfRows(numPoints)
    for i in range(numPoints):
      table.SetValue(i, 0, sliceOffsets[i] )
      table.SetValue(i, 1, normalizedSliceMeasurements[i])

    plotSeriesNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", "Axial uniformity")
    plotSeriesNode.SetAndObserveTableNodeID(tableNode.GetID())
    plotSeriesNode.SetXColumnName("offset")
    plotSeriesNode.SetYColumnName("intensity")
    plotSeriesNode.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
    plotSeriesNode.SetUniqueColor()

    # Create plot chart node
    plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")
    plotChartNode.AddAndObservePlotSeriesNodeID(plotSeriesNode.GetID())
    plotChartNode.SetTitle('Mean normalized intensity profile')
    plotChartNode.SetXAxisTitle('Offset from center in axial direction (mm)')
    plotChartNode.SetYAxisTitle('Mean normalized intensity')

    # Switch to a layout that contains a plot view to create a plot widget
    layoutManager = slicer.app.layoutManager()
    layoutWithPlot = slicer.modules.plots.logic().GetLayoutWithPlot(layoutManager.layout)
    layoutManager.setLayout(layoutWithPlot)

    # Select chart in plot view
    plotWidget = layoutManager.plotWidget(0)
    plotViewNode = plotWidget.mrmlPlotViewNode()
    plotViewNode.SetPlotChartNodeID(plotChartNode.GetID())
    
  def _createCylinderMesh(self, center, direction, radius=100.0, height=200.0):
    startPoint = [center[i]-direction[i]*height/2.0 for i in range(3)]
    cylinderSource = vtk.vtkCylinderSource()
    cylinderSource.SetResolution(365)
    cylinderSource.SetRadius(radius)
    cylinderSource.Update()
    normalizedX = direction
    normalizedY = [0,0,0]
    normalizedZ = [0,0,0]
    vtk.vtkMath.Normalize(normalizedX)
    arbitrary = [1,1,1]
    vtk.vtkMath.Cross(normalizedX, arbitrary, normalizedZ)
    vtk.vtkMath.Normalize(normalizedZ)
    vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)
    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    for i in range(3):
      matrix.SetElement(i, 0, normalizedX[i])
      matrix.SetElement(i, 1, normalizedY[i])
      matrix.SetElement(i, 2, normalizedZ[i])
    transform = vtk.vtkTransform()
    transform.Translate(startPoint)
    transform.Concatenate(matrix)
    transform.RotateZ(-90.0)
    transform.Scale(1.0, height, 1.0)
    transform.Translate(0, .5, 0)
    transformPD = vtk.vtkTransformPolyDataFilter()
    transformPD.SetTransform(transform)
    transformPD.SetInputConnection(cylinderSource.GetOutputPort())
    transformPD.Update()
    return transformPD.GetOutputDataObject(0)


#
# PETPhantomAnalysisLogic
#

class PETPhantomAnalysisLogic(ScriptedLoadableModuleLogic):

  def __init__(self):
    pass

  def __del__(self):
    pass

#
# PETPhantomAnalysisTest
#

from DICOMLib import DICOMUtils
class PETPhantomAnalysisTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_PETPhantomAnalysis()
    self.tearDown()

  def setUp(self):
    """ Open temporary DICOM database
    """
    slicer.mrmlScene.Clear(0)
    self.delayMs = 700
    self.tempDataDir = os.path.join(slicer.app.temporaryPath,'PETTest')
    self.tempDicomDatabaseDir = os.path.join(slicer.app.temporaryPath,'PETTestDicom')

  def doCleanups(self):
    """ cleanup temporary data in case an exception occurs
    """
    self.tearDown()

  def tearDown(self):
    """ Close temporary DICOM database and remove temporary data
    """
    try:
      import shutil
      if os.path.exists(self.tempDataDir):
        shutil.rmtree(self.tempDataDir)
    except Exception, e:
      import traceback
      traceback.print_exc()
      self.delayDisplay('Test caused exception!\n' + str(e),self.delayMs*2)

  def loadTestData(self):
    self.patienName = 'UNIFORMITY^Bio-mCT'
    #download data and add to dicom database
    zipFileUrl = 'https://github.com/QIICR/SlicerPETPhantomAnalysis/releases/download/test-data/PETCylinderPhantom.zip'
    zipFilePath = self.tempDataDir+'/dicom.zip'
    zipFileData = self.tempDataDir+'/dicom'
    expectedNumOfFiles = 171
    if not os.access(self.tempDataDir, os.F_OK):
      os.mkdir(self.tempDataDir)
    if not os.access(zipFileData, os.F_OK):
      os.mkdir(zipFileData)

    dicomWidget = slicer.modules.dicom.widgetRepresentation().self()
    dicomPluginCheckbox =  dicomWidget.detailsPopup.pluginSelector.checkBoxByPlugin
    dicomPluginStates = {(key,value.checked) for key,value in dicomPluginCheckbox.iteritems()}
    for cb in dicomPluginCheckbox.itervalues(): cb.checked=False
    dicomPluginCheckbox['DICOMScalarVolumePlugin'].checked = True

    # Download, unzip, import, and load data. Verify loaded nodes.
    loadedNodes = {'vtkMRMLScalarVolumeNode':1}
    with DICOMUtils.LoadDICOMFilesToDatabase(zipFileUrl, zipFilePath, zipFileData, expectedNumOfFiles, {}, loadedNodes) as success:
      self.assertTrue(success)

    self.assertEqual( len( slicer.util.getNodes('vtkMRMLSubjectHierarchyNode*') ), 1 )
    imageNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLScalarVolumeNode')

    for key,value in dicomPluginStates:
      dicomPluginCheckbox[key].checked=value

    # apply the SUVbw conversion factor and set units and quantity
    suvNormalizationFactor = 0.00012595161151
    quantity = slicer.vtkCodedEntry()
    quantity.SetFromString('CodeValue:126400|CodingSchemeDesignator:DCM|CodeMeaning:Standardized Uptake Value')
    units = slicer.vtkCodedEntry()
    units.SetFromString('CodeValue:{SUVbw}g/ml|CodingSchemeDesignator:UCUM|CodeMeaning:Standardized Uptake Value body weight')
    multiplier = vtk.vtkImageMathematics()
    multiplier.SetOperationToMultiplyByK()
    multiplier.SetConstantK(suvNormalizationFactor)
    multiplier.SetInput1Data(imageNode.GetImageData())
    multiplier.Update()
    imageNode.GetImageData().DeepCopy(multiplier.GetOutput())
    imageNode.GetVolumeDisplayNode().SetWindowLevel(6,3)
    imageNode.GetVolumeDisplayNode().SetAndObserveColorNodeID('vtkMRMLColorTableNodeInvertedGrey')
    imageNode.SetVoxelValueQuantity(quantity)
    imageNode.SetVoxelValueUnits(units)

    return imageNode

  def test_PETPhantomAnalysis(self):
    """ test standard processing
    """ 
    try:
      self.assertIsNotNone( slicer.modules.petphantomanalysis )
      with DICOMUtils.TemporaryDICOMDatabase(self.tempDicomDatabaseDir) as db:
        self.assertTrue(db.isOpen)
        self.assertEqual(slicer.dicomDatabase, db)

        self.delayDisplay('Loading PET DICOM dataset (including download if necessary)')
        petNode = self.loadTestData()

        self.delayDisplay('Running segmentation')
        m = slicer.util.mainWindow()
        m.moduleSelector().selectModule('PETPhantomAnalysis')
        qrWidget = slicer.modules.PETPhantomAnalysisWidget
        qrWidget.inputSelector.setCurrentNode(petNode)
        segmentationNode = qrWidget.segmentationSelector.addNode()
        qrWidget.inputVolumeSelected()
        qrWidget.segmentButton.click()

        # assert measurements are correct
        self.assertTrue(abs(float(qrWidget.meanValueLineEdit.text)-0.982502)<0.01)
        self.assertTrue(abs(float(qrWidget.stdValueLineEdit.text)-0.031612)<0.01)
        self.assertTrue(abs(float(qrWidget.maxRelDiffValueLineEdit.text)+0.0203663)<0.01)

        # clean up data from DICOM database
        import dicom
        patientUID = DICOMUtils.getDatabasePatientUIDByPatientName(self.patienName)
        db.removePatient(patientUID)

        self.delayDisplay('Test passed!')

    except Exception, e:
      import traceback
      traceback.print_exc()
      self.delayDisplay('Test caused exception!\n' + str(e),self.delayMs*2)
    
