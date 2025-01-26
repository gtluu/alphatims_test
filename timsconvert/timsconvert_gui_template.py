from PySide6.QtCore import QCoreApplication, QMetaObject, QRect, QSize
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (QCheckBox, QLabel, QLineEdit, QPushButton, QRadioButton, QSizePolicy, QSpinBox, QWidget,
                               QButtonGroup, QTableWidget, QHeaderView, QTableWidgetItem, QAbstractItemView, QMenu,
                               QMenuBar)


class Ui_TimsconvertGuiWindow(object):
    def setupUi(self, TimsconvertGuiWindow):
        # Main Window
        if not TimsconvertGuiWindow.objectName():
            TimsconvertGuiWindow.setObjectName(u"TimsconvertGuiWindow")
        TimsconvertGuiWindow.resize(800, 960)
        TimsconvertGuiWindow.setMinimumSize(QSize(0, 0))
        TimsconvertGuiWindow.setMaximumSize(QSize(800, 960))

        # Main Widget
        self.main = QWidget(TimsconvertGuiWindow)
        self.main.setObjectName(u"main")

        # Menu Bar
        self.MenuBar = QMenuBar(TimsconvertGuiWindow)
        self.MenuBar.setObjectName(u"MenuBar")
        self.MenuBar.setGeometry(QRect(0, 0, 800, 22))
        self.About = QMenu(self.MenuBar)
        self.About.setObjectName(u"About")
        self.AboutTimsconvert = QAction(TimsconvertGuiWindow)
        self.AboutTimsconvert.setObjectName(u"AboutTimsconvert")

        # Input
        # Button
        self.AddInputDialogueButton = QPushButton(self.main)
        self.AddInputDialogueButton.setObjectName(u"AddInputDialogueButton")
        self.AddInputDialogueButton.setGeometry(QRect(20, 810, 381, 24))
        # Button
        self.RemoveFromQueueButton = QPushButton(self.main)
        self.RemoveFromQueueButton.setObjectName(u"RemoveFromQueueButton")
        self.RemoveFromQueueButton.setGeometry(QRect(20, 850, 381, 24))
        # List
        self.InputList = QTableWidget(self.main)
        self.InputList.setObjectName(u"InputList")
        self.InputList.setGeometry(QRect(20, 20, 381, 771))
        self.InputList.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))
        # Set the resizing mode for the horizontal header to stretch
        self.InputList.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # Turn off row index display
        self.InputList.verticalHeader().setVisible(False)
        self.InputList.horizontalHeader().setVisible(True)
        self.InputList.setColumnCount(1)
        self.InputList.setHorizontalHeaderLabels(['Queue'])
        self.InputList.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.InputList.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        # Populate the table with some sample data
        self.InputList.setRowCount(0)
        old_row_count = self.InputList.rowCount()
        self.InputList.setRowCount(self.InputList.rowCount() + len([]))
        for row, input_filename in enumerate([], start=old_row_count):
            text_item = QTableWidgetItem(input_filename)
            self.InputList.setItem(row, 0, text_item)

        # Output Directory
        # Label
        self.OutputDirectoryLabel = QLabel(self.main)
        self.OutputDirectoryLabel.setObjectName(u"OutputDirectoryLabel")
        self.OutputDirectoryLabel.setGeometry(QRect(430, 20, 341, 20))
        # Line Edit
        self.OutputDirectoryLine = QLineEdit(self.main)
        self.OutputDirectoryLine.setObjectName(u"OutputDirectoryLine")
        self.OutputDirectoryLine.setGeometry(QRect(430, 50, 271, 22))
        self.OutputDirectoryLine.setMaxLength(32765)
        self.OutputDirectoryLine.setReadOnly(True)
        # Button
        self.OutputDirectoryBrowseButton = QPushButton(self.main)
        self.OutputDirectoryBrowseButton.setObjectName(u"OutputDirectoryBrowseButton")
        self.OutputDirectoryBrowseButton.setGeometry(QRect(710, 50, 61, 24))

        # Spectra Export Mode
        # Label
        self.ModeLabel = QLabel(self.main)
        self.ModeLabel.setObjectName(u"ModeLabel")
        self.ModeLabel.setGeometry(QRect(430, 89, 341, 21))
        # Profile
        self.ModeProfileRadio = QRadioButton(self.main)
        self.ModeProfileRadio.setObjectName(u"ModeProfileRadio")
        self.ModeProfileRadio.setGeometry(QRect(430, 120, 89, 20))
        # Centroid
        self.ModeCentroidRadio = QRadioButton(self.main)
        self.ModeCentroidRadio.setObjectName(u"ModeCentroidRadio")
        self.ModeCentroidRadio.setGeometry(QRect(430, 150, 89, 20))
        # Raw
        self.ModeRawRadio = QRadioButton(self.main)
        self.ModeRawRadio.setObjectName(u"ModeRawRadio")
        self.ModeRawRadio.setGeometry(QRect(430, 180, 89, 20))
        # Group
        self.ModeGroup = QButtonGroup()
        self.ModeGroup.addButton(self.ModeProfileRadio)
        self.ModeGroup.addButton(self.ModeCentroidRadio)
        self.ModeGroup.addButton(self.ModeRawRadio)
        self.ModeCentroidRadio.setChecked(True)

        # Encoding
        # Label
        self.EncodingLabel = QLabel(self.main)
        self.EncodingLabel.setObjectName(u"EncodingLabel")
        self.EncodingLabel.setGeometry(QRect(430, 220, 351, 21))
        # m/z Encoding
        # Label
        self.MZEncodingLabel = QLabel(self.main)
        self.MZEncodingLabel.setObjectName(u"MZEncodingLabel")
        self.MZEncodingLabel.setGeometry(QRect(430, 250, 351, 21))
        # 64
        self.MZEncoding64Radio = QRadioButton(self.main)
        self.MZEncoding64Radio.setObjectName(u"MZEncoding64Radio")
        self.MZEncoding64Radio.setGeometry(QRect(430, 280, 89, 21))
        # 32
        self.MZEncoding32Radio = QRadioButton(self.main)
        self.MZEncoding32Radio.setObjectName(u"MZEncoding32Radio")
        self.MZEncoding32Radio.setGeometry(QRect(550, 280, 89, 21))
        # Group
        self.MZEncodingGroup = QButtonGroup()
        self.MZEncodingGroup.addButton(self.MZEncoding64Radio)
        self.MZEncodingGroup.addButton(self.MZEncoding32Radio)
        self.MZEncoding64Radio.setChecked(True)
        # Intensity Encoding
        # Label
        self.IntensityEncodingLabel = QLabel(self.main)
        self.IntensityEncodingLabel.setObjectName(u"IntensityEncodingLabel")
        self.IntensityEncodingLabel.setGeometry(QRect(430, 310, 351, 21))
        # 64
        self.IntensityEncoding64Radio = QRadioButton(self.main)
        self.IntensityEncoding64Radio.setObjectName(u"IntensityEncoding64Radio")
        self.IntensityEncoding64Radio.setGeometry(QRect(430, 340, 89, 21))
        # 32
        self.IntensityEncoding32Radio = QRadioButton(self.main)
        self.IntensityEncoding32Radio.setObjectName(u"IntensityEncoding32Radio")
        self.IntensityEncoding32Radio.setGeometry(QRect(550, 340, 89, 21))
        # Group
        self.IntensityEncodingGroup = QButtonGroup()
        self.IntensityEncodingGroup.addButton(self.IntensityEncoding64Radio)
        self.IntensityEncodingGroup.addButton(self.IntensityEncoding32Radio)
        self.IntensityEncoding64Radio.setChecked(True)
        # Mobility Encoding
        # Label
        self.MobilityEncodingLabel = QLabel(self.main)
        self.MobilityEncodingLabel.setObjectName(u"MobilityEncodingLabel")
        self.MobilityEncodingLabel.setGeometry(QRect(430, 370, 351, 21))
        # 64
        self.MobilityEncoding64Radio = QRadioButton(self.main)
        self.MobilityEncoding64Radio.setObjectName(u"MobilityEncoding64Radio")
        self.MobilityEncoding64Radio.setGeometry(QRect(430, 400, 89, 21))
        # 32
        self.MobilityEncoding32Radio = QRadioButton(self.main)
        self.MobilityEncoding32Radio.setObjectName(u"MobilityEncoding32Radio")
        self.MobilityEncoding32Radio.setGeometry(QRect(550, 400, 89, 21))
        # Group
        self.MobilityEncodingGroup = QButtonGroup()
        self.MobilityEncodingGroup.addButton(self.MobilityEncoding64Radio)
        self.MobilityEncodingGroup.addButton(self.MobilityEncoding32Radio)
        self.MobilityEncoding64Radio.setChecked(True)

        # Mobility Compensation
        # Label
        self.MobilityCompensationLabel = QLabel(self.main)
        self.MobilityCompensationLabel.setObjectName(u"MobilityCompensationLabel")
        self.MobilityCompensationLabel.setGeometry(QRect(430, 430, 351, 21))
        # None
        self.MobilityCompensationNoneRadio = QRadioButton(self.main)
        self.MobilityCompensationNoneRadio.setObjectName(u"MobilityCompensationNoneRadio")
        self.MobilityCompensationNoneRadio.setGeometry(QRect(430, 459, 89, 21))
        # Global
        self.MobilityCompensationGlobalRadio = QRadioButton(self.main)
        self.MobilityCompensationGlobalRadio.setObjectName(u"MobilityCompensationGlobalRadio")
        self.MobilityCompensationGlobalRadio.setGeometry(QRect(550, 460, 89, 21))
        # Frame
        self.MobilityCompensationFrameRadio = QRadioButton(self.main)
        self.MobilityCompensationFrameRadio.setObjectName(u"MobilityCompensationFrameRadio")
        self.MobilityCompensationFrameRadio.setGeometry(QRect(670, 460, 89, 21))
        # Group
        self.MobilityCompensationGroup = QButtonGroup()
        self.MobilityCompensationGroup.addButton(self.MobilityCompensationNoneRadio)
        self.MobilityCompensationGroup.addButton(self.MobilityCompensationGlobalRadio)
        self.MobilityCompensationGroup.addButton(self.MobilityCompensationFrameRadio)
        self.MobilityCompensationGlobalRadio.setChecked(True)

        # Use Recalibrated Data
        self.RecalibratedDataCheckbox = QCheckBox(self.main)
        self.RecalibratedDataCheckbox.setObjectName(u"RecalibratedDataCheckbox")
        self.RecalibratedDataCheckbox.setGeometry(QRect(430, 500, 351, 21))
        self.RecalibratedDataCheckbox.setChecked(True)

        # Compression
        self.CompressionCheckbox = QCheckBox(self.main)
        self.CompressionCheckbox.setObjectName(u"CompressionCheckbox")
        self.CompressionCheckbox.setGeometry(QRect(430, 529, 351, 21))
        self.CompressionCheckbox.setChecked(True)

        # Exclude Mobility
        self.ExcludeMobilityCheckbox = QCheckBox(self.main)
        self.ExcludeMobilityCheckbox.setObjectName(u"ExcludeMobilityCheckbox")
        self.ExcludeMobilityCheckbox.setGeometry(QRect(430, 560, 351, 21))

        # MS2 Only
        self.Ms2OnlyCheckbox = QCheckBox(self.main)
        self.Ms2OnlyCheckbox.setObjectName(u"Ms2OnlyCheckbox")
        self.Ms2OnlyCheckbox.setGeometry(QRect(430, 590, 351, 21))

        # Barebones Metadata
        self.BarebonesMetadataCheckbox = QCheckBox(self.main)
        self.BarebonesMetadataCheckbox.setObjectName(u"BarebonesMetadataCheckbox")
        self.BarebonesMetadataCheckbox.setGeometry(QRect(430, 620, 351, 21))

        # Num Bins
        # Label
        self.NumBinsLabel = QLabel(self.main)
        self.NumBinsLabel.setObjectName(u"NumBinsLabel")
        self.NumBinsLabel.setGeometry(QRect(550, 150, 91, 21))
        self.NumBinsLabel.setVisible(False)
        # Checkbox
        self.BinSpectraCheckbox = QCheckBox(self.main)
        self.BinSpectraCheckbox.setObjectName(u"BinSpectraCheckbox")
        self.BinSpectraCheckbox.setGeometry(QRect(550, 120, 161, 21))
        self.BinSpectraCheckbox.setVisible(False)
        # Spin Box
        self.NumBinsSpinBox = QSpinBox(self.main)
        self.NumBinsSpinBox.setObjectName(u"NumBinsSpinBox")
        self.NumBinsSpinBox.setGeometry(QRect(550, 180, 91, 22))
        self.NumBinsSpinBox.setMaximum(1000000000)
        self.NumBinsSpinBox.setSingleStep(100)
        self.NumBinsSpinBox.setValue(0)
        self.NumBinsSpinBox.setVisible(False)

        # MALDI Dried Droplet Output Mode
        # Label
        self.MaldiOutputModeLabel = QLabel(self.main)
        self.MaldiOutputModeLabel.setObjectName(u"MaldiOutputModeLabel")
        self.MaldiOutputModeLabel.setGeometry(QRect(430, 660, 351, 21))
        # Combined
        self.MaldiOutputModeCombinedRadio = QRadioButton(self.main)
        self.MaldiOutputModeCombinedRadio.setObjectName(u"MaldiOutputModeCombinedRadio")
        self.MaldiOutputModeCombinedRadio.setGeometry(QRect(430, 690, 351, 20))
        # Individual
        self.MaldiOutputModeIndividualRadio = QRadioButton(self.main)
        self.MaldiOutputModeIndividualRadio.setObjectName(u"MaldiOutputModeIndividualRadio")
        self.MaldiOutputModeIndividualRadio.setGeometry(QRect(430, 720, 351, 20))
        # Group
        self.MaldiOutputModeGroupRadio = QRadioButton(self.main)
        self.MaldiOutputModeGroupRadio.setObjectName(u"MaldiOutputModeGroupRadio")
        self.MaldiOutputModeGroupRadio.setGeometry(QRect(430, 750, 351, 20))
        # Group
        self.MaldiOutputModeGroup = QButtonGroup()
        self.MaldiOutputModeGroup.addButton(self.MaldiOutputModeCombinedRadio)
        self.MaldiOutputModeGroup.addButton(self.MaldiOutputModeIndividualRadio)
        self.MaldiOutputModeGroup.addButton(self.MaldiOutputModeGroupRadio)
        self.MaldiOutputModeCombinedRadio.setChecked(True)

        # MALDI Plate Map
        # Line Edit
        self.MaldiPlateMapLine = QLineEdit(self.main)
        self.MaldiPlateMapLine.setObjectName(u"MaldiPlateMapLine")
        self.MaldiPlateMapLine.setGeometry(QRect(430, 820, 271, 22))
        self.MaldiPlateMapLine.setMaxLength(32765)
        self.MaldiPlateMapLine.setReadOnly(True)
        # Button
        self.MaldiPlateMapBrowseButton = QPushButton(self.main)
        self.MaldiPlateMapBrowseButton.setObjectName(u"MaldiPlateMapBrowseButton")
        self.MaldiPlateMapBrowseButton.setGeometry(QRect(710, 820, 61, 24))
        # Label
        self.MaldiPlateMapLabel = QLabel(self.main)
        self.MaldiPlateMapLabel.setObjectName(u"MaldiPlateMapLabel")
        self.MaldiPlateMapLabel.setGeometry(QRect(430, 790, 351, 20))

        # imzML Mode
        # Label
        self.MaldiImzmlModeLabel = QLabel(self.main)
        self.MaldiImzmlModeLabel.setObjectName(u"MaldiImzmlModeLabel")
        self.MaldiImzmlModeLabel.setGeometry(QRect(430, 860, 351, 21))
        # Continuous
        self.MaldiImzmlModeContinuousRadio = QRadioButton(self.main)
        self.MaldiImzmlModeContinuousRadio.setObjectName(u"MaldiImzmlModeContinuousRadio")
        self.MaldiImzmlModeContinuousRadio.setGeometry(QRect(550, 891, 89, 21))
        # Processed
        self.MaldiImzmlModeProcessedRadio = QRadioButton(self.main)
        self.MaldiImzmlModeProcessedRadio.setObjectName(u"MaldiImzmlModeProcessedRadio")
        self.MaldiImzmlModeProcessedRadio.setGeometry(QRect(430, 890, 89, 21))
        # Group
        self.MaldiImzmlModeGroup = QButtonGroup()
        self.MaldiImzmlModeGroup.addButton(self.MaldiImzmlModeProcessedRadio)
        self.MaldiImzmlModeGroup.addButton(self.MaldiImzmlModeContinuousRadio)
        self.MaldiImzmlModeProcessedRadio.setChecked(True)

        # Run
        self.RunButton = QPushButton(self.main)
        self.RunButton.setObjectName(u"RunButton")
        self.RunButton.setGeometry(QRect(20, 890, 381, 24))

        # retranslate
        TimsconvertGuiWindow.setCentralWidget(self.main)
        TimsconvertGuiWindow.setMenuBar(self.MenuBar)
        self.MenuBar.addAction(self.About.menuAction())
        self.About.addAction(self.AboutTimsconvert)
        self.retranslateUi(TimsconvertGuiWindow)

        QMetaObject.connectSlotsByName(TimsconvertGuiWindow)

    def retranslateUi(self, TimsconvertGuiWindow):
        TimsconvertGuiWindow.setWindowTitle(QCoreApplication.translate("TimsconvertGuiWindow", u"TIMSCONVERT", None))
        self.AboutTimsconvert.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"TIMSCONVERT Information", None))
        self.About.setTitle(QCoreApplication.translate("TimsconvertGuiWindow", u"About", None))
        self.AddInputDialogueButton.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Add Bruker timsTOF Data to Queue", None))
        self.RemoveFromQueueButton.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Remove from Queue", None))
        self.OutputDirectoryLabel.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Output Directory", None))
        self.OutputDirectoryBrowseButton.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Browse", None))
        self.ModeLabel.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Spectra Export Mode", None))
        self.ModeProfileRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Profile", None))
        self.ModeCentroidRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Centroid", None))
        self.ModeRawRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Raw", None))
        self.EncodingLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Binary Encoding Precision", None))
        self.MZEncodingLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"m/z Encoding Precision", None))
        self.MZEncoding64Radio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"64 bit", None))
        self.MZEncoding32Radio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"32 bit", None))
        self.IntensityEncodingLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Intensity Encoding Precision", None))
        self.IntensityEncoding32Radio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"32 bit", None))
        self.IntensityEncoding64Radio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"64 bit", None))
        self.MobilityEncodingLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Mobility Encoding Precision", None))
        self.MobilityEncoding64Radio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"64 bit", None))
        self.MobilityEncoding32Radio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"32 bit", None))
        self.MobilityCompensationLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Mobility Compensation Strategy", None))
        self.MobilityCompensationNoneRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"None", None))
        self.MobilityCompensationGlobalRadio.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Global", None))
        self.MobilityCompensationFrameRadio.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Per Frame", None))
        self.RecalibratedDataCheckbox.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Use Recalibrated Data", None))
        self.CompressionCheckbox.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Use zlib Compression", None))
        self.ExcludeMobilityCheckbox.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Exclude Trapped Ion Mobility (TIMS) Arrays", None))
        self.Ms2OnlyCheckbox.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Only Export MS/MS Spectra", None))
        self.NumBinsLabel.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Number of Bins", None))
        self.BinSpectraCheckbox.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Bin Profile Mode Spectra", None))
        self.MaldiOutputModeLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"MALDI Dried Droplet Output Mode", None))
        self.MaldiOutputModeCombinedRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow",
                                                                             u"Combine All Spectra In a *.d File to a Single mzML File",
                                                                             None))
        self.MaldiOutputModeIndividualRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow",
                                                                               u"Export Each Spectrum In a *.d File to Individual mzML Files",
                                                                               None))
        self.MaldiOutputModeGroupRadio.setText(QCoreApplication.translate("TimsconvertGuiWindow",
                                                                           u"Group Spectra with Matching Group Names into mzML Files",
                                                                           None))
        self.MaldiPlateMapBrowseButton.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Browse", None))
        self.MaldiPlateMapLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"MALDI Dried Droplet Plate Map", None))
        self.MaldiImzmlModeLabel.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"MALDI Imaging imzML Mode", None))
        self.MaldiImzmlModeContinuousRadio.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Continuous", None))
        self.MaldiImzmlModeProcessedRadio.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Processed", None))
        self.BarebonesMetadataCheckbox.setText(
            QCoreApplication.translate("TimsconvertGuiWindow", u"Only Include Barebones Metadata", None))
        self.RunButton.setText(QCoreApplication.translate("TimsconvertGuiWindow", u"Convert Data", None))

    def populate_table(self, data):
        old_row_count = self.InputList.rowCount()
        self.InputList.setRowCount(self.InputList.rowCount() + len(data))
        for row, input_filename in enumerate(data, start=old_row_count):
            text_item = QTableWidgetItem(input_filename)
            self.InputList.setItem(row, 0, text_item)
