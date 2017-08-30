import sys
import os
import numpy as np
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from visualization import StructureVisualization,BandStructureVisualization
import solid_state_tools as sst
from exciting_handler import Handler as Handler
import queue

esc_handler = Handler()
event_queue = queue.Queue()


class MayaviQWidget(QtGui.QWidget):
    def __init__(self, crystal_structure, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.visualization = StructureVisualization(crystal_structure)
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)

    def update_plot(self):
        self.visualization.update_plot()

    def update_crystal_structure(self, crystal_structure):
        self.visualization.crystal_structure = crystal_structure

class DftEngineWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        self.parent = parent
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QGridLayout(self)
        self.start_ground_state_calculation_button = QtGui.QPushButton('Start Ground\nState Calculation', self)
        self.start_ground_state_calculation_button.clicked.connect(self.start_ground_state_calculation)

    def start_ground_state_calculation(self):
        esc_handler.start_engine()

class BandStructureWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.bs_widget = BandStructureVisualization(parent=self)
        layout.addWidget(self.bs_widget)

class MainWindow(QtGui.QWidget):
    def __init__(self, *args, **kwargs):
        super(MainWindow,self).__init__(*args, **kwargs)
        self.crystal_structure = None
        self.saved_results = {}

        self.layout = QtGui.QGridLayout(self)
        self.mayavi_widget = MayaviQWidget(self.crystal_structure, parent=self)

        self.band_structure_window = BandStructureWindow(parent=self)
        self.dft_engine_window = DftEngineWindow(parent=self)

        self.tabWidget = QtGui.QTabWidget()
        self.layout.addWidget(self.tabWidget)

        self.tab_layout = QtGui.QVBoxLayout()
        self.tabWidget.setLayout(self.tab_layout)

        self.tabWidget.addTab(self.mayavi_widget,'Structure')
        self.tabWidget.addTab(self.dft_engine_window,'DFT-Engine')
        self.tabWidget.addTab(self.band_structure_window ,'Bandstructure')
        self.tabWidget.addTab(QtGui.QWidget(), 'Optical properties')


        # self.layout.addWidget(self.mayavi_widget, 0, 0)
        self.show()
        self.window = QtGui.QMainWindow()
        self.window.setWindowTitle("OpenDFT")
        self.window.setGeometry(50, 50, 1000, 1000)
        self.window.setCentralWidget(self)
        self.window.show()

        self.make_menu_bar()

    def make_new_project(self):
        folder_name = QtGui.QFileDialog().getExistingDirectory(parent = self)
        if len(folder_name)>1:
            self.project_directory = folder_name
            self.initiliaze_project()
            self.load_saved_results()

    def initialize_project(self):
        pass

    def load_project(self):
        folder_name = QtGui.QFileDialog().getExistingDirectory(parent = self)
        if len(folder_name)>1:
            self.project_directory = folder_name
            self.load_saved_results()

    def load_saved_results(self):
        os.chdir(self.project_directory)
        esc_handler.project_directory = self.project_directory

    def update_structure_plot(self):
        self.mayavi_widget.update_crystal_structure(self.crystal_structure)
        self.mayavi_widget.update_plot()

    def load_crystal_structure(self):
        file_dialog = QtGui.QFileDialog()
        file_dialog.setNameFilters(["Exciting (*.xml)", "All (*.*)"])

        if file_dialog.exec_():
            file_name = file_dialog.selectedFiles()
            if type(file_name) == list or type(file_name) is tuple:
                file_name = file_name[0]
            if len(file_name) == 0:
                return
            self.crystal_structure = esc_handler.parse_input_file(file_name)
            self.update_structure_plot()

    def make_menu_bar(self):
        self.menu_bar = self.window.menuBar()
        self.file_menu = self.menu_bar.addMenu('&File')

        new_project_action = QtGui.QAction("New project", self.window)
        new_project_action.setShortcut("Ctrl+n")
        new_project_action.setStatusTip('Start new project')
        new_project_action.triggered.connect(self.make_new_project)
        self.file_menu.addAction(new_project_action)

        load_project_action = QtGui.QAction("Load project", self.window)
        # load_project_action.setShortcut("Ctrl+n")
        load_project_action.setStatusTip('Load project')
        load_project_action.triggered.connect(self.load_project)
        self.file_menu.addAction(load_project_action)

        open_structure_action = QtGui.QAction("Load structure", self.window)
        open_structure_action.setShortcut("Ctrl+O")
        open_structure_action.setStatusTip('Load crystal structure')
        open_structure_action.triggered.connect(self.load_crystal_structure)
        self.file_menu.addAction(open_structure_action)

        close_app_action = QtGui.QAction("Exit", self.window)
        close_app_action.setShortcut("Ctrl+Q")
        close_app_action.setStatusTip('Leave The App')
        close_app_action.triggered.connect(self.close_application)
        self.file_menu.addAction(close_app_action)

    def close_application(self):
        sys.exit()


if __name__ == "__main__":
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    main = MainWindow()
    app.exec_()
