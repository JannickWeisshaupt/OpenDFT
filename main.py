from __future__ import division
import sys
import os
import numpy as np
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from visualization import StructureVisualization,BandStructureVisualization,ScfVisualization
import solid_state_tools as sst
from exciting_handler import Handler as Handler
import pickle
import time

try:
    import queue
except:
    import Queue as queue

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
    def __init__(self, parent):
        self.parent = parent
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QGridLayout(self)
        self.start_ground_state_calculation_button = QtGui.QPushButton('Start Ground\nState Calculation', self)
        self.start_ground_state_calculation_button.setFixedWidth(150)
        self.start_ground_state_calculation_button.setFixedHeight(50)

        self.start_ground_state_calculation_button.clicked.connect(self.start_ground_state_calculation)
        self.execute_error_dialog = QtGui.QErrorMessage(parent=self)
        self.execute_error_dialog.resize(500,200)

        trash_bs_points = np.array([[0,0,0],[0.750,0.500,0.250],[0.500,0.500,0.500]
                    ,[0.000,0.000,0.000],[0.500,0.500,0.000],[0.750,0.500,0.250],[0.750,0.375,0.375],[0.000,0.000,0.000]])
        trash_bs_labels = ['GAMMA','W','L','GAMMA','X','W','K','GAMMA']
        self.band_structure_points = zip(trash_bs_points,trash_bs_labels)

    def start_ground_state_calculation(self):
        try:
            tree = esc_handler.make_tree()
            esc_handler.add_scf_to_tree(tree,self.parent.crystal_structure)
            esc_handler.add_bs_to_tree(tree,self.band_structure_points)
            esc_handler.write_input_file(tree)
            time.sleep(0.1)
            esc_handler.start_ground_state_calculation()
            QtCore.QTimer.singleShot(1000, self.parent.check_engine)
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>'+repr(e)+'<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
             '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)

class ScfWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.scf_widget = ScfVisualization(parent=self)
        layout.addWidget(self.scf_widget)


class BandStructureWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.bs_widget = BandStructureVisualization(parent=self)
        layout.addWidget(self.bs_widget)

class MainWindow(QtGui.QMainWindow):
    def __init__(self,central_window, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.central_window = central_window

    def closeEvent(self, event):
        self.central_window.close_application()
        event.ignore()

class CentralWindow(QtGui.QWidget):
    def __init__(self, *args, **kwargs):
        super(CentralWindow, self).__init__(*args, **kwargs)
        self.project_loaded = False

        self.crystal_structure = None
        self.band_structure = None
        self.saved_results = {}
        self.project_properties = {'title':''}

        self.error_dialog = QtGui.QErrorMessage(parent=self)
        self.error_dialog.resize(700,400)

        self.layout = QtGui.QGridLayout(self)
        self.mayavi_widget = MayaviQWidget(self.crystal_structure, parent=self)

        self.band_structure_window = BandStructureWindow(parent=self)
        self.dft_engine_window = DftEngineWindow(self)
        self.scf_window = ScfWindow(parent=self)


        self.tabWidget = QtGui.QTabWidget()
        self.layout.addWidget(self.tabWidget)

        self.tab_layout = QtGui.QVBoxLayout()
        self.tabWidget.setLayout(self.tab_layout)

        self.tabWidget.addTab(self.mayavi_widget,'Structure')
        self.tabWidget.addTab(self.dft_engine_window,'DFT-Engine')
        self.tabWidget.addTab(self.band_structure_window ,'Bandstructure')
        self.tabWidget.addTab(QtGui.QWidget(), 'Optical properties')
        self.tabWidget.addTab(self.scf_window, 'Scf')


        self.show()
        self.window = MainWindow(self)
        self.window.setWindowTitle("OpenDFT")
        self.window.setGeometry(50, 50, 1000, 1000)
        self.window.setCentralWidget(self)
        self.make_menu_bar()

        self.window.show()

        if DEBUG:
            if sys.platform in ['linux','linux2']:
                self.project_directory = r"/home/jannick/OpenDFT_projects/test/"
            else:
                self.project_directory = r'D:\OpenDFT_projects\test\\'
            os.chdir(self.project_directory)
            self.load_saved_results()

    def make_new_project(self):
        folder_name = QtGui.QFileDialog().getExistingDirectory(parent = self)
        if len(folder_name)>1:
            if self.project_loaded:
                self.save_results()
                self.reset_results_and_plots()
            self.project_directory = folder_name
            self.initialize_project()

    def initialize_project(self):
        self.project_properties = {'title': ''}
        os.chdir(self.project_directory)

    def reset_results_and_plots(self):
        self.crystal_structure = None
        self.band_structure = None
        self.mayavi_widget.visualization.clear_plot()
        self.band_structure_window.bs_widget.clear_plot()
        self.scf_window.scf_widget.clear_plot()


    def load_project(self):
        folder_name = QtGui.QFileDialog().getExistingDirectory(parent = self)
        if len(folder_name)>1:
            os.chdir(self.project_directory)
            self.project_directory = folder_name
            self.load_saved_results()
            self.project_loaded = True

    def save_results(self):
        a = {'crystal structure': self.crystal_structure,'band structure':self.band_structure,'properties':self.project_properties}
        with open(self.project_directory+'save.pkl', 'wb') as handle:
            pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def load_saved_results(self):
        esc_handler.project_directory = self.project_directory
        try:
            with open('save.pkl', 'rb') as handle:
                b = pickle.load(handle)
                self.crystal_structure = b['crystal structure']
                if self.crystal_structure is not None:
                    self.mayavi_widget.update_crystal_structure(self.crystal_structure)
                    self.mayavi_widget.update_plot()

                self.band_structure = b['band structure']
                if self.band_structure is not None:
                    self.band_structure_window.bs_widget.plot(self.band_structure)

                self.project_properties = b['properties']

        except Exception as e:
            print(e)


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
        self.menu_bar.setNativeMenuBar(False)
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
        if not DEBUG:
            reply = QtGui.QMessageBox.question(self, 'Message',
                "Are you sure to quit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if DEBUG or reply == QtGui.QMessageBox.Yes:
            self.save_results()
            sys.exit()

    def check_engine(self):
        if esc_handler.is_engine_running():
            self.scf_data = esc_handler.read_scf_status()
            if self.scf_data is not None:
                self.scf_window.scf_widget.plot(self.scf_data)
            QtCore.QTimer.singleShot(500, self.check_engine)
        else:
            message,err = esc_handler.engine_process.communicate()
            if 'error' in message.lower():
                error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + message \
                                + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                         '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
                self.error_dialog.showMessage(error_message)

            self.band_structure = esc_handler.read_bandstructure()
            self.band_structure_window.bs_widget.plot(self.band_structure)

if __name__ == "__main__":
    DEBUG = True
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    main = CentralWindow()
    app.exec_()
