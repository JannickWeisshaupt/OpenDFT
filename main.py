import sys
import os
import numpy as np
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from visualization import StructureVisualization
import solid_state_tools as sst
from exciting_handler import Handler as Handler

esc_handler = Handler()

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


class MainWindow(QtGui.QWidget):
    def __init__(self, *args, **kwargs):
        super(MainWindow,self).__init__(*args, **kwargs)
        self.layout = QtGui.QGridLayout(self)
        self.crystal_structure = None
        self.mayavi_widget = MayaviQWidget(self.crystal_structure, parent=self)

        self.layout.addWidget(self.mayavi_widget, 0, 0)
        self.show()
        self.window = QtGui.QMainWindow()
        self.window.setWindowTitle("Embedding Mayavi in a PyQt4 Application")
        self.window.setGeometry(50, 50, 1000, 1000)
        self.window.setCentralWidget(self)
        self.window.show()

        self.make_menu_bar()

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
