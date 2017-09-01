from __future__ import division
import sys
import os
import numpy as np

os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from visualization import StructureVisualization, BandStructureVisualization, ScfVisualization
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

    def do_select_event(self):
        pass

class EntryWithLabel(QtGui.QWidget):
    def __init__(self, parent,label,value=None):
        QtGui.QWidget.__init__(self, parent)
        self.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.layout = QtGui.QHBoxLayout(self)
        self.textbox = QtGui.QLineEdit(self)
        self.textbox.setMaximumWidth(200)
        self.label_widget = QtGui.QLabel(label,parent=self)
        self.label_widget.setMaximumWidth(70)
        self.layout.setAlignment(QtCore.Qt.AlignLeft)
        self.layout.addWidget(self.label_widget)
        self.layout.addWidget(self.textbox)

        if value is not None:
            self.textbox.setText(value)

    def get_text(self):
        return self.textbox.text()

    def set_text(self,text):
        self.textbox.setText(text)

class OptionFrame(QtGui.QGroupBox):
    def __init__(self, parent,options,title='',tooltips={},checkbuttons=[]):
        QtGui.QGroupBox.__init__(self, parent)
        self.widgets_per_line = 4
        self.setTitle(title)
        self.tooltips = tooltips
        self.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.parent = parent
        self.options = options
        self.layout = QtGui.QGridLayout(self)
        # self.layout.setAlignment(QtCore.Qt.AlignTop)
        self.entry_dict = {}
        self.checkbuttons = []
        for text,state in checkbuttons:
            cb = QtGui.QCheckBox(text,parent=self)
            if state:
                cb.nextCheckState()
            self.checkbuttons.append(cb)
            self.layout.addWidget(cb)
        self.make_option_entries()

    def make_option_entries(self):
        counter = len(self.checkbuttons)//self.widgets_per_line+self.widgets_per_line
        for option_key,option_value in self.options.items():
            entry = EntryWithLabel(self,option_key,option_value)
            if option_key in self.tooltips.keys():
                entry.setToolTip(self.tooltips[option_key].replace('\n','<br>'))

            self.layout.addWidget(entry,counter//self.widgets_per_line,counter%self.widgets_per_line)
            self.entry_dict[option_key] = entry
            counter += 1

    def read_all_entries(self):
        for key,entry in self.entry_dict.items():
            self.options[key] = entry.get_text()

    def set_all_entries(self):
        for key,value in self.options.items():
            self.entry_dict[key].set_text(value)

    def read_checkbuttons(self):
        res = {}
        for cb in self.checkbuttons:
            res[cb.text()] = cb.checkState()
        return res

class DftEngineWindow(QtGui.QWidget):
    def __init__(self, parent):
        self.parent = parent
        QtGui.QWidget.__init__(self, parent)
        # self.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.layout = QtGui.QGridLayout(self)
        self.layout.setAlignment(QtCore.Qt.AlignTop)

        self.general_option_widget = OptionFrame(self,esc_handler.general_options,title='General options')
        self.layout.addWidget(self.general_option_widget)

        self.scf_option_widget = OptionFrame(self,esc_handler.scf_options,title='Groundstate options',tooltips=esc_handler.scf_options_tooltip)
        self.layout.addWidget(self.scf_option_widget)

        self.bs_option_widget = OptionFrame(self,esc_handler.bs_options,title='Bandstructure options',checkbuttons=[['Calculate',True]])
        self.layout.addWidget(self.bs_option_widget)

        self.relax_option_widget = OptionFrame(self,esc_handler.relax_options,title='Structure relaxation options')
        self.layout.addWidget(self.relax_option_widget)

        self.button_widget = QtGui.QWidget(self)
        self.button_layout = QtGui.QHBoxLayout(self.button_widget)
        self.button_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.start_ground_state_calculation_button = QtGui.QPushButton('Start Ground\nState Calculation', self.button_widget)
        self.start_ground_state_calculation_button.setFixedWidth(150)
        self.start_ground_state_calculation_button.setFixedHeight(50)
        self.start_ground_state_calculation_button.clicked.connect(self.start_ground_state_calculation)
        self.button_layout.addWidget(self.start_ground_state_calculation_button)

        self.start_relax_button = QtGui.QPushButton('Start Structure\nRelaxation', self.button_widget)
        self.start_relax_button.setFixedWidth(150)
        self.start_relax_button.setFixedHeight(50)
        self.start_relax_button.clicked.connect(self.start_relax)
        self.button_layout.addWidget(self.start_relax_button)

        self.abort_calculation_button = QtGui.QPushButton('Abort Calculation', self.button_widget)
        self.abort_calculation_button.setFixedWidth(150)
        self.abort_calculation_button.setFixedHeight(50)
        self.abort_calculation_button.clicked.connect(self.parent.abort_calculation)
        self.button_layout.addWidget(self.abort_calculation_button)

        self.execute_error_dialog = QtGui.QErrorMessage(parent=self)
        self.execute_error_dialog.resize(500, 200)

        self.layout.addWidget(self.button_widget)

        trash_bs_points = np.array([[0, 0, 0], [0.750, 0.500, 0.250], [0.500, 0.500, 0.500]
                                       , [0.000, 0.000, 0.000], [0.500, 0.500, 0.000], [0.750, 0.500, 0.250],
                                    [0.750, 0.375, 0.375], [0.000, 0.000, 0.000]])
        trash_bs_labels = ['GAMMA', 'W', 'L', 'GAMMA', 'X', 'W', 'K', 'GAMMA']
        self.band_structure_points = zip(trash_bs_points, trash_bs_labels)
        self.show()

    def update_all(self):
        self.scf_option_widget.set_all_entries()
        self.general_option_widget.set_all_entries()
        self.bs_option_widget.set_all_entries()

    def do_select_event(self):
        pass

    def check_if_engine_is_running_and_warn_if_so(self):
        pass

    def read_all_option_widgets(self):
        self.scf_option_widget.read_all_entries()
        self.general_option_widget.read_all_entries()
        self.bs_option_widget.read_all_entries()

    def start_ground_state_calculation(self):
        self.check_if_engine_is_running_and_warn_if_so()
        tasks = []
        self.read_all_option_widgets()
        bs_checkers = self.bs_option_widget.read_checkbuttons()
        if bs_checkers['Calculate']:
            bs_points = self.band_structure_points
            tasks.append('bandstructure')
        else:
            bs_points = None
        try:
            esc_handler.start_ground_state_calculation(self.parent.crystal_structure, band_structure_points=bs_points)
            QtCore.QTimer.singleShot(1000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_relax(self):
        self.check_if_engine_is_running_and_warn_if_so()
        tasks = ['relax']
        self.read_all_option_widgets()
        try:
            esc_handler.start_relax(self.parent.crystal_structure)
            QtCore.QTimer.singleShot(1000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

class ScfWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.scf_widget = ScfVisualization(parent=self)
        layout.addWidget(self.scf_widget)

    def do_select_event(self):
        pass

class BandStructureWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self)
        self.parent = parent
        self.layout = QtGui.QHBoxLayout(self)
        self.bs_widget = BandStructureVisualization(parent=self)
        self.treeview = QtGui.QTreeWidget(parent=self)
        self.treeview.setMaximumWidth(200)
        self.treeview.setHeaderHidden(True)
        self.treeview.itemClicked.connect(self.handle_item_changed)

        self.treeview.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeview.customContextMenuRequested.connect(self.openMenu)

        self.layout.addWidget(self.bs_widget)
        self.layout.addWidget(self.treeview)

        self.show()

    def delete_selected_item(self):
        index = self.treeview.selectedIndexes()[0]
        item = self.treeview.itemFromIndex(index)
        bs_name = item.text(0)
        del self.parent.band_structure[bs_name]
        self.update_tree()

    def openMenu(self, position):

        indexes = self.treeview.selectedIndexes()
        if len(indexes) > 0:

            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1

        menu = QtGui.QMenu()
        menu.addAction('Delete',self.delete_selected_item)
        menu.exec_(self.treeview.viewport().mapToGlobal(position))

    def handle_item_changed(self,item,column):
        bs_name = item.text(0)
        self.bs_widget.plot(self.parent.band_structure[bs_name])

    def add_bandstructure_key(self, title):
        item = QtGui.QTreeWidgetItem(self.treeview.invisibleRootItem(), [title])
        return item

    def update_tree(self):
        self.treeview.clear()
        for key,value in self.parent.band_structure.items():
            self.add_bandstructure_key(key)

    def do_select_event(self):
        self.update_tree()


class StatusBar(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self)
        # self.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.running_text = 'Engine is running'
        self.not_running_text = 'Engine is stopped'
        # self.setMaximumHeight(20)
        self.layout = QtGui.QVBoxLayout(self)
        self.layout.setAlignment(QtCore.Qt.AlignRight)
        self.status_label = QtGui.QLabel(self.not_running_text)
        # self.status_label.setMaximumHeight(20)
        self.layout.addWidget(self.status_label)
        self.show()

    def set_engine_status(self,status):
        if status:
            self.status_label.setText(self.running_text)
        else:
            self.status_label.setText(self.not_running_text)

class MainWindow(QtGui.QMainWindow):
    def __init__(self, central_window, *args, **kwargs):
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
        self.band_structure = {}
        self.saved_results = {}
        self.project_properties = {'title': ''}

        self.error_dialog = QtGui.QErrorMessage(parent=self)
        self.error_dialog.resize(700, 600)

        self.layout = QtGui.QGridLayout(self)
        self.mayavi_widget = MayaviQWidget(self.crystal_structure, parent=self)

        self.band_structure_window = BandStructureWindow(parent=self)
        self.dft_engine_window = DftEngineWindow(self)
        self.scf_window = ScfWindow(parent=self)

        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.currentChanged.connect(self.tab_is_changed)
        self.layout.addWidget(self.tabWidget)

        self.status_bar = StatusBar()
        self.layout.addWidget(self.status_bar)


        self.tab_layout = QtGui.QVBoxLayout()
        self.tabWidget.setLayout(self.tab_layout)

        self.list_of_tabs = [self.mayavi_widget,self.dft_engine_window,self.band_structure_window,self.scf_window]


        self.tabWidget.addTab(self.list_of_tabs[0], 'Structure')
        self.tabWidget.addTab(self.list_of_tabs[1], 'DFT-Engine')
        self.tabWidget.addTab(self.list_of_tabs[2],'Bandstructure')
        # self.tabWidget.addTab(QtGui.QWidget(), 'Optical properties')
        self.tabWidget.addTab(self.list_of_tabs[3], 'Scf')


        self.show()
        self.window = MainWindow(self)
        self.window.setWindowTitle("OpenDFT")
        self.window.setGeometry(50, 50, 1300, 900)
        self.window.setCentralWidget(self)
        self.make_menu_bar()

        self.window.show()

        if DEBUG:
            if sys.platform in ['linux', 'linux2']:
                project_directory = r"/home/jannick/OpenDFT_projects/diamond/"
            else:
                project_directory = r'D:\OpenDFT_projects\test\\'
            # self.load_saved_results()
            QtCore.QTimer.singleShot(500, lambda: self.load_project(folder_name=project_directory))


    def tab_is_changed(self,i):
        self.list_of_tabs[i].do_select_event()

    def make_new_project(self):
        folder_name = QtGui.QFileDialog().getExistingDirectory(parent=self)
        if len(folder_name) > 1:
            if self.project_loaded:
                self.save_results()
                self.reset_results_and_plots()
            self.project_directory = folder_name
            esc_handler.project_directory = self.project_directory
            self.initialize_project()

    def initialize_project(self):
        self.project_properties = {'title': ''}
        os.chdir(self.project_directory)

    def reset_results_and_plots(self):
        self.crystal_structure = None
        self.band_structure = {}
        self.mayavi_widget.visualization.clear_plot()
        self.band_structure_window.bs_widget.clear_plot()
        self.scf_window.scf_widget.clear_plot()

    def load_project(self,folder_name=None):
        if folder_name is None:
            folder_name = QtGui.QFileDialog().getExistingDirectory(parent=self)
        if len(folder_name) > 1:
            self.reset_results_and_plots()
            self.project_directory = folder_name
            os.chdir(self.project_directory)
            esc_handler.project_directory = self.project_directory
            self.load_saved_results()
            self.dft_engine_window.update_all()
            self.project_loaded = True

    def save_results(self):
        try:
            a = {'crystal structure': self.crystal_structure, 'band structure': self.band_structure,
                 'properties': self.project_properties,'scf_options':esc_handler.scf_options,
                 'dft engine':esc_handler.engine_name,'general options':esc_handler.general_options,'bs options':esc_handler.bs_options}
            with open(self.project_directory + '/save.pkl', 'wb') as handle:
                pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as e:
            print(e)

    def load_saved_results(self):
        esc_handler.project_directory = self.project_directory
        try:
            with open(self.project_directory + '/save.pkl', 'rb') as handle:
                b = pickle.load(handle)
                self.crystal_structure = b['crystal structure']
                if self.crystal_structure is not None:
                    self.mayavi_widget.update_crystal_structure(self.crystal_structure)
                    self.mayavi_widget.update_plot()
                loaded_bandstructure_dict = b['band structure']
                if type(loaded_bandstructure_dict) == dict:
                    self.band_structure = loaded_bandstructure_dict
                # if len(self.band_structure) !
                #     self.band_structure_window.bs_widget.plot(self.band_structure)
                load_scf_options = b['scf_options']
                if load_scf_options is not None and b['dft engine'] == esc_handler.engine_name:
                    for key,value in load_scf_options.items():
                        esc_handler.scf_options[key] = value

                load_general_options = b['general options']
                if load_general_options is not None and b['dft engine'] == esc_handler.engine_name:
                    for key,value in load_general_options.items():
                        esc_handler.general_options[key] = value

                load_bs_options = b['bs options']
                if load_bs_options is not None and b['dft engine'] == esc_handler.engine_name:
                    for key,value in load_bs_options.items():
                        esc_handler.bs_options[key] = value


                self.project_properties = b['properties']

        except Exception as e:
            print('Load failed with'+ repr(e))

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
        load_project_action.setShortcut("Ctrl+o")
        load_project_action.setStatusTip('Load project')
        load_project_action.triggered.connect(self.load_project)
        self.file_menu.addAction(load_project_action)

        save_project_action = QtGui.QAction("Save project", self.window)
        save_project_action.setShortcut("Ctrl+s")
        save_project_action.setStatusTip('Save project')
        save_project_action.triggered.connect(self.save_results)
        self.file_menu.addAction(save_project_action)

        open_structure_action = QtGui.QAction("Load structure", self.window)
        # open_structure_action.setShortcut("Ctrl+O")
        open_structure_action.setStatusTip('Load crystal structure')
        open_structure_action.triggered.connect(self.load_crystal_structure)
        self.file_menu.addAction(open_structure_action)

        close_app_action = QtGui.QAction("Exit", self.window)
        close_app_action.setShortcut("Ctrl+Q")
        close_app_action.setStatusTip('Leave The App')
        close_app_action.triggered.connect(self.close_application)
        self.file_menu.addAction(close_app_action)

    def abort_calculation(self):
        esc_handler.kill_engine()

    def close_application(self):
        if not DEBUG:
            reply = QtGui.QMessageBox.question(self, 'Message',
                                               "Are you sure to quit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if DEBUG or reply == QtGui.QMessageBox.Yes:
            self.save_results()
            sys.exit()

    def check_engine(self,tasks):

        def check_relax():
            new_struc = esc_handler.load_relax_structure()
            if new_struc is not None:
                self.crystal_structure = new_struc
                self.update_structure_plot()

        tasks = [x.lower() for x in tasks]
        if esc_handler.is_engine_running():
            self.scf_data = esc_handler.read_scf_status()
            if self.scf_data is not None:
                self.scf_window.scf_widget.plot(self.scf_data)
            QtCore.QTimer.singleShot(500,lambda: self.check_engine(tasks))
            self.status_bar.set_engine_status(True)
            if 'relax' in tasks:
                check_relax()
        else:
            self.status_bar.set_engine_status(False)
            message, err = esc_handler.engine_process.communicate()
            if ('error' in message.lower() or len(err)>0):
                error_message = 'DFT calculation finished with an error:<br><br>' + message+'<br>Error:<br>'+err \
                                + '<br><br>Try following:<br>1.Check if the selected dft engine is correctly installed<br>' \
                                  '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
                self.error_dialog.showMessage(error_message)

            if 'bandstructure' in tasks:
                read_bandstructure = esc_handler.read_bandstructure()
                if read_bandstructure is not None:
                    self.band_structure[esc_handler.general_options['title']] = read_bandstructure
                    self.band_structure_window.bs_widget.plot(self.band_structure[esc_handler.general_options['title']])
                    self.band_structure_window.update_tree()
            elif 'relax' in tasks:
                check_relax()


if __name__ == "__main__":
    DEBUG = True
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    main = CentralWindow()
    app.exec_()
