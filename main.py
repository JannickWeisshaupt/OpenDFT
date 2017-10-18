#!/usr/bin/python
from __future__ import division,absolute_import,print_function,unicode_literals
import sys
import os
import numpy as np

os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from visualization import StructureVisualization, BandStructureVisualization, ScfVisualization,OpticalSpectrumVisualization,colormap_list,BrillouinVisualization
import solid_state_tools as sst
from solid_state_tools import p_table,p_table_rev
from little_helpers import no_error_dictionary,CopySelectedCellsAction,PasteIntoTable
import pickle
import time
import threading
from collections import OrderedDict
import logging


try:
    import queue
except:
    import Queue as queue


general_handler = sst.GeneralHandler()
event_queue = queue.Queue()


class BrillouinWindow(QtGui.QDialog):

    def __init__(self, parent=None):
        super(BrillouinWindow, self).__init__(parent)

        self.k_path = None

        self.resize(900,600)
        layout = QtGui.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.mayavi_widget = BrillouinVisualization(self)
        self.ui = self.mayavi_widget.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)

        table_widget = QtGui.QWidget(parent=self)
        layout.addWidget(table_widget)

        table_layout = QtGui.QVBoxLayout(table_widget)

        self.table =  QtGui.QTableWidget(table_widget)
        self.table.setColumnCount(4)
        self.table.setRowCount(0)
        table_layout.addWidget(self.table)

        for i,label in enumerate(['x','y','z','Label']):
            item = QtGui.QTableWidgetItem()
            self.table.setHorizontalHeaderItem(i, item)
            item.setText(label)

        self.table.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        self.connect_tables()

        button_widget = QtGui.QWidget(parent=table_widget)
        table_layout.addWidget(button_widget)

        button_layout = QtGui.QHBoxLayout(button_widget)

        remove_all_button = QtGui.QPushButton('Remove path',parent=button_widget)
        button_layout.addWidget(remove_all_button)
        remove_all_button.clicked.connect(self.clear_path)

        remove_selected_button = QtGui.QPushButton('Remove point',parent=button_widget)
        button_layout.addWidget(remove_selected_button)
        remove_selected_button.clicked.connect(self.remove_point)

        add_atom_button = QtGui.QPushButton('Add point',parent=button_widget)
        button_layout.addWidget(add_atom_button)
        add_atom_button.clicked.connect(self.add_atom)

    def clear_path(self):
        for i in range(len(self.k_path)):
            del self.k_path[0]

        self.mayavi_widget.set_path(self.k_path)
        self.mayavi_widget.plot_path()
        self.update_table()

    def remove_point(self,*args,**kwargs):
        points = kwargs.pop('points',None)
        self.disconnect_tables()
        if points is None:
           points = sorted(set(index.row() for index in self.table.selectedIndexes()))
        for point in points[::-1]:
            self.table.removeRow(point)
        self.connect_tables()
        self.handle_change()

    def update_table(self):
        self.disconnect_tables()
        self.table.setRowCount(len(self.k_path))
        for i in range(len(self.k_path)):
            for j in range(4):
                item = QtGui.QTableWidgetItem()
                if j<3:
                    text = "{0:1.6f}".format(self.k_path[i][0][j])
                else:
                    text = self.k_path[i][1]
                item.setText(text)
                self.table.setItem(i,j,item)
        self.connect_tables()

    def set_path(self,k_path):

        self.k_path = k_path
        self.mayavi_widget.set_path(k_path)
        self.update_table()

    def add_atom(self):
        self.disconnect_tables()
        n = self.table.rowCount()
        self.table.setRowCount(n+1)
        for i in range(4):
            item = QtGui.QTableWidgetItem()
            self.table.setItem(n+1,i,item)
        self.connect_tables()

    def read_table(self):
        n_k = self.table.rowCount()
        k_points = []
        for i in range(n_k):
            try:
                coords = np.array([float(self.table.item(i,j).text()) for j in range(3)])
                k_point = [coords,str(self.table.item(i,3).text())]
                k_points.append(k_point)
            except Exception:
                pass

        return k_points

    def disconnect_tables(self):
        try:
            self.table.itemChanged.disconnect()
        except Exception as e:
            logging.exception(e)

    def connect_tables(self):
        self.table.itemChanged.connect(self.handle_change)

    def handle_change(self):
        k_path_read = self.read_table()
        for i in range(len(self.k_path)):
            del self.k_path[0]
        for el in k_path_read:
            self.k_path.append(el)

        self.mayavi_widget.set_path(self.k_path)
        self.mayavi_widget.plot_path()


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

    def update_plot(self,keep_view = False):
        self.visualization.update_plot(keep_view=keep_view)

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
        self.label_widget.setMaximumWidth(90)
        self.layout.setAlignment(QtCore.Qt.AlignLeft)
        self.layout.addWidget(self.label_widget)
        self.layout.addWidget(self.textbox)
        self.editFinished_command = None

        if value is not None:
            self.textbox.setText(value)

    def get_text(self):
        return self.textbox.text()

    def set_text(self,text):
        self.textbox.setText(text)

    def connect_editFinished(self,command):
        self.editFinished_command = command
        self.textbox.editingFinished.connect(self.handleEditingFinished)

    def handleEditingFinished(self):
        if self.textbox.isModified():
            self.editFinished_command()
        self.textbox.setModified(False)


class OptionFrame(QtGui.QGroupBox):
    def __init__(self, parent,options,title='',tooltips={},checkbuttons=[],buttons=[]):
        QtGui.QGroupBox.__init__(self, parent)
        self.widgets_per_line = 4
        self.setTitle(title)
        self.tooltips = tooltips
        self.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.parent = parent
        self.options = options
        self.layout = QtGui.QGridLayout(self)
        self.layout.setAlignment(QtCore.Qt.AlignLeft)
        self.entry_dict = {}
        self.checkbuttons = []
        self.buttons = []
        counter = 0

        for text,state in checkbuttons:
            cb = QtGui.QCheckBox(text,parent=self)
            if state:
                cb.nextCheckState()
            self.checkbuttons.append(cb)
            self.layout.addWidget(cb, counter // self.widgets_per_line, counter % self.widgets_per_line)
            counter += 1
        for text,function in buttons:
            button = QtGui.QPushButton(text)
            button.clicked.connect(function)
            button.setFixedHeight(30)
            self.buttons.append(button)
            self.layout.addWidget(button, counter // self.widgets_per_line, counter % self.widgets_per_line)
            counter += 1
        self.make_option_entries()

    def make_option_entries(self):
        counter = (len(self.checkbuttons)+len(self.buttons))//self.widgets_per_line+self.widgets_per_line
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
        self.abort_bool = False

        QtGui.QWidget.__init__(self, parent)
        # self.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.layout = QtGui.QGridLayout(self)
        self.layout.setAlignment(QtCore.Qt.AlignTop)

        mygroupbox = QtGui.QWidget()
        myform = QtGui.QFormLayout()

        self.scroll_area = QtGui.QScrollArea(parent=self)
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setFixedHeight(800)
        self.scroll_area.setWidget(mygroupbox)

        self.layout.addWidget(self.scroll_area)

        self.general_option_widget = OptionFrame(self,esc_handler.general_options,title='General options')
        myform.addRow(self.general_option_widget)

        self.scf_option_widget = OptionFrame(self,esc_handler.scf_options,title='Groundstate options',tooltips=esc_handler.scf_options_tooltip)
        myform.addRow(self.scf_option_widget)

        self.bs_option_widget = OptionFrame(self,esc_handler.bs_options,title='Bandstructure options',checkbuttons=[['Calculate',True]],buttons=[['Choose k-path',self.parent.open_brillouin_window]])
        myform.addRow(self.bs_option_widget)

        self.relax_option_widget = OptionFrame(self,esc_handler.relax_options,title='Structure relaxation options')
        myform.addRow(self.relax_option_widget)

        self.gw_option_widget = OptionFrame(self,esc_handler.gw_options,title='GW options',tooltips=esc_handler.gw_options_tooltip)
        myform.addRow(self.gw_option_widget)

        self.phonons_option_widget = OptionFrame(self,esc_handler.phonons_options,title='Phonon options')
        myform.addRow(self.phonons_option_widget)

        self.optical_spectrum_option_widget = OptionFrame(self,esc_handler.optical_spectrum_options,title='Excited states options')
        myform.addRow(self.optical_spectrum_option_widget)

        mygroupbox.setLayout(myform)

        self.button_widget = QtGui.QWidget(self)
        self.button_widget.show()
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

        self.start_gw_button = QtGui.QPushButton('Start GW', self.button_widget)
        self.start_gw_button.setFixedWidth(150)
        self.start_gw_button.setFixedHeight(50)
        self.start_gw_button.clicked.connect(self.start_gw)
        self.button_layout.addWidget(self.start_gw_button)

        self.start_phonon_button = QtGui.QPushButton('Start Phonon\nBandstructure', self.button_widget)
        self.start_phonon_button.setFixedWidth(150)
        self.start_phonon_button.setFixedHeight(50)
        self.start_phonon_button.clicked.connect(self.start_phonons)
        self.button_layout.addWidget(self.start_phonon_button)

        self.start_optical_spectrum_button = QtGui.QPushButton('Calculate optical\nspectrum', self.button_widget)
        self.start_optical_spectrum_button.setFixedWidth(150)
        self.start_optical_spectrum_button.setFixedHeight(50)
        self.start_optical_spectrum_button.clicked.connect(self.start_optical_spectrum_calculation)
        self.button_layout.addWidget(self.start_optical_spectrum_button)

        self.abort_calculation_button = QtGui.QPushButton('Abort Calculation', self.button_widget)
        self.abort_calculation_button.setFixedWidth(150)
        self.abort_calculation_button.setFixedHeight(50)
        self.abort_calculation_button.clicked.connect(self.abort_calculation)
        self.button_layout.addWidget(self.abort_calculation_button)

        self.execute_error_dialog = QtGui.QErrorMessage(parent=self)
        self.execute_error_dialog.resize(500, 200)

        self.layout.addWidget(self.button_widget)

        trash_bs_points = np.array([[0, 0, 0], [0.750, 0.500, 0.250], [0.500, 0.500, 0.500]
                                       , [0.000, 0.000, 0.000], [0.500, 0.500, 0.000], [0.750, 0.500, 0.250],
                                    [0.750, 0.375, 0.375], [0.000, 0.000, 0.000]])
        trash_bs_labels = ['GAMMA', 'W', 'L', 'GAMMA', 'X', 'W', 'K', 'GAMMA']
        self.band_structure_points = list(zip(trash_bs_points, trash_bs_labels))
        self.show()

    def update_all(self):
        self.scf_option_widget.set_all_entries()
        self.general_option_widget.set_all_entries()
        self.gw_option_widget.set_all_entries()
        self.bs_option_widget.set_all_entries()
        self.optical_spectrum_option_widget.set_all_entries()

    def do_select_event(self):
        pass

    def check_engine_for_compatibility(self,tasks_in):
        tasks = [x for x in tasks_in]
        if 'g0w0 bands' in tasks:
            tasks.remove('g0w0 bands')
        struc_type = type(self.parent.crystal_structure)
        if struc_type == sst.MolecularStructure:
            sym_type = 'non-periodic'
        elif struc_type == sst.CrystalStructure:
            sym_type = 'periodic'
        else:
            raise ValueError('bad type: '+ str(struc_type)+' for structure object')
        if sym_type not in esc_handler.supported_methods:
            raise Exception(sym_type+' structures are not supported in the selected dft engine')
        for task in tasks:
            if task not in esc_handler.supported_methods:
                raise Exception(task + ' is not supported by the selected dft engine')

    def check_if_engine_is_running_and_warn_if_so(self):
        if esc_handler.custom_command_active:
            return
        if esc_handler.is_engine_running():
            raise Exception('Engine is already running')

    def read_all_option_widgets(self):
        self.scf_option_widget.read_all_entries()
        self.general_option_widget.read_all_entries()
        self.gw_option_widget.read_all_entries()
        self.bs_option_widget.read_all_entries()
        self.optical_spectrum_option_widget.read_all_entries()

    def prepare_start(self,tasks):
        self.abort_bool = False
        self.check_if_engine_is_running_and_warn_if_so()
        self.check_engine_for_compatibility(tasks)
        self.read_all_option_widgets()

    def start_ground_state_calculation(self):
        tasks = []
        if esc_handler.will_scf_run():
            tasks.append('scf')

        bs_checkers = self.bs_option_widget.read_checkbuttons()
        if bs_checkers['Calculate'] and type(self.parent.crystal_structure) is sst.CrystalStructure:
            bs_points = self.band_structure_points
            tasks.append('bandstructure')
        else:
            bs_points = None
        try:
            self.prepare_start(tasks)
            esc_handler.start_ground_state(self.parent.crystal_structure, band_structure_points=bs_points)
            QtCore.QTimer.singleShot(1000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following:<br> 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_relax(self):
        tasks = ['relax']
        try:
            self.prepare_start(tasks)
            esc_handler.start_relax(self.parent.crystal_structure)
            QtCore.QTimer.singleShot(1000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_gw(self):
        tasks = []
        if esc_handler.will_scf_run():
            tasks.append('scf')
        bs_checkers = self.bs_option_widget.read_checkbuttons()
        if bs_checkers['Calculate']:
            tasks.append('bandstructure')
        tasks.extend(['g0w0','g0w0 bands'])
        try:
            self.prepare_start(tasks)
            esc_handler.start_gw(self.parent.crystal_structure,self.band_structure_points)
            QtCore.QTimer.singleShot(1000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_phonons(self):
        tasks = ['phonons']
        try:
            self.prepare_start(tasks)
            esc_handler.start_phonon(self.parent.crystal_structure, self.band_structure_points)
            QtCore.QTimer.singleShot(2000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_optical_spectrum_calculation(self):
        tasks = ['optical spectrum']
        try:
            self.prepare_start(tasks)
            esc_handler.start_optical_spectrum(self.parent.crystal_structure)
            QtCore.QTimer.singleShot(2000,lambda: self.parent.check_engine(tasks))
        except Exception as e:
            error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + repr(
                e) + '<br><br>Try following<br>: 1.Check if the selected dft engine is correctly installed<br>' \
                     '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
            self.execute_error_dialog.showMessage(error_message)
        else:
            self.parent.status_bar.set_engine_status(True)

    def abort_calculation(self):
        self.abort_bool = True
        esc_handler.kill_engine()


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


class InfoWindow(QtGui.QWidget):
    def __init__(self,parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.text_widget = QtGui.QTextBrowser(parent=self)
        self.text_widget.setFontFamily('monospace')
        layout.addWidget(self.text_widget)
        self.vertical_scrollbar = self.text_widget.verticalScrollBar()
        self.last_text = ''

        self.combobox = QtGui.QComboBox(self)
        layout.addWidget(self.combobox)
        self.combobox.addItem('Output')
        self.combobox.addItem('Input')

        self.combobox.setCurrentIndex(0)
        self.combobox.currentIndexChanged.connect(self.do_select_event)
        self.combobox.setMaximumWidth(150)


    def do_select_event(self):
        try:
            files = [esc_handler.info_file,esc_handler.input_filename]
            file = files[self.combobox.currentIndex()]
            self.update_text(esc_handler.project_directory+esc_handler.working_dirctory+file)
        except IOError:
            self.text_widget.setHtml('')
            self.last_text = ''

    def update_text(self,filename):
        cur_pos = self.vertical_scrollbar.value()
        with open(filename,'r') as f:
            text = f.read()
        if text == self.last_text:
            return

        self.last_text = text

        # text = text.replace('\n','<br>')
        self.text_widget.setPlainText(text)
        self.vertical_scrollbar.setValue(cur_pos)


class PlotWithTreeview(QtGui.QWidget):
    def __init__(self,Visualizer,data_dictionary,parent=None):
        QtGui.QWidget.__init__(self)
        self.parent = parent
        self.data_dictionary = data_dictionary
        self.layout = QtGui.QHBoxLayout(self)
        self.plot_widget = Visualizer(parent=self)
        self.treeview = QtGui.QTreeWidget(parent=self)
        self.treeview.setMaximumWidth(200)
        self.treeview.setHeaderHidden(True)
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed)

        self.treeview.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeview.customContextMenuRequested.connect(self.openMenu)

        self.layout.addWidget(self.plot_widget)
        self.layout.addWidget(self.treeview)

        self.show()

    def delete_selected_item(self):
        index = self.treeview.selectedIndexes()[0]
        item = self.treeview.itemFromIndex(index)
        bs_name = item.text(0)
        del self.data_dictionary[bs_name]
        self.update_tree()

    # def rename_selected_item(self):
    #     index = self.treeview.selectedIndexes()[0]
    #     item = self.treeview.itemFromIndex(index)
    #     bs_name = item.text(0)


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
        # menu.addAction('Rename', self.rename_selected_item)
        menu.exec_(self.treeview.viewport().mapToGlobal(position))


    def handle_item_changed(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return
        item = self.treeview.itemFromIndex(indexes[0])
        bs_name = item.text(0)
        self.plot_widget.plot(self.data_dictionary[bs_name])

    def add_result_key(self, title):
        item = QtGui.QTreeWidgetItem(self.treeview.invisibleRootItem(), [title])
        return item

    def clear_treeview(self):
        self.treeview.clear()

    def update_tree(self):
        self.treeview.clear()
        for key,value in OrderedDict(sorted(self.data_dictionary.items())).items():
            self.add_result_key(key)

    def do_select_event(self):
        self.update_tree()


class ChooseEngineWindow(QtGui.QDialog):

    def __init__(self, parent, defaults):
        super(ChooseEngineWindow, self).__init__(parent=parent)
        self.setWindowTitle('Please choose a DFT engine')
        self.defaults = defaults

        self.parent = parent
        self.resize(900,700)
        self.handlers = general_handler.handlers
        self.selected_handler = None
        self.success = False

        main_layout = QtGui.QVBoxLayout(self)

        info_widget = QtGui.QWidget(parent=self)
        main_layout.addWidget(info_widget)
        layout = QtGui.QHBoxLayout(info_widget)

        self.treeview = QtGui.QTreeWidget(parent=self)
        self.treeview.setMaximumWidth(200)
        self.treeview.setHeaderHidden(True)
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed)

        self.update_tree()
        layout.addWidget(self.treeview)

        self.text_widget = QtGui.QTextBrowser(parent=self)
        self.text_widget.setOpenExternalLinks(True)
        # self.text_widget.setReadOnly(True)
        layout.addWidget(self.text_widget)
        self.vertical_scrollbar = self.text_widget.verticalScrollBar()

        options_widget = QtGui.QWidget(parent=self)
        main_layout.addWidget(options_widget)

        options_layout = QtGui.QHBoxLayout(options_widget)
        self.remember_checkbox = QtGui.QCheckBox('Remember choice', parent=options_widget)
        options_layout.addWidget(self.remember_checkbox)

        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        main_layout.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.accept_own)
        self.buttonBox.rejected.connect(self.reject_own)

    def add_result_key(self, title):
        item = QtGui.QTreeWidgetItem(self.treeview.invisibleRootItem(), [title])
        return item

    def update_tree(self):
        self.treeview.clear()
        for key,value in self.handlers.items():
            item = self.add_result_key(key)
            if general_handler.is_handler_available(key):
                color = QtGui.QColor("green")
            else:
                color = QtGui.QColor("red")
            item.setForeground(0, QtGui.QBrush(color))

    def handle_item_changed(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return
        item = self.treeview.itemFromIndex(indexes[0])
        bs_name = item.text(0)

        self.selected_handler = self.handlers[bs_name]()
        self.selected_handler_class = self.handlers[bs_name]

        if general_handler.is_handler_available(bs_name):
            install_text = '<p style="color:Green;font-weight:bold">{} installation found in: {}</p>'.format(bs_name,self.selected_handler.find_engine_folder())
        else:
            install_text = '<p style="color:Red;font-weight:bold">No {} installation found</p>'.format(bs_name)

        method_descriptions = [self.selected_handler.supported_methods.get_description(x) for x in self.selected_handler.supported_methods]
        text = install_text+'Supported methods:\n\n - '+'\n - '.join(method_descriptions)+'\n\n'+self.selected_handler.info_text

        text = text.replace('\n', '<br>')
        self.text_widget.setHtml(text)

    def accept_own(self):
        if self.selected_handler is None:
            return

        if not general_handler.is_handler_available(self.selected_handler.engine_name):
            reply = QtGui.QMessageBox.question(self, 'Engine not installed', "The selected dft engine seems not to be installed. "
                                  "The program will not be able to calculate any electronic properties. "
                                  "You can however still visualize structures. Are you sure to proceed?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return

        global esc_handler
        esc_handler = self.selected_handler
        global Handler
        Handler = self.selected_handler_class


        if self.remember_checkbox.checkState():
            self.defaults['default engine'] = esc_handler.engine_name
        self.success = True
        self.close()

    def reject_own(self):
        sys.exit()
        # self.reject()

    def closeEvent(self, event):
        event.accept()

    def link(self, linkStr):
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))


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

    def set_engine_status(self,status,tasks=None):
        if status:
            if tasks:
                tasks_string = ', '.join(tasks)
                tasks_string2 = ' with Tasks: '+tasks_string
            else:
                tasks_string2 = ''
            self.status_label.setText(self.running_text+tasks_string2)
        else:
            self.status_label.setText(self.not_running_text)


class EngineOptionsDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(EngineOptionsDialog, self).__init__(parent)

        self.parent = parent
        self.command_filename = ''

        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Apply)

        self.buttonBox.accepted.connect(self.accept_own)
        self.buttonBox.rejected.connect(self.reject_own)
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(self.apply)


        self.grid_layout_widget =QtGui.QWidget(self)
        self.grid_layout = QtGui.QGridLayout(self.grid_layout_widget)

        self.custom_command_checkbox = QtGui.QCheckBox('Use custom command', parent=self)
        self.grid_layout.addWidget(self.custom_command_checkbox, 0, 0, 2, 1)


        self.load_custom_command_button = QtGui.QPushButton('Select command file',self)
        self.load_custom_command_button.setFixedWidth(150)
        self.load_custom_command_button.setFixedHeight(30)
        self.load_custom_command_button.clicked.connect(self.load_custom_command)
        self.grid_layout.addWidget(self.load_custom_command_button,0, 1, 2, 1)

        self.filename_label = QtGui.QLabel(self.grid_layout_widget)
        self.grid_layout.addWidget(self.filename_label, 2, 0, 1, 2)

        self.species_path_entry = EntryWithLabel(self,'Dft engine path')
        self.grid_layout.addWidget(self.species_path_entry, 3, 0, 1, 2)

        self.verticalLayout = QtGui.QVBoxLayout(self)
        self.verticalLayout.addWidget(self.grid_layout_widget)
        self.verticalLayout.addWidget(self.buttonBox)

    def apply(self):
        self.parent.project_properties['custom command'] = self.filename_label.text()
        self.parent.project_properties['custom command active'] = bool(self.custom_command_checkbox.checkState())
        esc_handler.custom_command = self.filename_label.text()
        esc_handler.custom_command_active = bool(self.custom_command_checkbox.checkState())
        species_path = self.species_path_entry.get_text()
        if len(species_path) > 0:
            self.parent.project_properties['custom dft folder'] = species_path
            esc_handler.dft_installation_folder = species_path

    def accept_own(self):
        self.apply()
        self.close()

    def reject_own(self):
        self.reject()

    def load_custom_command(self):
        self.custom_command_checkbox.setEnabled(True)
        file_dialog = QtGui.QFileDialog()
        file_dialog.setNameFilters(["sh script (*.sh)", "All (*.*)"])

        if file_dialog.exec_():
            file_name = file_dialog.selectedFiles()
            if type(file_name) == list or type(file_name) is tuple:
                file_name = file_name[0]
            if len(file_name) == 0:
                return
            self.filename_label.setText(file_name)

    def update_all(self):

        if not self.parent.project_properties['custom command']:
            self.custom_command_checkbox.setEnabled(False)
        else:
            if self.parent.project_properties['custom command active']:
                if not self.custom_command_checkbox.checkState():
                    self.custom_command_checkbox.toggle()
        self.species_path_entry.set_text(esc_handler.dft_installation_folder)
        self.filename_label.setText(self.parent.project_properties['custom command'])


class OptionWithTreeview(PlotWithTreeview):
    def __init__(self,side_panel,data_dictionary,parent=None):
        super(OptionWithTreeview, self).__init__(side_panel,data_dictionary,parent)
        self.add_result_key('None')

    def handle_item_changed(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return
        item = self.treeview.itemFromIndex(indexes[0])
        bs_name = item.text(0)

        plot_options = self.plot_widget.get_options()

        if main.mayavi_widget.visualization.cp is not None:
            main.mayavi_widget.update_plot(keep_view=True)
        if bs_name != 'None':
            main.mayavi_widget.visualization.plot_density((self.data_dictionary[bs_name]),**plot_options)

    def update_tree(self):
        self.treeview.clear()
        for key,value in OrderedDict(sorted(self.data_dictionary.items())).items():
            self.add_result_key(key)
        self.add_result_key('None')


class SliderWithEntry(QtGui.QWidget):
    def __init__(self,parent=None,label=None,limits=[0,1],value=None):
        super(SliderWithEntry, self).__init__(parent)
        # self.horizontalLayoutWidget.setGeometry(QtCore.QRect(90, 150, 160, 31))
        self.limits = limits
        if value is None:
            self.value = limits[0]
        if value < limits[0] or value > limits[1]:
            raise ValueError('Value must be within bounds')

        self.value = value
        self.horizontalLayout = QtGui.QGridLayout(self)

        if label is not None:
            self.label = QtGui.QLabel(label)
            self.horizontalLayout.addWidget(self.label,0,0)
            counter = 1
        else:
            counter = 0

        limit_range = limits[1]-limits[0]

        self.horizontalSlider = QtGui.QSlider(self)
        self.horizontalSlider.setMinimumWidth(200)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setValue((self.value-limits[0])/limit_range*100)
        QtCore.QObject.connect(self.horizontalSlider, QtCore.SIGNAL('valueChanged(int)'), self.change_text)

        self.horizontalLayout.addWidget(self.horizontalSlider,0,counter)

        self.lineEdit = QtGui.QLineEdit(self)
        self.lineEdit.setText("{0:1.1f}".format(self.value))
        self.horizontalLayout.addWidget(self.lineEdit,0,1+counter)
        self.horizontalLayout.setColumnStretch(0+counter, 3)
        self.horizontalLayout.setColumnStretch(1+counter, 1)


    def change_text(self):
        val = self.horizontalSlider.value()*(self.limits[1]-self.limits[0])/100+self.limits[0]
        self.lineEdit.setText("{0:1.1f}".format(val))
        self.value = val

    # def change_slider(self):
    #     try:
    #         val = float(self.lineEdit.text())
    #     except:
    #         return
    #     self.value = val
    #     self.horizontalSlider.setValue(val)

    def get_value(self):
        try:
            val = float(self.lineEdit.text())
        except:
            val = self.value
        return val


class KsStatePlotOptionWidget(QtGui.QWidget):

    def __init__(self,parent):
        super(KsStatePlotOptionWidget, self).__init__(parent)
        self.parent = parent
        self.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        self.verticalLayoutWidget = QtGui.QWidget(self)
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)

        self.opacity_slider = SliderWithEntry(self.verticalLayoutWidget,label='Opacity',limits=[0,1],value=0.5)
        self.verticalLayout.addWidget(self.opacity_slider)

        self.contours_entry = EntryWithLabel(self,'Contours:','10')
        self.verticalLayout.addWidget(self.contours_entry)

        self.transparent_checkbox = QtGui.QCheckBox('Transparent')
        self.transparent_checkbox.toggle()
        self.verticalLayout.addWidget(self.transparent_checkbox)


        self.colormap_combobox = QtGui.QComboBox(self)
        self.verticalLayout.addWidget(self.colormap_combobox)
        for el in colormap_list:
            self.colormap_combobox.addItem(el)
        index = self.colormap_combobox.findText('hot', QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.colormap_combobox.setCurrentIndex(index)
        self.colormap_combobox.currentIndexChanged.connect(self.parent.handle_item_changed)


        button_frame = QtGui.QWidget(self)
        self.button_layout = QtGui.QHBoxLayout(button_frame)
        self.verticalLayout.addWidget(button_frame)
        self.button_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.apply_button = QtGui.QPushButton('Apply')
        self.apply_button.setFixedSize(100,50)
        self.apply_button.clicked.connect(self.parent.handle_item_changed)
        self.button_layout.addWidget(self.apply_button)

    def get_options(self):
        opacity = self.opacity_slider.get_value()

        contour_str = self.contours_entry.get_text()
        if ',' in contour_str or '.' in contour_str:
            contours = contour_str.split(',')
            contours = [float(x) for x in contours]
        else:
            contours = int(contour_str)

        transparent = bool(self.transparent_checkbox.checkState())
        colormap = colormap_list[self.colormap_combobox.currentIndex()]
        out_dic = {'opacity':opacity,'contours':contours,'transparent':transparent,'colormap':colormap}


        return out_dic
        # self.test_label = QtGui.QLabel('asdas')
        # self.verticalLayout.addWidget(self.test_label)

        # self.verticalLayoutWidget.show()


class KsStateWindow(QtGui.QDialog):
    def __init__(self,parent):
        super(KsStateWindow, self).__init__(parent)
        self.calc_queue = queue.Queue()
        self.current_calc_properties = {}
        self.setFixedSize(700,500)
        self.parent = parent
        self.main_widget = QtGui.QWidget(parent=self)
        self.layout = QtGui.QVBoxLayout(self)
        self.calc_ks_group = QtGui.QGroupBox(parent=self.main_widget)
        self.calc_ks_group.setTitle('Calculation Options')
        self.layout.addWidget(self.calc_ks_group)

        self.sub_layout = QtGui.QGridLayout(self.calc_ks_group)

        self.k_point_entry = EntryWithLabel(self.calc_ks_group,'k point')
        self.sub_layout.addWidget(self.k_point_entry,0,0)

        self.n_band_entry = EntryWithLabel(self.calc_ks_group,'Band index')
        self.sub_layout.addWidget(self.n_band_entry,0,1)

        self.label_entry = EntryWithLabel(self.calc_ks_group,'Label')
        self.sub_layout.addWidget(self.label_entry,0,2)

        button_frame = QtGui.QWidget(self.calc_ks_group)
        self.sub_layout.addWidget(button_frame,2,0,1,0)

        button_layout = QtGui.QHBoxLayout(button_frame)

        self.calculate_button = QtGui.QPushButton('Calculate KS State',button_frame)
        self.calculate_button.setFixedWidth(150)
        self.calculate_button.setFixedHeight(50)
        self.calculate_button.clicked.connect(self.calculate_ks_state)
        button_layout.addWidget(self.calculate_button)

        self.calculate_density_button = QtGui.QPushButton('Calculate\nelectron density',button_frame)
        self.calculate_density_button.setFixedWidth(150)
        self.calculate_density_button.setFixedHeight(50)
        self.calculate_density_button.clicked.connect(self.calculate_electron_density)
        button_layout.addWidget(self.calculate_density_button)

        self.choose_nk_button = QtGui.QPushButton(button_frame)
        self.choose_nk_button.setFixedWidth(150)
        self.choose_nk_button.setFixedHeight(50)
        self.choose_nk_button.clicked.connect(self.choose_nk)
        button_layout.addWidget(self.choose_nk_button)

        self.plot_group = QtGui.QGroupBox(parent=self.main_widget)
        self.plot_group.setTitle('Plot Options')
        self.layout.addWidget(self.plot_group)

        self.plot_widget = OptionWithTreeview(KsStatePlotOptionWidget,self.parent.ks_densities,parent=self)

        self.sub_layout2 = QtGui.QVBoxLayout(self.plot_group)
        self.sub_layout2.addWidget(self.plot_widget)

        self.plot_widget.update_tree()

    def calculate_electron_density(self):
        esc_handler.calculate_electron_density(self.parent.crystal_structure)
        self.current_calc_properties['type'] = 'density'
        self.current_calc_properties['label'] = self.label_entry.get_text()
        QtCore.QTimer.singleShot(100, self.check_engine)


    def calculate_ks_state(self):
        n_band_str = self.n_band_entry.get_text()
        if ',' in n_band_str:
            tr = n_band_str.split(',')
            n_band = [int(x) for x in tr]
        elif '-' in n_band_str:
            n_band_split = n_band_str.split('-')
            n_band = list(range(int(n_band_split[0]),1+int(n_band_split[1])))
        else:
            n_band = [int(n_band_str)]

        k_band_str = self.k_point_entry.get_text()
        if ',' in k_band_str:
            tr = k_band_str.split(',')
            k = [int(x) for x in tr]
        else:
            k = [int(k_band_str)]

        if len(k) == 1:
            k = k*len(n_band)
        if len(n_band) == 1:
            n_band = n_band*len(k)

        label = self.label_entry.get_text()
        label_list = [label]*len(k)

        to_do_list = zip(k,n_band,label_list)
        map(self.calc_queue.put,to_do_list)
        self.start_ks_calculation()


    def choose_nk(self):
        pass

    def start_ks_calculation(self):
        k,n_band,label = self.calc_queue.get()
        self.current_calc_properties = {'type':'ks density','k':k,'n_band':n_band,'label':label}
        esc_handler.calculate_ks_density(self.parent.crystal_structure,[k,n_band])
        QtCore.QTimer.singleShot(100, self.check_engine)

    def check_engine(self):
        tasks=['ks density']
        if esc_handler.is_engine_running(tasks=tasks):
            self.parent.status_bar.set_engine_status(True,tasks=tasks)
            QtCore.QTimer.singleShot(100, self.check_engine)
        else:
            self.parent.status_bar.set_engine_status(False)
            message, err = esc_handler.engine_process.communicate()
            if ('error' in message.lower() or len(err)>0):
                error_message = 'DFT calculation finished with an error:<br><br>' + message.replace('\n','<br>')+'<br>Error:<br>'+err.replace('\n','<br>') \
                                + '<br><br>Try following:<br>1.Check if the selected dft engine is correctly installed<br>' \
                                  '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
                self.parent.error_dialog.showMessage(error_message)

            ks_dens = esc_handler.read_ks_state()
            label = self.current_calc_properties['label']
            if self.current_calc_properties['type'] == 'ks density':
                n_band = self.current_calc_properties['n_band']
                k = self.current_calc_properties['k']
                key = "{} k{} n{}".format(label,k,n_band)
            elif self.current_calc_properties['type'] == 'density':
                key = label +' density'

            if ks_dens is not None:
                self.parent.ks_densities[key] = ks_dens
                self.plot_widget.update_tree()

            if not self.calc_queue.empty():
                self.start_ks_calculation()


class MainWindow(QtGui.QMainWindow):
    def __init__(self, central_window, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.central_window = central_window

    def closeEvent(self, event):
        self.central_window.close_application()
        event.ignore()


class EditStructureWindow(QtGui.QDialog):
    "TODO: Clean unit cell correctly for new structure. Make the buttons to really set the structure"
    def __init__(self,parent):
        super(EditStructureWindow, self).__init__(parent)
        self.setWindowTitle('Edit Structure')
        self.setFixedSize(650, 700)
        self.parent = parent
        self.anything_changed = False

        self.crystal_structure = None
        self.number_of_atoms = 1

        self.main_layout = QtGui.QHBoxLayout(self)

        self.structure_widget = QtGui.QWidget(self)
        self.main_layout.addWidget(self.structure_widget)

        self.verticalLayout = QtGui.QVBoxLayout(self.structure_widget)
        self.unit_cell_box = QtGui.QGroupBox(self.structure_widget)
        self.unit_cell_box.setTitle('Unit Cell')
        self.verticalLayout.addWidget(self.unit_cell_box)


        self.unit_cell_layout = QtGui.QVBoxLayout(self.unit_cell_box)
        self.unit_cell_layout.setAlignment(QtCore.Qt.AlignTop)

        unit_cell_option_widget = QtGui.QWidget(self.unit_cell_box)
        # unit_cell_option_widget.setFixedHeight(10)
        self.unit_cell_option_layout = QtGui.QHBoxLayout(unit_cell_option_widget)
        self.unit_cell_option_layout.setAlignment(QtCore.Qt.AlignLeft)
        self.unit_cell_layout.addWidget(unit_cell_option_widget)

        self.scale_entry = EntryWithLabel(self,'Scale')
        self.scale_entry.setFixedHeight(50)
        self.unit_cell_option_layout.addWidget(self.scale_entry)
        self.scale_entry.set_text('1.0')
        self.scale_entry.connect_editFinished(self.handle_change)

        self.periodic_checkbox = QtGui.QCheckBox('Periodic',parent=self)
        self.periodic_checkbox.toggle()
        self.periodic_checkbox.stateChanged.connect(self.handle_change)
        self.unit_cell_option_layout.addWidget(self.periodic_checkbox)

        self.unit_cell_table =  QtGui.QTableWidget(self.unit_cell_box)
        self.unit_cell_table.setColumnCount(3)
        self.unit_cell_table.setRowCount(3)
        self.unit_cell_table.setFixedWidth(328)
        self.unit_cell_table.setFixedHeight(128)

        # copy_action_unit = CopySelectedCellsAction(self.unit_cell_table)
        # self.unit_cell_table.addAction(copy_action_unit)

        self.unit_cell_layout.addWidget(self.unit_cell_table)

        item = QtGui.QTableWidgetItem()
        self.unit_cell_table.setHorizontalHeaderItem(0, item)
        item.setText('x')

        item = QtGui.QTableWidgetItem()
        self.unit_cell_table.setHorizontalHeaderItem(1, item)
        item.setText('y')

        item = QtGui.QTableWidgetItem()
        self.unit_cell_table.setHorizontalHeaderItem(2, item)
        item.setText('z')

        for i in range(3):
            for j in range(3):
                item = QtGui.QTableWidgetItem()
                self.unit_cell_table.setItem(i,j,item)

        self.atom_box = QtGui.QGroupBox(self.structure_widget)
        self.atom_box.setTitle('Atoms')
        self.verticalLayout.addWidget(self.atom_box)

        self.atom_layout = QtGui.QVBoxLayout(self.atom_box)

        self.atom_table = QtGui.QTableWidget(self.atom_box)
        copy_action_atoms = CopySelectedCellsAction(self.atom_table)
        self.atom_table.addAction(copy_action_atoms)
        paste_action = PasteIntoTable(self.atom_table,self)
        self.atom_table.addAction(paste_action)

        self.atom_table.setColumnCount(4)
        self.atom_table.setRowCount(1)
        self.atom_table.setFixedWidth(450)
        self.atom_table.setFixedHeight(300)
        self.make_header()
        self.atom_layout.addWidget(self.atom_table)

        self.atom_table_buttons_widget = QtGui.QWidget(self)
        self.atom_layout.addWidget(self.atom_table_buttons_widget)

        self.atom_table_buttons_layout = QtGui.QHBoxLayout(self.atom_table_buttons_widget)
        self.atom_table_buttons_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.add_atom_button = QtGui.QPushButton('Add atom',self)
        self.add_atom_button.setFixedWidth(150)
        self.add_atom_button.clicked.connect(self.add_atom)
        self.atom_table_buttons_layout.addWidget(self.add_atom_button)

        self.remove_atom_button = QtGui.QPushButton('Remove atoms',self)
        self.remove_atom_button.setFixedWidth(150)
        self.remove_atom_button.clicked.connect(self.remove_atoms)
        self.atom_table_buttons_layout.addWidget(self.remove_atom_button)

        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Apply)
        self.verticalLayout.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(self.apply)

        header = self.atom_table.horizontalHeader()
        header.setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        header.setResizeMode(1, QtGui.QHeaderView.Stretch)
        header.setResizeMode(2, QtGui.QHeaderView.Stretch)
        header.setResizeMode(3, QtGui.QHeaderView.Stretch)

        self.unit_cell_table.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        self.unit_cell_table.verticalHeader().setResizeMode(QtGui.QHeaderView.Stretch)

        self.atom_table.itemChanged.connect(self.handle_change)
        self.unit_cell_table.itemChanged.connect(self.handle_change)

    def apply(self):
        if self.anything_changed:
            crystal_structure = self.read_tables()
            if crystal_structure is not None:
                main.crystal_structure = crystal_structure


    def accept(self):
        self.apply()
        super(EditStructureWindow, self).accept()

    def reject(self):
        if self.anything_changed:
            main.mayavi_widget.update_crystal_structure(main.crystal_structure)
            main.mayavi_widget.update_plot()
        super(EditStructureWindow, self).reject()

    def make_header(self):
        item = QtGui.QTableWidgetItem()
        self.atom_table.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.atom_table.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.atom_table.setHorizontalHeaderItem(2, item)
        item = QtGui.QTableWidgetItem()
        self.atom_table.setHorizontalHeaderItem(3, item)
        item = self.atom_table.horizontalHeaderItem(0)
        item.setText("Species")
        item = self.atom_table.horizontalHeaderItem(1)
        item.setText("x")
        item = self.atom_table.horizontalHeaderItem(2)
        item.setText("y")
        item = self.atom_table.horizontalHeaderItem(3)
        item.setText("z")

    def set_structure(self,structure):
        self.crystal_structure = structure

    def clear_unit_cell_table(self):
        self.unit_cell_table.clearContents()

    def clear_atom_table(self):
        self.atom_table.clearContents()
        # self.make_header()

    def disconnect_tables(self):
        try:
            self.unit_cell_table.itemChanged.disconnect()
            self.atom_table.itemChanged.disconnect()
        except Exception as e:
            logging.exception(e)

    def connect_tables(self):
        self.unit_cell_table.itemChanged.connect(self.handle_change)
        self.atom_table.itemChanged.connect(self.handle_change)

    def add_atom(self):
        self.disconnect_tables()

        n_rows = self.atom_table.rowCount()
        self.atom_table.setRowCount(n_rows+1)
        for j in range(4):
            item = QtGui.QTableWidgetItem()
            self.atom_table.setItem(n_rows,j,item)
        self.connect_tables()

    def remove_atoms(self,*args,**kwargs):
        atoms = kwargs.pop('atoms',None)
        self.disconnect_tables()
        if atoms is None:
           atoms = sorted(set(index.row() for index in self.atom_table.selectedIndexes()))
        for atom in atoms[::-1]:
            self.atom_table.removeRow(atom)
        self.connect_tables()
        self.handle_change()

    def update_fields(self):
        self.disconnect_tables()
        self.scale_entry.set_text('1.0')
        try:
            if self.crystal_structure is None:
                self.clear_atom_table()
                self.clear_unit_cell_table()
                self.set_number_of_atoms(6)
            else:
                scale = self.crystal_structure.scale
                self.scale_entry.set_text('{0:1.6f}'.format(scale))
                if type(self.crystal_structure) is sst.CrystalStructure:
                    unit_cell = self.crystal_structure.lattice_vectors
                    for i in range(3):
                        for j in range(3):
                            self.unit_cell_table.item(i,j).setText("{0:1.6f}".format(unit_cell[i,j]/scale))
                    if not self.periodic_checkbox.checkState():
                        self.periodic_checkbox.toggle()
                elif type(self.crystal_structure) is sst.MolecularStructure:
                    if self.periodic_checkbox.checkState():
                        self.periodic_checkbox.toggle()

                n_atoms = self.crystal_structure.atoms.shape[0]
                self.set_number_of_atoms(n_atoms)
                for i,atom in enumerate(self.crystal_structure.atoms):
                    coords = atom[0:3]
                    for j,coord in enumerate(coords):
                        item = self.atom_table.item(i,j+1)
                        item.setText('{0:1.6f}'.format(coord))
                    item = self.atom_table.item(i, 0)
                    item.setText(p_table[atom[3]])
        except Exception as e:
            logging.exception(e)

        self.connect_tables()

    def set_number_of_atoms(self,N):
        self.atom_table.setRowCount(N)
        self.number_of_atoms = N
        for i in range(N):
            for j in range(4):
                item = QtGui.QTableWidgetItem()
                self.atom_table.setItem(i,j,item)

    def read_tables(self):
        try:
            scale_string = self.scale_entry.get_text()
            scale = float(scale_string)
            if scale == 0:
                raise Exception('Scale cannot be zero')
        except Exception as e:
            logging.exception(e)
            scale = 1.0
        if self.periodic_checkbox.checkState():
            try:
                unit_cell = np.zeros((3,3))
                for i in range(3):
                    for j in range(3):
                        item = self.unit_cell_table.item(i,j)
                        unit_cell[i,j] = float(item.text())

                unit_cell = unit_cell*scale
            except Exception:
                return None
        else:
            unit_cell = None

        n_rows = self.atom_table.rowCount()
        atoms = np.zeros((n_rows,4))
        for i in range(n_rows):
            a_type = self.atom_table.item(i,0).text()
            try:
                a_type = int(a_type)
                a_type_is_number = True
            except:
                a_type_is_number = False
            if a_type not in p_table_rev.keys():
                continue
            coord = np.zeros((1,3))
            skip_bool = False
            for j in range(1,4):
                try:
                    coord[0,j-1] = float(self.atom_table.item(i,j).text())
                except:
                    skip_bool = True
                    break
            if skip_bool:
                continue
            atoms[i,:3] = coord
            if not a_type_is_number:
                a_type = p_table_rev[a_type]
            atoms[i,3] = a_type

        atoms_clean = atoms[atoms[:,3]!=0,:]

        if self.periodic_checkbox.checkState():
            out_struc = sst.CrystalStructure(unit_cell,atoms_clean, scale=scale)
        else:
            out_struc = sst.MolecularStructure(atoms_clean,scale=scale)

        return out_struc

    def handle_change(self):
        self.anything_changed = True
        crystal_structure = self.read_tables()
        if crystal_structure is not None:
            main.mayavi_widget.update_crystal_structure(crystal_structure)
            main.mayavi_widget.update_plot()


class CentralWindow(QtGui.QWidget):
    def __init__(self,parent=None, *args, **kwargs):
        super(CentralWindow, self).__init__(*args, **kwargs)
        self.project_loaded = False
        self.project_directory = None
        self.parent=parent
        self.crystal_structure = None
        self.band_structures = {}
        self.optical_spectra = {}
        self.ks_densities = {}
        self.project_properties = {'title': '','dft engine':'','custom command':'','custom command active':False,'custom dft folder':''}
        self.esc_handler_options = {}

        self.installation_folder = os.path.dirname(__file__)

        self.load_defaults()

        if self.defaults['default engine'] is None:
            choose_engine_window = ChooseEngineWindow(self,self.defaults)
            choose_engine_window.exec_()
            if not choose_engine_window.success:
                sys.exit()
        else:
            global esc_handler
            global Handler
            bs_name = self.defaults['default engine']
            esc_handler = general_handler.handlers[bs_name]()
            Handler = general_handler.handlers[bs_name]

        self.error_dialog = QtGui.QErrorMessage(parent=self)
        self.error_dialog.resize(700, 600)

        self.layout = QtGui.QGridLayout(self)

        self.mayavi_widget = MayaviQWidget(self.crystal_structure, parent=self)
        self.band_structure_window = PlotWithTreeview(Visualizer=BandStructureVisualization, data_dictionary=self.band_structures, parent=self)
        self.optical_spectra_window = PlotWithTreeview(Visualizer=OpticalSpectrumVisualization, data_dictionary=self.optical_spectra, parent=self)
        self.dft_engine_window = DftEngineWindow(self)
        self.scf_window = ScfWindow(parent=self)
        self.info_window = InfoWindow(parent=self)

        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.currentChanged.connect(self.tab_is_changed)
        self.layout.addWidget(self.tabWidget)

        self.status_bar = StatusBar()
        self.layout.addWidget(self.status_bar)

        self.engine_option_window = EngineOptionsDialog(self)
        self.ks_state_window = KsStateWindow(self)
        self.structure_window = EditStructureWindow(self)
        self.brillouin_window = BrillouinWindow(self)

        self.tab_layout = QtGui.QVBoxLayout()
        self.tabWidget.setLayout(self.tab_layout)

        self.list_of_tabs = [self.mayavi_widget,self.dft_engine_window,self.band_structure_window,self.optical_spectra_window,self.scf_window,self.info_window]


        self.tabWidget.addTab(self.list_of_tabs[0], 'Structure')
        self.tabWidget.addTab(self.list_of_tabs[1], 'DFT-Engine')
        self.tabWidget.addTab(self.list_of_tabs[2],'Bandstructure')
        self.tabWidget.addTab(self.list_of_tabs[3],'Optical Spectrum')
        self.tabWidget.addTab(self.list_of_tabs[4], 'Scf')
        self.tabWidget.addTab(self.list_of_tabs[5], 'Info')

        self.tabWidget.show()

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
                # project_directory = r"/home/jannick/exciting_cluster/GaN"
            else:
                project_directory = r'D:\OpenDFT_projects\test'
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
        self.project_properties.update({'title': '','dft engine':'','custom command':'','custom command active':False,'custom dft folder':''})
        self.window.setWindowTitle("OpenDFT - " + self.project_directory)
        os.chdir(self.project_directory)
        if (esc_handler.pseudo_directory is not None) and (not os.path.isdir(self.project_directory + esc_handler.pseudo_directory)):
            os.mkdir(self.project_directory + esc_handler.pseudo_directory)
        self.project_loaded = True

    def reset_results_and_plots(self):
        self.crystal_structure = None
        esc_handler.reset_to_defaults()
        self.project_properties.clear()
        for key, value in self.band_structures.items():
            del self.band_structures[key]
        for key, value in self.optical_spectra.items():
            del self.optical_spectra[key]
        for key, value in self.ks_densities.items():
            del self.ks_densities[key]
        self.mayavi_widget.visualization.clear_plot()
        self.band_structure_window.plot_widget.clear_plot()
        self.band_structure_window.clear_treeview()
        self.scf_window.scf_widget.clear_plot()
        self.optical_spectra_window.plot_widget.clear_plot()
        self.optical_spectra_window.clear_treeview()
        self.dft_engine_window.update_all()

    def load_project(self,*args,**kwargs):
        folder_name = kwargs.pop('folder_name', None)
        if folder_name is None:
            folder_name = QtGui.QFileDialog().getExistingDirectory(parent=self)
        if len(folder_name) > 1:
            self.reset_results_and_plots()
            self.project_directory = folder_name
            os.chdir(self.project_directory)
            esc_handler.project_directory = self.project_directory
            self.load_saved_results()
            self.dft_engine_window.update_all()
            self.window.setWindowTitle("OpenDFT - "+self.project_directory)
            self.project_loaded = True

    def save_results(self):
        # TODO save options for each handler seperately
        try:
            self.dft_engine_window.read_all_option_widgets()

            option_dic_specific_handler = {'scf_options':esc_handler.scf_options,'general options':esc_handler.general_options,'bs options':esc_handler.bs_options,
                 'phonon options':esc_handler.phonons_options,'optical spectrum options':esc_handler.optical_spectrum_options,
                 'gw options':esc_handler.gw_options}
            self.esc_handler_options[esc_handler.engine_name] = option_dic_specific_handler
            a = {'crystal structure': self.crystal_structure, 'band structure': self.band_structures, 'optical spectra':self.optical_spectra,'esc handler options': self.esc_handler_options,
                 'properties': self.project_properties,'dft engine':esc_handler.engine_name,'ks densities':self.ks_densities,'k path':self.dft_engine_window.band_structure_points}
            with open(self.project_directory + '/save.pkl', 'wb') as handle:
                pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as e:
            logging.exception(e)

    def load_saved_results(self):
        esc_handler.project_directory = self.project_directory
        try:
            with open(self.project_directory + '/save.pkl', 'rb') as handle:
                b = pickle.load(handle)

                self.crystal_structure = b.pop('crystal structure',None)
                if self.crystal_structure is not None:
                    self.mayavi_widget.update_crystal_structure(self.crystal_structure)
                    self.mayavi_widget.update_plot()

                loaded_bandstructure_dict = b.pop('band structure',None)
                if type(loaded_bandstructure_dict) == dict:
                    for key,value in loaded_bandstructure_dict.items():
                        self.band_structures[key] = value

                loaded_optical_spectra_dict = b.pop('optical spectra',None)
                if type(loaded_optical_spectra_dict) == dict:
                    for key,value in loaded_optical_spectra_dict.items():
                        self.optical_spectra[key] = value

                loaded_ksdens_dict = b.pop('ks densities',None)
                if type(loaded_ksdens_dict ) == dict:
                    for key,value in loaded_ksdens_dict.items():
                        self.ks_densities[key] = value

                self.esc_handler_options = b.pop('esc handler options',None)

                option_dic_specific_handler = self.esc_handler_options.pop(esc_handler.engine_name,None)
                if option_dic_specific_handler is not None:

                    load_scf_options = option_dic_specific_handler.pop('scf_options', None)
                    if load_scf_options is not None:
                        for key,value in load_scf_options.items():
                            esc_handler.scf_options[key] = value

                    load_general_options = option_dic_specific_handler.pop('general options',None)
                    if load_general_options is not None:
                        for key,value in load_general_options.items():
                            esc_handler.general_options[key] = value

                    load_bs_options = option_dic_specific_handler.pop('bs options',None)
                    if load_bs_options is not None:
                        for key,value in load_bs_options.items():
                            esc_handler.bs_options[key] = value

                    load_phonon_options = option_dic_specific_handler.pop('phonon options',None)
                    if load_phonon_options is not None:
                        for key,value in load_phonon_options.items():
                            esc_handler.phonons_options[key] = value

                    load_gw_options = option_dic_specific_handler.pop('gw options',None)
                    if load_gw_options is not None:
                        for key,value in load_gw_options.items():
                            esc_handler.gw_options[key] = value

                    load_optical_spectrum_options = option_dic_specific_handler.pop('optical spectrum options',None)
                    if load_optical_spectrum_options is not None:
                        for key,value in load_optical_spectrum_options.items():
                            esc_handler.optical_spectrum_options[key] = value

                k_path = b.pop('k path', None)
                if k_path is not None:
                    self.dft_engine_window.band_structure_points = k_path

                self.project_properties.update(b['properties'])
                ## Update esc_handler ! DANGER ZONE !
                try:
                    esc_handler.custom_command_active = self.project_properties['custom command active']
                    esc_handler.custom_command = self.project_properties['custom command']
                    if self.project_properties['custom dft folder']:
                        esc_handler.dft_installation_folder = self.project_properties['custom dft folder']
                except:
                    self.project_properties['custom command active'] = False
                    self.project_properties['custom command'] = ''
                    self.project_properties['custom dft folder'] = ''

        except IOError:
            print('file not found')

    def load_defaults(self):
        self.defaults = {}

        try:
            with open(self.installation_folder + '/defaults.pkl', 'rb') as handle:
                b = pickle.load(handle)
        except IOError:
            logging.info('Default file not found')
            b = {}

        default_engine = b.pop('default engine',None)

        self.defaults['default engine'] = default_engine

    def save_defaults(self):
        with open(self.installation_folder + '/defaults.pkl', 'wb') as handle:
            pickle.dump(self.defaults, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def update_structure_plot(self):
        self.mayavi_widget.update_crystal_structure(self.crystal_structure)
        self.mayavi_widget.update_plot()

        # t = MyQThread(self.mayavi_widget.update_plot)
        # t.start()

    def load_crystal_structure(self,filetype):
        file_dialog = QtGui.QFileDialog()

        if filetype == 'exciting':
            file_dialog.setNameFilters(["Exciting (*.xml)", "All (*.*)"])
        elif filetype == 'quantum espresso':
            file_dialog.setNameFilters(["Quantum Espresso (*.in)", "All (*.*)"])
        elif filetype == 'cif':
            file_dialog.setNameFilters(["Cif (*.cif)", "All (*.*)"])
        else:
            raise Exception('bad filetype')

        if file_dialog.exec_():
            file_name = file_dialog.selectedFiles()
            if type(file_name) == list or type(file_name) is tuple:
                file_name = file_name[0]
            if len(file_name) == 0:
                return

            if filetype in ['exciting','quantum espresso']:
                self.crystal_structure = general_handler.parse_input_file(filetype,file_name)
            elif filetype == 'cif':
                parser = sst.StructureParser()
                self.crystal_structure = parser.parse_cif_file(file_name)

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

        self.file_menu.addSeparator()

        new_structure_action = QtGui.QAction("New structure", self.window)
        new_structure_action.setShortcut('Ctrl+Shift+n')
        new_structure_action.setStatusTip('Make new structure by hand')
        new_structure_action.triggered.connect(lambda: self.open_structure_window(new=True))
        self.file_menu.addAction(new_structure_action)

        edit_structure_action = QtGui.QAction("Edit structure", self.window)
        edit_structure_action.setShortcut('Ctrl+Shift+e')
        edit_structure_action.setStatusTip('Edit existing structure by hand')
        edit_structure_action.triggered.connect(lambda: self.open_structure_window(new=False))
        self.file_menu.addAction(edit_structure_action)

        import_structure_menu = self.file_menu.addMenu('Import structure from')

        open_structure_action_exciting = QtGui.QAction("exciting input file", self.window)
        open_structure_action_exciting.setStatusTip('Load crystal structure from exciting xml')
        open_structure_action_exciting.triggered.connect(lambda: self.load_crystal_structure('exciting'))
        import_structure_menu.addAction(open_structure_action_exciting)

        open_structure_action_quantum_espresso = QtGui.QAction("quantum_espresso input file", self.window)
        open_structure_action_quantum_espresso.setStatusTip('Load crystal structure from quantum espresso input file')
        open_structure_action_quantum_espresso.triggered.connect(lambda: self.load_crystal_structure('quantum espresso'))
        import_structure_menu.addAction(open_structure_action_quantum_espresso)

        open_structure_action_cif = QtGui.QAction("cif", self.window)
        open_structure_action_cif.setShortcut('Ctrl+Shift+c')
        open_structure_action_cif.setStatusTip('Load crystal structure from cif file')
        open_structure_action_cif.triggered.connect(lambda: self.load_crystal_structure('cif'))
        import_structure_menu.addAction(open_structure_action_cif)

        self.file_menu.addSeparator()

        close_app_action = QtGui.QAction("Exit", self.window)
        close_app_action.setShortcut("Ctrl+Q")
        close_app_action.setStatusTip('Leave The App')
        close_app_action.triggered.connect(self.close_application)
        self.file_menu.addAction(close_app_action)

        self.vis_menu = self.menu_bar.addMenu('&Visualize')

        ks_vis_action = QtGui.QAction("Visualize KS state", self.window)
        ks_vis_action.setStatusTip('Visualize a Kohn-Sham state in the structure window')
        ks_vis_action.triggered.connect(self.open_state_vis_window)
        self.vis_menu.addAction(ks_vis_action)

        self.dft_menu = self.menu_bar.addMenu('&DFT Engine')

        dft_options_action = QtGui.QAction("Options", self.window)
        dft_options_action.setStatusTip('Options for dft engine')
        dft_options_action.triggered.connect(self.open_engine_option_window)
        self.dft_menu.addAction(dft_options_action)

    def close_application(self):
        if not DEBUG:
            reply = QtGui.QMessageBox.question(self, 'Message',
                                               "Are you sure to quit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if DEBUG or reply == QtGui.QMessageBox.Yes:
            self.save_defaults()
            self.save_results()
            self.parent.quit()
            logging.info('Program stopped normally')
            sys.exit()

    def check_engine(self,tasks):
        def check_relax():
            new_struc = esc_handler.load_relax_structure()
            if new_struc is not None:
                self.crystal_structure = new_struc
                self.update_structure_plot()

        tasks = [x.lower() for x in tasks]
        if self.dft_engine_window.abort_bool:
            self.status_bar.set_engine_status(False)
            self.dft_engine_window.abort_bool = False
            return
        elif esc_handler.is_engine_running(tasks=tasks):
            selected_tab_index = self.tabWidget.currentIndex()
            if selected_tab_index == 4:
                self.scf_data = esc_handler.read_scf_status()
                if self.scf_data is not None:
                    self.scf_window.scf_widget.plot(self.scf_data)
            elif selected_tab_index == 5:
                self.info_window.do_select_event()
            QtCore.QTimer.singleShot(500,lambda: self.check_engine(tasks))
            self.status_bar.set_engine_status(True,tasks=tasks)
            if 'relax' in tasks:
                check_relax()
        else:
            self.scf_data = esc_handler.read_scf_status()
            if self.scf_data is not None:
                self.scf_window.scf_widget.plot(self.scf_data)
            self.status_bar.set_engine_status(False)
            message, err = esc_handler.engine_process.communicate()
            if ('error' in message.lower() or len(err)>0):
                error_message = 'DFT calculation finished with an error:<br><br>' + message.replace('\n',"<br>")+'<br>Error:<br>'+err.replace('\n','<br>') \
                                + '<br><br>Try following:<br>1.Check if the selected dft engine is correctly installed<br>' \
                                  '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
                self.error_dialog.showMessage(error_message)

            if 'scf' in tasks and type(self.crystal_structure) is sst.MolecularStructure:
                energy_diagram = esc_handler.read_energy_diagram()
                if energy_diagram is not None:
                    self.band_structures[esc_handler.general_options['title']] = energy_diagram

            if 'bandstructure' in tasks or 'g0w0' in tasks:
                read_bandstructures = []
                titles = [esc_handler.general_options['title']]
                if 'bandstructure' in tasks:
                    if esc_handler.engine_name == 'quantum espresso':
                        coords = [x[0] for x in self.dft_engine_window.band_structure_points]
                        labels = [x[1] for x in self.dft_engine_window.band_structure_points]
                        new_coords = self.crystal_structure.convert_to_tpiba(coords)
                        band_structure_points = zip(new_coords,labels)
                        read_bandstructure = esc_handler.read_bandstructure(special_k_points=band_structure_points)
                    else:
                        read_bandstructure = esc_handler.read_bandstructure()
                    read_bandstructures.append(read_bandstructure)

                if 'g0w0' in tasks:
                    read_bandstructures.append(esc_handler.read_gw_bandstructure())
                if 'bandstructure' in tasks and 'g0w0' in tasks:
                    titles.append(esc_handler.general_options['title']+'_gw')


                for read_bandstructure,title in zip(read_bandstructures,titles):
                    if read_bandstructure is None:
                        continue
                    self.band_structures[title] = read_bandstructure
                    self.band_structure_window.update_tree()
                if len(read_bandstructures)!=0 and self.band_structure_window.plot_widget.first_plot_bool:
                    self.band_structure_window.plot_widget.plot(self.band_structures[esc_handler.general_options['title']])
            if 'relax' in tasks:
                check_relax()
            if 'phonons' in tasks:
                read_bandstructure = esc_handler.read_phonon_bandstructure()
                self.band_structures[esc_handler.general_options['title'] + '_phonon'] = read_bandstructure
                self.band_structure_window.update_tree()
            if 'optical spectrum' in tasks:
                read_spectrum = esc_handler.read_optical_spectrum()
                self.optical_spectra[esc_handler.general_options['title']] = read_spectrum
                self.optical_spectra_window.update_tree()

    def open_engine_option_window(self):
        self.engine_option_window.update_all()
        self.engine_option_window.show()

    def open_state_vis_window(self):
        self.ks_state_window.plot_widget.update_tree()
        self.ks_state_window.show( )

    def open_structure_window(self,new=False):
        if new:
            self.structure_window.set_structure(None)
        else:
            self.structure_window.set_structure(self.crystal_structure)

        self.structure_window.anything_changed = False
        self.structure_window.update_fields()
        self.structure_window.show()

    def open_brillouin_window(self):
        if type(self.crystal_structure) is not sst.CrystalStructure:
            return

        if self.crystal_structure is not None and self.brillouin_window.mayavi_widget.crystal_structure is not self.crystal_structure and type(self.crystal_structure):
            self.brillouin_window.close() # This is a hack because for some reason the picker is broken when you update the plot
            self.brillouin_window = BrillouinWindow(self)
            self.brillouin_window.mayavi_widget.set_crystal_structure(self.crystal_structure)
            self.brillouin_window.set_path(self.dft_engine_window.band_structure_points)
            self.brillouin_window.mayavi_widget.update_plot()

        self.brillouin_window.show()


if __name__ == "__main__":
    DEBUG = True

    current_time = time.localtime()
    current_time_string = [str(x) for x in current_time[:3]]

    installation_folder = os.path.dirname(__file__)

    if not os.path.exists(installation_folder + "/logfiles/"):
        os.makedirs(installation_folder + "/logfiles/")

    logging.basicConfig(level=logging.DEBUG,
                        filename=installation_folder + "/logfiles/" + "_".join(current_time_string) + ".log")
    logging.info('Program started')

    app = QtGui.QApplication.instance()
    main = CentralWindow(parent=app)

    app.exec_()
