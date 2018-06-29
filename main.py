#!/usr/bin/env python3
from __future__ import division, absolute_import, print_function, unicode_literals
import six
import sys
from pathlib import Path
import psutil

if not sys.getfilesystemencoding():
    sys.getfilesystemencoding = lambda: 'UTF-8'
import os
import numpy as np
import warnings
from collections import OrderedDict

DEBUG = False
if DEBUG:
    warnings.simplefilter('always', UserWarning)
else:
    warnings.simplefilter('ignore', UserWarning)

os.environ['ETS_TOOLKIT'] = 'qt4'

import imp

try:
    imp.find_module('PySide')  # test if PySide if available
except ImportError:
    os.environ['QT_API'] = 'pyqt'  # signal to pyface that PyQt4 should be used

from pyface.qt import QtGui, QtCore
from visualization import StructureVisualization, BandStructureVisualization, ScfVisualization, \
    OpticalSpectrumVisualization, colormap_list, BrillouinVisualization, VolumeSlicer,DosVisualization,PhononVisualization
import solid_state_tools as sst
from solid_state_tools import p_table, p_table_rev
from little_helpers import CopySelectedCellsAction, PasteIntoTable, set_procname, get_proc_name, \
    find_data_file, get_stacktrace_as_string,eval_expr,find_fraction
from TerminalClass import PythonTerminal
import pickle
import time
import threading
from collections import OrderedDict
import logging
import syntax
import re
import copy

try:
    import queue
except:
    import Queue as queue

general_handler = sst.GeneralHandler()


class BrillouinWindow(QtGui.QDialog):
    def __init__(self, parent=None):
        super(BrillouinWindow, self).__init__(parent)

        self.k_path = None

        self.resize(900, 600)
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

        self.table = QtGui.QTableWidget(table_widget)

        copy_action = CopySelectedCellsAction(self.table)
        self.table.addAction(copy_action)
        paste_action = PasteIntoTable(self.table, self)
        self.table.addAction(paste_action)

        self.table.setColumnCount(4)
        self.table.setRowCount(0)
        table_layout.addWidget(self.table)

        for i, label in enumerate(['x', 'y', 'z', 'Label']):
            item = QtGui.QTableWidgetItem()
            self.table.setHorizontalHeaderItem(i, item)
            item.setText(label)

        self.table.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        self.connect_tables()

        button_widget = QtGui.QWidget(parent=table_widget)
        table_layout.addWidget(button_widget)

        button_layout = QtGui.QHBoxLayout(button_widget)

        remove_all_button = QtGui.QPushButton('Remove path', parent=button_widget)
        button_layout.addWidget(remove_all_button)
        remove_all_button.clicked.connect(self.clear_path)

        remove_selected_button = QtGui.QPushButton('Remove point', parent=button_widget)
        button_layout.addWidget(remove_selected_button)
        remove_selected_button.clicked.connect(self.remove_point)

        add_atom_button = QtGui.QPushButton('Add point', parent=button_widget)
        button_layout.addWidget(add_atom_button)
        add_atom_button.clicked.connect(self.add_atom)

        load_standard_kpath_button = QtGui.QPushButton('Standard path', parent=button_widget)
        button_layout.addWidget(load_standard_kpath_button)
        load_standard_kpath_button.clicked.connect(self.load_standard_path)

    def clear_path(self):
        for i in range(len(self.k_path)):
            del self.k_path[0]

        self.mayavi_widget.set_path(self.k_path)
        self.mayavi_widget.plot_path()
        self.update_table()

    def remove_point(self, *args, **kwargs):
        points = kwargs.pop('points', None)
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
                if j < 3:
                    text = find_fraction(self.k_path[i][0][j])
                else:
                    text = self.k_path[i][1]
                item.setText(text)
                self.table.setItem(i, j, item)
        self.connect_tables()

    def set_path(self, k_path):

        self.k_path = k_path
        self.mayavi_widget.set_path(k_path)
        self.update_table()

    def add_atom(self):
        self.disconnect_tables()
        n = self.table.rowCount()
        self.table.setRowCount(n + 1)
        for i in range(4):
            item = QtGui.QTableWidgetItem()
            self.table.setItem(n, i, item)
        self.connect_tables()

    def read_table(self):
        n_k = self.table.rowCount()
        k_points = []
        for i in range(n_k):
            try:
                coords = np.array([float(eval_expr(self.table.item(i, j).text())) for j in range(3)])
                k_point = [coords, str(self.table.item(i, 3).text())]
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

    def load_standard_path(self):

        try:
            conv_path = sst.calculate_standard_path(main.crystal_structure)
        except Exception as e:
            conv_path = sst.get_emergency_path()
            main.error_dialog.showMessage(str(e)+'<br> A standard non structure specific path will be loaded. Use with care!')

        self.clear_path()
        for el in conv_path:
            self.k_path.append(el)
        self.mayavi_widget.plot_path()
        self.update_table()


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

    def update_plot(self, keep_view=False):
        self.visualization.update_plot(keep_view=keep_view)

    def update_crystal_structure(self, crystal_structure):
        self.visualization.crystal_structure = crystal_structure

    def do_select_event(self):
        pass


class MayaviPhononWindow(QtGui.QMainWindow):
    def __init__(self, crystal_structure,data_dictionary, parent=None):
        super(MayaviPhononWindow, self).__init__(parent)

        self.resize(1000, 600)

        self.data_dictionary = data_dictionary
        self.setWindowTitle('Phonon visualization')

        self.main_widget = QtGui.QWidget(self)
        self.setCentralWidget(self.main_widget)

        main_layout = QtGui.QVBoxLayout(self.main_widget)

        self.sub_widget = QtGui.QWidget(self.main_widget)
        main_layout.addWidget(self.sub_widget)

        layout = QtGui.QHBoxLayout(self.sub_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.visualization = PhononVisualization(crystal_structure)
        self.ui = self.visualization.edit_traits(parent=self.sub_widget,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self.sub_widget)

        self.treeview = QtGui.QTreeWidget(parent=self)
        self.treeview.setMaximumWidth(200)
        self.treeview.setHeaderHidden(True)
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed_main)
        layout.addWidget(self.treeview)

        self.mode_treeview = QtGui.QTreeWidget(parent=self)
        self.mode_treeview.setMaximumWidth(60)
        self.mode_treeview.setHeaderHidden(True)
        self.mode_treeview.itemSelectionChanged.connect(self.handle_item_changed)
        layout.addWidget(self.mode_treeview)

        self.k_treeview = QtGui.QTreeWidget(parent=self)
        self.k_treeview.setMaximumWidth(150)
        self.k_treeview.setHeaderHidden(True)
        self.k_treeview.itemSelectionChanged.connect(self.handle_item_changed)
        layout.addWidget(self.k_treeview)

        self.info_widget = QtGui.QWidget(self.main_widget)
        main_layout.addWidget(self.info_widget)

        info_layout = QtGui.QHBoxLayout(self.info_widget)
        info_layout.setAlignment(QtCore.Qt.AlignLeft)

        infos = ['frequency']
        self.info_labels = {}

        for info in infos:
            new_label = LabeledLabel(info.title()+':')
            info_layout.addWidget(new_label)
            self.info_labels[info] = new_label


    def update_plot(self, keep_view=False):
        self.visualization.update_plot(keep_view=keep_view)

    def update_crystal_structure(self, crystal_structure):
        self.visualization.crystal_structure = crystal_structure

    def handle_item_changed_main(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return

        item = self.treeview.itemFromIndex(indexes[0])
        bs_name = item.text(0)
        if bs_name == 'None':
            self.k_treeview.clear()
            self.mode_treeview.clear()
            return


        phonon_eigenvectors = self.get_phonon_eigenvectors()
        self.update_mode_and_k_tree(phonon_eigenvectors)
        self.handle_item_changed()

    def get_phonon_eigenvectors(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return

        item = self.treeview.itemFromIndex(indexes[0])
        bs_name = item.text(0)

        if bs_name == 'None':
            self.visualization.remove_phonons()
            return

        data = self.data_dictionary[bs_name]
        return data.phonon_eigenvectors

    def handle_item_changed(self):

        phonon_eigenvectors = self.get_phonon_eigenvectors()
        self.visualization.phonon_eigenvectors = phonon_eigenvectors

        mode,k = self.get_mode_and_k()
        self.visualization.plot_phonons(mode,k)

        self.update_infos(phonon_eigenvectors,mode,k)


    def add_result_key(self, title):
        item = QtGui.QTreeWidgetItem(self.treeview.invisibleRootItem(), [title])
        return item

    def clear_treeview(self):
        self.treeview.itemSelectionChanged.disconnect()

        self.treeview.clear()
        self.add_result_key('None')
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed_main)

    def update_tree(self):
        self.treeview.itemSelectionChanged.disconnect()
        self.treeview.clear()
        self.add_result_key('None')
        for key, value in OrderedDict(sorted(self.data_dictionary.items())).items():
            if hasattr(value,'phonon_eigenvectors'):
                self.add_result_key(key)
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed_main)

    def update_mode_and_k_tree(self,phonon_eigenvectors):
        self.mode_treeview.clear()
        self.k_treeview.clear()

        try:
            self.mode_treeview.itemSelectionChanged.disconnect()
            self.k_treeview.itemSelectionChanged.disconnect()
        except Exception:
            pass

        modes = range(phonon_eigenvectors.modes.shape[0])
        ks = range(phonon_eigenvectors.modes.shape[1])
        k_points = phonon_eigenvectors.k_vectors

        for mode in modes:
            QtGui.QTreeWidgetItem(self.mode_treeview.invisibleRootItem(), ['{0}'.format(mode)])
        for k in ks:
            QtGui.QTreeWidgetItem(self.k_treeview.invisibleRootItem(), ['{0:1.2f} {1:1.2f} {2:1.2f}'.format(*k_points[k,:])])

        self.mode_treeview.itemSelectionChanged.connect(self.handle_item_changed)
        self.k_treeview.itemSelectionChanged.connect(self.handle_item_changed)

    def get_mode_and_k(self):
        indexes_mode = self.mode_treeview.selectedIndexes()
        if len(indexes_mode) == 0:
            mode = 0
        else:
            item = self.mode_treeview.itemFromIndex(indexes_mode[0])
            mode = int(item.text(0))

        indexes_k = self.k_treeview.selectedIndexes()
        if len(indexes_k) == 0:
            k = 0
        else:
            k = indexes_k[0].row()

        return mode,k

    def update_infos(self,eigs,mode,k):
        freqs = eigs.frequencies
        freq = freqs[mode,k]
        data = {'frequency':'{0:1.1f} cm^-1'.format(freq)}

        for label,widget in self.info_labels.items():
            if label in data.keys():
                widget.text(data[label])



class EntryWithLabel(QtGui.QWidget):
    def __init__(self, parent, label, value=None, width_text=200, width_label=90):
        QtGui.QWidget.__init__(self, parent)
        self.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.layout = QtGui.QHBoxLayout(self)
        self.textbox = QtGui.QLineEdit(self)
        if width_text:
            self.textbox.setMaximumWidth(width_text)
        self.label_widget = QtGui.QLabel(label, parent=self)
        if width_label:
            self.label_widget.setMaximumWidth(width_label)
        self.layout.setAlignment(QtCore.Qt.AlignLeft)
        self.layout.addWidget(self.label_widget)
        self.layout.addWidget(self.textbox)
        self.editFinished_command = None

        if value is not None:
            self.textbox.setText(value)

    def get_text(self):
        return self.textbox.text()

    def set_text(self, text):
        self.textbox.setText(text)

    def connect_editFinished(self, command):
        self.editFinished_command = command
        self.textbox.editingFinished.connect(self.handleEditingFinished)

    def handleEditingFinished(self):
        if self.textbox.isModified():
            self.editFinished_command()
        self.textbox.setModified(False)


class LabeledLabel(QtGui.QWidget):
    def __init__(self,label,value='',width_label=None):
        super(LabeledLabel,self).__init__()
        self.layout = QtGui.QHBoxLayout(self)
        self.label_label = QtGui.QLabel(text=label)
        self.label_label.setAlignment(QtCore.Qt.AlignLeft)
        if width_label:
            self.label_label.setFixedWidth(width_label)
        self.layout.addWidget(self.label_label)
        self.value_label = QtGui.QLabel(text=value)
        self.value_label.setAlignment(QtCore.Qt.AlignLeft)
        self.value_label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.layout.addWidget(self.value_label)

    def text(self,value):
        self.value_label.setText(value)

    def change_label(self,label):
        self.label_label.setText(label)


class ConsoleWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(ConsoleWindow, self).__init__(parent)

        self.main_widget = QtGui.QWidget(self)
        self.setCentralWidget(self.main_widget)

        self.setMinimumSize(1000, 800)
        self.resize(1200, 800)
        self.parent = parent
        self.main_layout = QtGui.QVBoxLayout(self.main_widget)
        self.make_menubar()
        self.custom_window_title = 'OpenDFT Python scripting '
        self.setWindowTitle(self.custom_window_title)
        app_icon = QtGui.QIcon()
        app_icon.addFile(find_data_file('icon.ico'), QtCore.QSize(238, 238))
        self.setWindowIcon(app_icon)
        self.error_widget = QtGui.QErrorMessage(parent=self)

        self.update_fields_timer = QtCore.QTimer()
        self.update_fields_timer.timeout.connect(self.check_queue_and_update)

        self.code_thread = threading.Thread(target=self.run_code)
        self.queue = queue.Queue()

        self.splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.main_layout.addWidget(self.splitter)

        class CodingTextEdit(QtGui.QTextEdit):
            # TODO replace tab by some nice functionality (4 spaces + multiselect)
            def __init__(self, *args, **kwargs):
                super(CodingTextEdit, self).__init__(*args, **kwargs)
                self.setStyleSheet('QTextEdit { font-size: 10pt; font-family: monospace; }')
                self.setTabStopWidth(8)

            def keyPressEvent(self, event):
                # Shift + Tab is not the same as trying to catch a Shift modifier and a tab Key.
                # Shift + Tab is a Backtab!!
                if event.key() == QtCore.Qt.Key_Backtab:
                    cur = self.textCursor()
                    # Copy the current selection
                    pos = cur.position()  # Where a selection ends
                    anchor = cur.anchor()  # Where a selection starts (can be the same as above)

                    # Can put QtGui.QTextCursor.MoveAnchor as the 2nd arg, but this is the default
                    cur.setPosition(pos)

                    # Move the position back one, selection the character prior to the original position
                    cur.setPosition(pos - 1, QtGui.QTextCursor.KeepAnchor)

                    if str(cur.selectedText()) == "\t":
                        # The prior character is a tab, so delete the selection
                        cur.removeSelectedText()
                        # Reposition the cursor with the one character offset
                        cur.setPosition(anchor - 1)
                        cur.setPosition(pos - 1, QtGui.QTextCursor.KeepAnchor)
                    else:
                        # Try all of the above, looking before the anchor (This helps if the achor is before a tab)
                        cur.setPosition(anchor)
                        cur.setPosition(anchor - 1, QtGui.QTextCursor.KeepAnchor)
                        if str(cur.selectedText()) == "\t":
                            cur.removeSelectedText()
                            cur.setPosition(anchor - 1)
                            cur.setPosition(pos - 1, QtGui.QTextCursor.KeepAnchor)
                        else:

                            # Its not a tab, so reset the selection to what it was
                            cur.setPosition(anchor)
                            cur.setPosition(pos, QtGui.QTextCursor.KeepAnchor)
                # elif event.key() == QtCore.Qt.Key_Tab:
                #     cur = self.textCursor()
                #     # Copy the current selection
                #     pos = cur.position()  # Where a selection ends
                #     anchor = cur.anchor()  # Where a selection starts (can be the same as above)
                #     selected_text = cur.selectedText()
                #     if not selected_text:
                #         selected_text = cur.block().text()
                #         cur.select(QtGui.QTextCursor.BlockUnderCursor)
                #         cur.beginEditBlock()
                #         cur.removeSelectedText()
                #         cur.insertText('  '+selected_text)
                #         cur.endEditBlock()
                #         # cur.deleteChar()
                #     else:
                #         select_list = selected_text.splitlines()
                #         new_sel = ['  '+x for x in select_list]
                #         new_text = '\n'.join(new_sel)
                #         cur.beginEditBlock()
                #         cur.removeSelectedText()
                #         cur.insertText(new_text)
                #         cur.endEditBlock()
                #
                #     cur.setPosition(anchor)
                #     cur.setPosition(pos, QtGui.QTextCursor.KeepAnchor)

                else:
                    return QtGui.QTextEdit.keyPressEvent(self, event)

        self.input_text_widget = CodingTextEdit(self)
        highlight = syntax.PythonHighlighter(self.input_text_widget.document())
        self.input_scrollbar = self.input_text_widget.verticalScrollBar()
        self.splitter.addWidget(self.input_text_widget)

        self.output_text_widget = QtGui.QTextBrowser(self)
        self.output_text_widget.setStyleSheet('QTextBrowser { font-size: 10pt; font-family: monospace; }')
        self.output_scrollbar = self.output_text_widget.verticalScrollBar()
        self.splitter.addWidget(self.output_text_widget)

        sub_frame = QtGui.QWidget(self)
        sub_layout = QtGui.QHBoxLayout(sub_frame)

        class InteractiveText(EntryWithLabel):
            def __init__(self, *args, **kwargs):
                super(InteractiveText, self).__init__(*args, **kwargs)
                self.parent = args[0]

            def keyPressEvent(self, event):
                if event.key() in [QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return]:
                    self.parent.handle_interactive_text()
                elif event.key() == QtCore.Qt.Key_Up:
                    self.parent.history_move(1)
                elif event.key() == QtCore.Qt.Key_Down:
                    self.parent.history_move(-1)
                else:
                    event.accept()

        self.interactive_text = InteractiveText(self, '>>', width_label=50, width_text=None)
        # self.interactive_text.textbox.returnPressed.connect(self.handle_interactive_text)
        self.interactive_text.textbox.setStyleSheet('font-size: 10pt; font-family: monospace; ')
        sub_layout.addWidget(self.interactive_text)

        self.status_bar = StatusBar(parent=self, running_text='Code is running', not_running_text='awaiting input')
        sub_layout.addWidget(self.status_bar)

        self.main_layout.addWidget(sub_frame)

        self.python_interpreter = PythonTerminal({})

        self.welcome_text = """# Welcome to the OpenDFT scripting console
#
# You can normally(*) use python to script here in addition to OpenDFT objects that can be
# used to calculate electronic properties and visualize them. The structure, engine options etc. are loaded, when the
# scripting window is opened and can be changed within the script. 
#
# Predefined variables:
# 
# engine:           Dft engine object that can be used to calculate electronic properties. 
#
# structure:        crystal or molecular structure from the main application. 
#                   If you defined a structure with the main window you can directly use it here.
#
# plot_structure:   Function that updates the plot in the main window.
#
# plot_scf:         Function that plots the current scf convergence in the main window.
#
# Use the help function to learn more about the variables, 
# e.g. help(engine) and help(engine.start_ground_state) should be quite helpful
#
# (*) For technical reasons matplotlib can be used but has to be put at the end of the script after a special seperator,
# namely [dollarsign]matplotlib (see example below)
#
# Following is an runnable (press f5) example of how to find the optimal cell scale (with very few steps for faster run-time)
#
import numpy as np

atoms = structure.atoms # Save the current atom information for future use
lattice_vectors = structure.lattice_vectors # save the lattice vectors

energies = []
scales = np.linspace(0.8,1.2,5) # scale factor for unit cell
scales_succes = []

for scale in scales:
  structure_temp = CrystalStructure(scale*lattice_vectors,atoms)
  plot_structure(structure_temp) # plot the structure in the main window
  engine.start_ground_state(structure_temp,blocking=True) # start the calculation
  try:
    scf_list = engine.read_scf_status() # Read the full list of scf energies
    plot_scf(scf_list)
    energies.append(scf_list[-1,1]) # append only the last to the energies list
    scales_succes.append(scale)
  except Exception as e:
    print(repr(e))

print('Emin = {0:1.2f}'.format(min(energies)))
optimal_scale = scales_succes[energies.index(min(energies))]
print('Optimal cell scale = {0:1.2f}'.format(optimal_scale))

plot_structure(structure) 
 
## Plot with matplotlib

$matplotlib

import matplotlib.pyplot as plt

plt.figure(1)
plt.clf()
plt.plot(scales_succes,energies,'k',linewidth=2)
plt.show()
        """
        self.input_text_widget.setPlainText(self.welcome_text)

        self.last_history = None

        self.saved_code = self.welcome_text
        self.saved_code_filename = None
        self.interactive_history = []
        self.current_history_element = -1

    def make_menubar(self):
        self.menu_bar = QtGui.QMenuBar(self)
        self.main_layout.addWidget(self.menu_bar)
        self.menu_bar.setNativeMenuBar(False)
        self.file_menu = self.menu_bar.addMenu('&File')

        new_file_action = QtGui.QAction("New", self)
        new_file_action.setShortcut("Ctrl+n")
        new_file_action.triggered.connect(self.new_file)
        self.file_menu.addAction(new_file_action)

        load_file_action = QtGui.QAction("Load", self)
        load_file_action.setShortcut("Ctrl+o")
        load_file_action.triggered.connect(self.load_file)
        self.file_menu.addAction(load_file_action)

        save_file_action = QtGui.QAction("Save", self)
        save_file_action.setShortcut("Ctrl+s")
        save_file_action.triggered.connect(self.save_code)
        self.file_menu.addAction(save_file_action)

        save_as_file_action = QtGui.QAction("Save as", self)
        save_as_file_action.setShortcut("Ctrl+shift+s")
        save_as_file_action.triggered.connect(lambda: self.save_code(ask_filename=True))
        self.file_menu.addAction(save_as_file_action)

        reload_vars_action = QtGui.QAction("Reload variables from main", self)
        reload_vars_action.setShortcut('Ctrl+F5')
        reload_vars_action.triggered.connect(self.parent.open_scripting_console)
        self.file_menu.addAction(reload_vars_action)
        # Todo change this to multiprocessing process and use queue or pipe to move around information

        self.run_menu = self.menu_bar.addMenu('&Run')
        run_code_action = QtGui.QAction("Run code", self)
        run_code_action.setShortcut("F5")
        run_code_action.triggered.connect(self.start_code_execution)
        self.run_menu.addAction(run_code_action)

        run_cell_action = QtGui.QAction("Run cell", self)
        run_cell_action.setToolTip('Runs the selected cell. Cells can be seperated with ##')
        run_cell_action.setShortcut("F7")
        run_cell_action.triggered.connect(self.run_cell)
        self.run_menu.addAction(run_cell_action)

        run_cell_and_jump_action = QtGui.QAction("Run cell and jump", self)
        run_cell_and_jump_action.setToolTip('Runs the selected cell. Cells can be seperated with ##')
        run_cell_and_jump_action.setShortcut("F8")
        run_cell_and_jump_action.triggered.connect(lambda: self.run_cell(jump=True))
        self.run_menu.addAction(run_cell_and_jump_action)

        run_selection_action = QtGui.QAction("Run selection", self)
        run_selection_action.setShortcut("F9")
        run_selection_action.triggered.connect(self.run_selection)
        self.run_menu.addAction(run_selection_action)

        # terminate_execution_action = QtGui.QAction("Terminate execution", self)
        # terminate_execution_action.setShortcut("F12")
        # terminate_execution_action.triggered.connect(self.terminate_execution)
        # self.run_menu.addAction(terminate_execution_action)

    def update_output(self):
        if self.code_thread.is_alive():
            self.status_bar.set_engine_status(True)
        else:
            self.status_bar.set_engine_status(False)

        if self.last_history == self.python_interpreter.out_history:
            return
        else:
            self.last_history = copy.deepcopy(self.python_interpreter.out_history)

        cur_pos = self.output_scrollbar.value()
        if cur_pos == self.output_scrollbar.maximum():
            final_pos = None
        else:
            final_pos = cur_pos

        history = [item.replace(u"\u2029", '\n') for sublist in self.python_interpreter.out_history for item in sublist
                   if len(item) > 0]  # collapse list of list to one list
        out_text = u'\n'.join(history)
        self.output_text_widget.setPlainText(out_text)
        if final_pos is None:
            final_pos = self.output_scrollbar.maximum()
        self.output_scrollbar.setValue(final_pos)
        QtGui.QApplication.processEvents()

    def history_move(self, x):
        self.current_history_element += x
        if self.current_history_element < 0:
            self.current_history_element = -1
        elif self.current_history_element >= len(self.interactive_history):
            self.current_history_element = len(self.interactive_history) - 1
        if self.current_history_element < 0:
            self.interactive_text.set_text('')
        else:
            self.interactive_text.set_text(self.interactive_history[self.current_history_element])

    def run_selection(self):
        cursor = self.input_text_widget.textCursor()
        code_text = cursor.selectedText().replace(u"\u2029", '\n')
        out = self.python_interpreter.run_code(code_text)
        # self.update_output()

    def run_cell(self, jump=False):
        cursor = self.input_text_widget.textCursor()
        code_text = self.input_text_widget.toPlainText()
        cells = re.split('##[^#\n]*\n', code_text)
        y = cursor.blockNumber()
        code_split = code_text.split('\n')
        del_pos = [i for i, x in enumerate(code_split) if '##' in x]
        tr = [i + 1 for i, x in enumerate(del_pos) if x <= y]
        if len(tr) == 0:
            sel_cell = 0
        else:
            sel_cell = max(tr)
        out = self.python_interpreter.run_code(cells[sel_cell])
        # self.update_output()
        if jump:
            index = 0
            regex = QtCore.QRegExp('##')
            for i in range(sel_cell + 1):
                res = regex.indexIn(code_text, index)
                if res < 0:
                    break
                else:
                    index = res + 2
            cursor.setPosition(index)
            self.input_text_widget.setTextCursor(cursor)

    def run_code(self, code_text=None):
        if code_text is None:
            code_text = self.input_text_widget.toPlainText()
        s_time = time.time()
        out = self.python_interpreter.run_code(code_text)
        run_time = time.time() - s_time
        if run_time < 600:
            unit = 's'
            show_time = run_time
        elif run_time < 6000:
            show_time = run_time / 60
            unit = 'min'
        elif run_time < 60 ** 2 * 50:
            show_time = run_time / 60 ** 2
            unit = 'h'
        else:
            show_time = run_time / 60 ** 2 / 24
            unit = 'days'
        self.python_interpreter.out_history.append(
            ['------- Code execution finished after {0:1.1f} {1} -------\n'.format(show_time, unit)])

    def handle_interactive_text(self):
        code_text = self.interactive_text.get_text()
        self.python_interpreter.out_history.append(['>> ' + code_text])
        self.interactive_text.set_text('')
        out = self.python_interpreter.run_code(code_text)
        # self.update_output()
        self.interactive_history.insert(0, code_text)
        self.current_history_element = -1

    def terminate_execution(self):
        if self.code_thread.is_alive():
            # self.code_thread.__stop()
            pass

    def new_file(self):
        if self.check_saved_progress():
            self.input_text_widget.setPlainText('')
        self.saved_code = ''
        self.saved_code_filename = None

    def check_saved_progress(self):
        if self.saved_code != self.input_text_widget.toPlainText():
            msg = "There is unsaved progress. Do you want to save before starting a new file?"
            reply = QtGui.QMessageBox.question(self, 'Unsaved progress',
                                               msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

            if reply == QtGui.QMessageBox.No:
                return False
            elif reply == QtGui.QMessageBox.Yes:
                self.save_code()
        return True

    def load_file(self):
        if not self.check_saved_progress():
            return

        file_name = QtGui.QFileDialog.getOpenFileName(self, 'Load File')[0]
        if not file_name:
            return
        with open(file_name, 'r') as f:
            text = f.read()
        self.input_text_widget.setPlainText(text)
        self.saved_code = text
        self.saved_code_filename = file_name
        self.setWindowTitle(self.custom_window_title + ' - ' + file_name)

    def save_code(self, ask_filename=False):
        if self.saved_code_filename is None or ask_filename:
            file_name = QtGui.QFileDialog.getSaveFileName(self, 'Save File')[0]
        else:
            file_name = self.saved_code_filename

        if not file_name:
            return
        code = self.input_text_widget.toPlainText()
        with open(file_name, 'w') as f:
            f.write(code)
        self.saved_code_filename = file_name
        self.saved_code = code
        self.setWindowTitle(self.custom_window_title + ' - ' + file_name)

    def start_code_execution(self):
        code_text = self.input_text_widget.toPlainText()
        splitted_code = re.split(r'[\s\t\n]+\$matplotlib[\s\t]*\n', code_text, 1)
        main_code = splitted_code[0]
        if len(splitted_code) > 1:
            matplotlib_code = splitted_code[1]
        else:
            matplotlib_code = None

        if self.contains_matplotlib(main_code):
            self.error_widget.showMessage(
                'There seems to be matplotlib code in the main body. Please seperate them by using a #$matplotlib seperator')
            return

        self.start_code_thread(lambda: self.run_code(code_text=main_code))
        if matplotlib_code is not None:
            self.queue.put({'task': 'matplotlib', 'code': matplotlib_code})

    def start_code_thread(self, target):
        if self.code_thread.is_alive():
            return
        self.code_thread = threading.Thread(target=target)
        self.code_thread.start()

    def show(self):
        self.update_fields_timer.start(100)
        super(ConsoleWindow, self).show()

    def closeEvent(self, event):
        self.update_fields_timer.stop()
        event.accept()

    def handle_queue_item(self, item):
        task = item['task']

        if task == 'matplotlib':
            if not self.code_thread.is_alive():
                self.python_interpreter.run_code(item['code'])
            else:
                self.queue.put({'task': 'matplotlib', 'code': item['code']})

    def check_queue_and_update(self):
        if not self.queue.empty():
            q_item = self.queue.get()
            self.handle_queue_item(q_item)
        self.update_output()

    def contains_matplotlib(self, code):
        code_lines = code.splitlines()
        matplotlib_found = False
        for line in code_lines:
            line_wo_comments = line.split('#')
            if 'matplotlib' in line_wo_comments:
                matplotlib_found = True
                break
        return matplotlib_found


class CodeInformationWindow(QtGui.QDialog):
    def __init__(self, parent=None):
        super(CodeInformationWindow, self).__init__(parent)

        self.parent = self
        self.resize(300, 500)
        self.setWindowTitle('Run details')
        layout = QtGui.QHBoxLayout(self)
        self.text_view = QtGui.QTextBrowser(self)
        layout.addWidget(self.text_view)

    def show_information(self, information, name=None):
        if information is None:
            text = 'No details available'
            self.text_view.setPlainText(text)
            return

        try:
            if name is None:
                sub_text = 'selected result'
            else:
                sub_text = name
            text = 'Run details for ' + sub_text + '\n'

            def add_dic_to_text(text, indic, increment=0):
                dic = OrderedDict(sorted(indic.items(), key=lambda t: t[0]))
                for key, value in dic.items():
                    text += ' ' * increment + key + ': ' + value + '\n'
                return text

            if 'scf' in information.keys():
                text += '\nSCF:\n'
                text = add_dic_to_text(text, information['scf'], increment=2)

            if 'bandstructure' in information.keys():
                text += '\nBandstructure:\n'
                k_path = information['bandstructure']['k path']
                text += '  K path:\n'
                for el, label in k_path:
                    text += '    {0:1.3f} {1:1.3f} {2:1.3f}'.format(*el) + '  ' + label + '\n'
                for key, value in information['bandstructure'].items():
                    if key != 'k path':
                        text += key + ': ' + value + '\n'

            rest_dic = {key: value for key, value in information.items() if key not in ('scf', 'bandstructure')}
            for key, value in rest_dic.items():
                text += '\n' + key.title() + ':\n'
                text = add_dic_to_text(text, value, increment=2)
        except Exception as e:
            logging.error(e)
            text = 'An error occured while reading the information'

        self.text_view.setPlainText(text)


class LoadResultsWindow(QtGui.QDialog):
    def __init__(self, parent, tasks):
        super(LoadResultsWindow, self).__init__(parent)
        self.setWindowTitle('Edit Structure')
        self.setFixedSize(300, 100)
        self.parent = parent
        self.tasks = tasks

        main_layout = QtGui.QVBoxLayout(self)

        self.result_name_entry = EntryWithLabel(self, 'name')
        main_layout.addWidget(self.result_name_entry)

        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel | QtGui.QDialogButtonBox.Ok)
        main_layout.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.accept_own)
        self.buttonBox.rejected.connect(self.reject_own)

    def accept_own(self):
        self.parent.load_results_from_engine(self.tasks, title=self.result_name_entry.get_text())
        self.close()

    def reject_own(self):
        self.reject()


class OptionFrame(QtGui.QGroupBox):
    def __init__(self, parent, options, title='', tooltips={}, checkbuttons=[], buttons=[]):
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

        for text, state in checkbuttons:
            cb = QtGui.QCheckBox(text, parent=self)
            if state:
                cb.nextCheckState()
            self.checkbuttons.append(cb)
            self.layout.addWidget(cb, counter // self.widgets_per_line, counter % self.widgets_per_line)
            counter += 1
        for text, function in buttons:
            button = QtGui.QPushButton(text)
            button.clicked.connect(function)
            button.setFixedHeight(30)
            self.buttons.append(button)
            self.layout.addWidget(button, counter // self.widgets_per_line, counter % self.widgets_per_line)
            counter += 1
        self.make_option_entries()

    def make_option_entries(self):
        counter = (len(self.checkbuttons) + len(self.buttons)) // self.widgets_per_line + self.widgets_per_line
        for option_key, option_value in self.options.items():
            entry = EntryWithLabel(self, option_key, option_value)
            if option_key in self.tooltips.keys():
                entry.setToolTip(self.tooltips[option_key].replace('\n', '<br>'))

            self.layout.addWidget(entry, counter // self.widgets_per_line, counter % self.widgets_per_line)
            self.entry_dict[option_key] = entry
            counter += 1

    def read_all_entries(self):
        for key, entry in self.entry_dict.items():
            self.options[key] = entry.get_text()

    def set_all_entries(self):
        for key, value in self.options.items():
            self.entry_dict[key].set_text(value)

    def read_checkbuttons(self):
        res = {}
        for cb in self.checkbuttons:
            res[cb.text()] = bool(cb.checkState())
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
        # self.scroll_area.setFixedHeight(800)
        self.scroll_area.setWidget(mygroupbox)

        self.layout.addWidget(self.scroll_area)

        self.general_option_widget = OptionFrame(self, esc_handler.general_options, title='General options')
        myform.addRow(self.general_option_widget)

        self.scf_option_widget = OptionFrame(self, esc_handler.scf_options, title='Groundstate options',
                                             tooltips=esc_handler.scf_options_tooltip)
        myform.addRow(self.scf_option_widget)

        self.bs_option_widget = OptionFrame(self, esc_handler.bs_options, title='Bandstructure options',
                                            checkbuttons=[['Calculate BS', True],['Calculate DOS',False]],
                                            buttons=[['Choose k-path', self.parent.open_brillouin_window]])
        myform.addRow(self.bs_option_widget)

        self.relax_option_widget = OptionFrame(self, esc_handler.relax_options, title='Structure relaxation options',
                                               tooltips=esc_handler.relax_options_tooltip)
        myform.addRow(self.relax_option_widget)

        self.gw_option_widget = OptionFrame(self, esc_handler.gw_options, title='GW options',
                                            tooltips=esc_handler.gw_options_tooltip)
        myform.addRow(self.gw_option_widget)

        self.phonons_option_widget = OptionFrame(self, esc_handler.phonons_options, title='Phonon options',
                                                 tooltips=esc_handler.phonons_options_tooltip)
        myform.addRow(self.phonons_option_widget)

        self.optical_spectrum_option_widget = OptionFrame(self, esc_handler.optical_spectrum_options,
                                                          title='Excited states options',
                                                          tooltips=esc_handler.optical_spectrum_options_tooltip)
        myform.addRow(self.optical_spectrum_option_widget)

        mygroupbox.setLayout(myform)

        self.button_widget = QtGui.QWidget(self)
        self.button_widget.show()
        self.button_layout = QtGui.QHBoxLayout(self.button_widget)
        self.button_layout.setAlignment(QtCore.Qt.AlignLeft)

        button_size = (150, 50)

        self.start_ground_state_calculation_button = QtGui.QPushButton('Start Ground\nState Calculation',
                                                                       self.button_widget)
        self.start_ground_state_calculation_button.setFixedSize(*button_size)
        self.start_ground_state_calculation_button.clicked.connect(self.start_ground_state_calculation)
        self.button_layout.addWidget(self.start_ground_state_calculation_button)

        self.start_relax_button = QtGui.QPushButton('Start Structure\nRelaxation', self.button_widget)
        self.start_relax_button.setFixedSize(*button_size)
        self.start_relax_button.clicked.connect(self.start_relax)
        self.button_layout.addWidget(self.start_relax_button)

        self.start_gw_button = QtGui.QPushButton('Start GW', self.button_widget)
        self.start_gw_button.setFixedSize(*button_size)
        self.start_gw_button.clicked.connect(self.start_gw)
        self.button_layout.addWidget(self.start_gw_button)

        self.start_phonon_button = QtGui.QPushButton('Start Phonon\nBandstructure', self.button_widget)
        self.start_phonon_button.setFixedSize(*button_size)
        self.start_phonon_button.clicked.connect(self.start_phonons)
        self.button_layout.addWidget(self.start_phonon_button)

        self.start_optical_spectrum_button = QtGui.QPushButton('Calculate optical\nspectrum', self.button_widget)
        self.start_optical_spectrum_button.setFixedSize(*button_size)
        self.start_optical_spectrum_button.clicked.connect(self.start_optical_spectrum_calculation)
        self.button_layout.addWidget(self.start_optical_spectrum_button)

        self.abort_calculation_button = QtGui.QPushButton('Abort Calculation', self.button_widget)
        self.abort_calculation_button.setFixedSize(*button_size)
        self.abort_calculation_button.clicked.connect(self.abort_calculation)
        self.button_layout.addWidget(self.abort_calculation_button)

        self.execute_error_dialog = QtGui.QErrorMessage(parent=self)
        self.execute_error_dialog.resize(500, 200)

        self.layout.addWidget(self.button_widget)

        self.band_structure_points = None
        self.show()

    def update_all(self):
        self.scf_option_widget.set_all_entries()
        self.general_option_widget.set_all_entries()
        self.gw_option_widget.set_all_entries()
        self.bs_option_widget.set_all_entries()
        self.optical_spectrum_option_widget.set_all_entries()
        self.relax_option_widget.set_all_entries()

    def do_select_event(self):
        pass

    def check_engine_for_compatibility(self, tasks_in):
        tasks = [x for x in tasks_in]
        if 'g0w0 bands' in tasks:
            tasks.remove('g0w0 bands')
        struc_type = type(self.parent.crystal_structure)
        if struc_type == sst.MolecularStructure:
            sym_type = 'non-periodic'
        elif struc_type == sst.CrystalStructure:
            sym_type = 'periodic'
        else:
            raise ValueError('bad type: ' + str(struc_type) + ' for structure object')
        if sym_type not in esc_handler.supported_methods:
            raise Exception(sym_type + ' structures are not supported in the selected dft engine')
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
        self.relax_option_widget.read_all_entries()
        self.phonons_option_widget.read_all_entries()

    def prepare_start(self, tasks):
        self.abort_bool = False
        self.check_if_engine_is_running_and_warn_if_so()
        self.check_engine_for_compatibility(tasks)
        self.read_all_option_widgets()
        overwrite_bool = self.check_data_overwrite(tasks)

    def start_ground_state_calculation(self):
        tasks = []
        if esc_handler.will_scf_run():
            tasks.append('scf')

        bs_checkers = self.bs_option_widget.read_checkbuttons()
        if bs_checkers['Calculate BS'] and type(self.parent.crystal_structure) is sst.CrystalStructure:
            bs_points = self.band_structure_points
            tasks.append('bandstructure')
        else:
            bs_points = None
        if bs_checkers['Calculate DOS']:
            tasks.append('dos')
        try:
            self.prepare_start(tasks)
            esc_handler.start_ground_state(self.parent.crystal_structure, band_structure_points=bs_points,dos=bs_checkers['Calculate DOS'])
            QtCore.QTimer.singleShot(1000, lambda: self.parent.check_engine(tasks))
        except Exception as e:
            self.show_stacktrace_in_error_dialog()
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_relax(self):
        tasks = ['relax']
        try:
            self.prepare_start(tasks)
            esc_handler.start_relax(self.parent.crystal_structure)
            QtCore.QTimer.singleShot(1000, lambda: self.parent.check_engine(tasks))
        except Exception as e:
            self.show_stacktrace_in_error_dialog()
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_gw(self):
        tasks = []

        # bs_checkers = self.bs_option_widget.read_checkbuttons()
        # if bs_checkers['Calculate BS']:
        #     tasks.append('bandstructure')
        tasks.extend(['g0w0', 'g0w0 bands'])
        try:
            self.prepare_start(tasks)
            esc_handler.start_gw(self.parent.crystal_structure, self.band_structure_points)
            QtCore.QTimer.singleShot(1000, lambda: self.parent.check_engine(tasks))
        except Exception as e:
            self.show_stacktrace_in_error_dialog()
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_phonons(self):
        tasks = ['phonons']
        try:
            self.prepare_start(tasks)
            esc_handler.start_phonon(self.parent.crystal_structure, self.band_structure_points)
            QtCore.QTimer.singleShot(2000, lambda: self.parent.check_engine(tasks))
        except Exception as e:
            self.show_stacktrace_in_error_dialog()
        else:
            self.parent.status_bar.set_engine_status(True)

    def start_optical_spectrum_calculation(self):
        tasks = ['optical spectrum']
        try:
            self.prepare_start(tasks)
            esc_handler.start_optical_spectrum(self.parent.crystal_structure)
            QtCore.QTimer.singleShot(2000, lambda: self.parent.check_engine(tasks))
        except Exception as e:
            self.show_stacktrace_in_error_dialog()
        else:
            self.parent.status_bar.set_engine_status(True)

    def abort_calculation(self):
        self.abort_bool = True
        esc_handler.kill_engine()

    def configure_buttons(self, disable_all=False):
        button_list = [self.start_ground_state_calculation_button, self.start_gw_button,
                       self.start_optical_spectrum_button, self.start_phonon_button, self.start_relax_button,
                       self.abort_calculation_button]
        if disable_all:
            for button in button_list:
                button.setEnabled(False)
        else:
            for button in button_list:
                button.setEnabled(True)

    def check_data_overwrite(self, tasks):
        title = esc_handler.general_options['title']
        if 'bandstructure' in tasks or 'g0w0' in tasks or 'phonons' in tasks or (
                        'scf' in tasks and type(self.parent.crystal_structure) is sst.MolecularStructure):
            data_dic = self.parent.band_structures
        elif 'optical spectrum' in tasks:
            data_dic = self.parent.optical_spectra
        else:
            data_dic = {}

        if title in data_dic.keys():
            del_msg = title + ' is already saved. Do you want to overwrite it?'
            reply = QtGui.QMessageBox.question(self, 'Sure to overwrite?', del_msg, QtGui.QMessageBox.Yes,
                                               QtGui.QMessageBox.No)

            if reply == QtGui.QMessageBox.No:
                raise Exception('Data already exist')

    def show_stacktrace_in_error_dialog(self):
        stacktrace = get_stacktrace_as_string()
        error_message = 'Could not perform Dft Calculation. Task failed with message:<br><br>' + stacktrace + '<br><br>Try following:<br> 1.Check if the selected dft engine is correctly installed<br>' \
                                                                                                              '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
        self.execute_error_dialog.showMessage(error_message)


class ScfWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        # layout.setContentsMargins(0, 0, 0, 0)
        # layout.setSpacing(0)
        self.scf_widget = ScfVisualization(parent=self)
        layout.addWidget(self.scf_widget)

    def do_select_event(self):
        if self.scf_widget.ax is not None:
            self.scf_widget.figure.tight_layout()


class InfoWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.text_widget = QtGui.QTextBrowser(parent=self)
        self.text_widget.setFontFamily('monospace')
        layout.addWidget(self.text_widget)
        self.vertical_scrollbar = self.text_widget.verticalScrollBar()
        self.last_hash = 0

        self.combobox = QtGui.QComboBox(self)
        layout.addWidget(self.combobox)
        self.combobox.addItem('Output')
        self.combobox.addItem('Input')

        self.combobox.setCurrentIndex(0)
        self.combobox.currentIndexChanged.connect(self.do_select_event)
        self.combobox.setMaximumWidth(150)

    def do_select_event(self):
        if esc_handler.project_directory is None:
            return
        try:
            files = [esc_handler.current_output_file, esc_handler.current_input_file]
            file = files[self.combobox.currentIndex()]
            self.update_text(esc_handler.project_directory + esc_handler.working_dirctory + file)
        except IOError:
            self.text_widget.setHtml('')
            self.last_hash = 0

    def update_text(self, filename):
        cur_pos = self.vertical_scrollbar.value()
        with open(filename, 'r') as f:
            text = f.read()
        if hash(text) == self.last_hash:
            return

        self.last_hash = hash(text)

        # text = text.replace('\n','<br>')
        self.text_widget.setPlainText(text)
        self.vertical_scrollbar.setValue(cur_pos)


class PlotWithTreeview(QtGui.QWidget):
    def __init__(self, Visualizer, data_dictionary, parent=None, selection_mode=None):
        QtGui.QWidget.__init__(self)
        self.parent = parent
        self.data_dictionary = data_dictionary
        self.layout = QtGui.QHBoxLayout(self)
        self.plot_widget = Visualizer(parent=self)
        self.treeview = QtGui.QTreeWidget(parent=self)
        self.treeview.setMaximumWidth(200)
        self.treeview.setHeaderHidden(True)
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed)
        if selection_mode == 'multiple':
            self.treeview.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.treeview.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeview.customContextMenuRequested.connect(self.openMenu)

        self.layout.addWidget(self.plot_widget)
        self.layout.addWidget(self.treeview)

        self.show()

    def delete_selected_item(self):
        indexes = self.treeview.selectedIndexes()
        items = [self.treeview.itemFromIndex(index) for index in indexes]
        bs_names = [item.text(0) for item in items]

        del_msg = "Are you sure you want to delete " + ', '.join(bs_names) + ' permanently?'
        reply = QtGui.QMessageBox.question(self, 'Sure to delete?', del_msg, QtGui.QMessageBox.Yes,
                                           QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            for bs_name in bs_names:
                del self.data_dictionary[bs_name]
            self.update_tree()

    def export_selected_item(self, code=False):
        index = self.treeview.selectedIndexes()[0]
        item = self.treeview.itemFromIndex(index)
        bs_name = item.text(0)
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Select filename')

        if filename:
            if type(filename) in [list, tuple]:
                filename = filename[0]
            self.plot_widget.export(filename, self.data_dictionary[bs_name], code=code)

    def show_info_selected_item(self):
        index = self.treeview.selectedIndexes()[0]
        item = self.treeview.itemFromIndex(index)
        bs_name = item.text(0)
        main.information_window.show()
        main.information_window.show_information(self.data_dictionary[bs_name].engine_information, name=bs_name)

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
        menu.addAction('Details', self.show_info_selected_item)
        menu.addAction('Export', self.export_selected_item)
        menu.addAction('Delete', self.delete_selected_item)

        # menu.addAction('Export with code',lambda: self.export_selected_item(code=True))
        # menu.addAction('Rename', self.rename_selected_item)
        menu.exec_(self.treeview.viewport().mapToGlobal(position))

    def handle_item_changed(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return

        data_list = []
        name_list = []
        for index in indexes:
            item = self.treeview.itemFromIndex(index)
            bs_name = item.text(0)
            data = self.data_dictionary[bs_name]
            data_list.append(data)
            name_list.append(bs_name)
        self.plot_widget.plot(data_list, name_list=name_list)
        main.information_window.show_information(self.data_dictionary[bs_name].engine_information, name=bs_name)

    def add_result_key(self, title):
        item = QtGui.QTreeWidgetItem(self.treeview.invisibleRootItem(), [title])
        return item

    def clear_treeview(self):
        self.treeview.clear()

    def update_tree(self):
        self.treeview.itemSelectionChanged.disconnect()
        self.treeview.clear()
        for key, value in OrderedDict(sorted(self.data_dictionary.items())).items():
            self.add_result_key(key)
        self.treeview.itemSelectionChanged.connect(self.handle_item_changed)

    def do_select_event(self):
        self.update_tree()


class ChooseEngineWindow(QtGui.QDialog):
    def __init__(self, parent, defaults):
        super(ChooseEngineWindow, self).__init__(parent=parent)
        self.setWindowTitle('Please choose a DFT engine')
        self.defaults = defaults

        self.parent = parent
        self.resize(900, 700)
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
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel | QtGui.QDialogButtonBox.Ok)
        main_layout.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.accept_own)
        self.buttonBox.button(QtGui.QDialogButtonBox.Ok).setEnabled(False)
        self.buttonBox.rejected.connect(self.reject_own)

    def add_result_key(self, title):
        item = QtGui.QTreeWidgetItem(self.treeview.invisibleRootItem(), [title])
        return item

    def update_tree(self):
        self.treeview.clear()
        for key, value in self.handlers.items():
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
            install_text = '<p style="color:Green;font-weight:bold">{} installation found in: {}</p>'.format(bs_name,
                                                                                                             self.selected_handler.find_engine_folder())
        else:
            install_text = '<p style="color:Red;font-weight:bold">No {} installation found</p>'.format(bs_name)

        method_descriptions = [self.selected_handler.supported_methods.get_description(x) for x in
                               self.selected_handler.supported_methods]
        text = install_text + 'Supported methods:\n\n - ' + '\n - '.join(
            method_descriptions) + '\n\n' + self.selected_handler.info_text

        text = text.replace('\n', '<br>')
        self.text_widget.setHtml(text)
        self.buttonBox.button(QtGui.QDialogButtonBox.Ok).setEnabled(True)

    def accept_own(self):
        if self.selected_handler is None:
            return

        if not general_handler.is_handler_available(self.selected_handler.engine_name):
            reply = QtGui.QMessageBox.question(self, 'Engine not installed',
                                               "The selected dft engine seems not to be installed. "
                                               "The program will not be able to calculate any electronic properties. "
                                               "You can however still visualize structures. Are you sure to proceed?",
                                               QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
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


class VolumeSlicerWidget(QtGui.QDialog):
    def __init__(self, parent=None):
        super(VolumeSlicerWidget, self).__init__(parent)

        self.resize(900, 600)
        layout = QtGui.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.mayavi_widget = VolumeSlicer()
        # self.mayavi_widget.configure_traits()
        self.ui = self.mayavi_widget.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)


class StatusBar(QtGui.QWidget):
    def __init__(self, parent=None, running_text='Engine is running', not_running_text='Engine inactive'):
        QtGui.QWidget.__init__(self)
        # self.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.running_text = running_text
        self.not_running_text = not_running_text
        # self.setMaximumHeight(20)
        self.layout = QtGui.QHBoxLayout(self)
        self.layout.setAlignment(QtCore.Qt.AlignRight)
        self.status_label = QtGui.QLabel(self.not_running_text)
        self.status_label.setStyleSheet("color: red;font:bold 14px")
        self.cpu_label = QtGui.QLabel('')
        self.cpu_label.setStyleSheet("color: green;font:bold 14px")
        self.memory_label = QtGui.QLabel('')
        self.memory_label.setStyleSheet("color: green;font:bold 14px")

        self.layout.addWidget(self.status_label)
        self.layout.addWidget(self.cpu_label)
        self.layout.addWidget(self.memory_label)
        self.show()

    def set_engine_status(self, status, tasks=None):
        if status:
            if tasks:
                tasks_string = ', '.join(tasks)
                tasks_string2 = ' with tasks: ' + tasks_string
            else:
                tasks_string2 = ''
            if not esc_handler.custom_command_active:
                cpu_usage = psutil.cpu_percent()
                mem_usage = psutil.virtual_memory().percent
                cpu_string = 'CPU: {0:1.1f}%'.format(cpu_usage)
                mem_string = 'Memory: {0:1.1f}%'.format(mem_usage)
            else:
                cpu_string = 'Running on cluster'
                mem_string = ''
                mem_usage = 0.0

            self.status_label.setText(self.running_text + tasks_string2+'.')
            self.status_label.setStyleSheet("color: green;font:bold 14px")
            self.cpu_label.setText(cpu_string)
            self.memory_label.setText(mem_string)

            if mem_usage > 90:
                self.memory_label.setStyleSheet("color: red;font:bold 14px")
            elif mem_usage > 70:
                self.memory_label.setStyleSheet("color: orange;font:bold 14px")
            else:
                self.memory_label.setStyleSheet("color: green;font:bold 14px")
        else:
            self.status_label.setText(self.not_running_text)
            self.status_label.setStyleSheet("color: red;font:bold 14px")
            self.cpu_label.setText('')
            self.memory_label.setText('')


class EngineOptionsDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(EngineOptionsDialog, self).__init__(parent)

        self.parent = parent
        self.command_filename = ''

        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtGui.QDialogButtonBox.Cancel | QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Apply)

        self.buttonBox.accepted.connect(self.accept_own)
        self.buttonBox.rejected.connect(self.reject_own)
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(self.apply)

        self.grid_layout_widget = QtGui.QWidget(self)
        self.grid_layout = QtGui.QGridLayout(self.grid_layout_widget)

        self.custom_command_checkbox = QtGui.QCheckBox('Use custom command', parent=self)
        self.grid_layout.addWidget(self.custom_command_checkbox, 0, 0, 1, 1)

        self.load_custom_command_button = QtGui.QPushButton('Select command file', self)
        self.load_custom_command_button.setFixedWidth(150)
        self.load_custom_command_button.setFixedHeight(30)
        self.load_custom_command_button.clicked.connect(self.load_custom_command)
        self.grid_layout.addWidget(self.load_custom_command_button, 0, 1, 1, 1)

        self.filename_label = QtGui.QLabel(self.grid_layout_widget)
        self.grid_layout.addWidget(self.filename_label, 1, 0, 1, 2)

        self.species_path_entry = EntryWithLabel(self, 'Dft engine path')
        self.grid_layout.addWidget(self.species_path_entry, 2, 0, 1, 2)

        engine_label = QtGui.QLabel(self)
        engine_label.setText('Current engine')
        self.grid_layout.addWidget(engine_label, 3, 0, 1, 1)

        current_engine_label = QtGui.QLabel(self)
        current_engine_label.setText(esc_handler.engine_name)
        self.grid_layout.addWidget(current_engine_label, 3, 1, 1, 1)

        combo_label = QtGui.QLabel(self)
        combo_label.setText('Engine at statup')
        self.grid_layout.addWidget(combo_label, 4, 0, 1, 1)

        self.ask_engine_combobox = QtGui.QComboBox(self)
        self.grid_layout.addWidget(self.ask_engine_combobox, 4, 1, 1, 1)

        self.startup_text = 'Ask at startup'
        self.ask_engine_combobox.addItem(self.startup_text)
        for handler in general_handler.handlers.keys():
            self.ask_engine_combobox.addItem(handler)

        def set_combo_by_text(combo, text):
            index = combo.findText(text, QtCore.Qt.MatchFixedString)
            if index >= 0:
                combo.setCurrentIndex(index)

        if self.parent.defaults['default engine'] is None:
            set_combo_by_text(self.ask_engine_combobox, self.startup_text)
        else:
            set_combo_by_text(self.ask_engine_combobox, self.parent.defaults['default engine'])

        self.verticalLayout = QtGui.QVBoxLayout(self)
        self.verticalLayout.addWidget(self.grid_layout_widget)
        self.verticalLayout.addWidget(self.buttonBox)

    def apply(self):
        self.parent.project_properties['custom command'] = self.filename_label.text()
        self.parent.project_properties['custom command active'] = bool(self.custom_command_checkbox.checkState())
        esc_handler.custom_command = self.filename_label.text()
        esc_handler.custom_command_active = bool(self.custom_command_checkbox.checkState())
        species_path = self.species_path_entry.get_text()
        if len(species_path) > 0 and species_path != esc_handler.dft_installation_folder:
            self.parent.project_properties['custom dft folder'] = species_path
            esc_handler.dft_installation_folder = species_path

        startup_text = self.ask_engine_combobox.currentText()
        if startup_text == self.startup_text:
            self.parent.defaults['default engine'] = None
        else:
            self.parent.defaults['default engine'] = startup_text

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
    def __init__(self, side_panel, data_dictionary, parent=None):
        super(OptionWithTreeview, self).__init__(side_panel, data_dictionary, parent)
        self.add_result_key('None')

    def open_slice_widget(self):
        indexes = self.treeview.selectedIndexes()
        if len(indexes) == 0:
            return
        item = self.treeview.itemFromIndex(indexes[0])
        bs_name = item.text(0)

        plot_options = self.plot_widget.get_options()

        main.volume_slicer_window.mayavi_widget.set_data(self.data_dictionary[bs_name].density, main.crystal_structure)
        main.volume_slicer_window.mayavi_widget.display_scene3d(colormap=plot_options['colormap'])
        main.volume_slicer_window.show()

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
            main.mayavi_widget.visualization.plot_density((self.data_dictionary[bs_name]), **plot_options)
            main.information_window.show_information(self.data_dictionary[bs_name].engine_information, bs_name)
            if main.volume_slicer_window.isVisible():
                self.open_slice_widget()  # this forces a replot and does not do any other harm

    def update_tree(self):
        self.treeview.clear()
        for key, value in OrderedDict(sorted(self.data_dictionary.items())).items():
            self.add_result_key(key)
        self.add_result_key('None')


class SliderWithEntry(QtGui.QWidget):
    def __init__(self, parent=None, label=None, limits=(0, 1), value=None):
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
            self.horizontalLayout.addWidget(self.label, 0, 0)
            counter = 1
        else:
            counter = 0

        limit_range = limits[1] - limits[0]

        self.horizontalSlider = QtGui.QSlider(self)
        self.horizontalSlider.setMinimumWidth(200)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setValue((self.value - limits[0]) / limit_range * 100)
        self.horizontalSlider.valueChanged.connect(self.change_text)

        self.horizontalLayout.addWidget(self.horizontalSlider, 0, counter)

        self.lineEdit = QtGui.QLineEdit(self)
        self.lineEdit.setText("{0:1.1f}".format(self.value))
        self.lineEdit.textChanged.connect(self.change_slider)
        self.horizontalLayout.addWidget(self.lineEdit, 0, 1 + counter)
        self.horizontalLayout.setColumnStretch(0 + counter, 3)
        self.horizontalLayout.setColumnStretch(1 + counter, 1)

    def change_text(self):
        val = self.horizontalSlider.value() * (self.limits[1] - self.limits[0]) / 100 + self.limits[0]
        self.lineEdit.textChanged.disconnect()
        self.lineEdit.setText("{0:1.1f}".format(val))
        self.lineEdit.textChanged.connect(self.change_slider)

        self.value = val

    def change_slider(self):
        try:
            val = float(self.lineEdit.text())
        except:
            return
        self.value = (val - self.limits[0]) / (self.limits[1] - self.limits[0]) * 100
        self.horizontalSlider.valueChanged.disconnect()
        self.horizontalSlider.setValue(self.value)
        self.horizontalSlider.valueChanged.connect(self.change_text)

    def get_value(self):
        try:
            val = float(self.lineEdit.text())
        except:
            val = self.value
        return val


class KsStatePlotOptionWidget(QtGui.QWidget):
    def __init__(self, parent):
        super(KsStatePlotOptionWidget, self).__init__(parent)
        self.parent = parent
        self.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.verticalLayoutWidget = QtGui.QWidget(self)
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)

        self.opacity_slider = SliderWithEntry(self.verticalLayoutWidget, label='Opacity', limits=[0, 1], value=0.5)
        self.verticalLayout.addWidget(self.opacity_slider)

        self.contours_entry = EntryWithLabel(self, 'Contours:', '10')
        self.contours_entry.connect_editFinished(self.parent.handle_item_changed)
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
        self.apply_button.setFixedSize(100, 50)
        self.apply_button.clicked.connect(self.parent.handle_item_changed)
        self.button_layout.addWidget(self.apply_button)

        self.slice_button = QtGui.QPushButton('Slice')
        self.slice_button.setFixedSize(100, 50)
        self.slice_button.clicked.connect(self.parent.open_slice_widget)
        self.button_layout.addWidget(self.slice_button)

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
        out_dic = {'opacity': opacity, 'contours': contours, 'transparent': transparent, 'colormap': colormap}

        return out_dic
        # self.test_label = QtGui.QLabel('asdas')
        # self.verticalLayout.addWidget(self.test_label)

        # self.verticalLayoutWidget.show()


class KsStateWindow(QtGui.QDialog):
    def __init__(self, parent):
        super(KsStateWindow, self).__init__(parent)
        self.calc_queue = queue.Queue()
        self.current_calc_properties = {}
        self.setFixedSize(700, 500)
        self.parent = parent
        self.main_widget = QtGui.QWidget(parent=self)
        self.layout = QtGui.QVBoxLayout(self)
        self.calc_ks_group = QtGui.QGroupBox(parent=self.main_widget)
        self.calc_ks_group.setTitle('Calculation Options')
        self.layout.addWidget(self.calc_ks_group)

        self.sub_layout = QtGui.QGridLayout(self.calc_ks_group)

        self.k_point_entry = EntryWithLabel(self.calc_ks_group, 'k point')
        self.sub_layout.addWidget(self.k_point_entry, 0, 0)

        self.n_band_entry = EntryWithLabel(self.calc_ks_group, 'Band index')
        self.sub_layout.addWidget(self.n_band_entry, 0, 1)

        self.label_entry = EntryWithLabel(self.calc_ks_group, 'Label')
        self.sub_layout.addWidget(self.label_entry, 0, 2)

        button_frame = QtGui.QWidget(self.calc_ks_group)
        self.sub_layout.addWidget(button_frame, 2, 0, 1, 0)

        button_layout = QtGui.QHBoxLayout(button_frame)

        self.calculate_button = QtGui.QPushButton('Calculate KS State', button_frame)
        self.calculate_button.setFixedWidth(150)
        self.calculate_button.setFixedHeight(50)
        self.calculate_button.clicked.connect(self.calculate_ks_state)
        button_layout.addWidget(self.calculate_button)

        self.calculate_density_button = QtGui.QPushButton('Calculate\nelectron density', button_frame)
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

        self.plot_widget = OptionWithTreeview(KsStatePlotOptionWidget, self.parent.ks_densities, parent=self)

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
            n_band = list(range(int(n_band_split[0]), 1 + int(n_band_split[1])))
        else:
            n_band = [int(n_band_str)]

        k_band_str = self.k_point_entry.get_text()
        if ',' in k_band_str:
            tr = k_band_str.split(',')
            k = [int(x) for x in tr]
        elif '-' in k_band_str:
            k_split = k_band_str.split('-')
            k = list(range(int(k_split[0]), 1 + int(k_split[1])))
        else:
            k = [int(k_band_str)]

        if len(k) == 1:
            k = k * len(n_band)
        if len(n_band) == 1:
            n_band = n_band * len(k)

        label = self.label_entry.get_text()
        label_list = [label] * len(k)

        to_do_list = zip(k, n_band, label_list)
        for el in to_do_list:
            self.calc_queue.put(el)
        self.start_ks_calculation()

    def choose_nk(self):
        pass

    def start_ks_calculation(self):
        k, n_band, label = self.calc_queue.get()
        self.current_calc_properties = {'type': 'ks density', 'k': k, 'n_band': n_band, 'label': label}
        esc_handler.calculate_ks_density(self.parent.crystal_structure, [k, n_band])
        QtCore.QTimer.singleShot(20, self.check_engine)

    def check_engine(self):
        tasks = ['ks density']
        if esc_handler.is_engine_running(tasks=tasks):
            self.parent.status_bar.set_engine_status(True, tasks=tasks)
            QtCore.QTimer.singleShot(100, self.check_engine)
        else:
            self.parent.status_bar.set_engine_status(False)
            message, err = esc_handler.engine_process.communicate()
            if type(message) == bytes:
                message, err = message.decode(), err.decode()
            if ('error' in message.lower() or len(err) > 0):
                error_message = 'DFT calculation finished with an error:<br><br>' + message.replace('\n',
                                                                                                    '<br>') + '<br>Error:<br>' + err.replace(
                    '\n', '<br>') \
                                + '<br><br>Try following:<br>1.Check if the selected dft engine is correctly installed<br>' \
                                  '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
                self.parent.error_dialog.showMessage(error_message)

            ks_dens = esc_handler.read_ks_state()
            ks_dens.engine_information = {'scf': copy.deepcopy(self.parent.last_run_information['scf'])}
            label = self.current_calc_properties['label']
            if self.current_calc_properties['type'] == 'ks density':
                n_band = self.current_calc_properties['n_band']
                k = self.current_calc_properties['k']
                key = "{} k{} n{}".format(label, k, n_band)
            elif self.current_calc_properties['type'] == 'density':
                key = label + ' density'

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


class EditStructureWindow(QtGui.QMainWindow):
    def __init__(self, parent):
        super(EditStructureWindow, self).__init__(parent)

        self.main_widget = QtGui.QWidget(self)
        self.setCentralWidget(self.main_widget)

        self.setWindowTitle('Edit Structure')
        self.resize(1300, 800)
        self.setFixedHeight(800)
        self.parent = parent
        self.anything_changed = False

        self.error_dialog = QtGui.QErrorMessage(parent=self)
        self.error_dialog.resize(700, 600)

        self.crystal_structure = None

        self.main_layout = QtGui.QVBoxLayout(self.main_widget)
        # self.make_menubar()

        self.sub_main_widget = QtGui.QWidget(self)
        self.main_layout.addWidget(self.sub_main_widget)
        self.sub_main_layout = QtGui.QHBoxLayout(self.sub_main_widget)

        self.structure_visualization = MayaviQWidget(self.crystal_structure, parent=self)
        self.sub_main_layout.addWidget(self.structure_visualization)

        self.structure_widget = QtGui.QWidget(self)
        self.sub_main_layout.addWidget(self.structure_widget)

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

        self.scale_entry = EntryWithLabel(self, 'Scale')
        self.scale_entry.setFixedHeight(50)
        self.unit_cell_option_layout.addWidget(self.scale_entry)
        self.scale_entry.set_text('1.0')
        self.scale_entry.connect_editFinished(self.handle_change)

        self.periodic_checkbox = QtGui.QCheckBox('Periodic', parent=self)
        self.periodic_checkbox.toggle()
        self.periodic_checkbox.stateChanged.connect(self.handle_change)
        self.unit_cell_option_layout.addWidget(self.periodic_checkbox)

        self.unit_cell_table = QtGui.QTableWidget(self.unit_cell_box)
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
                self.unit_cell_table.setItem(i, j, item)

        self.atom_box = QtGui.QGroupBox(self.structure_widget)
        self.atom_box.setTitle('Atoms')
        self.verticalLayout.addWidget(self.atom_box)

        self.atom_layout = QtGui.QVBoxLayout(self.atom_box)

        atom_options_widget = QtGui.QWidget(self)
        self.atom_layout.addWidget(atom_options_widget)

        atom_options_layout = QtGui.QHBoxLayout(atom_options_widget)
        atom_options_layout.setAlignment(QtCore.Qt.AlignLeft)

        unit_label = QtGui.QLabel(self)
        unit_label.setText('Unit')
        atom_options_layout.addWidget(unit_label)

        self.unit_combobox = QtGui.QComboBox(self)
        atom_options_layout.addWidget(self.unit_combobox)
        self.unit_combobox.addItem('Crystal')
        self.unit_combobox.addItem('Cartesian (Bohr)')

        self.unit_combobox.setCurrentIndex(0)
        self.unit_combobox.currentIndexChanged.connect(self.switch_units)
        self.unit_combobox.setMaximumWidth(150)

        self.atom_table = QtGui.QTableWidget(self.atom_box)
        copy_action_atoms = CopySelectedCellsAction(self.atom_table)
        self.atom_table.addAction(copy_action_atoms)
        paste_action = PasteIntoTable(self.atom_table, self)
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

        self.add_atom_button = QtGui.QPushButton('Add atom', self)
        self.add_atom_button.setFixedWidth(150)
        self.add_atom_button.clicked.connect(self.add_atom)
        self.atom_table_buttons_layout.addWidget(self.add_atom_button)

        self.remove_atom_button = QtGui.QPushButton('Remove atoms', self)
        self.remove_atom_button.setFixedWidth(150)
        self.remove_atom_button.clicked.connect(self.remove_atoms)
        self.atom_table_buttons_layout.addWidget(self.remove_atom_button)

        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtGui.QDialogButtonBox.Cancel | QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Apply)
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

        self.info_widget = QtGui.QWidget(self)
        self.info_sub_layout = QtGui.QHBoxLayout(self.info_widget)
        self.info_box = QtGui.QGroupBox(self.info_widget)
        self.info_box.setFixedWidth(240)
        self.info_box.setTitle('Structure Properties')
        self.info_sub_layout.addWidget(self.info_box)

        self.info_layout = QtGui.QVBoxLayout(self.info_box)
        self.info_layout.setAlignment(QtCore.Qt.AlignTop)
        self.sub_main_layout.addWidget(self.info_widget)

        self.info_items = OrderedDict([('volume',{'format':'{0:1.1f} Angstrom<sup>3</sup>','widget':None}),('density',{'format':'{0:1.2f} g/cm<sup>3</sup>','widget':None}),
                                       ('crystal system', {'format': '{0}', 'widget': None}),
                                       ('space group',{'format':'{0}','widget':None}),('point group',{'format':'{0}','widget':None}),
                                       ])

        for name,data_dic in self.info_items.items():
            label = LabeledLabel(name.title()+':',width_label=90)
            self.info_layout.addWidget(label)
            data_dic['widget'] = label

    def apply(self):
        if self.anything_changed:
            crystal_structure = self.read_tables()
            if crystal_structure is not None:
                main.crystal_structure = crystal_structure
                self.update_plot(mode='main')
            else:
                raise ValueError('Bad structure')

    def accept(self):
        try:
            self.apply()
        except ValueError:
            msg = "The entered structure is not valid. Possible cause: Invalid unit cell. Do you want to close this window anyway? All entries will be lost. "
            reply = QtGui.QMessageBox.question(self, 'Invalid structure',
                                               msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

            if reply == QtGui.QMessageBox.Yes:
                self.close()
            else:
                return
        else:
            self.close()

    def reject(self):
        if self.anything_changed:
            main.mayavi_widget.update_crystal_structure(main.crystal_structure)
            main.mayavi_widget.update_plot()
        self.close()

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

    def make_menubar(self):
        self.menu_bar = QtGui.QMenuBar(self)
        self.main_layout.addWidget(self.menu_bar)
        self.menu_bar.setNativeMenuBar(False)
        self.transformation_menu = self.menu_bar.addMenu('&Transformations')

        spatial_transform_action = QtGui.QAction("Spatial", self)
        spatial_transform_action.triggered.connect(self.open_spatial_transform_window)
        self.transformation_menu.addAction(spatial_transform_action)

    def set_structure(self, structure):
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
        self.atom_table.setRowCount(n_rows + 1)
        for j in range(4):
            item = QtGui.QTableWidgetItem()
            self.atom_table.setItem(n_rows, j, item)
        self.connect_tables()

    def remove_atoms(self, *args, **kwargs):
        atoms = kwargs.pop('atoms', None)
        self.disconnect_tables()
        if atoms is None:
            atoms = sorted(set(index.row() for index in self.atom_table.selectedIndexes()))
        for atom in atoms[::-1]:
            self.atom_table.removeRow(atom)
        self.connect_tables()
        self.handle_change()

    def update_fields(self):
        self.disconnect_tables()
        try:
            if self.crystal_structure is None:
                self.clear_atom_table()
                self.clear_unit_cell_table()
                self.set_number_of_atoms(6)
            else:
                if type(self.crystal_structure) is sst.CrystalStructure:
                    scale = self.crystal_structure.scale
                    self.scale_entry.set_text('{0:1.6f}'.format(scale))
                    unit_cell = self.crystal_structure.lattice_vectors
                    for i in range(3):
                        for j in range(3):
                            self.unit_cell_table.item(i, j).setText(find_fraction(unit_cell[i, j] / scale))
                    if not self.periodic_checkbox.checkState():
                        self.periodic_checkbox.stateChanged.disconnect()
                        self.periodic_checkbox.toggle()
                        self.periodic_checkbox.stateChanged.connect(self.handle_change)
                elif type(self.crystal_structure) is sst.MolecularStructure:
                    if self.periodic_checkbox.checkState():
                        self.periodic_checkbox.stateChanged.disconnect()
                        self.periodic_checkbox.toggle()
                        self.periodic_checkbox.stateChanged.connect(self.handle_change)

                if self.unit_combobox.currentIndex() == 0:
                    atoms = self.crystal_structure.atoms
                elif self.unit_combobox.currentIndex() == 1:
                    atoms = self.crystal_structure.calc_absolute_coordinates()

                n_atoms = atoms.shape[0]
                self.set_number_of_atoms(n_atoms)

                for i, atom in enumerate(atoms):
                    coords = atom[0:3]
                    for j, coord in enumerate(coords):
                        item = self.atom_table.item(i, j + 1)
                        item.setText(find_fraction(coord))
                    item = self.atom_table.item(i, 0)
                    item.setText(p_table[atom[3]])
        except Exception as e:
            print(e)
            logging.exception(e)

        self.connect_tables()

    def set_number_of_atoms(self, N):
        self.atom_table.setRowCount(N)
        for i in range(N):
            for j in range(4):
                item = QtGui.QTableWidgetItem()
                self.atom_table.setItem(i, j, item)

    def read_tables(self):
        unit_cell, scale = self.read_unit_cell_table()
        atoms = self.read_atom_table()

        try:
            if self.periodic_checkbox.checkState():
                out_struc = sst.CrystalStructure(unit_cell, atoms, scale=scale,
                                                 relative_coords=not bool(self.unit_combobox.currentIndex()))
            else:
                out_struc = sst.MolecularStructure(atoms)
        except Exception:
            return None
        else:
            return out_struc

    def handle_change(self):
        self.anything_changed = True

        if not self.periodic_checkbox.checkState():
            if self.unit_combobox.currentIndex() == 0:
                self.unit_combobox.currentIndexChanged.disconnect()
                self.unit_combobox.setCurrentIndex(1)
                self.unit_combobox.currentIndexChanged.connect(self.switch_units)
                self.switch_units(convert_periodic=True, periodic=True)

        self.switch_enabled_fields()
        self.update_plot()
        self.update_info()

    def update_plot(self, mode='local'):
        crystal_structure = self.read_tables()
        if crystal_structure is not None:
            if mode == 'local':
                self.structure_visualization.update_crystal_structure(crystal_structure)
                self.structure_visualization.update_plot()
            elif mode == 'main':
                main.mayavi_widget.update_crystal_structure(crystal_structure)
                main.mayavi_widget.update_plot()

    def switch_units(self, event=None, periodic=None, unit=None, convert_periodic=False):
        if periodic is None:
            periodic = self.periodic_checkbox.checkState()
        if unit is None:
            unit = self.unit_combobox.currentIndex()  # 0 for crystal 1 for cartesian

        unit_cell, scale = self.read_unit_cell_table(periodic=periodic)
        try:
            atoms = self.read_atom_table()
        except Exception:
            atoms = None

        if atoms is None or len(atoms) == 0:
            return

        if unit_cell is None:
            msg = 'Could not convert units. Unit cell is not valid. Do you want to proceed without conversion?'
            reply = QtGui.QMessageBox.question(self, 'Invalid unit cell',
                                               msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

            if reply == QtGui.QMessageBox.Yes:
                if not self.periodic_checkbox.checkState():
                    self.update_plot()
                return
            elif reply == QtGui.QMessageBox.No:
                if convert_periodic:
                    self.periodic_checkbox.stateChanged.disconnect()
                    self.periodic_checkbox.toggle()
                    self.periodic_checkbox.stateChanged.connect(self.handle_change)
                else:
                    self.unit_combobox.currentIndexChanged.disconnect()
                    self.unit_combobox.setCurrentIndex(not self.unit_combobox.currentIndex())
                    self.unit_combobox.currentIndexChanged.connect(self.switch_units)


        elif convert_periodic:
            if not periodic:
                raise ValueError('If convert_periodic is true, periodic must be True as well')
            out_struc = sst.CrystalStructure(unit_cell, atoms, scale=scale, relative_coords=unit)
            atoms = out_struc.calc_absolute_coordinates()
            conv_struc = sst.MolecularStructure(atoms)
            self.crystal_structure = conv_struc
            self.update_fields()
        else:
            out_struc = sst.CrystalStructure(unit_cell, atoms, scale=scale, relative_coords=unit)
            self.crystal_structure = out_struc
            self.update_fields()

    def read_atom_table(self):
        n_rows = self.atom_table.rowCount()
        atoms = np.zeros((n_rows, 4))
        for i in range(n_rows):
            try:
                a_type = self.atom_table.item(i, 0).text()
                a_type = int(a_type)
                a_type_is_number = True
            except ValueError:
                a_type_is_number = False
            except Exception:
                continue
            if a_type not in p_table_rev.keys():
                continue
            coord = np.zeros((1, 3))
            skip_bool = False
            for j in range(1, 4):
                try:
                    coord[0, j - 1] = float(eval_expr(self.atom_table.item(i, j).text()))
                except:
                    skip_bool = True
                    break
            if skip_bool:
                continue
            atoms[i, :3] = coord
            if not a_type_is_number:
                a_type = p_table_rev[a_type]
            atoms[i, 3] = a_type

        atoms_clean = atoms[atoms[:, 3] != 0, :]
        return atoms_clean

    def read_unit_cell_table(self, periodic=None):
        if periodic is None:
            periodic = self.periodic_checkbox.checkState()

        try:
            scale_string = self.scale_entry.get_text()
            scale = float(scale_string)
            if scale == 0:
                raise ValueError('Scale cannot be zero')
        except Exception as e:
            logging.exception(e)
            scale = 1.0
        if periodic:
            try:
                unit_cell = np.zeros((3, 3))
                for i in range(3):
                    for j in range(3):
                        item = self.unit_cell_table.item(i, j)
                        unit_cell[i, j] = float(eval_expr(item.text()))

                unit_cell = unit_cell * scale
            except Exception:
                return None, scale
        else:
            unit_cell = None
        return unit_cell, scale

    def switch_enabled_fields(self):
        bool_list = [True, True, True]
        if self.periodic_checkbox.checkState():
            bool_list_out = bool_list
        else:
            bool_list_out = [not (x) for x in bool_list]

        self.scale_entry.setEnabled(bool_list_out[0])
        self.unit_cell_table.setEnabled(bool_list_out[1])
        self.unit_combobox.setEnabled(bool_list_out[2])

    def open_spatial_transform_window(self):
        pass

    def update_info(self):
        crystal_structure = self.read_tables()

        if crystal_structure is None:
            self.clear_info()
            return

        data = {}
        if isinstance(crystal_structure,sst.CrystalStructure):
            from solid_state_tools import bohr
            lattice_vectors = crystal_structure.lattice_vectors
            volume = np.dot(np.cross(lattice_vectors[0, :], lattice_vectors[1, :]), lattice_vectors[2, :])
            data['volume'] = volume
            density = crystal_structure.density(unit='g/cm^3')
            data['density'] = density
            data.update(crystal_structure.lattice_information())

        elif isinstance(crystal_structure,sst.MolecularStructure):
            data.update(crystal_structure.symmetry_information())

        else:
            raise ValueError('Structure most be either CrystalStructure of MolecularStructure type but got {}'.format(type(crystal_structure))+' instead')

        for key, info_dic in self.info_items.items():
            if not key in data.keys():
                info_dic['widget'].text('')
            else:
                info_dic['widget'].text(info_dic['format'].format(data[key]))

    def clear_info(self):
        for name,dic in self.info_items.items():
            dic['widget'].text('')


class CentralWindow(QtGui.QWidget):
    def __init__(self, parent=None, *args, **kwargs):
        super(CentralWindow, self).__init__(*args, **kwargs)
        self.project_loaded = False
        self.project_directory = None
        self.parent = parent
        self._crystal_structure = None
        self.band_structures = {}
        self.dos = {}
        self.optical_spectra = {}
        self.ks_densities = {}
        self.project_properties = {'title': '', 'dft engine': '', 'custom command': '', 'custom command active': False,
                                   'custom dft folder': ''}
        self.esc_handler_options = {}
        self.last_run_information = {'scf': {}, 'bandstructure': {}, 'gw': {}, 'optical spectrum': {}, 'relax': {},
                                     'phonon': {}}

        self.temp_folder = os.path.expanduser("~") + "/.OpenDFT"

        self.queue = queue.Queue()

        self.load_defaults()

        if self.defaults['default engine'] is None:
            choose_engine_window = ChooseEngineWindow(self, self.defaults)
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
        self.band_structure_window = PlotWithTreeview(Visualizer=BandStructureVisualization,
                                                      data_dictionary=self.band_structures, parent=self)

        self.dos_window = PlotWithTreeview(Visualizer=DosVisualization,
                                                      data_dictionary=self.dos, parent=self,selection_mode='multiple')

        self.optical_spectra_window = PlotWithTreeview(Visualizer=OpticalSpectrumVisualization,
                                                       data_dictionary=self.optical_spectra, parent=self,
                                                       selection_mode='multiple')
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
        self.console_window = ConsoleWindow(self)
        self.information_window = CodeInformationWindow(self)
        self.volume_slicer_window = VolumeSlicerWidget(self)
        self.phonon_window = MayaviPhononWindow(self.crystal_structure,self.band_structures,parent=self)


        self.tab_layout = QtGui.QVBoxLayout()
        self.tabWidget.setLayout(self.tab_layout)

        self.list_of_tabs = [self.mayavi_widget, self.dft_engine_window, self.band_structure_window,self.dos_window,
                             self.optical_spectra_window, self.scf_window, self.info_window]

        self.tabWidget.addTab(self.list_of_tabs[0], 'Structure')
        self.tabWidget.addTab(self.list_of_tabs[1], 'DFT-Engine')
        self.tabWidget.addTab(self.list_of_tabs[2], 'Band Structure')
        self.tabWidget.addTab(self.list_of_tabs[3], 'DoS')
        self.tabWidget.addTab(self.list_of_tabs[4], 'Optical Spectrum')
        self.tabWidget.addTab(self.list_of_tabs[5], 'Scf')
        self.tabWidget.addTab(self.list_of_tabs[6], 'Info')

        self.tabWidget.show()

        self.show()
        self.window = MainWindow(self)
        self.window.setWindowTitle("OpenDFT")
        self.window.setWindowIcon(QtGui.QIcon('icon.ico'))
        self.window.setMinimumSize(700, 700)
        self.window.resize(1300, 1100)
        self.window.setCentralWidget(self)
        self.make_menu_bar()

        self.window.show()

        self.configure_buttons(disable_all=True)

        self.queue_timer = QtCore.QTimer()
        self.queue_timer.timeout.connect(self.update_program)
        self.queue_timer.start(200)

        self.save_timer = QtCore.QTimer()
        self.save_timer.timeout.connect(self.save_results)
        self.save_timer.start(300000)

        if DEBUG:
            if sys.platform in ['linux', 'linux2']:
                project_directory = r"/home/jannick/OpenDFT_projects/CH4_exciting/"
                # project_directory = r"/home/jannick/exciting_cluster/GaN"
            else:
                project_directory = r'D:\OpenDFT_projects\test'
            # self.load_saved_results()
            QtCore.QTimer.singleShot(500, lambda: self.load_project(folder_name=project_directory))
            # QtCore.QTimer.singleShot(510, self.open_scripting_console)

    @property
    def crystal_structure(self):
        return self._crystal_structure

    @crystal_structure.setter
    def crystal_structure(self, value):
        self._crystal_structure = value

        if self.dft_engine_window.band_structure_points is None and isinstance(value,sst.CrystalStructure):
            try:
                self.dft_engine_window.band_structure_points = sst.calculate_standard_path(value)
            except Exception as e:
                self.dft_engine_window.band_structure_points = sst.get_emergency_path()
                logging.error('Could not calculate standard path')
                self.error_dialog.showMessage(str(e)+'<br>A standard non structure specific path will be loaded. Use with care!')

            self.brillouin_window.set_path(self.dft_engine_window.band_structure_points)

    def tab_is_changed(self, i):
        self.list_of_tabs[i].do_select_event()

    def make_new_project(self):
        folder_name = QtGui.QFileDialog().getExistingDirectory(parent=self,options=QtGui.QFileDialog.ShowDirsOnly)

        if folder_name != self.project_directory:
            try:
                self.check_and_set_lock(folder_name)
            except Exception:
                return

        if os.path.isfile(folder_name+'/save.pkl'):
            msg = "This folder contains an existing OpenDFT project. If you proceed this project will be overwritten. Proceed?"
            reply = QtGui.QMessageBox.question(self, 'Overwrite project?',msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

            if reply == QtGui.QMessageBox.No:
                return

        if len(folder_name) > 1:
            if self.project_loaded:
                self.save_results()
                self.reset_results_and_plots()
            self.project_directory = folder_name
            esc_handler.project_directory = self.project_directory
            self.initialize_project()
            self.configure_buttons()

    def initialize_project(self):
        self.project_properties.update(
            {'title': '', 'dft engine': '', 'custom command': '', 'custom command active': False,
             'custom dft folder': ''})
        self.window.setWindowTitle("OpenDFT - " + self.project_directory)
        os.chdir(self.project_directory)
        if (esc_handler.pseudo_directory is not None) and (
                not os.path.isdir(self.project_directory + esc_handler.pseudo_directory)):
            os.mkdir(self.project_directory + esc_handler.pseudo_directory)
        self.project_loaded = True

    def reset_results_and_plots(self):
        self.release_lock()
        self.crystal_structure = None
        self.project_directory = None
        self.project_loaded = False
        esc_handler.reset_to_defaults()
        self.project_properties.clear()
        self.dft_engine_window.band_structure_points = None
        self.band_structures.clear()
        self.dos.clear()
        self.optical_spectra.clear()
        self.ks_densities.clear()
        self.mayavi_widget.visualization.clear_plot()
        self.band_structure_window.plot_widget.clear_plot()
        self.band_structure_window.clear_treeview()
        self.dos_window.plot_widget.clear_plot()
        self.scf_window.scf_widget.clear_plot()
        self.optical_spectra_window.plot_widget.clear_plot()
        self.optical_spectra_window.clear_treeview()
        self.dft_engine_window.update_all()
        self.window.setWindowTitle("OpenDFT")

    def load_project(self, *args, **kwargs):

        folder_name = kwargs.pop('folder_name', None)
        if folder_name is None:
            if self.project_directory:
                p = str(Path(self.project_directory).parents[0])
            else:
                p = os.path.expanduser("~")

            folder_name = QtGui.QFileDialog().getExistingDirectory(parent=self,directory=p,options=QtGui.QFileDialog.ShowDirsOnly)

        try:
            self.check_and_set_lock(folder_name)
        except Exception:
            return

        save_filename = folder_name + '/save.pkl'
        if folder_name and os.path.isfile(save_filename):

            if self.project_loaded:
                self.save_results()
                self.reset_results_and_plots()

            self.project_directory = folder_name
            os.chdir(self.project_directory)
            esc_handler.project_directory = self.project_directory
            try:
                self.load_saved_results()
            except Exception as e:
                self.reset_results_and_plots()
                stack_trace = get_stacktrace_as_string()
                err_msg = 'Project seems to be corrupt and could not be loaded. Loading failed with error message:<br>' + stack_trace
                self.error_dialog.showMessage(err_msg)
            else:
                self.dft_engine_window.update_all()
                self.tab_is_changed(self.tabWidget.currentIndex())
                self.window.setWindowTitle("OpenDFT - " + self.project_directory)
                self.project_loaded = True
                self.configure_buttons()
        elif not folder_name:
            return
        else:
            error_message = 'Could not load project. ' + str(folder_name) + ' is not a valid OpenDFT project'
            self.error_dialog.showMessage(error_message)

    def save_results(self):
        # Todo move engine folder and custom command into engine specific dic
        if not self.project_loaded:
            return

        try:
            self.dft_engine_window.read_all_option_widgets()

            option_dic_specific_handler = {'scf_options': esc_handler.scf_options,
                                           'general options': esc_handler.general_options,
                                           'bs options': esc_handler.bs_options,
                                           'phonon options': esc_handler.phonons_options,
                                           'optical spectrum options': esc_handler.optical_spectrum_options,
                                           'gw options': esc_handler.gw_options,
                                           'relax options': esc_handler.relax_options,
                                           'last run information': self.last_run_information}
            self.esc_handler_options[esc_handler.engine_name] = option_dic_specific_handler
            a = {'crystal structure': self.crystal_structure, 'band structure': self.band_structures,'dos':self.dos,
                 'optical spectra': self.optical_spectra, 'esc handler options': self.esc_handler_options,
                 'properties': self.project_properties, 'dft engine': esc_handler.engine_name,
                 'ks densities': self.ks_densities, 'k path': self.dft_engine_window.band_structure_points}
            with open(self.project_directory + '/save.pkl', 'wb') as handle:
                pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as e:
            logging.exception(e)

    def load_saved_results(self):
        esc_handler.project_directory = self.project_directory
        try:
            with open(self.project_directory + '/save.pkl', 'rb') as handle:
                b = pickle.load(handle, encoding='latin1')

                k_path = b.pop('k path', None)
                if k_path is not None:
                    self.dft_engine_window.band_structure_points = k_path

                self.crystal_structure = b.pop('crystal structure', None)
                if self.crystal_structure is not None:
                    self.mayavi_widget.update_crystal_structure(self.crystal_structure)
                    self.mayavi_widget.update_plot()

                loaded_bandstructure_dict = b.pop('band structure', None)
                if type(loaded_bandstructure_dict) == dict:
                    for key, value in loaded_bandstructure_dict.items():
                        self.band_structures[key] = value

                loaded_dos_dict = b.pop('dos',None)
                if type(loaded_dos_dict) == dict:
                    for key, value in loaded_dos_dict.items():
                        self.dos[key] = value


                loaded_optical_spectra_dict = b.pop('optical spectra', None)
                if type(loaded_optical_spectra_dict) == dict:
                    for key, value in loaded_optical_spectra_dict.items():
                        self.optical_spectra[key] = value

                loaded_ksdens_dict = b.pop('ks densities', None)
                if type(loaded_ksdens_dict) == dict:
                    for key, value in loaded_ksdens_dict.items():
                        self.ks_densities[key] = value

                self.esc_handler_options = b.pop('esc handler options', {})

                option_dic_specific_handler = self.esc_handler_options.pop(esc_handler.engine_name, None)

                if option_dic_specific_handler is not None:

                    last_run_information = option_dic_specific_handler.pop('last run information', None)
                    if last_run_information is not None:
                        self.last_run_information = last_run_information

                    def set_esc_handler_dic_to_loaded_dic(esc_dic, loaded_dic):
                        if loaded_dic is not None:
                            for key, value in loaded_dic.items():
                                if key in esc_dic.keys():
                                    esc_dic[key] = value

                    load_scf_options = option_dic_specific_handler.pop('scf_options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.scf_options, load_scf_options)

                    load_general_options = option_dic_specific_handler.pop('general options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.general_options, load_general_options)

                    load_bs_options = option_dic_specific_handler.pop('bs options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.bs_options, load_bs_options)

                    load_phonon_options = option_dic_specific_handler.pop('phonon options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.phonons_options, load_phonon_options)

                    load_gw_options = option_dic_specific_handler.pop('gw options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.gw_options, load_gw_options)

                    load_optical_spectrum_options = option_dic_specific_handler.pop('optical spectrum options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.optical_spectrum_options,
                                                      load_optical_spectrum_options)

                    load_relax_options = option_dic_specific_handler.pop('relax options', None)
                    set_esc_handler_dic_to_loaded_dic(esc_handler.relax_options, load_relax_options)

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
            with open(self.temp_folder + '/defaults.pkl', 'rb') as handle:
                b = pickle.load(handle)
        except Exception as e:
            logging.info('Default file not found or invalid')
            b = {}

        default_engine = b.pop('default engine', None)

        self.defaults['default engine'] = default_engine

    def save_defaults(self):
        with open(self.temp_folder + '/defaults.pkl', 'wb') as handle:
            pickle.dump(self.defaults, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def update_structure_plot(self):
        self.mayavi_widget.update_crystal_structure(self.crystal_structure)
        self.mayavi_widget.update_plot()

        # t = MyQThread(self.mayavi_widget.update_plot)
        # t.start()

    def load_crystal_structure(self, filetype):
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

            if filetype in ['exciting', 'quantum espresso']:
                self.crystal_structure = general_handler.parse_input_file(filetype, file_name)
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

        self.save_project_action = QtGui.QAction("Save project", self.window)
        self.save_project_action.setShortcut("Ctrl+s")
        self.save_project_action.setStatusTip('Save project')
        self.save_project_action.triggered.connect(self.save_results)
        self.file_menu.addAction(self.save_project_action)

        self.file_menu.addSeparator()

        self.new_structure_action = QtGui.QAction("New structure", self.window)
        self.new_structure_action.setShortcut('Ctrl+Shift+n')
        self.new_structure_action.setStatusTip('Make new structure by hand')
        self.new_structure_action.triggered.connect(lambda: self.open_structure_window(new=True))
        self.file_menu.addAction(self.new_structure_action)

        self.edit_structure_action = QtGui.QAction("Edit structure", self.window)
        self.edit_structure_action.setShortcut('Ctrl+Shift+e')
        self.edit_structure_action.setStatusTip('Edit existing structure by hand')
        self.edit_structure_action.triggered.connect(lambda: self.open_structure_window(new=False))
        self.file_menu.addAction(self.edit_structure_action)

        self.import_structure_menu = self.file_menu.addMenu('Import structure from')

        open_structure_action_exciting = QtGui.QAction("exciting input file", self.window)
        open_structure_action_exciting.setStatusTip('Load crystal structure from exciting xml')
        open_structure_action_exciting.triggered.connect(lambda: self.load_crystal_structure('exciting'))
        self.import_structure_menu.addAction(open_structure_action_exciting)

        open_structure_action_quantum_espresso = QtGui.QAction("quantum_espresso input file", self.window)
        open_structure_action_quantum_espresso.setStatusTip('Load crystal structure from quantum espresso input file')
        open_structure_action_quantum_espresso.triggered.connect(
            lambda: self.load_crystal_structure('quantum espresso'))
        self.import_structure_menu.addAction(open_structure_action_quantum_espresso)

        open_structure_action_cif = QtGui.QAction("cif", self.window)
        open_structure_action_cif.setShortcut('Ctrl+Shift+c')
        open_structure_action_cif.setStatusTip('Load crystal structure from cif file')
        open_structure_action_cif.triggered.connect(lambda: self.load_crystal_structure('cif'))
        self.import_structure_menu.addAction(open_structure_action_cif)

        self.file_menu.addSeparator()

        self.import_results_menu = self.file_menu.addMenu('Import results')

        import_result_action_bandstructure = QtGui.QAction("bandstructure", self.window)
        import_result_action_bandstructure.triggered.connect(lambda: self.open_load_result_window(['bandstructure']))
        self.import_results_menu.addAction(import_result_action_bandstructure)

        import_result_action_gw_bandstructure = QtGui.QAction("gw bandstructure", self.window)
        import_result_action_gw_bandstructure.triggered.connect(lambda: self.open_load_result_window(['g0w0']))
        self.import_results_menu.addAction(import_result_action_gw_bandstructure)

        import_result_action_optical_spectrum = QtGui.QAction("optical spectrum", self.window)
        import_result_action_optical_spectrum.triggered.connect(
            lambda: self.open_load_result_window(['optical spectrum']))
        self.import_results_menu.addAction(import_result_action_optical_spectrum)

        import_result_action_phonon_bandstructure = QtGui.QAction("phonon bandstructure", self.window)
        import_result_action_phonon_bandstructure.triggered.connect(lambda: self.open_load_result_window(['phonons']))
        self.import_results_menu.addAction(import_result_action_phonon_bandstructure)

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

        phonon_action = QtGui.QAction("Phonons", self.window)
        phonon_action.triggered.connect(self.open_phonon_window)
        self.vis_menu.addAction(phonon_action)

        self.dft_menu = self.menu_bar.addMenu('&DFT Engine')

        dft_options_action = QtGui.QAction("Options", self.window)
        dft_options_action.setStatusTip('Options for dft engine')
        dft_options_action.triggered.connect(self.open_engine_option_window)
        self.dft_menu.addAction(dft_options_action)

        self.scripting_menu = self.menu_bar.addMenu('&Scripting')

        open_console_action = QtGui.QAction("Open console", self.window)
        open_console_action.triggered.connect(self.open_scripting_console)
        self.scripting_menu.addAction(open_console_action)

    def close_application(self):
        if not DEBUG:
            reply = QtGui.QMessageBox.question(self, 'Message',
                                               "Are you sure to quit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if DEBUG or reply == QtGui.QMessageBox.Yes:
            self.save_defaults()
            if self.project_loaded:
                self.save_results()
                self.release_lock()
            self.parent.quit()
            logging.info('Program stopped normally')
            sys.exit()

    def check_relax(self):
        new_struc = esc_handler.load_relax_structure()
        if new_struc is not None:
            self.crystal_structure = new_struc
            self.update_structure_plot()

    def check_engine(self, tasks):

        tasks = [x.lower() for x in tasks]
        if self.dft_engine_window.abort_bool:
            self.status_bar.set_engine_status(False)
            self.dft_engine_window.abort_bool = False
            return
        elif esc_handler.is_engine_running(tasks=tasks):
            selected_tab_index = self.tabWidget.currentIndex()
            if self.list_of_tabs[selected_tab_index] == self.scf_window:
                self.scf_data = esc_handler.read_scf_status()
                if self.scf_data is not None:
                    self.scf_window.scf_widget.plot(self.scf_data)
            elif self.list_of_tabs[selected_tab_index] == self.info_window:
                self.info_window.do_select_event()
            QtCore.QTimer.singleShot(500, lambda: self.check_engine(tasks))
            self.status_bar.set_engine_status(True, tasks=tasks)
            if 'relax' in tasks:
                self.check_relax()
        else:
            self.scf_data = esc_handler.read_scf_status()
            if self.scf_data is not None:
                self.scf_window.scf_widget.plot(self.scf_data)
            self.status_bar.set_engine_status(False)
            message, err = esc_handler.engine_process.communicate()
            if type(message) == bytes:
                message, err = message.decode(), err.decode()

            actual_error_bool = len(err)>0 and not err.strip() == 'Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL'
            if ('error' in message.lower() or actual_error_bool):
                error_message = 'DFT calculation finished with an error:<br><br>' + message.replace('\n',
                                                                                                    "<br>") + '<br>Error:<br>' + err.replace(
                    '\n', '<br>') \
                                + '<br><br>Try following:<br>1.Check if the selected dft engine is correctly installed<br>' \
                                  '2. Check if the input file was correctly parsed into the respective folder (e.g. input.xml in exciting_files for exciting)'
                self.error_dialog.showMessage(error_message)

            try:
                self.load_results_from_engine(tasks)
                self.update_run_information(tasks)
            except Exception as e:
                stacktrace = get_stacktrace_as_string()
                error_message_load = 'Reading of the results of the calculation failed with error<br>' + stacktrace
                self.error_dialog.showMessage(error_message_load)

    def update_run_information(self, tasks):
        if 'scf' in tasks:
            self.last_run_information['scf'].update(esc_handler.scf_options)
            self.last_run_information['scf']['timestamp'] = time.ctime()
        if 'bandstructure' in tasks or 'g0w0' in tasks:
            bandstructure_options = {'k path': copy.deepcopy(self.dft_engine_window.band_structure_points)}
            self.last_run_information['bandstructure'].update(bandstructure_options)
        if 'g0w0' in tasks:
            self.last_run_information['gw'].update(esc_handler.gw_options)
        if 'relax' in tasks:
            self.last_run_information['relax'].update(esc_handler.relax_options)
        if 'phonons' in tasks:
            self.last_run_information.update(esc_handler.phonons_options)
        if 'optical spectrum' in tasks:
            self.last_run_information['optical spectrum'].update(esc_handler.optical_spectrum_options)

    def load_results_from_engine(self, tasks, title=None):
        if title is None:
            title = esc_handler.general_options['title']

        if 'scf' in tasks:
            scf_info = copy.deepcopy(esc_handler.scf_options)
            scf_info['timestamp'] = time.ctime()
        else:
            scf_info = copy.deepcopy(self.last_run_information['scf'])

        if 'scf' in tasks and type(self.crystal_structure) is sst.MolecularStructure:
            energy_diagram = esc_handler.read_energy_diagram()
            if energy_diagram is not None:
                energy_diagram.engine_information = {'scf': scf_info}
                self.band_structures[title] = energy_diagram

        if 'g0w0' in tasks:
            gw_info = copy.deepcopy(esc_handler.gw_options)
        else:
            gw_info = copy.deepcopy(self.last_run_information['gw'])

        if 'bandstructure' in tasks or 'g0w0' in tasks:
            read_bandstructures = []
            engine_informations = []
            titles = [title]

            if 'bandstructure' in tasks:
                if esc_handler.engine_name == 'quantum espresso':
                    coords = [x[0] for x in self.dft_engine_window.band_structure_points]
                    labels = [x[1] for x in self.dft_engine_window.band_structure_points]
                    new_coords = self.crystal_structure.convert_to_tpiba(coords)
                    band_structure_points = list(zip(new_coords, labels))
                    read_bandstructure = esc_handler.read_bandstructure(special_k_points=band_structure_points)
                elif esc_handler.engine_name == 'abinit':
                    read_bandstructure = esc_handler.read_bandstructure(
                        special_k_points=self.dft_engine_window.band_structure_points,
                        crystal_structure=self.crystal_structure)
                else:
                    read_bandstructure = esc_handler.read_bandstructure()
                read_bandstructures.append(read_bandstructure)
                engine_information = {'scf': scf_info, 'bandstructure': {
                    'k path': copy.deepcopy(self.dft_engine_window.band_structure_points)}}
                engine_informations.append(engine_information)

            if 'g0w0' in tasks:
                read_bandstructures.append(
                    esc_handler.read_gw_bandstructure(special_k_points=self.dft_engine_window.band_structure_points,
                                                      structure=self.crystal_structure))
                engine_information = {'scf': scf_info, 'bandstructure': {
                    'k path': copy.deepcopy(self.dft_engine_window.band_structure_points)},
                                      'gw': copy.deepcopy(esc_handler.gw_options)}
                engine_informations.append(engine_information)
            if 'bandstructure' in tasks and 'g0w0' in tasks:
                titles.append(title + '_gw')

            for read_bandstructure, title, engine_information in zip(read_bandstructures, titles, engine_informations):
                if read_bandstructure is None:
                    continue
                read_bandstructure.engine_information = engine_information
                self.band_structures[title] = read_bandstructure
                self.band_structure_window.update_tree()

        if 'dos' in tasks:
            dos = esc_handler.read_dos()
            engine_information = {'scf': scf_info, 'bandstructure': {
                'k path': copy.deepcopy(self.dft_engine_window.band_structure_points)}}
            dos.engine_information = engine_information
            self.dos[title] = dos
            self.dos_window.update_tree()

        if 'relax' in tasks:
            self.check_relax()
        if 'phonons' in tasks:
            engine_information = {'scf': scf_info, 'bandstructure': {
                'k path': copy.deepcopy(self.dft_engine_window.band_structure_points)},
                                  'phonon': copy.deepcopy(esc_handler.phonons_options)}
            read_bandstructure = esc_handler.read_phonon_bandstructure(
                special_k_points=self.dft_engine_window.band_structure_points, structure=self.crystal_structure)
            read_bandstructure.engine_information = engine_information
            self.band_structures[title] = read_bandstructure
            self.band_structure_window.update_tree()
        if 'optical spectrum' in tasks:
            engine_information = {'scf': scf_info,
                                  'optical spectrum': copy.deepcopy(esc_handler.optical_spectrum_options),
                                  'gw': gw_info}
            read_spectrum = esc_handler.read_optical_spectrum()
            read_spectrum.engine_information = engine_information
            self.optical_spectra[title] = read_spectrum
            self.optical_spectra_window.update_tree()

    def open_engine_option_window(self):
        self.engine_option_window.update_all()
        self.engine_option_window.show()

    def open_state_vis_window(self):
        self.ks_state_window.plot_widget.update_tree()
        self.ks_state_window.show()

    def open_structure_window(self, new=False):
        if new:
            self.structure_window.set_structure(None)
        else:
            self.structure_window.set_structure(self.crystal_structure)

        self.structure_window.anything_changed = False
        self.structure_window.update_fields()
        self.structure_window.switch_enabled_fields()
        self.structure_window.update_plot()
        self.structure_window.update_info()
        self.structure_window.show()

    def open_brillouin_window(self):
        if type(self.crystal_structure) is not sst.CrystalStructure:
            return

        if self.crystal_structure is not None and self.brillouin_window.mayavi_widget.crystal_structure is not self.crystal_structure and type(
                self.crystal_structure):
            self.brillouin_window.close()  # This is a hack because for some reason the picker is broken when you update the plot
            self.brillouin_window = BrillouinWindow(self)
            self.brillouin_window.mayavi_widget.set_crystal_structure(self.crystal_structure)
            self.brillouin_window.set_path(self.dft_engine_window.band_structure_points)
            self.brillouin_window.mayavi_widget.update_plot()

        self.brillouin_window.show()

    def open_load_result_window(self, task):
        result_window = LoadResultsWindow(self, task)
        result_window.show()

    def open_scripting_console(self):
        self.dft_engine_window.read_all_option_widgets()
        self.console_window.show()

        def add_plot_to_queue(structure):
            q_item = {'task': 'plot structure', 'structure': structure}
            self.queue.put(q_item)

        def add_scf_to_queue(scf_data):
            q_item = {'task': 'plot scf', 'scf data': scf_data}
            self.queue.put(q_item)
            time.sleep(0.3)

        shared_vars = {'structure': self.crystal_structure, 'engine': esc_handler, 'plot_structure': add_plot_to_queue,
                       'CrystalStructure': sst.CrystalStructure,
                       'MolecularStructure': sst.MolecularStructure, 'OpticalSpectrum': sst.OpticalSpectrum,
                       'BandStructure': sst.BandStructure, 'EnergyDiagram': sst.EnergyDiagram,
                       'KohnShamDensity': sst.KohnShamDensity, 'MolecularDensity': sst.MolecularDensity,
                       'plot_scf': add_scf_to_queue}
        self.console_window.python_interpreter.update_vars(shared_vars)

    def open_phonon_window(self):
        self.phonon_window.update_crystal_structure(self.crystal_structure)
        self.phonon_window.update_plot()
        self.phonon_window.update_tree()
        self.phonon_window.show()

    def configure_buttons(self, disable_all=False):
        self.dft_engine_window.configure_buttons(disable_all=disable_all)
        if self.project_directory is None:
            self.save_project_action.setEnabled(False)
            self.edit_structure_action.setEnabled(False)
            self.new_structure_action.setEnabled(False)
            self.import_structure_menu.setEnabled(False)
            self.vis_menu.setEnabled(False)
        else:
            self.save_project_action.setEnabled(True)
            self.edit_structure_action.setEnabled(True)
            self.new_structure_action.setEnabled(True)
            self.import_structure_menu.setEnabled(True)
            self.vis_menu.setEnabled(True)

    def handle_queue(self):
        if self.queue.empty():
            return

        queue_item = self.queue.get()
        taskname = queue_item['task']

        if taskname == 'plot structure':
            structure = queue_item['structure']
            self.mayavi_widget.update_crystal_structure(structure)
            self.mayavi_widget.update_plot()
            QtGui.QApplication.processEvents()
        elif taskname == 'plot scf':
            scf_data = queue_item['scf data']
            if scf_data is not None:
                self.scf_window.scf_widget.plot(scf_data)
            QtGui.QApplication.processEvents()

    def check_integrety(self):
        scf_check = self.dft_engine_window.scf_option_widget.options == esc_handler.scf_options
        gw_check = self.dft_engine_window.gw_option_widget.options == esc_handler.gw_options
        phonon_check = self.dft_engine_window.phonons_option_widget.options == esc_handler.phonons_options
        optical_check = self.dft_engine_window.optical_spectrum_option_widget.options == esc_handler.optical_spectrum_options

        checks = [scf_check, gw_check, phonon_check, optical_check]

        if not all(checks):
            warnings.warn('Option dictionaries were not connected anymore')
            logging.warning('Option dictionaries were not connected')
            self.dft_engine_window.scf_option_widget.options = esc_handler.scf_options
            self.dft_engine_window.gw_option_widget.options = esc_handler.gw_options
            self.dft_engine_window.phonons_option_widget.options = esc_handler.phonons_options
            self.dft_engine_window.optical_spectrum_option_widget.options = esc_handler.optical_spectrum_options

    def update_program(self):
        self.check_integrety()
        self.handle_queue()

    def check_and_set_lock(self,folder_name):
        lock_path = folder_name+'/.lock'
        if os.path.isfile(lock_path):
            msg = "This project seems to be opened by another instance of OpenDFT. Opening one project by two instances is not safe. Proceed?"
            reply = QtGui.QMessageBox.question(self, 'Project locked',
                                               msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

            if reply == QtGui.QMessageBox.No:
                raise Exception('Project is locked')

        open(lock_path, 'a').close()

    def release_lock(self):
        lock_path = self.project_directory +'/.lock'
        if os.path.isfile(lock_path):
            os.remove(lock_path)


if __name__ == "__main__":

    current_time = time.localtime()
    current_time_string = [str(x) for x in current_time[:3]]

    # installation_folder = os.path.dirname(find_data_file(''))
    temp_folder = os.path.expanduser("~") + "/.OpenDFT"

    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    if not os.path.exists(temp_folder + "/logfiles/"):
        os.makedirs(temp_folder + "/logfiles/")

    logging.basicConfig(level=logging.DEBUG,
                        filename=temp_folder + "/logfiles/" + "_".join(current_time_string) + ".log")
    logging.info('Program started')

    import vtk

    vtk_output = vtk.vtkFileOutputWindow()  # redirects vtk errors to log file
    vtk_output.SetFileName(temp_folder + "/logfiles/vtk_log.txt")
    vtk.vtkOutputWindow().SetInstance(vtk_output)

    set_procname(b'OpenDFT')

    app = QtGui.QApplication.instance()
    main = CentralWindow(parent=app)

    app.exec_()
