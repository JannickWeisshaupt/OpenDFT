from pyface.qt import QtGui, QtCore
import numpy as np
import sys
import os
import traceback
from collections import OrderedDict


class no_error_dictionary:
    def __init__(self,dic_in):
        self.dic = dic_in

    def __getitem__(self, key):
        try:
            return self.dic[key]
        except KeyError:
            return None


def flatten_dictionary(dictionary):
    new_dic = {}
    for key,value in dictionary.items():
        for key2,value in value.items():
            new_dic[key+'_'+key2] = value
    return new_dic


class CopySelectedCellsAction(QtGui.QAction):
    def __init__(self, table_widget):
        if not isinstance(table_widget, QtGui.QTableWidget):
            raise ValueError(str('CopySelectedCellsAction must be initialised with a QTableWidget. A %s was given.' % type(table_widget)))
        super(CopySelectedCellsAction, self).__init__("Copy", table_widget)
        self.setShortcut('Ctrl+c')
        self.triggered.connect(self.copy_cells_to_clipboard)
        self.table_widget = table_widget

    def copy_cells_to_clipboard(self):
        if len(self.table_widget.selectionModel().selectedIndexes()) > 0:
            # sort select indexes into rows and columns
            previous = self.table_widget.selectionModel().selectedIndexes()[0]
            columns = []
            rows = []
            for index in self.table_widget.selectionModel().selectedIndexes():
                if previous.column() != index.column():
                    columns.append(rows)
                    rows = []
                rows.append(index.data())
                previous = index
            columns.append(rows)

            # add rows and columns to clipboard
            clipboard = ""
            nrows = len(columns[0])
            ncols = len(columns)
            for r in range(nrows):
                for c in range(ncols):
                    clipboard += columns[c][r]
                    if c != (ncols-1):
                        clipboard += '\t'
                clipboard += '\n'

            # copy to the system clipboard
            sys_clip = QtGui.QApplication.clipboard()
            sys_clip.setText(clipboard)

class PasteIntoTable(QtGui.QAction):
    def __init__(self, table_widget,parent):
        if not isinstance(table_widget, QtGui.QTableWidget):
            raise ValueError(str('CopySelectedCellsAction must be initialised with a QTableWidget. A %s was given.' % type(table_widget)))
        super(PasteIntoTable, self).__init__("Copy", table_widget)
        self.setShortcut('Ctrl+v')
        self.parent = parent
        self.triggered.connect(self.paste_cells_from_clipboard)
        self.table_widget = table_widget

    def paste_cells_from_clipboard(self):
        self.parent.disconnect_tables()
        index = self.table_widget.selectedIndexes()[0]
        i0 = index.row()
        j0 = index.column()

        sys_clip = QtGui.QApplication.clipboard()
        clipboard = sys_clip.text()


        rows = clipboard.split('\n')
        n_rows = len(rows)
        col_0 = rows[0].split()
        n_col = len(col_0)

        # data = np.zeros((n_rows,n_col))

        for i,row in enumerate(rows):
            col = row.split()
            for j,el in enumerate(col):
                self.table_widget.item(i+i0,j+j0).setText(el)
        self.parent.connect_tables()
        self.parent.handle_change()

def convert_to_ordered(d):
    return OrderedDict(sorted(d.items(), key=lambda t: t[0]))

def set_procname(newname):
    from ctypes import cdll, byref, create_string_buffer
    libc = cdll.LoadLibrary('libc.so.6')    #Loading a 3rd party library C
    buff = create_string_buffer(len(newname)+1) #Note: One larger than the name (man prctl says that)
    buff.value = newname                 #Null terminated string as it should be
    libc.prctl(15, byref(buff), 0, 0, 0)

def get_proc_name():
    from ctypes import cdll, byref, create_string_buffer
    libc = cdll.LoadLibrary('libc.so.6')
    buff = create_string_buffer(128)
    # 16 == PR_GET_NAME from <linux/prctl.h>
    libc.prctl(16, byref(buff), 0, 0, 0)
    return buff.value

def find_data_file(filename):
    if getattr(sys, 'frozen', False):
        # The application is frozen
        datadir = os.path.dirname(sys.executable)
    else:
        # The application is not frozen
        # Change this bit to match where you store your data files:
        datadir = os.path.dirname(__file__)

    return datadir + filename

def get_stacktrace_as_string():
    exc_type, exc_value, exc_traceback = sys.exc_info()
    error = traceback.format_exception(exc_type, exc_value, exc_traceback)
    joined_error = '<br>'.join(error)
    joined_error = joined_error.replace(' ','&#160;')
    return joined_error