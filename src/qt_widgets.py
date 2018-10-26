from pyface.qt import QtGui, QtCore
from little_helpers import find_data_file

class CopySelectedCellsAction(QtGui.QAction):
    def __init__(self, table_widget):
        if not isinstance(table_widget, QtGui.QTableWidget):
            raise ValueError(str(
                'CopySelectedCellsAction must be initialised with a QTableWidget. A %s was given.' % type(
                    table_widget)))
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
                    if c != (ncols - 1):
                        clipboard += '\t'
                clipboard += '\n'

            # copy to the system clipboard
            sys_clip = QtGui.QApplication.clipboard()
            sys_clip.setText(clipboard)


class PasteIntoTable(QtGui.QAction):
    def __init__(self, table_widget, parent):
        if not isinstance(table_widget, QtGui.QTableWidget):
            raise ValueError(str(
                'CopySelectedCellsAction must be initialised with a QTableWidget. A %s was given.' % type(
                    table_widget)))
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

        for i, row in enumerate(rows):
            col = row.split()
            for j, el in enumerate(col):
                self.table_widget.item(i + i0, j + j0).setText(el)
        self.parent.connect_tables()
        self.parent.handle_change()


class EntryWithLabel(QtGui.QWidget):
    def __init__(self, parent, label, value=None, width_text=200, width_label=90,tooltip=None):
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

        if tooltip is not None:
            self.textbox.setToolTip(tooltip)

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

    def setValidator(self,validator):
        self.textbox.setValidator(validator)


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


class MySearchLineEdit(QtGui.QLineEdit):
    def __init__(self, parent=None, delay=500, callback=None):
        super(MySearchLineEdit, self).__init__(parent)
        self.delay = delay
        self.callback = callback
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.setTextMargins(20, 0, 0, 0)
        self.connect(self, QtCore.SIGNAL("textEdited(QString)"),
                self.eventDelay)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                self.startCallback)

    def eventDelay(self, text):
        self.timer.stop()
        self.timer.start(self.delay)

    def startCallback(self):
        text = self.text()
        if self.callback:
            self.callback(text)

    def paintEvent(self, event=None):
        super(MySearchLineEdit, self).paintEvent(event)
        pixmap = QtGui.QPixmap(find_data_file("/data/icons/search-icon.png")).scaled(16, 16, QtCore.Qt.KeepAspectRatio,QtCore.Qt.SmoothTransformation)
        painter = QtGui.QPainter(self)
        pixmap_width = pixmap.width()
        pixmap_height = pixmap.height()
        left_border = 6
        painter.drawPixmap((left_border),
                           (self.height() - pixmap.height()) / 2,
                           pixmap)


class QFloatTableWidgetItem(QtGui.QTableWidgetItem):
    def __init__ (self, value):
        super(QFloatTableWidgetItem, self).__init__()
        self.setText('{0:1.2f}'.format(value))

    def __lt__ (self, other):
        if (isinstance(other, QFloatTableWidgetItem)):
            selfDataValue  = float(self.text())
            otherDataValue = float(other.text())
            return selfDataValue < otherDataValue
        else:
            return QtGui.QTableWidgetItem.__lt__(self, other)

def make_splash_screen():
    splash_pix = QtGui.QPixmap(find_data_file('/data/artwork/gaas_cubic_with_header.png')).scaled(500, 500, QtCore.Qt.KeepAspectRatio,QtCore.Qt.SmoothTransformation)
    splash = QtGui.QSplashScreen(splash_pix, QtCore.Qt.WindowStaysOnTopHint)
    splash.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint | QtCore.Qt.FramelessWindowHint)
    splash.setEnabled(False)
    splash.setWindowIcon(QtGui.QIcon(find_data_file('/data/icons/icon.ico')))
    splash.setWindowTitle('OpenDFT')
    # splash = QSplashScreen(splash_pix)
    # adding progress bar
    progressBar = QtGui.QProgressBar(splash)
    progressBar.setTextVisible(False)
    progressBar.setMaximum(10)
    progressBar.setGeometry(0, splash_pix.height() - 50, splash_pix.width(), 20)

    # splash.setMask(splash_pix.mask())

    splash.show()
    # splash.showMessage("<h1><font color='white',font-size: 350%;>OpenDFT</font></h1>", QtCore.Qt.AlignTop | QtCore.Qt.AlignCenter, QtCore.Qt.black)
    return splash,progressBar