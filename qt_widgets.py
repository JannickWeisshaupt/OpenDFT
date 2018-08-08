from pyface.qt import QtGui, QtCore


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