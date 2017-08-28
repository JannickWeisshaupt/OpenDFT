
import solid_state_tools as sst
import xml.etree.ElementTree as ET


class Handler:

    def parse_input_file(self,filename):
        tree = ET.parse(filename)
        root = tree.getroot()

    def __init__(self):
        self.default_extension = '.xml'

