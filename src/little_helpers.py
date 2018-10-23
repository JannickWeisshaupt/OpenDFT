import numpy as np
import sys
import os
import traceback
from collections import OrderedDict
import ast
import operator as op
from fractions import Fraction
from collections import deque



class DequeSet(deque):
    def __init__(self,iterable=(), maxlen=None):
        iterable_new = []
        for x in iterable:
            if x not in iterable_new:
                iterable_new.append(x)

        super().__init__(iterable_new,maxlen)

    def append(self, x):
        if x in self:
          self.remove(x)
        super().append(x)

    def appendleft(self, x):
        if x in self:
          self.remove(x)
        super().appendleft(x)

    def extend(self, x):
        for el in x:
            self.append(el)

    def extendleft(self, x):
        for el in x[::-1]:
            self.appendleft(el)

def find_possible_sums(repeat):
    possible_sums = set()

    for i in range(repeat[0] + 1):
        for j in range(repeat[1] + 1):
            for l in range(repeat[2] + 1):
                pos_sum = (i, j, l)
                possible_sums.add(pos_sum)
    return possible_sums


def find_grid_connections(possible_sums):
    connections = set()
    for i, pos_sum in enumerate(possible_sums):
        for j, pos_sum_1 in enumerate(possible_sums):
            diff_sum = (abs(x - y) for x, y in zip(pos_sum, pos_sum_1))
            if sum(diff_sum) == 1:
                connections.add(tuple(sorted([i, j])))
    return connections


def convert_to_ordered(d):
    return OrderedDict(sorted(d.items(), key=lambda t: t[0]))


def set_procname(newname):
    if sys.platform in ['linux', 'linux2']:
        from ctypes import cdll, byref, create_string_buffer
        libc = cdll.LoadLibrary('libc.so.6')  # Loading a 3rd party library C
        buff = create_string_buffer(len(newname) + 1)  # Note: One larger than the name (man prctl says that)
        buff.value = newname  # Null terminated string as it should be
        libc.prctl(15, byref(buff), 0, 0, 0)
    else:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(newname.decode())


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


def get_stacktrace_as_string(html=True):
    exc_type, exc_value, exc_traceback = sys.exc_info()
    error = traceback.format_exception(exc_type, exc_value, exc_traceback)
    if html:
        page_break = '<br>'
    else:
        page_break = '\n'
    joined_error = page_break.join(error)
    if html:
        joined_error = joined_error.replace(' ', '&#160;')
    return joined_error

# supported operators
operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg,'sin':np.sin,'cos':np.cos,'tan':np.tan,'sqrt':np.sqrt,'pi':np.pi,
             'sind': lambda x: np.sin(x*np.pi/180), 'cosd':lambda x: np.cos(x*np.pi/180), 'tand': lambda x: np.tan(x*np.pi/180)}


def eval_expr(expr):
    return eval_(ast.parse(expr, mode='eval').body)


def eval_(node):
    if isinstance(node, ast.Num): # <number>
        return node.n
    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
        return operators[type(node.op)](eval_(node.left), eval_(node.right))
    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return operators[type(node.op)](eval_(node.operand))
    elif isinstance(node, ast.Call): # <func>(<arg>)
        return operators[node.func.id](eval_(node.args[0]))
    elif isinstance(node, ast.Name):
        return operators[node.id]
    else:
        raise TypeError(node)


def find_fraction(num):
    frac = Fraction(num).limit_denominator(16)
    if frac.numerator == 0:
        return '0'
    elif frac.denominator == 1:
        return '{}'.format(frac.numerator)
    elif abs(frac-num) > 1e-14:
        return '{0:1.11f}'.format(num)
    else:
        return '{0}/{1}'.format(frac.numerator, frac.denominator)


def convert_to_bool(x,match_case = False):
    if not match_case:
        x = x.lower()
    x = x.strip()

    if x in ['true', '1']:
        return True
    elif x in ['false', '0']:
        return False
    else:
        raise ValueError('Bad input to convert to bool with' + str(x))


def fix_requests_if_frozen():
    def override_where():
        """ overrides certifi.core.where to return actual location of cacert.pem"""
        # change this to match the location of cacert.pem
        return find_data_file('/data/cacert.pem')

    # is the program compiled?
    if hasattr(sys, "frozen"):
        import certifi.core

        os.environ["REQUESTS_CA_BUNDLE"] = override_where()
        certifi.core.where = override_where

        # delay importing until after where() has been replaced
        import requests.utils
        import requests.adapters
        # replace these variables in case these modules were
        # imported before we replaced certifi.core.where
        requests.utils.DEFAULT_CA_BUNDLE_PATH = override_where()
        requests.adapters.DEFAULT_CA_BUNDLE_PATH = override_where()

if __name__ == "__main__":
    # print(eval_expr('sind(45)'))
    print(find_fraction(6/32))