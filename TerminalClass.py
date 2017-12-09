import code
import contextlib
import sys
if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from io import BytesIO as StringIO
import copy


@contextlib.contextmanager
def capture():
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out=[StringIO(), StringIO()]
        sys.stdout,sys.stderr = out
        yield out
    finally:
        sys.stdout,sys.stderr = oldout, olderr
        out[0] = out[0].getvalue()
        out[1] = out[1].getvalue()


# class Parser(threading.Thread):
#     output_lock = threading.RLock()
#
#     def __init__(self,target,arguments, *args,**kwargs):
#         super().__init__(*args, **kwargs)
#         self.target=target
#         self.arguments=arguments
#
#
#     def run(self):
#         with self.output_lock:
#             self.target(*self.arguments)

class PythonTerminal(code.InteractiveConsole,object):

    def __init__(self, shared_vars):
        self.shared_vars_start = copy.deepcopy(shared_vars)
        self.shared_vars = shared_vars
        super(PythonTerminal,self).__init__(shared_vars)
        self.out_history = []

    def run_code(self,code_string):
        with capture() as out:
            # for line in code_string.split('\n'):
            #     self.push(line)

            # t = Parser(self.runcode,(code_string,))
            # t.start()
            if type(code_string) == str:
                code_string = unicode(code_string)
            self.runcode(code_string)

        self.out_history.append(out)
        return out

    def restart_interpreter(self):
        self.__init__(self.shared_vars_start)

    def stop(self):
        raise NotImplementedError

    def update_vars(self,vars):
        self.shared_vars.update(vars)

if __name__ == '__main__':
    import numpy as np
    a = np.array(range(10))
    PyTerm = PythonTerminal({'Betrag': a})
    test_code = u"""
print(Betrag.sum())
print(Betrag.mean())
print(Betrag.std())
for i in range(10):
    print(i)
    print('This is working')


print("Code ends here")
import matplotlib.pyplot as plt

# plt.figure()
# plt.plot([1,2,3],[1,2,3])
# plt.show()
#
# plt.figure()
# plt.plot([1,2,3],[1,4,9])
# plt.show()

test_awesomeness = 'This is awesome'

"""
    test2 = """
test_awesomeness = 'This is awesome'
print(test_awesomeness)
print(test_awesomeness)
print('Yeah')

"""

    out = PyTerm.run_code(test_code)
    print(out)

    out2 = PyTerm.run_code(test2)
    print(out2)

