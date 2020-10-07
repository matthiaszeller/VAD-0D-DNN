import os
import time


class Logger(object):
    """Convenient object to manage logging (writing in file, print on stdout).
    Implements a buffer to reduce disk access and reduce latency.
    WARNING: you must call `flush()` before deleting a `Logger` instance, otherwise data
    retained in the buffer will be lost.

    Note: `Logger` can be used to redirect stdout:
        >>> import sys
        >>> old_stdout = sys.stdout
        >>> sys.stdout = Logger(...)
        >>> # stuff that outputs on stdout...
        >>> sys.stdout = old_stdout
    """
    # TODO: make this class a context manager !

    def __init__(self, filepath, print=False, buffer_size=50):
        """
        :param str filepath: the file to write in
        :param bool print: whether to print on stdout in addition to writing in `filepath`
        :param int buffer_size: the number of lines to store before writing in file
        """
        self.buffer = []
        self.pid = os.getpid()
        self.filepath = filepath
        self.__k = 0
        self.__logfun = self.__log_and_print if print else self.__log
        self.__buffer_size = buffer_size

    def log(self, msg):
        self.__logfun(msg)

    def write(self, msg):
        """Equivalent to `log(msg)`, but does not log if msg == '\n'.
        This is useful to redirect the sys.stdout to a custom logger."""
        if msg == '\n':
            return

        self.__logfun(msg)

    def __log(self, msg):
        t = time.time()
        self.buffer.append(
            (self.pid, t, msg)
        )
        self.__k += 1
        if self.__k > self.__buffer_size:
            self.flush()

    def __log_and_print(self, msg):
        t = time.time()
        self.buffer.append(
            (self.pid, t, msg)
        )
        print(os.getpid(), msg)
        self.__k += 1
        if self.__k > self.__buffer_size:
            self.flush()

    def flush(self):
        with open(self.filepath, 'a') as f:
            f.write('\n'.join(f'{e[0]};{e[1]};{e[2]}' for e in self.buffer) + '\n')
        self.buffer = []
        self.__k = 0

