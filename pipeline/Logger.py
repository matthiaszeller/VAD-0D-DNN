import os
import time


class Logger(object):

    def __init__(self, filepath, print=False):
        self.buffer = []
        self.pid = os.getpid()
        self.filepath = filepath
        self.k = 0
        self._logfun = self.__log_and_print if print else self.__log

    def log(self, msg):
        self._logfun(msg)

    def write(self, msg):
        """Equivalent to `log(msg)`, but does not log if msg == '\n'.
        This is useful to redirect the sys.stdout to a custom logger."""
        if msg == '\n':
            return

        self._logfun(msg)

    def __log(self, msg):
        t = time.time()
        self.buffer.append(
            (self.pid, t, msg)
        )
        self.k += 1
        if self.k > 50:
            self.flush()

    def __log_and_print(self, msg):
        t = time.time()
        self.buffer.append(
            (self.pid, t, msg)
        )
        print(os.getpid(), msg)
        self.k += 1
        if self.k > 50:
            self.flush()

    def flush(self):
        with open(self.filepath, 'a') as f:
            f.write('\n'.join(f'{e[0]},{e[1]},{e[2]}' for e in self.buffer) + '\n')
        self.buffer = []
        self.k = 0

