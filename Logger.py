import os
import time


class Logger:

    def __init__(self, filepath):
        self.buffer = []
        self.pid = os.getpid()
        self.filepath = filepath
        self.k = 0

    def log(self, msg):
        t = time.time()
        self.buffer.append(
            (self.pid, t, msg)
        )
        self.k += 1
        if self.k > 50:
            self.write()

    def write(self):
        with open(self.filepath, 'a') as f:
            f.write('\n'.join(f'{e[0]},{e[1]},{e[2]}' for e in self.buffer) + '\n')
        self.buffer = []
        self.k = 0

