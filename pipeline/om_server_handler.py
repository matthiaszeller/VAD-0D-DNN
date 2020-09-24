
import threading
import psutil
from time import sleep, time
from typing import List


class OMServerHanlder(object):
    """Manage the OpenModelica servers. The OM sessions launched with multiprocessing
    do not terminate once the variable bound to the session is deleted, which wastes
    some RAM. Note that those sessions are only used to compile the 0D model, which
    only takes a few seconds.

    Monitoring old omc processes and kill them after a timeout is the only way through
    the issue so far."""

    monitored_omc_ps: List[psutil.Process]

    def __init__(self, timeout=15, logfun=None):
        """
        :param float timeout: time after which to kill the OM session in seconds
        """
        super().__init__()
        self.monitored_omc_ps = []
        self._sleep_time = 4
        self._omc_session_timeout = timeout
        self._logfun = (lambda _: None) if logfun is None else logfun
        self.__maxiter = 1000

        thread = threading.Thread(target=self.__loop, args=())
        thread.daemon = True
        thread.start()

    def kill_sessions(self):
        self.__update_processes()
        N = len(self.monitored_omc_ps)
        for p in self.monitored_omc_ps:
            p.kill()
        self._logfun(f'OMServerHanlder: Killed {N} remaining OM sessions')

    def __loop(self):
        while True:
            # ---- Monitor
            self.__update_processes()
            # --- Kill if necessary
            self.__kill_outdated_sessions()
            # --- Sleep
            print('here i am')
            sleep(self._sleep_time)

    def __kill_outdated_sessions(self):
        t = time()
        for ps in self.monitored_omc_ps:
            if t - ps.create_time() > self._omc_session_timeout:
                pid = ps.pid
                ps.kill()
                self._logfun(f'OMServerHanlder: Killed OM session [pid {pid}]')

    def __update_processes(self):
        self.monitored_omc_ps = [
            ps for ps in psutil.process_iter() if ps.name() == 'omc'
        ]
