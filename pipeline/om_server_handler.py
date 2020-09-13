
import threading
import psutil
from time import sleep, time
from typing import List


class OMServerHanlder(object):
    """Manage the OpenModelica servers. The OM sessions launched with multiprocessing
    do not terminate once the variable bound to the session is deleted, which wastes
    some RAM.

    Monitoring old omc processes and kill them after a timeout is the only way through
    the issue so far."""

    monitored_omc_ps: List[psutil.Process]

    def __init__(self):
        super().__init__()
        self.monitored_omc_ps = []
        self._sleep_time = 4
        self._omc_session_timeout = 20

        thread = threading.Thread(target=self.__loop, args=())
        thread.daemon = True
        thread.start()

    def kill_sessions(self):
        self.__update_processes()
        for p in self.monitored_omc_ps:
            p.kill()

    def __loop(self):
        # ---- Monitor
        self.__update_processes()
        # --- Kill if necessary
        self.__kill_outdated_sessions()
        # --- Loop
        sleep(self._sleep_time)
        self.__loop()

    def __kill_outdated_sessions(self):
        t = time()
        for ps in self.monitored_omc_ps:
            if t - ps.create_time() > self._omc_session_timeout:
                ps.kill()

    def __update_processes(self):
        self.monitored_omc_ps = [
            ps for ps in psutil.process_iter() if ps.name() == 'omc'
        ]
