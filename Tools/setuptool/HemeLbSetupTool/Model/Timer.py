import time


class Timer(object):

    def __init__(self):
        self._running = False

    def Start(self):
        assert not self._running
        self._running = True
        self._startTime = time.clock()
        return

    def Stop(self):
        assert self._running
        self._stopTime = time.clock()
        self._running = False
        return

    def GetTime(self):
        if self._running:
            return time.clock() - self._startTime
        return self._stopTime - self._startTime

    pass
