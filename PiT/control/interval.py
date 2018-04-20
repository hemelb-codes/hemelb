import numpy as np

class Interval(object):
    """Simple time interval representation.
    """
    def __init__(self, ti, tf, nsteps=None, dt=None):
        if nsteps is None:
            if dt is None:
                raise ValueError("nsteps and dt can't both be None")
            else:
                nsteps = int((tf-ti)/dt)
                pass
        else:
            if dt is None:
                dt = (tf - ti)/nsteps
            else:
                raise ValueError("Can't specify both nsteps and dt")
            pass

        self.start = ti
        self.end = tf
        self.dt = dt
        self.nsteps = nsteps
        return

    def generate(self):
        return np.arange(self.start, self.end, self.dt)
    
    def subdivide(self, num_time_slices):
        slice_steps = self.nsteps / num_time_slices
        bounds = np.arange(num_time_slices + 1, dtype=float) * ((self.end-self.start) / float(num_time_slices)) + self.start
        starts = bounds[:-1]
        ends = bounds[1:]
        return [Interval(float(starts[i]), float(ends[i]), nsteps=slice_steps) for i in range(num_time_slices)]

    def to_dict(self):
        return {
            'start': self.start,
            'stop': self.end,
            'nstep': self.nsteps
            }
    pass
