import numpy as np

class SortedVector(object):
    """Simple vector that maintains its elements in sorted order.
    """
    
    def __init__(self, capacity=16, dtype=int):
        self.capacity = capacity
        self.array = np.zeros(capacity, dtype=dtype)
        self.size = 0

    def add(self, val):
        """Insert a value.
        Find - log N
        Insert - N
        """
        idx = self.array[:self.size].searchsorted(val)
        if self.size + 1 > self.capacity:
            # enlarge if needed
            new = np.zeros(self.capacity*2, dtype=self.array.dtype)
            new[:self.size] = self.array[:]
            self.array = new
            
        self.array[idx+1:self.size+1] = self.array[idx:self.size]
        self.array[idx] = val
        return

    def __contains__(self, val):
        """Return True if val is in the array.
        Log N complexity.
        """
        idx = self.array[:self.size].searchsorted(val)
        return self.array[idx] == val
    pass
