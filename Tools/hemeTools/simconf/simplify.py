
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
from functools import wraps
import quantities as pq

def simplify(func):
    """Decorator to simplify the units of the return value of a function.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        ans = func(*args, **kwargs)
        if isinstance(ans, pq.Quantity):
            return ans.simplified

        return ans
    
    return wrapper
