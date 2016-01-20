
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
def IterPairs(iterable):
    """Return adjacent pairs of items from iterable. I.e.:
    IterPairs(lst) -> (lst[0], lst[1]), (lst[1], lst[2]), (lst[2], lst[3]), ..., (lst[i], lst[i+1])
    """
    first = True
    for item in iterable:
        if first:
            new = item
            first = False
            continue
        else:
            old = new
            new = item
            yield old, new
        continue
    return
