# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from ..Util.Observer import Observable

class EmptySelection(Observable):
    """Represent an empty selection.
    """
    def __getattr__(self, attr):
        return self
    
    pass

EmptySelection = EmptySelection()


def isNone(value):
    return (value is None) or (value is EmptySelection)
