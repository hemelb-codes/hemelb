from HemeLbSetupTool.Util.Observer import Observable

class EmptySelection(Observable):
    """Represent an empty selection.
    """
    def __getattr__(self, attr):
        return self
    
    pass

EmptySelection = EmptySelection()


def isNone(value):
    return (value is None) or (value is EmptySelection)
