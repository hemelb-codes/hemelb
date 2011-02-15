from HemeLbSetupTool.Util.Observer import Observable

class Vector(Observable):
    def __init__(self, *args):
        if len(args) == 0:
            nan = float("nan")
            self.x = nan
            self.y = nan
            self.z = nan
        elif len(args) == 1:
            self.x, self.y, self.z = args[0]
        elif len(args) == 3:
            self.x, self.y, self.z = args
        else:
            raise ValueError

        return
    
    pass


