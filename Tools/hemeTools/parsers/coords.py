"""Parse the coordinates file that results from the HemeLB segtool.

"""
def parse(filename):
    parDict = {}
    execfile(filename, {}, parDict)
    return parDict

