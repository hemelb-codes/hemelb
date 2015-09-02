"""_fixpaths.py

Sort out imports of support code from the PV scripts.
"""
import sys
import os.path
import inspect
fullPathToSelf = os.path.abspath(inspect.getfile(inspect.currentframe()))
# Walk up the tree starting from this file to get the scripts directory
current = fullPathToSelf
last = None
while last != "scripts":
    current, last = os.path.split(current)
scripts = os.path.join(current, last)
# If scripts isn't on the $PYTHONPATH, add it.
if scripts not in sys.path:
    sys.path.insert(0, scripts)
