import contextlib
import os

def IndentXml(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            IndentXml(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
    return

@contextlib.contextmanager
def cd(path):
    start = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(start)
