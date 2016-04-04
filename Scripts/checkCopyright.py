#!/usr/bin/env python
#
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""
This script checks all files in the repo for the correct copyright statement.
"""

import os.path
import subprocess
import re
import sys

py_coding_re = re.compile('#.*coding[=:]\s*([-\w.]+)')

def info(filename, msg):
    if verbose:
        print filename, msg

def error(filename, msg):
    error.occurred = True
    print filename, msg
error.occured = False

comment_chars = {
    '.applescript': '--',
    '.cc' : '//',
    '.cmake': '#',
    '.cpp' : '//',
    '.sh':'#',
    '.i' : '//',
    '.hpp' : '//',
    '.h' : '//',
    '.py' : '#',
    '.pyx' : '#',
    '.pxd' : '#',
    '.yml' : '#',
    'README': '',
    }

excluded_files = set((
    'Tools/setuptool/HemeLbSetupTool/Model/Generation/Point_inside_polyhedron_3.h',
    'Tools/setuptool/HemeLbSetupTool/Model/Generation/Triangle_3_Ray_3_do_intersect.h',
    ))

def comment_char(filename):
    base = os.path.basename(filename)
    root, ext = os.path.splitext(base)
    if base == 'CMakeLists.txt':
        ext = '.cmake'
    
    return comment_chars[ext]
    
def is_source(filename):
    try:
        cc = comment_char(filename)
        return True
    except KeyError:
        return False    

def is_in_git(filename):
    git_ls_cmd = ['git', 'ls-tree', '--name-only', '--full-name', 'HEAD', filename]
    output = subprocess.check_output(git_ls_cmd).strip()
    return output == filename

def check_cr(filename):
    if not is_in_git(filename):
        info(filename, "Skipping non-VCed file")
        return
    if not is_source(filename):
        info(filename, "Skipping non-source file")
        return

    if filename in excluded_files:
        info(filename, "Skipping file as not checked")
        return
    
    info(filename, "Checking file")
        
    cc = comment_char(filename)
    ncc = len(cc)
    root, ext = os.path.splitext(filename)
    
    f = open(filename)
    
    line = f.readline()

    if ext in ('.py', '.sh'):
        # .sh and .py files can have a #!
        if line.startswith('#!'):
            # discard
            line = f.readline()

    if ext == '.py':
        # .py may have an encoding spec in lines 1 or 2
        if py_coding_re.match(line):
            # discard
            line = f.readline()
    
    # Can start with blank lines
    while line.strip() in ('', cc):
        line = f.readline()
    
    # CR statement is 4 lines and must(?) appear now
    file_msg = []
    for i in xrange(4):
        if not line.startswith(cc):
            warn(filename + ": line that should be part of the CR statement doesn't start with a comment character")
            
        file_msg.append(line[ncc:].strip())
        line = f.readline()
        
    file_msg = '\n'.join(file_msg)
    if file_msg != cr_msg:
        error(filename, "copyright statement doesn't read as expected")
    else:
        info(filename, "OK")
    return

cr_msg =  """This file is part of HemeLB and is Copyright (C)
the HemeLB team and/or their institutions, as detailed in the
file AUTHORS. This software is provided under the terms of the
license in the file LICENSE."""

if __name__ == "__main__":
    try:
        if sys.argv[1] in ('-v', '--verbose'):
            verbose = True
        else:
            verbose = False
    except IndexError:
        verbose = False
        pass

    self_path = os.path.abspath(__file__)
    scripts_dir = os.path.dirname(self_path)
    repo_root = os.path.dirname(scripts_dir)
    
    for dirpath, dirnames, filenames in os.walk(repo_root):
        if '.git' in dirnames:
            dirnames.remove('.git')
        
        for fn in filenames:
            check_cr(os.path.relpath(os.path.join(dirpath, fn), repo_root))
    
    sys.exit(1 if error.occurred else 0)
