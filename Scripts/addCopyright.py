#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 


import os
import re
import sys
import shutil
from contextlib import contextmanager
 
copyblob="""
Copyright (C) University College London, 2007-2012, all rights reserved.

This file is part of HemeLB and is CONFIDENTIAL. You may not work 
with, install, use, duplicate, modify, redistribute or share this
file, or any part thereof, other than as allowed by any agreement
specifically made by you with University College London.
"""
 
@contextmanager
def WorkingDirectory(dst_dir):
    current_dir = os.getcwd()
    os.chdir(dst_dir)
    try:
        yield dst_dir
    finally:
        os.chdir(current_dir)

commentchars={'.cc':'//','.hpp':'//','.h':'//','.py':'#','CMakeLists.txt':'#','.yml':'#'}

def Walk(codeDir):
    for dirpath, dirnames, filenames in os.walk(codeDir):
        if 'build' in dirnames:
            dirnames.remove('build')
            pass
        
        for name in filenames:
            base, ext = os.path.splitext(name)
            if ext=='.in':
              ext=os.path.splitext(base)[1]
            if ext in commentchars.keys():
                yield os.path.normpath(os.path.join(dirpath, name)),ext

def SplitAll(path):
    head, tail = os.path.split(path)
    if head == '':
        ans = []
    else:
        ans = SplitAll(head)
        
    ans.append(tail)
    return ans


if __name__ == '__main__':
    script = os.path.abspath(sys.argv[0])
    scriptDir = os.path.dirname(script)
    codeDir = os.path.normpath(os.path.join(scriptDir, os.pardir))

    with WorkingDirectory(codeDir):
        
        for sourceFile,ext in Walk('.'):
          print sourceFile
          backupFile=sourceFile+'.bak'
          shutil.copy(sourceFile,backupFile)
          to_prepend="\n".join([commentchars[ext]+' '+line for line in copyblob.split('\n')])
          oldtext=open(backupFile).read()
          if oldtext[0:2]=='#!':
            envline=oldtext.split('\n')[0]+'\n'
            oldtext='\n'.join(oldtext.split('\n')[1:-1])
          else:
            envline=''
          with open(sourceFile,'w') as replacement:
            replacement.write(envline)
            replacement.write(to_prepend)
            replacement.write('\n'*2)
            replacement.write(oldtext)
    
        