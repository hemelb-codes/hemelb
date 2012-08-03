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
from contextlib import contextmanager
 
@contextmanager
def WorkingDirectory(dst_dir):
    current_dir = os.getcwd()
    os.chdir(dst_dir)
    try:
        yield dst_dir
    finally:
        os.chdir(current_dir)

# Want to match: start of line; '#include'; one or more space;
# double quote; anything (but non-greedily); double quote.
nonSystemIncludeFinder = re.compile(r'\A#include\s+"(.+?)"')
# As above, but with <...> 
systemIncludeFinder = re.compile(r'\A#include\s+<(.+?)>')

def IncludedFileGenerator(filename, includeFinder):
    """Given a filename, return a generator that will yield the line
    numbers and paths included by include statements in the
    file.
    """    
    codefile = file(filename)
    for iLine, line in enumerate(codefile):
        match = includeFinder.match(line)
        if match:
            yield (iLine + 1), match.group(1)

FileNonSystemIncludedFileGenerator = lambda filename: IncludedFileGenerator(filename, nonSystemIncludeFinder)
FileSystemIncludedFileGenerator = lambda filename: IncludedFileGenerator(filename, systemIncludeFinder)

extensions = set(('.cc', '.hpp', '.h'))
def Walk(codeDir):
    for dirpath, dirnames, filenames in os.walk(codeDir):
        if 'build' in dirnames:
            dirnames.remove('build')
            pass
        
        for name in filenames:
            base, ext = os.path.splitext(name)
            if ext in extensions or (ext == '.in' and os.path.splitext(base)[1] in extensions):
                yield os.path.normpath(os.path.join(dirpath, name))

def CheckIncludePaths(sourceFile):
    errors = False
    for lineNumber, include in FileNonSystemIncludedFileGenerator(sourceFile):
        if include in ignoredIncludes:
            continue
        
        if not os.path.exists(include) and not os.path.exists(include + '.in'):
            sys.stderr.write(
                '{file}:{line} Bad include path "{dodgy}"\n'.format(file=sourceFile,
                                                                    line=lineNumber,
                                                                    dodgy=include)
                )
            errors = True
    return errors

cLibraryHeaders = set((
    'assert.h',
    'ctype.h',
    'errno.h',
    'float.h',
    'iso646.h',
    'limits.h',
    'locale.h',
    'math.h',
    'setjmp.h',
    'signal.h',
    'stdarg.h',
    'stddef.h',
    'stdio.h',
    'stdlib.h',
    'string.h',
    'time.h',
    ))

def CheckSystemIncludePaths(sourceFile):
    errors = False
    for lineNumber, include in FileSystemIncludedFileGenerator(sourceFile):
        if include in cLibraryHeaders:
            sys.stderr.write(
                '{file}:{line} Use of deprecated C library header <{include}>\n'.format(
                    file=sourceFile,
                    line=lineNumber,
                    include=include
                    )
                )
            errors = True
    return errors

def ignoreCopyright(f):
  line=f.readline()
  while line.strip()=='' or line[0:2]=='//':
    line=f.readline()
  return line

def GetGuardLines(filename):
    f = file(filename)
    lines=[ignoreCopyright(f)]
    lines.append(f.readline())
    
    for line in f:
        continue
    lines.append(line)
    return lines

def SplitAll(path):
    head, tail = os.path.split(path)
    if head == '':
        ans = []
    else:
        ans = SplitAll(head)
        
    ans.append(tail)
    return ans

def CheckGuardErrors(sourceFile):
    lines= GetGuardLines(sourceFile)
    
    parts = SplitAll(sourceFile)
    parts[-1] = parts[-1].replace('.', '_')
    define = ('HEMELB_' + '_'.join(parts)).upper()

    error = False
    
    line0 = '#ifndef {define}\n'.format(define=define)
    if lines[0] != line0:
        sys.stderr.write(
                '{file}:1 Bad include guard; must be {required!r} but was {actual!r}\n'.format(file=sourceFile,
                                                                            required=line0,actual=lines[0])
            )
        error = True

    line1 = '#define {define}\n'.format(define=define)
    if lines[1] != line1:
        sys.stderr.write(
                '{file}:2 Bad include guard; must be {required!r} but was {actual!r}\n'.format(file=sourceFile,
                                                                            required=line1,actual=lines[1])
            )
        error = True
        
    return error

if __name__ == '__main__':
    script = os.path.abspath(sys.argv[0])
    scriptDir = os.path.dirname(script)
    codeDir = os.path.normpath(os.path.join(scriptDir, os.pardir, 'Code'))

    ignoredIncludes = set((
        'MPWide.h',
	'tinyxml.h',
        'parmetis.h',
        'ctemplate/template.h',
        ))

    errors = False
    with WorkingDirectory(codeDir):
        
        for sourceFile in Walk('.'):
            includeErrors = CheckIncludePaths(sourceFile)
            sysIncludeErrors = CheckSystemIncludePaths(sourceFile)
            
            base, ext = os.path.splitext(sourceFile)
            if ext == '.h' or ext == '.hpp' or ext == '.in':
                guardErrors = CheckGuardErrors(sourceFile)
            else:
                guardErrors = False
            errors = errors or (includeErrors or sysIncludeErrors or guardErrors)

    if errors:
        raise SystemExit(1)
    
        
