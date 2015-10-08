#!/usr/bin/env python

## Program:   PyNS
## Module:    Profiling.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

'''
This module describes the run time performance of pyNS, providing a variety of statistics.
This module use Hotshot, which is a replacement for the existing default profile module.
As it is written mostly in C, it should result in a much smaller performance impact than the existing profile module.
This module automatically prints a simple profiling report, sorted by time and calls.
ncalls: for the number of calls,
tottime: for the total time spent in the given function (and excluding time made in calls to sub-functions),
percall :is the quotient of tottime divided by ncalls
cumtime :is the total time spent in this and all subfunctions (from invocation till exit). This figure is accurate even for recursive functions.
percall: is the quotient of cumtime divided by primitive calls
filename:lineno(function): provides the respective data of each function 
'''

from pyNS import runSimulation

def main():  
    
    runSimulation()
   

import hotshot, hotshot.stats
prof = hotshot.Profile("pyNS.profile")
command = """main()"""
prof.runctx( command, globals(), locals())
prof.close()
stats = hotshot.stats.load("pyNS.profile")
stats.strip_dirs()
stats.sort_stats('cumulative')
stats.print_stats(30)