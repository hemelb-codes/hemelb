#!/usr/bin/env python

## Program:   PyNS
## Module:    export.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from json import load
import sys, argparse

def export(fileName):
    '''
    Retrieving information about time, pressure and flow of given mesh
    from json file and writing them in .txt files.
    time [sec];pressure[Pa];flow[m3/s]
    '''
    f = open(fileName)
    infos = load(f)
    f.close()
    name = fileName.split('.')[0]
    text_file = open(name+'.txt', "w")
    text_file.write('time[s];pressure[Pa];flow[m3/s]\n')
    data = infos['items'][0]
    time = []
    flow = []
    pressure = []
    flow_data = data['flow']
    pressure_data = data['pressure']
    for values in pressure_data:
        pressure.append(values[1]*133.322)
    for values in flow_data:
        time.append(values[0])
        flow.append(values[1]/6e7)
    i = 0
    while i <len(time):
        text_file.write(str('{:.4e}'.format(time[i]))+';'+str('{:.4e}'.format(pressure[i]))+';'+str('{:.4e}'.format(flow[i])+'\n'))
        i+=1
    text_file.close()
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Export module for pyNS')
    parser.add_argument("-f", "--file", action="store", dest='fileName', default=None, help="Specify pyNS json input file path")
    args = parser.parse_args()
    
    fileName = args.fileName
    if fileName is None:
        sys.exit("Please specify json input file path")
    export(fileName)