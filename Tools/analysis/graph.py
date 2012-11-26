# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import copy
import itertools
import csv

class Graph(object):
    def __init__(self,config):
        self.select={}
        for prop in ['name','select','curves','dependent','independent']:
            if prop in config:
                setattr(self,prop,config[prop])
            
    def specialise(self,*specialisations):
        """Return a copy of this graph, with the configuration modified by some additional dicts of dicts.
        """
        result=copy.deepcopy(self)
        for specialisation in specialisations:
            if not specialisation:
                continue
            if 'select' in specialisation:
                result.select.update(specialisation['select'])
            for prop in ['curves','dependent','independent']:
                if prop in specialisation:
                    if not hasattr(result,prop):
                        setattr(result,prop,[])
                    getattr(result,prop).extend(specialisation[prop])
        return result
        
    def prepare(self,results):
        self.filtered_results=results.filter(self.select)
        def tuplify_lists(arg):
            return tuple(arg) if hasattr(arg, '__iter__') else arg
        def curve_key(result):
            return tuple([tuplify_lists(result.datum(prop)) for prop in self.curves])
        def sort_key(result):
            return tuple([tuplify_lists(result.datum(prop)) for prop in self.independent])
        self.curve_results={key:list(group) for key,group in itertools.groupby(sorted(self.filtered_results,key=curve_key),curve_key)}
        # for now, support only a single dependent and a single independent variable on each curve
        self.curve_data={
            curve:[ 
                [result.datum(independent) for independent in self.independent]+[result.datum(dependent) for dependent in self.dependent]+list(curve)
                    for result in sorted(results,key=sort_key)
                        ] for curve,results in self.curve_results.iteritems()}
                        
    def write_data(self,csv_file):
        writer=csv.writer(csv_file,lineterminator="\n",quoting=csv.QUOTE_NONNUMERIC)
        fieldnames=self.independent+self.dependent+self.curves
        csv_file.write("#CSV File produced by hemelb reporter\n")
        
        #Write a comment explaining what the groups between which blank lines are inserted are.
        #This makes gnuplot not join-the-dots across these breaks.
        csv_file.write("#Groups are separated by "+str(tuple(self.curves))+"\n#")
        writer.writerow(fieldnames)
        for curve,results in sorted(self.curve_data.iteritems()):
            writer.writerow([])
            writer.writerow([])
            csv_file.write("#"+str(curve)+"\n")
            writer.writerows(results)

        
