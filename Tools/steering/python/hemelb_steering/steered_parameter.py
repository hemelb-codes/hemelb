# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

#Jiggery-pokery so that client classes of RemoteHemeLB can write remote.Latitude=50
# We define a python "descriptor protocol" class
# See the descriptors discussion in http://docs.python.org/reference/datamodel.html#new-style-and-classic-classes if you never tried these.

class SteeredParameter(object):
    """
    Descriptor implementation for a HemeLB steered parameter
    """
    def __init__(self, index, name):
        self.name = name
        self.index = index
        
    def initialise_in_instance(self, instance, value):
        if not hasattr(instance, 'properties'):
            instance.properties={}
            instance.property_old_values={}
        if not self.name in instance.properties:
            instance.properties[self.name]=value
            instance.property_old_values[self.name]=value
            
    def changed(self, instance):
        return not self.old_value(instance) == self.value(instance)
        
    def old_value(self, instance):
        return instance.property_old_values[self.name]
        
    def value(self, instance):
        return instance.properties[self.name]
        
    def __get__(self, instance, owner):
        return self.value(instance)
        
    def __set__(self, instance, value):
        instance.property_old_values[self.name] = instance.properties[self.name]
        instance.properties[self.name] = value

    
