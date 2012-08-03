# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from paged_socket import PagedSocket
from steered_parameter import SteeredParameter
from image import Image
from config import config
import xdrlib

class RemoteHemeLB(object):
    """
    Represent a remote HemeLB as a python object
    Steerable parameters become properties of the object:
    e.g. myheme.Latitude=50
    """
   
    def __init__(self, address, port, steering_id):
        self.address = address
        self.port = port
        self.steering_id = steering_id
        self.socket = PagedSocket(address=self.address,
            port=self.port,
            receive_length= 3 * RemoteHemeLB.xdr_int_bytes,
            additional_receive_length_function=RemoteHemeLB._calculate_receive_length)
        self.latitude = 0
        self.image = None
        for steered_parameter in self.steered_parameters:
            steered_parameter.initialise_in_instance(self, config['steered_parameter_defaults'][steered_parameter.name])

    xdr_int_bytes = 4
    xdr_double_bytes = 8

    def step(self):
        self.receive()
        self.send()
        
    def send(self):
        if any([par.changed(self) for par in self.steered_parameters]):
            vals = xdrlib.Packer()
            for par in self.steered_parameters:
                vals.pack_float(float(par.value(self)))
            self.socket.send(vals.get_buffer())
    
    @staticmethod
    def _calculate_receive_length(header):
        unpacker = xdrlib.Unpacker(header)
        width = unpacker.unpack_int()
        height = unpacker.unpack_int()
        frame = unpacker.unpack_int()
        return frame + 3*RemoteHemeLB.xdr_double_bytes + 3*RemoteHemeLB.xdr_int_bytes
        
    def receive(self):
        page = self.socket.receive()
        unpacker = xdrlib.Unpacker(page)
        self.width = unpacker.unpack_int()
        self.height = unpacker.unpack_int()
        self.frame = unpacker.unpack_int()
        self.image = Image(self.width, self.height, self.frame / Image.bytes_per_pixel, unpacker)
        self.time_step = unpacker.unpack_int()
        self.time = unpacker.unpack_double()
        unpacker.unpack_int() # throw away cycle
        self.inlets = unpacker.unpack_int()
        self.mouse_pressure = unpacker.unpack_double()
        self.mouse_stress = unpacker.unpack_double()
        
    def __str__(self):
        return "HemeLB at %s: Step %s" % (self.address, self.time_step)

# Metaprogram property access
# We use metaprogramming techniques to define the properties which can be steered,
# based on the contents of config.yml
# We use the SteeredParameter python "descriptor" class to define these attributes.
# See the descriptors discussion in http://docs.python.org/reference/datamodel.html#new-style-and-classic-classes if you never tried these.
RemoteHemeLB.steered_parameters=[]        

for index,name in enumerate(config['steered_parameters']):
    steered_parameter = SteeredParameter(index,name)
    setattr(RemoteHemeLB, name, steered_parameter)
    RemoteHemeLB.steered_parameters.append(steered_parameter)
