from paged_socket import PagedSocket
from steered_parameter import SteeredParameter
from image import Image
from config import config
import xdrlib

class RemoteHemeLB(object):
    # Metaprogram property access       
    def __init__(self,address,port,steering_id):
        self.address=address
        self.port=port
        self.steering_id=steering_id
        self.socket=PagedSocket(address=self.address,
            port=self.port,
            receive_length=3*RemoteHemeLB.xdr_int_bytes,
            additional_receive_length_function=RemoteHemeLB._calculate_receive_length)
        self.latitude=0
        self.image=None
        for steered_parameter in self.steered_parameters:
            steered_parameter.initialise_in_instance(self,config['steered_parameter_defaults'][steered_parameter.name])

    xdr_int_bytes=4
    xdr_double_bytes=8

    def step(self):
        if any([par.changed(self) for par in self.steered_parameters]):
            self.socket.send([par.value(self) for par in self.steered_parameters])
        self.receive()
        self.image=self.socket.receive() or self.image
    
    @staticmethod
    def _calculate_receive_length(header):
        unpacker=xdrlib.Unpacker(header)
        print header
        width=unpacker.unpack_int()
        height=unpacker.unpack_int()
        frame=unpacker.unpack_int()
        return frame+3*RemoteHemeLB.xdr_double_bytes+3*RemoteHemeLB.xdr_int_bytes
        
    def receive(self):
        page=self.socket.receive()
        unpacker=xdrlib.Unpacker(page)
        self.width=unpacker.unpack_int()
        self.height=unpacker.unpack_int()
        self.frame=unpacker.unpack_int()
        self.image=Image(self.frame/Image.bytes_per_pixel,unpacker)
        self.step=unpacker.unpack_int()
        self.time=unpacker.unpack_double()
        unpacker.unpack_int() # throw away cycle
        self.inlets=unpacker.unpack_int()
        self.mouse_pressure=unpacker.unpack_double()
        self.mouse_stress=unpacker.unpack_double()

#Metaprogram property access
RemoteHemeLB.steered_parameters=[]        
for index,name in enumerate(config['steered_parameters']):
    steered_parameter=SteeredParameter(index,name)
    setattr(RemoteHemeLB,name,steered_parameter)
    RemoteHemeLB.steered_parameters.append(steered_parameter)