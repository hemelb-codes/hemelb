from paged_socket import PagedSocket
from steered_parameter import SteeredParameter
from config import config

class RemoteHemeLB(object):
    # Metaprogram property access       
    def __init__(self,address,port,steering_id):
        self.address=address
        self.port=port
        self.steering_id=steering_id
        self.socket=PagedSocket(address=self.address,port=self.port)
        self.latitude=0
        self.image=None
        for steered_parameter in self.steered_parameters:
            steered_parameter.initialise_in_instance(self)

    def step(self):
        self.image=self.socket.receive() or self.image
        if any([par.changed(self) for par in self.steered_parameters]):
            self.socket.send([par.value(self) for par in self.steered_parameters])

#Metaprogram property access
RemoteHemeLB.steered_parameters=[]        
for index,name in enumerate(config['steered_parameters']):
    steered_parameter=SteeredParameter(index,name)
    setattr(RemoteHemeLB,name,steered_parameter)
    RemoteHemeLB.steered_parameters.append(steered_parameter)