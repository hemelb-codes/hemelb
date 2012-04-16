from paged_socket import PagedSocket

class RemoteHemeLB(object):
    def __init__(self,address,port,steering_id):
        self.address=address
        self.port=port
        self.steering_id=steering_id
        self.socket=PagedSocket(address=self.address,port=self.port)
        self.latitude=0
        self.image=None
    def step(self):
        self.image=self.socket.receive() or self.image
        if not self.latitude == 0:
            self.socket.send(self.latitude)
        