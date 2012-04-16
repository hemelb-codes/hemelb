import socket

class PagedSocket(object):
    def __init__(self,address,port):
        self.address=address
        self.port=port
        self.socket=socket.socket()
        
        self.socket.connect(address,port)
    def send(self,page):
        self.socket.send(page)
    def receive(self):
        contents=self.socket.recv(1024)
        if contents=='':
            return None
        return contents