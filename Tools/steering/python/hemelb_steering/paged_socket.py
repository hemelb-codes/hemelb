import socket

#Model of the way HemeLB steering treats socket-comms.
#On send, is a known-length buffer of float values.
class PagedSocket(object):
    def __init__(self,address,port):
        self.address=address
        self.port=port
        self.socket=socket.socket()
        self.socket.connect(address,port)
    def send(self,vals):
        totalsent = 0
        while totalsent < len(vals):
            sent = self.socket.send(str(bytearray(vals[totalsent:])))
            if sent == 0:
                raise RuntimeError("socket connection broken")
            totalsent += sent
    def receive(self):
        contents=self.socket.recv(1024)
        if contents=='':
            return None
        return contents
