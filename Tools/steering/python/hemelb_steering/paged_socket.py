import socket

#Model of the way HemeLB steering treats socket-comms.
#On send, is a known-length buffer of values.
#On receive, there is a header, which is used to determine the length of the rest of the message
class PagedSocket(object):
    def __init__(self,address,port,receive_length,additional_receive_length_function=lambda header:0):
        self.address=address
        self.port=port
        self.socket=socket.socket()
        self.socket.connect((address,port))
        self.receive_length=receive_length
        self.additional_receive_length_function=additional_receive_length_function
    def send(self,vals):
        totalsent = 0
        while totalsent < len(vals):
            sent = self.socket.send(vals[totalsent:])
            if sent == 0:
                raise RuntimeError("socket connection broken")
            totalsent += sent
            print totalsent,len(vals)
    def receive(self):
        header=self._receive_fixed_length(self.receive_length)
        additional_length=self.additional_receive_length_function(header)
        return header+self._receive_fixed_length(additional_length)
    def _receive_fixed_length(self,length):
        msg = ''
        while len(msg) < length:
            chunk = self.socket.recv(length-len(msg))
            if chunk == '':
                raise RuntimeError("socket connection broken")
            msg = msg + chunk
        return msg
