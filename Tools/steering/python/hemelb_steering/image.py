class Image(object):
    class Pixel(object):
        def __init__(self,unpacker):
            self.index=unpacker.unpack_int()
            self.x=self.index&((1<<17)-1)
            self.y=self.index>>16
            self.velocity=[]
            self.stress=[]
            self.pressure=[]
            self.stress2=[]
            self.components=[self.velocity,self.stress,self.pressure,self.stress2]
            for _ in ['red','green','blue']:
                blob=unpacker.unpack_fopaque(4)
                for index,component in enumerate(self.components):
                    component.append(blob[index])
    def __init__(self,pixel_count,unpacker):
        pixels=[]
        for _ in xrange(pixel_count):
            pixels.append(Image.Pixel(unpacker))
    bytes_per_pixel=2*2+3*4 #each of three colors with four sub-images per color and two two-byte coordinates
