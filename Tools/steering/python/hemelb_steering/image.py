# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 


try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
    
import numpy as N


"""
Class describing the image data streamed by HemeLB
"""
class Image(object):

    def __init__(self, width, height, pixel_count, unpacker):
        self.pixels = []
        self.width = width
        self.height = height
        self.given_pixel_count = pixel_count
        self.full_pixel_count = width * height # Might be sparse
        self.data = N.frombuffer(unpacker.unpack_fopaque(pixel_count * Image.bytes_per_pixel), dtype=Image.pixel)
    
        
    def pil(self, component='velocity'):
        """ 
        Transform the data to python image library format
        """
        try:
            from PIL import Image as PILImage
        except ImportError:
            import Image as PILImage
        pil_string_data = bytearray([255] * 3 * self.full_pixel_count)
        fields_wanted = ["%s_%s" % (component,color) for color in Image.colors]
        for pixel in self.data:
            offset = pixel['y'] * self.width + pixel['x']
            pil_string_data[3 * offset : 3 * offset + 3] = [pixel[field] for field in fields_wanted]
        return PILImage.fromstring("RGB", (self.width,self.height), str(pil_string_data))
        
    # Constants useful for parsing the image
    subimages=['velocity', 'stress', 'pressure', 'stress2']
    colors=['red', 'green', 'blue']
    fields=["%s_%s" % (subimage, color) for subimage in subimages for color in colors ]
    bytes_per_pixel=2*2 + 3*4 #each of three colors with four sub-images per color and two two-byte coordinates
    pixel=N.dtype({'names': ['x', 'y'] + fields, 'formats': [N.dtype('>H')] * 2 + [N.uint8] * len(subimages) * len(colors)})
