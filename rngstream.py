#load RngStream dll library
from ctypes import *
import os

#def corresponding structure
class RngStream_InfoState(Structure):
    _fields_ =[('Cg',c_double),
               ('Bg',c_double),
               ('Ig',c_double),
               ('Anti',c_int),
               ('IncPrec',c_int),
               ('name',c_char_p)]

class RngStream():
    def __init__(self):
        self.lib = CDLL(os.path.dirname(__file__)+"/RngStream.so")

        #initialize with pointer
        self.lib.RngStream_CreateStream.restype = POINTER(RngStream_InfoState)
        # g = RngStream.RngStream_CreateStream()

        #uniform [0,1]
        self.lib.RngStream_RandU01.restype = c_double

        #RandInt [20,30]
        self.lib.RngStream_RandInt.restype = c_int
        # i = c_int(20)
        # j = c_int(30)

        self.lib.RngStream_SetPackageSeed.restype = c_int
        # seed = c_ulong * 6
        # err = RngStream.RngStream_SetPackageSeed(seed)

    def set_package_seed(self,s):
        return self.lib.RngStream_SetPackageSeed(s)
    def create_stream(self):
        return self.lib.RngStream_CreateStream(None)
    def float_uniform_01(self,g):
        return self.lib.RngStream_RandU01(g)
    def int_uniform(self,g,l,h):
        return self.lib.RngStream_RandInt(g, l, h)
    def get_state(self, g, s):
        self.lib.RngStream_GetState(g, s)
