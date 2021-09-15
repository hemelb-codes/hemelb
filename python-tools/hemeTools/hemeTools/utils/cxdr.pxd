# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from libc.stdint cimport int32_t, uint32_t

cdef extern from "XdrSerialisation.h":
    void xdr_deserialise_int32(int32_t* val, const char* src_buf)
    void xdr_deserialise_uint32(uint32_t* val, const char* src_buf)
    void xdr_deserialise_float(float* val, const char* src_buf)
    void xdr_deserialise_double(double* val, const char* src_buf)
