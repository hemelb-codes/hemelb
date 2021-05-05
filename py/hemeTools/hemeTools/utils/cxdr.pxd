
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cdef extern from "rpc/types.h":
    ctypedef bint bool_t
    pass

cdef extern from "rpc/xdr.h":
    cdef enum xdr_op:
        XDR_ENCODE=0
        XDR_DECODE=1
        XDR_FREE=2
    ctypedef struct XDR:
        pass
    void xdrmem_create(XDR* x, char* buf, unsigned int len, xdr_op op)
    bool_t xdr_double(XDR *x, double* v)
    bool_t xdr_float(XDR *x, float* v)
    bool_t xdr_int(XDR *x, int* v)
    bool_t xdr_u_int(XDR *x, unsigned int* v)
    # These aren't actually functions but macros that call f pointers in teh XDR struct
    # xdr_getpos(XDR *x)
    # xdr_setpos(XDR *x, iPosition)
    void xdr_destroy(XDR* c)
    # xdr_destroy(XDR *x)
