# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#

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
