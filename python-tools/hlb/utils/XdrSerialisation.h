// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_SERIALISATION_H
#define HEMELB_IO_WRITERS_XDR_SERIALISATION_H

#include <stdint.h>

// Signed and unsigned integers of 32 bits or less
void xdr_deserialise_int32(int32_t* val, const char* src_buf);
void xdr_deserialise_uint32(uint32_t* val, const char* src_buf);

// Signed and unsigned 64 bit ints
void xdr_deserialise_int64(int64_t* val, const char* src_buf);
void xdr_deserialise_uint64(uint64_t* val, const char* src_buf);

// 32 bit floats
void xdr_deserialise_float(float* val, const char* src_buf);

// 64 bit doubles
void xdr_deserialise_double(double* val, const char* src_buf);

#endif
