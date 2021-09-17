// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <arpa/inet.h>

// Signed and unsigned integers of 32 bits or less
void xdr_deserialise_int32(int32_t* val, const char* src_buf) {
  const uint32_t* in = (const uint32_t*)src_buf;
  *val = ntohl(*in);
}

void xdr_deserialise_uint32(uint32_t* val, const char* src_buf) {
  const uint32_t* in = (const uint32_t*)src_buf;
  *val = ntohl(*in);
}

// Signed and unsigned 64 bit ints
void xdr_deserialise_int64(int64_t* val, const char* src_buf) {
  const uint32_t* in = (const uint32_t*)src_buf;
  uint32_t msw = ntohl(in[0]);
  uint32_t lsw = ntohl(in[1]);
  *val = ((int64_t)msw << 32) ^ (int64_t)lsw;
}
void xdr_deserialise_uint64(uint64_t* val, const char* src_buf) {
  const uint32_t* in = (const uint32_t*)src_buf;
  uint32_t msw = ntohl(in[0]);
  uint32_t lsw = ntohl(in[1]);
  *val = ((uint64_t)msw << 32) ^ (uint64_t)lsw;
}

// 32 bit floats
void xdr_deserialise_float(float* val, const char* src_buf) {
  const uint32_t* in = (const uint32_t*)src_buf;
  uint32_t* out = (uint32_t*)(val);
  *out = ntohl(*in);
}

// 64 bit doubles
void xdr_deserialise_double(double* val, const char* src_buf) {
  uint64_t* out = (uint64_t*)(val);
  xdr_deserialise_uint64(out, src_buf);
}
