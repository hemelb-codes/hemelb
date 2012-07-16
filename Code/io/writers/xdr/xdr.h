// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_IO_WRITERS_XDR_XDR_H
#define HEMELB_IO_WRITERS_XDR_XDR_H

/*
 * This include file should be used to include XDR rather than the standard
 * rpc/xdr.h so we can ensure cross platform compatibility.
 *
 */

#include <rpc/types.h>
#include <rpc/xdr.h>

/* The XDR implementation on BSD based machines uses
 * xdr_u_int??_t() to pack/unpack these types, while Linux
 * based machines have xdr_uint??_t(). Use the config macro to
 * create aliases on BSD.
 */
#ifdef HEMELB_CFG_ON_BSD
#define xdr_uint16_t xdr_u_int16_t
#define xdr_uint32_t xdr_u_int32_t
#define xdr_uint64_t xdr_u_int64_t
#endif //  HEMELB_CFG_ON_BSD


#endif // HEMELB_IO_WRITERS_XDR_XDR_H
