# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

if(NOT HAVE_XDRUINTXX_T)
  #Without the standard names for the xdr functions, create aliases for the existing ones.
  add_definitions(-Dxdr_uint16_t=xdr_u_int16_t)
  add_definitions(-Dxdr_uint32_t=xdr_u_int32_t)
  add_definitions(-Dxdr_uint64_t=xdr_u_int64_t)
endif()

add_library(hemelb_io
  PathManager.cc 
  writers/ascii/AsciiFileWriter.cc writers/ascii/AsciiStreamWriter.cc
  writers/xdr/XdrVectorWriter.cc
  writers/xdr/XdrFileReader.cc writers/xdr/XdrFileWriter.cc
  writers/xdr/XdrMemReader.cc writers/xdr/XdrMemWriter.cc
  writers/xdr/XdrWriter.cc
  writers/null/NullWriter.cc
  writers/Writer.cc
  formats/geometry.cc
  xml/XmlAbstractionLayer.cc
  )

target_link_libraries(hemelb_io PUBLIC
  hemelb_util
  )
target_link_libraries(hemelb_io PRIVATE
  TinyXML::TinyXML
  )
