# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_io
  PathManager.cc 
  writers/ascii/AsciiFileWriter.cc writers/ascii/AsciiStreamWriter.cc
  writers/xdr/XdrFileReader.cc writers/xdr/XdrFileWriter.cc
  writers/xdr/XdrMemReader.cc writers/xdr/XdrMemWriter.cc
  writers/xdr/XdrReader.cc writers/xdr/XdrWriter.cc
  writers/null/NullWriter.cc
  writers/Writer.cc
  formats/geometry.cc
  xml/XmlAbstractionLayer.cc
  )
target_link_libraries(hemelb_io
  hemelb_util)
