add_library(hemelb_io
	PathManager.cc 
	writers/ascii/AsciiFileWriter.cc writers/ascii/AsciiStreamWriter.cc
	writers/xdr/XdrFileReader.cc writers/xdr/XdrFileWriter.cc
	writers/xdr/XdrMemReader.cc writers/xdr/XdrMemWriter.cc
	writers/xdr/XdrReader.cc writers/xdr/XdrWriter.cc
	writers/Writer.cc
	formats/geometry.cc
	xml/XmlAbstractionLayer.cc
	)
