HemelbDir = ../../../../../Code
VtkDir = $(EPDROOT)/include/vtk-5.6

HemelbIoSrc = XdrFileWriter.cc XdrMemWriter.cc XdrWriter.cc Writer.cc

LocalSrc = ConfigWriter.cpp BlockWriter.cpp ConfigGenerator.cpp Domain.cpp Block.cpp

CXX = g++
CXXFLAGS = -Wall -fno-strict-aliasing -fno-common -dynamic -arch i386 -DNDEBUG -g -O3 -arch i386 -Wno-deprecated -I$(HemelbDir) -I$(VtkDir)
LDFLAGS = -bundle -undefined dynamic_lookup -arch i386 -g


HemelbIoObj = $(subst .cc,.o, $(HemelbIoSrc))
LocalObj = $(subst .cpp,.o, $(LocalSrc))

Generation.so : $(LocalObj) $(HemelbIoObj)
	$(CXX) $(LDFLAGS) $^ -o $@

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o : $(HemelbDir)/io/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY : clean
clean :
	-rm Generation.so $(LocalObj) $(HemelbIoObj)
