# numpy, vtk
generationDir = HemeLbSetupTool/Model/Generation
epdInclude = $(EPDROOT)/include
pythonDir = $(epdInclude)/python2.7
vtkDir = $(epdInclude)/vtk-5.6
HemeLbDir = ../../Code
include_dirs = $(vtkDir) $(HemeLbDir)
libraries = 
library_dirs = 
extra_compile_args =
extra_link_args = 

generation_cpp = $(addprefix $(generationDir)/, Neighbours.cpp Block.cpp BlockWriter.cpp ConfigWriter.cpp ConfigGenerator.cpp Domain.cpp Site.cpp Debug.cpp)

hemelb_cpp = $(addprefix $(HemeLbDir)/, D3Q15.cc io/XdrFileWriter.cc io/XdrMemWriter.cc io/XdrWriter.cc io/Writer.cc)

# SWIG wrapper
swig_cpp = $(generationDir)/Wrap.cpp

CXX = g++
CXXFLAGS = -Wall -MMD -fno-strict-aliasing -fno-common -dynamic -arch i386 -DNDEBUG -DHEMELB_CFG_ON_BSD -g -O3 -arch i386  -Wno-deprecated -no-cpp-precomp -Werror
INCLUDES = -isystem$(pythonDir) -isystem$(vtkDir) -iquote$(HemeLbDir)

LDFLAGS = -g -bundle -undefined dynamic_lookup -arch i386 # -L/usr/local/lib -lmpi
objects = $(addprefix build/, $(notdir $(generation_cpp:.cpp=.o) $(hemelb_cpp:.cc=.o) $(swig_cpp:.cpp=.o)))

main = HemeLbSetupTool/Model/_Generation.so
$(main) : $(objects)
	$(CXX) $(LDFLAGS) $^ -o $@

-include $(addprefix build/, $(notdir $(generation_cpp:.cpp=.d) $(swig_cpp:.cpp=.d)))

build = $(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ 
build/%.o : $(generationDir)/%.cpp
	$(build)
build/%.o : $(HemeLbDir)/%.cc
	$(build)
build/%.o : $(HemeLbDir)/io/%.cc
	$(build)

$(swig_cpp) : $(swig_cpp:.cpp=.i) $(generationDir)/ConfigGenerator.h $(generationDir)/Iolet.h
	swig -c++ -python -o $@ -outdir HemeLbSetupTool/Model $<

.PHONY : test
test : $(main)
	env PYTHONPATH=$(PWD) scripts/setuptool --profile ~/tmp/setuptooltest/cyl_048.pro

debug : $(main)
	env PYTHONPATH=$(PWD) gdb -x break.gdb -x run.gdb python

.PHONY : clean
clean :
	-rm $(main) $(objects)

.PHONY : cleanall
cleanall : clean
	-rm $(swig_cpp) HemeLbSetupTool/Model/Generation.py
