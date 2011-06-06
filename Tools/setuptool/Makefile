# numpy, vtk
epdInclude = $(EPDROOT)/include
pythonDir = $(epdInclude)/python2.7
vtkDir = $(epdInclude)/vtk-5.6
HemeLbDir = ../../Code
include_dirs = $(vtkDir) $(HemeLbDir)
libraries = 
library_dirs = 
extra_compile_args =
extra_link_args = 

generation_cpp = $(addprefix HemeLbSetupTool/Model/Generation/, Block.cpp BlockWriter.cpp ConfigWriter.cpp ConfigGenerator.cpp Domain.cpp Site.cpp)

hemelb_cpp = $(addprefix $(HemeLbDir)/, io/XdrFileWriter.cc io/XdrMemWriter.cc io/XdrWriter.cc io/Writer.cc)

# SWIG wrapper
swig_cpp = HemeLbSetupTool/Model/Generation/Wrap.cpp
cmd = 'swig -c++ -python -o HemeLbSetupTool/Model/Generation/Wrap.cpp -outdir HemeLbSetupTool/Model HemeLbSetupTool/Model/Generation/Wrap.i'

CXX = g++
CXXFLAGS = -MD -fno-strict-aliasing -fno-common -dynamic -arch i386 -DNDEBUG -DHEMELB_CFG_ON_BSD -g -O0 -arch i386  -Wno-deprecated -no-cpp-precomp
INCLUDES = -isystem$(pythonDir) -isystem$(vtkDir) -iquote$(HemeLbDir)

LDFLAGS = -g -bundle -undefined dynamic_lookup -arch i386 -L/usr/local/lib -lmpi
objects = $(addprefix build/, $(notdir $(generation_cpp:.cpp=.o) $(hemelb_cpp:.cc=.o) $(swig_cpp:.cpp=.o)))

main = HemeLbSetupTool/Model/_Generation.so
$(main) : $(objects)
	$(CXX) $(LDFLAGS) $^ -o $@

-include $(addprefix build/, $(notdir $(generation_cpp:.cpp=.d) $(swig_cpp:.cpp=.d)))

build = $(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ 
build/%.o : HemeLbSetupTool/Model/Generation/%.cpp
	$(build)
build/%.o : $(HemeLbDir)/io/%.cc
	$(build)

$(swig_cpp) : $(swig_cpp:.cpp=.i)
	swig -c++ -python -o $@ -outdir HemeLbSetupTool/Model $<

.PHONY : test
test : $(main)
	env PYTHONPATH=$(PWD) python -m pdb scripts/setuptool --profile ~/tmp/setuptooltest/cyl_048.pro

debug : $(main)
	env PYTHONPATH=$(PWD) gdb -x break.gdb -x run.gdb python

.PHONY : clean
clean :
	-rm $(main) $(objects)

.PHONY : cleanall
cleanall : clean
	-rm $(swig_cpp) HemeLbSetupTool/Model/Generation.py
