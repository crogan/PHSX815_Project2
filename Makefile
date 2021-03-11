# C++ compiler
CXX = g++ 

# necessary compiler flags for using ROOT (root.cern.ch) - remove these if you're not using root
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

# ROOT shared library flags
GLIBS = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))

# some compiler flags
CXXFLAGS = -std=c++11
# ROOT flags
CXXFLAGS += -fPIC $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))

# location of source code
SRCDIR = ./src/

#location of header files
INCLUDEDIR = ./include/

CXXFLAGS += -I$(INCLUDEDIR)

# location of object files (from compiled library files)
OUTOBJ = ./obj/

CC_FILES := $(wildcard src/*.cc)
HH_FILES := $(wildcard include/*.hh)
OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))

# targets to make
all: PlotProb.x ModelSim.x ModelHypoTest.x

# recipe for building PlotProb.x
PlotProb.x:  $(SRCDIR)PlotProb.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o PlotProb.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch PlotProb.x

# recipe for building ModelSim.x
ModelSim.x:  $(SRCDIR)ModelSim.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o ModelSim.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch ModelSim.x

# recipe for building ModelHypoTest.x
ModelHypoTest.x:  $(SRCDIR)ModelHypoTest.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o ModelHypoTest.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch ModelHypoTest.x

$(OUTOBJ)%.o: src/%.cc include/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean-up target (make clean)
clean:
	rm -f *.x
	rm -rf *.dSYM
	rm -f $(OUTOBJ)*.o
