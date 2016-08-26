#
# the makefile that will be used by the 'make' utility
#

#a variable which holds the name of the executable file to be created
name := EdepNSRL2014

#a variable that puts this into a variable that is used in the binmake.gmk file
#(see below)
G4TARGET := $(name)

#variable, what does it do?
G4EXLIB := true

#define the target 'all' to be 'phony'; ie, isn't a file
.PHONY: all

#target is all, dependencies are 'lib' and 'bin', ie, will check to see if anything was changed in those sub-directories
all: lib bin

#now, include the "geant4" makefile, builds the source code in a standard way
include $(G4INSTALL)/config/binmake.gmk

#Another G4 forum suggestion, analysis forum, posted 5/4/4
#ifdef ROOTSYS
#ifdef G4UI_USE_ROOT
#NB: The -DG4ANALYSIS_USE_ROOT part below is CRITICAL!
CPPFLAGS += -DG4ANALYSIS_USE_ROOT $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS := $(filter-out -lNew, $(ROOTLIBS))
ROOTLIBS := $(filter-out -lThread, $(ROOTLIBS))
ROOTLIBS := $(filter-out -lpThread, $(ROOTLIBS))
LDLIBS += $(ROOTLIBS)
#try some extra stuff, seems to only require 4th one to make...
#EXTRALIBS += $(shell$(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
#INCFLAGS += -I$(ROOTSYS)/include/root
#LDFLAGS += -L$(ROOTSYS)/lib/root
LDLIBS_PREFINAL += $(shell $(ROOTSYS)/bin/root-config --libs)
#endif
#endif
