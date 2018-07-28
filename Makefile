#
#      --> all main programs in the root directory
#      --> a directory named $(INCDIR) containing header files
#      --> a directory named $(SRCDIR) containing source files
#
##############################################

# Get the build date
######################
BUILD_DATE := "\"`date`\""
DATEFLAGS := -DBUILD_DATE=$(BUILD_DATE)

# Directories settings
########################
INCDIR := include
SRCDIR := src
TOOLDIR := toolkits
BUILDDIR := build
BINDIR := bin
# List of sub-directories contained into INCDIR and SRCDIR source directories
DIRS = algorithm analytic_simulator datafile dynamic image management optimizer projector scanner

# Extensions settings
#######################
SRCEXT := .cc
HDREXT := .hh
DEPEXT := .d
OBJEXT := .o
SRCBUILDPREF  := __src
MAINBUILDPREF := __main
TOOLBUILDPREF := __tool
# We create a sub-build-directory for each sub-source-directory
BUILDDIRS = $(foreach DIR, $(DIRS), $(BUILDDIR)/$(SRCBUILDPREF)_$(DIR))

# Input and output files
##########################
# Includes
INCLUDE = $(foreach DIR, $(DIRS), -I$(INCDIR)/$(DIR))
# Class sources
SRC    = $(wildcard $(SRCDIR)/*/*$(SRCEXT))
HDR    = $(wildcard $(SRCDIR)/*/*$(HDREXT))
# Main sources
MAIN   = $(wildcard *$(SRCEXT))
# Tool sources
TOOL   = $(wildcard $(TOOLDIR)/*$(SRCEXT))
# Executables (main programs)
EXE_MAIN    = $(patsubst %$(SRCEXT),$(BINDIR)/%,$(MAIN))
# Executables (tool programs)
EXE_TOOL    = $(patsubst $(TOOLDIR)/%$(SRCEXT),$(BINDIR)/%,$(TOOL))
# Class objects
OBJ    = $(foreach DIR, $(DIRS), \
           $(patsubst $(SRCDIR)/$(DIR)/%$(SRCEXT),$(BUILDDIR)/$(SRCBUILDPREF)_$(DIR)/%$(OBJEXT),$(wildcard $(SRCDIR)/$(DIR)/*$(SRCEXT)) \
            ) \
          )
# Toolkit objects
OBJ   += $(patsubst $(TOOLDIR)/%$(SRCEXT),$(BUILDDIR)/$(TOOLBUILDPREF)_%$(OBJEXT),$(TOOL))
# Main objects
OBJ   += $(patsubst %$(SRCEXT),$(BUILDDIR)/$(MAINBUILDPREF)_%$(OBJEXT),$(MAIN))
# Class dependencies
DEP    = $(foreach DIR, $(DIRS), \
           $(patsubst $(SRCDIR)/$(DIR)/%$(SRCEXT),$(BUILDDIR)/$(SRCBUILDPREF)_$(DIR)/%$(DEPEXT),$(wildcard $(SRCDIR)/$(DIR)/*$(SRCEXT)) \
            ) \
          )
# Toolkit dependencies
DEP   += $(patsubst $(TOOLDIR)/%$(SRCEXT),$(BUILDDIR)/$(TOOLBUILDPREF)_%$(DEPEXT),$(TOOL))
# Main dependencies
DEP   += $(patsubst %$(SRCEXT),$(BUILDDIR)/$(MAINBUILDPREF)_%$(DEPEXT),$(MAIN))

# Architecture
#########################
ARCH  := $(shell uname -p)

# Main variables
###############################
ifndef CXX
CXX = g++
endif
# C++11 is used
CXXFLAGS = $(INCLUDE) -std=c++11
# We want the cleanest code possible
CXXFLAGS += -Wall
# To avoid warnings if FLTNB is set to long double (warnings are thrown when a cast from long double to double occurs)
CXXFLAGS += -Wno-narrowing
# To avoid warnings that says some functions are unused (mandatory when being generic)
CXXFLAGS += -Wno-unused-function
# To avoid warnings that says some variables are unused (mandatory when being generic)
CXXFLAGS += -Wno-unused-parameter
# Linker
LDFLAGS = $(INCLUDE)
LINKER = $(CXX)
# For 64-bits architectures
ifeq ($(ARCH),x86_64)
CXXFLAGS += -D_64 -m64
LDFLAGS += -D_64
endif

# Static link of CASToR
#########################
ifeq ($(CASTOR_STATIC), 1)
  LDFLAGS += -static-libgcc -static-libstdc++ -static -lpthread
endif

# Cross-compilation with MINGW 32/64 bits
###########################################
ifeq ($(CASTOR_MINGW), 32)
  # Change compiler
  CXX = i686-w64-mingw32-g++
  # Add static linking flags to the linker
  LDFLAGS += -static-libgcc -static-libstdc++ -static -lpthread -m32
  # Specify to CASToR that we cross-compile
  CXXFLAGS += -DCASTOR_USE_MINGW -m32
  # We also need the CASTOR_CONFIG variable to hard-link it
  ifdef CASTOR_CONFIG
  CXXFLAGS += -DCASTOR_CONFIG=$(CASTOR_CONFIG)
  endif
endif
ifeq ($(CASTOR_MINGW), 64)
  # Change compiler
  CXX = x86_64-w64-mingw32-g++
  # Add static linking flags to the linker
  LDFLAGS += -static-libgcc -static-libstdc++ -static -lpthread
  # Specify to CASToR that we cross-compile
  CXXFLAGS += -DCASTOR_USE_MINGW
  # We also need the CASTOR_CONFIG variable to hard-link it
  ifdef CASTOR_CONFIG
  CXXFLAGS += -DCASTOR_CONFIG=$(CASTOR_CONFIG)
  endif
endif

# CPU optimization options
############################
ifeq ($(CASTOR_DEBUG), 1)
CXXFLAGS += -g3 -fno-inline -O0
BINDIR := $(BINDIR)_debug
CXXFLAGS += -DCASTOR_DEBUG
else
  CXXFLAGS += -O3
  # OpenMP multi-threading
  ifeq ($(CASTOR_OMP), 1)
    CXXFLAGS += -fopenmp -DCASTOR_OMP
    LDFLAGS += -fopenmp
  else
    CXXFLAGS += -Wno-unknown-pragmas
  endif
  # MPI multi-computer
  ifeq ($(CASTOR_MPI), 1)
    CXX = mpic++
    CXXFLAGS += -DCASTOR_MPI
  endif
  # SIMD auto
  ifeq ($(CASTOR_SIMD), 1)
    CXXFLAGS += -ftree-vectorize -msse2 -ffast-math -fassociative-math
    CXXFLAGS += -ftree-vectorizer-verbose=2
  endif
endif

# Verbosity of compilation and execution
##########################################
ifeq ($(CASTOR_VERBOSE), 1)
#  ifeq ($(CASTOR_SIMD), 1)
#  CXXFLAGS += -ftree-vectorizer-verbose=2
#  endif
  CXXFLAGS += -Wall -DCASTOR_VERBOSE
endif

# ROOT support for datafile conversion
##########################################
ifeq ($(CASTOR_ROOT), 1)
  CXXFLAGS += $(shell root-config --cflags) -DCASTOR_ROOT
  LDLIBS   = $(shell root-config --libs)
endif

#4D Deformation Ledesma
############################
ifeq ($(CASTOR_DEFORMATION_LEDESMA), 1)
  LDLIBS += -L./include/image/lib -lElasticRegistration
  INCLUDE += -I./include/image/lib
  CXXFLAGS += -DCASTOR_DEFORMATION_LEDESMA
endif

#Deformation Elastix // TODO
############################
#ifeq ($(CASTOR_DEFORMATION_ELASTIX), 1)
#  LDLIBS += -L./include/image/lib -lElasticRegistration
#  INCLUDE += -I./include/image/lib
#  CXXFLAGS += -DCASTOR_DEFORMATION_ELASTIX
#endif

#######################################################################
##                               Rules                               ##
#######################################################################

################
##  all Part  ##
################
all : $(BUILDDIR)/exe.last

################
##  EXE Part  ##
################
# Linking $(MAIN)
$(BINDIR)/% : $(BUILDDIR)/$(MAINBUILDPREF)_%$(OBJEXT) $(BUILDDIR)/obj.last
	@if [ ! -d $(BINDIR) ] ; then mkdir -p $(BINDIR) ; fi
	@echo "Linking" $(subst $(BINDIR)/,,$@) "..."
	@$(LINKER) $(LDFLAGS) -o $@ $(BUILDDIR)/$(SRCBUILDPREF)_*/*$(OBJEXT) $< $(LDLIBS)
# Linking $(TOOL)
$(BINDIR)/% : $(BUILDDIR)/$(TOOLBUILDPREF)_%$(OBJEXT) $(BUILDDIR)/obj.last
	@if [ ! -d $(BINDIR) ] ; then mkdir -p $(BINDIR) ; fi
	@echo "Linking" $(subst $(BINDIR)/,,$@) "..."
	@$(LINKER) $(LDFLAGS) -o $@ $(BUILDDIR)/$(SRCBUILDPREF)_*/*$(OBJEXT) $< $(LDLIBS)

################
##  OBJ Part  ##
################
# Compiling $(SRCDIR) for CPU
$(BUILDDIR)/$(SRCBUILDPREF)_%.o : $(SRCDIR)/%$(SRCEXT) $(BUILDDIR)/$(SRCBUILDPREF)_%$(DEPEXT)
	@if [ ! -d $(BUILDDIR) ] ; then mkdir -p $(BUILDDIR) ; fi
	@echo "Compiling" $(subst $(SRCDIR)/,,$<) "..."
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -c -o $@ $<

# Compiling $(TOOLDIR)
$(BUILDDIR)/$(TOOLBUILDPREF)_%.o : $(TOOLDIR)/%$(SRCEXT) $(BUILDDIR)/$(TOOLBUILDPREF)_%$(DEPEXT)
	@if [ ! -d $(BUILDDIR) ] ; then mkdir -p $(BUILDDIR) ; fi
	@echo "Compiling" $< "..."
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) $(DATEFLAGS) -c -o $@ $<

# Compiling $(SRCEXT)
$(BUILDDIR)/$(MAINBUILDPREF)_%.o : %$(SRCEXT) $(BUILDDIR)/$(MAINBUILDPREF)_%$(DEPEXT)
	@if [ ! -d $(BUILDDIR) ] ; then mkdir -p $(BUILDDIR) ; fi
	@echo "Compiling" $< "..."
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) $(DATEFLAGS) -c -o $@ $<

######################
##  .PHONY targets  ##
######################
.PHONY : all obj clean

obj : $(BUILDDIR)/obj.last

$(BUILDDIR)/obj.last : $(OBJ)
	@touch $@

$(BUILDDIR)/exe.last : $(EXE_MAIN) $(EXE_TOOL)
	@echo "All executables are in" $(BINDIR)
	@touch $@

################
##  DEP Part  ##
################
# Making dependencies for $(MAIN)
$(BUILDDIR)/$(MAINBUILDPREF)_%$(DEPEXT) : %$(SRCEXT)
	@mkdir -p $(BUILDDIRS)
	@echo "Making dependencies for" $< "..."
	@set -e; rm -f $@; \
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
ifneq ($(DEP),)
ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif
endif

# Making dependencies for $(TOOL)
$(BUILDDIR)/$(TOOLBUILDPREF)_%$(DEPEXT) : $(TOOLDIR)/%$(SRCEXT)
	@mkdir -p $(BUILDDIRS)
	@echo "Making dependencies for" $< "..."
	@set -e; rm -f $@; \
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
ifneq ($(DEP),)
ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif
endif

# Making dependencies for $(SRCDIR) for CPU
$(BUILDDIR)/$(SRCBUILDPREF)_%$(DEPEXT) : $(SRCDIR)/%$(SRCEXT) $(INCDIR)/%$(HDREXT)
	@mkdir -p $(BUILDDIRS)
	@echo "Making dependencies for" $(subst $(SRCDIR)/,,$<) "..."
	@set -e; rm -f $@; \
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
ifneq ($(DEP),)
ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif
endif

##################
##  clean Part  ##
##################
clean ::
	@echo ""
	@echo "Cleaning up ..."
	@rm -rf $(BUILDDIR)/
	@echo "Done."
	@echo ""
