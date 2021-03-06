#
# File:  GNUmakefile
#
#----------------------------------------------------------------------------
#
#  This is the main driver makefile for compiling SCRIP.  It relies on
#  sub-makefiles to perform preprocessing, dependency analysis and compiling.
#
#  Several variables must be set either as environment variables or on the 
#  command line (e.g. make SCRIPDIR=/your/path).  These variables are:
#
#  SCRIPEXEDIR = the directory where you want the executable to end up.  A
#              common value might be . (the current working directory)
#  SCRIPDIR    = the path to the SCRIP distribution from which you want
#              to grab source files, etc.  This should be the path to the 
#              top level SCRIP directory, not the individual source,mpi 
#              directories.
#  SCRIPARCH   = the base name of a file in the build directory where 
#              architecture and location-specific information is defined
#
#  The optimize variable has a default value of yes and can not be set
#  as an environment variable.  It can only be changed to no via the
#  command line.
#  OPTIMIZE  = yes/no  whether you want full optimization or no optimization
#
#----------------------------------------------------------------------------
#
#  The bourne shell is safest for make
SHELL = /bin/sh
#
#
#  By default, you build optimized code.  To override, add "OPTIMIZE=no" 
#  to the make command line.
OPTIMIZE   = yes
#
#----------------------------------------------------------------------------
#
#  Check that path variables SCRIPEXEDIR, SCRIPDIR have been set.
#
#----------------------------------------------------------------------------

VARSDEFINED =

ifneq (,$(SCRIPEXEDIR))
  VARSDEFINED := $(VARSDEFINED)SCRIPEXEDIR
  export SCRIPEXEDIR
else
  bogusexedir:
	@echo "  Please setenv SCRIPEXEDIR"
endif

ifneq (,$(SCRIPDIR))
  VARSDEFINED := $(VARSDEFINED)SCRIPDIR
  export SCRIPDIR
else
  bogusSCRIPdir:
	@echo "  Please setenv SCRIPDIR"
endif

#  Set the directory where dependency info will reside.
DepDir = $(SCRIPEXEDIR)/compile/Depends

#----------------------------------------------------------------------------
#
#  Include architecture-specific flags and options. 
#
#----------------------------------------------------------------------------

ifneq (,$(SCRIPARCH))
  export SCRIPARCH
  VARSDEFINED := $(VARSDEFINED)SCRIPARCH
  ifneq (,$(SCRIPDIR))
     include $(SCRIPDIR)/build/$(SCRIPARCH).gnu
  endif
else
  bogusSCRIParch:
	@echo "  Please setenv SCRIPARCH"
endif

#----------------------------------------------------------------------------
#
#  Define a symbol (TARGETX) for the executable name (SCRIP)
#
#----------------------------------------------------------------------------

TARGETX = scrip

ifeq ($(OPTIMIZE),no)
  #  If building a debug version, append "_db" to the executable name
  TARGETX = scrip_db
endif

#----------------------------------------------------------------------------
#
#  Make the SCRIP executable
#
#----------------------------------------------------------------------------

.PHONY: depends preprocess compile dirs

$(TARGETX): depends
	@echo "  GNUmakefile is making target '$(TARGETX)'"
	@make -f $(SCRIPDIR)/build/compile.mk SCRIPEXEDIR=$(SCRIPEXEDIR) SCRIPDIR=$(SCRIPDIR) OPTIMIZE=$(OPTIMIZE) TARGETX=$(TARGETX) 

#----------------------------------------------------------------------------
#
#  Create the file dependencies for compiling.
#
#----------------------------------------------------------------------------

depends: preprocess
	@echo "$(SCRIPEXEDIR) $(SCRIPDIR)"
	@echo "$(SCRIPARCH) automatically generating dependencies for compiling"
	@make -f $(SCRIPDIR)/build/depends.mk SCRIPEXEDIR=$(SCRIPEXEDIR) OPTIMIZE=$(OPTIMIZE)

#----------------------------------------------------------------------------
#
#  Perform the preprocessing
#
#----------------------------------------------------------------------------

preprocess: dirs
	@echo "$(SCRIPARCH) preprocess SCRIP source"
	@make -f $(SCRIPDIR)/build/preprocess.mk

#----------------------------------------------------------------------------
#
#  Create the directory structure if it does not exist.
#
#----------------------------------------------------------------------------

dirs:
	@echo "$(SCRIPARCH) is creating the required directory structure"
	@cd $(SCRIPEXEDIR) &&  if [ ! -d compile ]; then mkdir compile; fi
	@cd $(SCRIPEXEDIR)/compile &&  \
            if [ ! -d Depends ]; then mkdir Depends; fi 
	@cd $(SCRIPEXEDIR)

#----------------------------------------------------------------------------
#
#  Include utilities makefile with rules for clean, clobber, etc.
#
#----------------------------------------------------------------------------

clean:
ifeq ($(VARSDEFINED),SCRIPEXEDIRSCRIPDIRSCRIPARCH)
	@cd $(SCRIPEXEDIR) && $(RM) core scrip scrip_db 
	@cd $(SCRIPEXEDIR) && $(RM) compile/Depends/*.* 
	@cd $(SCRIPEXEDIR) && $(RM) compile/*.f compile/*.f90 compile/*.c
	@cd $(SCRIPEXEDIR) && $(RM) compile/*.o compile/*.$(MODSUF) 
else
	@echo "  Please setenv SCRIPEXEDIR SCRIPDIR and SCRIPARCH"
	@echo "  Only $(VARSDEFINED) have been defined."
endif

#----------------------------------------------------------------------------
