##
## NOTE: When compiling python, this Makefile is only used to call the
## python build process.  No variables set in this file affect python.
## If the python build is having trouble, look in setup.py to fix it.
##



# By default, just call the python build process
all :
	@echo "The default build process has changed, and this Makefile is just for useful hints in unusual cases."
	@echo ""
	@echo "You probably want to run the following to build and install the code:"
	@echo "  python setup.py install --user"

# Optionally, exclude GSL and the dependent functions
no-GSL :
	$(shell python setup.py install --user --no-GSL)

# If needed, we can also make object files to use in other C++ programs
cpp : Quaternions.o QuaternionUtilities.o IntegrateAngularVelocity.o



# Mike likes to use python's virtual environments to build with
# various versions and arrangements of python.  As such, he finds that
# the following is, unfortunately, necessary for him.  This makes him
# happy.  It is probably not necessary for most users, however.
MikeHappy :
	$(VIRTUAL_ENV)/bin/python setup.py install --prefix=$(VIRTUAL_ENV)

















############################################
## Flags for building only the c++ files: ##
############################################
# The compiler needs to be able to find the GSL (GNU Scientific
# Library) headers and libraries.  The following paths are the most
# common places for these to be installed.  If compilation doesn't
# work, correct these paths.
INCFLAGS = -I/opt/local/include -I/usr/local/include
LIBFLAGS = -L/opt/local/lib -L/usr/local/lib
ifdef GSL_HOME
	INCFLAGS = -I${GSL_HOME}/include ${INCFLAGS}
	LIBFLAGS = -L${GSL_HOME}/lib ${LIBFLAGS}
endif
# Set compiler name and optimization flags here, if desired
C++ = g++
OPT = -O3 -Wall -Wno-deprecated
## DON'T USE -ffast-math in OPT


#############################################################################
## The following are pretty standard and probably won't need to be changed ##
#############################################################################

# Tell 'make' not to look for files with the following names
.PHONY : all cpp clean allclean realclean swig doc

# This rebuilds the documentation, assuming doxygen is working
doc :
	make -C docs

# If needed, we can also make object files to use in other C++ programs
cpp : Quaternions.o IntegrateAngularVelocity.o

# This is how to build those object files
%.o : %.cpp %.hpp Errors.hpp
	$(C++) $(OPT) -c $(INCFLAGS) -DUSE_GSL $< -o $@

# The following are just handy targets for removing compiled stuff
clean :
	-/bin/rm -f *.o
allclean : clean
	-/bin/rm -rf build
realclean : allclean

# This just runs SWIG, which can be handy for committing pre-swigged files
swig :
	swig -python -globals constants -c++ -o Quaternions_wrap.cpp Quaternions.i
