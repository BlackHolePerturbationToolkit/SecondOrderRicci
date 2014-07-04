import sys

env = DefaultEnvironment()

# pretty output
if ARGUMENTS.get('VERBOSE') != '1':
  if sys.stdout.isatty():
    env['CXXCOMSTR'] = "\033[92mCompiling\033[0m $TARGET"
    env['LINKCOMSTR'] = "\033[94mLinking\033[0m $TARGET"
  else:
    env['CXXCOMSTR'] = "Compiling $TARGET"
    env['LINKCOMSTR'] = "Linking $TARGET"

# Build options
env['LIBS']     = ['gsl', 'm', 'gslcblas', 'mpi_cxx', 'mpi', 'hdf5', 'hdf5_hl', 'boost_program_options-mt']
env['LIBPATH']  = ['/usr/local/lib/']
env['CPPPATH']  = ['/usr/local/include']
env['CXXFLAGS'] = ['-O2', '-std=c++11', '-Wall', '-g', '-O2', '-D__GIT_VERSION=\"$(GIT_VERSION)\"]']
env['GIT_VERSION'] = "$(shell sh -c 'git describe --dirty --always 2>/dev/null || git describe --always')"
env['CXX']      = 'clang++'


sources = ['Coupling.cc', 'h1.cc', 'h5wrapper.cc', 'R2.cc', 'Ricci.cc', 'utils.cc']
executable = 'Ricci'

Program(executable, sources)
Decider('MD5-timestamp')
