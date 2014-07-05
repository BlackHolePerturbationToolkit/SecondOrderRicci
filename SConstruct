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
env['LIBS']     = ['gsl', 'm', 'hdf5']
env['LIBPATH']  = ['/usr/local/lib/']
env['LINKFLAGS']= ['-fopenmp']
env['CPPPATH']  = ['/usr/local/include']
env['CXXFLAGS'] = ['-O3', '-DBOOST_DISABLE_ASSERTS', '-fopenmp', '-std=c++11', '-Wall', '-g', '-D__GIT_VERSION=\"$(GIT_VERSION)\"]']
env['GIT_VERSION'] = "$(shell sh -c 'git describe --dirty --always 2>/dev/null || git describe --always')"
env['CXX']      = 'g++-4.9'


sources = ['Coupling.cc', 'h1.cc', 'h5wrapper.cc', 'Ricci.cc', 'utils.cc',
           'R2_1.cc', 'R2_2.cc', 'R2_3.cc', 'R2_4.cc', 'R2_5.cc',
           'R2_6.cc', 'R2_7.cc', 'R2_8.cc', 'R2_9.cc', 'R2_10.cc']
executable = 'Ricci'

Program(executable, sources)
Decider('MD5-timestamp')
