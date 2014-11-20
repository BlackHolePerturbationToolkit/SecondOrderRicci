import sys
from subprocess import Popen, CalledProcessError, PIPE

def git_version():
    try:
        version = Popen(['git', 'describe', '--dirty', '--always'], stdout=PIPE).communicate()[0]
    except CalledProcessError:
        version = Popen(['git', 'describe', '--always'], stdout=PIPE).communicate()[0]
    return version.decode('utf_8').rstrip('\n')

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
env['LIBS']     = ['gsl', 'gslcblas', 'm', 'hdf5']
env['LIBPATH']  = ['/usr/local/lib/']
env['LINKFLAGS']= ['-fopenmp']
env['RPATH']    = ['/opt/gcc/4.9.2/lib64']
env['CPPPATH']  = ['/usr/local/include', './Quaternions', './SphericalFunctions']
env['CXXFLAGS'] = ['-O3', '-DBOOST_DISABLE_ASSERTS', '-fopenmp', '-std=c++11', '-Wall', '-g', '$(-D__GIT_VERSION="\\"' + git_version() + '\\""$)']
env['CXX']      = '/opt/gcc/4.9.2/bin/g++'


sources = ['Coupling.cc', 'h1.cc', 'h5wrapper.cc', 'Ricci.cc', 'utils.cc',
           'R2_1.cc', 'R2_2.cc', 'R2_3.cc', 'R2_4.cc', 'R2_5.cc',
           'R2_6.cc', 'R2_7.cc', 'R2_8.cc', 'R2_9.cc', 'R2_10.cc',
           'hh_1.cc', 'hh_2.cc', 'hh_3.cc', 'hh_4.cc', 'hh_5.cc',
           'hh_6.cc', 'hh_7.cc', 'hh_8.cc', 'hh_9.cc', 'hh_10.cc',
           'Quaternions/Quaternions.cpp', 'Quaternions/QuaternionUtilities.cpp',
           'SphericalFunctions/Combinatorics.cpp', 'SphericalFunctions/WignerDMatrices.cpp']
executable = 'Ricci'

Program(executable, sources)
Decider('MD5-timestamp')
