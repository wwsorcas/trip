CC = 'clang', 
CCFLAGS = '-Ofast -pedantic -Wunused -Wno-dangling-else -Wno-deprecated-register -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE'
CFLAGS = '-std=c99 -Wimplicit'
LINKFLAGS = '-pthread'
CXXFLAGS = '-std=c++11'
CXX = 'clang++'
