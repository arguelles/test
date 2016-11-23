# the compiler		
CC = clang++

# compiler flags
CFLAGS = -Wall -g -std=c++11

# define any directories containing header files

INCLUDES = -I/Users/marjon/Dropbox\ \(MIT\)/work/IceCube/decay/neutrino_decay/verosimilitud/inc -I/usr/include/python2.7 -I/usr/local/Cellar/hdf5/1.8.16_1/include -I/Users/marjon/local/dlib-18.18
INCLUDES += -I/Users/carguelles/DropboxMIT/NeutrinoDecay/verosimilitud/inc

# define any library paths in addition to /usr/lib
LFLAGS = -L/usr/local/Cellar/hdf5/1.8.16_1/lib
LFLAGS += -L/Users/marjon/Dropbox\ \(MIT\)/work/IceCube/decay/neutrino_decay/verosimilitud/
LFLAGS += -L/Users/carguelles/DropboxMIT/NeutrinoDecay/verosimilitud/

# define any libraries to link into executable:
LIBS = -lpython2.7 -lm -lhdf5 -lhdf5_cpp 
# adding reference to llh library
LIBS += -lverosimilitud

SRCS := test.cpp

#SRCS := test.cpp $(wildcard /Users/marjon/Dropbox\ #\(MIT\)/work/IceCube/decay/neutrino_decay/verosimilitud/src/*)

#SRCS := test.cpp /Users/marjon/Dropbox\ \(MIT\)/work/IceCube/decay/neutrino_decay/verosimilitud/src/* 

#EXCLUDE=$(subst /Users/marjon/Dropbox\ \(MIT\)/work/IceCube/decay/neutrino_decay/verosimilitud/src/Verosimilitud.cpp,,${SRCS})

OBJS = $(SRCS:.c=.o)

# the build target executable:
TARGET = test

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

clean:
	$(RM) $(TARGET)
