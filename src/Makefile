
CFLAGS = $(shell gsl-config --cflags)
##CFLAGS += -ggdb -gdwarf-2 -Wall -I include -O
CFLAGS +=-O3 -I include
MPCFLAG = -I /cm/shared/mpi/openmpi/2.1/intel/17.0/include
LIBS   = $(shell gsl-config --libs)

CC     = g++
MPCC 	= mpicxx

BDIR   = bin
ODIR   = obj
SDIR   = subroutines
EDIR   = EXEs
OS    := $(shell uname)

ifeq ($(OS),Linux)
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -r 's/(\.cc|.cpp)/.o/')
else
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -E 's/(\.cc|.cpp)/.o/')
endif

_EXECUTABLES = rd_gen_reweight_NOPBC_complex \




OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

EXECUTABLES = $(patsubst %,$(BDIR)/%,$(_EXECUTABLES))

PAREXECUTABLES = $(patsubst %,$(BDIR)/%,$(_PAREXECUTABLES))

_SOURCES = ${_EXECUTABLES:=.cpp}
SOURCES = $(patsubst %,$(EDIR)/%,$(_SOURCES))

_PARSOURCES = ${_PAREXECUTABLES:=.cpp}
PARSOURCES = $(patsubst %,$(EDIR)/%,$(_PARSOURCES))

all: dirs $(EXECUTABLES) $(PAREXECUTABLES)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) $(CFLAGS2) -c $< -o $@  
	@echo "------------"

$(EXECUTABLES): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS)
	@echo "------------"

dirs: 
	mkdir -p bin
	mkdir -p obj

clean:
	rm -f $(ODIR)/*.o $(EXECUTABLES)


