###############################
#      Compilers & Tools      #
###############################
GCC ?= gcc
GFORTRAN ?= gfortran
DOXYGEN ?= doxygen

###############################
#   OS & Architecture Flags   #
###############################

OS_ARCH ?= x86-64
STDC ?= c99

###############################
#       Compiler Flags        #
###############################

# C & related flags, with debug switch
CFLAGS := -march=$(OS_ARCH) -std=$(STDC) -Wall -Wextra -Wpedantic

ifdef DEBUG
CFLAGS += -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wuninitialized
CFLAGS += -Og -g3
else
CFLAGS += -O2 -flto -fstack-protector -Wno-unused
endif

###############################
#        Libraries            #
###############################

MERCURY_OBJS := forbel.o
CMERCURY_LIB_OBJS := orbel.o danby.o mxx.o mio.o mco.o
CMERCURY_LIB := libcmercury.so

###############################
#         Binaries            #
###############################

CMERCURY_TEST_BIN := main.bin
CMERCURY_TEST_OBJ := main.o

###############################
#      Linker (LD) Flags      #
###############################

LDFLAGS := -Wl,-rpath=`pwd`,--enable-new-dtags -L. -lm

###############################
#       Definition Flags      #
###############################

DEFINES ?=
ifdef DEBUG
DEFINES += -DDEBUG
else
DEFINES += -DNDEBUG
endif

###############################
#    Include Path & Flags     #
###############################

INCLUDES ?= -I. -I/usr/include

###############################
#        Make Targets         #
###############################

OBJECTS		:= $(wildcard *.o)
LIBRARIES	:= $(wildcard *.so)
BINARIES	:= $(wildcard *.bin)
TEMPFILES   := $(wildcard *.out)

###############################
#         Make Rules          #
###############################

.PHONY: clean all build doc

build: all

all: $(CMERCURY_LIB)

test: $(CMERCURY_TEST_BIN)
	./$<

clean:
	@rm -f $(OBJECTS) $(LIBRARIES) $(BINARIES) $(TEMPFILES)

###############################
# Mercury Regression Testing  #
###############################

$(MERCURY_OBJS): %.o : %.for
	$(GFORTRAN) -o $@ -c $<

#####################
#  CMercury Library #
#####################

$(CMERCURY_LIB_OBJS): %.o : %.c
	$(GCC) -fPIC $(CFLAGS) $(INCLUDES) $(DEFINES) -o $@ -c $<

$(CMERCURY_LIB): $(CMERCURY_LIB_OBJS)
	$(GCC) -shared -Wl,-soname,$@ -o $@ $^ $(LDFLAGS)

########################
# CMercury Test Binary #
########################

$(CMERCURY_TEST_OBJ): %.o : %.c
	$(GCC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

$(CMERCURY_TEST_BIN): $(CMERCURY_LIB) $(MERCURY_OBJS) $(CMERCURY_TEST_OBJ)
	$(GCC) $(MERCURY_OBJS) $(CMERCURY_TEST_OBJ) -o $@ $(LDFLAGS) -lcmercury -lgfortran

#######################
# Documentation Build #
#######################

doc: 
	$(DOXYGEN) Doxygen.conf
