###############################
#      Compilers & Tools      #
###############################
GCC := gcc

###############################
#       COMPILER FLAGS        #
###############################

# C & related flags, with debug switch
CFLAGS := -std=c99 -Wall -Wextra -Wpedantic

ifdef DEBUG
CFLAGS += -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wuninitialized
CFLAGS += -Og -g3 -DDEBUG
else
CFLAGS += -march=x86-64 -O2 -DNDEBUG
endif

###############################
#        Libraries            #
###############################

CMERCURY_LIB_OBJS := orbel.o danby.o
CMERCURY_LIB := libcmercury.so

###############################
#         Binaries            #
###############################

CMERCURY_TEST_BIN := main.bin
CMERCURY_TEST_SRC := main.c

###############################
#      Linker (LD) Flags      #
###############################

LDFLAGS := -L. -lm

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

INCLUDES := -I.

###############################
#        Make Targets         #
###############################

OBJECTS		:= $(wildcard *.o)
LIBRARIES	:= $(wildcard *.so)
BINARIES	:= $(wildcard *.bin)

###############################
#         Make Rules          #
###############################

.PHONY: clean all build

build: all

all: $(CMERCURY_LIB)

test: $(CMERCURY_TEST_BIN)
	./$<

clean:
	@rm -f $(OBJECTS) $(LIBRARIES) $(BINARIES) *.s *.i *.ii

# ------- CMercury Library -------

$(CMERCURY_LIB_OBJS): %.o : %.c
	$(GCC) -fPIC $(CFLAGS) $(INCLUDES) $(DEFINES) -o $@ -c $<

$(CMERCURY_LIB): $(CMERCURY_LIB_OBJS)
	$(GCC) -shared -Wl,-soname,$@ -o $@ $^ $(LDFLAGS)

# ------- CMercury Test Binary -------

$(CMERCURY_TEST_BIN): $(CMERCURY_LIB)
	$(GCC) $(CFLAGS) $(INCLUDES) $(CMERCURY_TEST_SRC) -o $@ $(LDFLAGS) -lcmercury
