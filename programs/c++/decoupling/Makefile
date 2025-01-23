# default optimization level, change for debugging
LEVEL = 3

CXX = c++
FC  = gfortran

LHADIR = $(shell lhapdf-config --prefix)

ifeq ($(FC),ifort)
CXXFLAGS = -O$(LEVEL) -std=c++17 -I.. -I$(LHADIR)/include $(OPT)
else
CXXFLAGS = -O$(LEVEL) -std=c++17 -I.. -I$(LHADIR)/include -DGFORTRAN $(OPT)
endif
MAKE_SHR = true

ifeq ($(strip $(MAKE_SHR)),true)
CXXFLAGS += -fPIC -shared
endif

ifeq ($(shell uname), Darwin)
  LD      = $(CXX)
  LDFLAGS =
  LIBS    = -L$(LHADIR)/lib -lLHAPDF
  ifeq ($(FC),gfortran)
  LIBS += -lgfortran -lquadmath
  endif
else
  ifeq ($(FC),ifort)
  LDFLAGS = -nofor-main
  endif
  LD      = $(FC)
  LIBS    = -L.. -L$(LHADIR)/lib -lLHAPDF -lstdc++ -ldl
endif

LIBS += -L/home/tom/Documents/software/software/REvolver-master/code/lib -lREvolver

SRCS =main.cpp
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
EXEC = $(patsubst %.cpp,%,$(SRCS))

.PHONY: clean

target: $(OBJS) $(EXEC)

%: %.o
	$(LD) CRunDec.o $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *.grd *.wgt $(EXEC)
