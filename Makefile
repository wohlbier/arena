HDRS =
SRCS = main.cc
EMU_ASM = $(subst .cc,.emu.s,$(SRCS))
EMU_OBJS = $(subst .cc,.emu.o,$(SRCS))

EMU_PATH = /home/jgwohlbier/devel/packages/emu-19.10
EMU_CXX = $(EMU_PATH)/bin/emu-cc
EMU_SIM = $(EMU_PATH)/bin/emusim.x

CPPFLAGS =
LDFLAGS =

EXE  = arena
EMU_EXE = $(EXE).mwx

all : $(EMU_EXE) main.emu.s

$(EMU_EXE) : $(EMU_OBJS)
	$(EMU_CXX) -o $(EMU_EXE) $(EMU_OBJS) $(LDFLAGS)

%.emu.o: %.cc $(HDRS)
	$(EMU_CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.emu.s: %.cc $(HDRS)
	$(EMU_CXX) $(CPPFLAGS) $(CXXFLAGS) -S -c $< -o $@

.PHONY : clean

clean :
	-$(RM) *~ $(EMU_ASM) $(EMU_OBJS) $(EMU_EXE)
