ROOTCONFIG = root-config

CXX     =$(shell $(ROOTCONFIG) --cxx)
GCC	=$(shell $(ROOTCONFIG) --cc)
CXXFLAGS	= -O3 -Wall -fPIC -Wextra
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs) -lMinuit
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)
ROOTINC :=$(shell $(ROOTCONFIG) --incdir)

CUSTOMDEFINES = -DNOISY -DDEBUG 
CXXFLAGS+= -g -std=c++11 $(CUSTOMDEFINES) $(ROOTCFLAGS)
MAKEDEPEND =$(CXX)
INCLUDES = -I$(ROOTINC)

LIBS	:= $(ROOTLIBS)
LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LD	:= $(shell $(ROOTCONFIG) --ld)
OBJS	:= TaConfig.o TaInput.o TaOutput.o TaChannel.o TaAccumulator.o VAnalysisModule.o TaRegression.o
HDR	:= $(OBJS:.o=.hh) 
DEPS 	:= $(OBJS:.o=.d)

all:  $(OBJS) main.o dance_Dict dance libdance.so

dance_Dict: $(HDR) danceLinkDef.hh
	rootcint -f $@.cc -c $(INCLUDES) $^;
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $@.cc -o $@.o ;
$(OBJS):
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(@:.o=.cc) -o  $@ ;

main.o: main.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^ -o $@;

dance:$(OBJS) dance_Dict.o main.o
	  $(LD) $(CXXFLAGS) $(ROOTLIBS) $^ -o $@;
libdance.so: $(OBJS) dance_Dict.o
	$(LD) $(LDFLAGS) $(ROOTLIBS) -shared $^ -o $@;
clean:
	rm -f *.o;
	rm -f *Dict*;	
install:all
	mv dance ../;
%.d:	%.cc
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -std=c++11 -MM $(INCLUDES) -c $< \
	      | sed '\''s%^.*\.o%$*\.o%g'\'' \
	      | sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
	      [ -s $@ ] || rm -f $@'

-include $(DEPS)
