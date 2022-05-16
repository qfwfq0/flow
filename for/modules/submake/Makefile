FLAGS=
#DEBUG=-g -traceback -warn interface
DEBUG=
SRC := $(wildcard *for)
MAIN = specyl.for cyl1.inc.for files.inc.for specyl.deb.for
SRCS = $(filter-out $(MAIN),$(SRC))
PWD=`pwd`
NSEC=`pwd | awk '{sub(/for\//," "); print $$2}'`
PARDIR=$(shell dirname $(PWD)/../../)
BASENAME="/ricercatori/ft/specyl/veranda/flow/"
SPDIR=$(shell readlink -f $(BASENAME))
BTCHFILE=$(BATCHFILE)_$(NSEC)
LGFILE=$(LOGFILE)_$(NSEC)

subcompile : obj
	@echo "compilo specyl, sec $(NSEC)"
	@echo $(FC) $(DEBUG) $(OMP) $(FLAGS) $(SRCS:.for=.o) -I$(PWD)/ -o $(SIMDIR)bin/flwspc_$(NSEC).bin $(PWD)/specyl.for
	@$(FC) $(DEBUG) $(OMP) $(FLAGS) $(SRCS:.for=.o) -I$(PWD)/ -o $(SIMDIR)bin/flwspc_$(NSEC).bin $(PWD)/specyl.for

obj : $(SRCS:.for=.o) 

%.o : %.for
	@echo "Compiling" $< "sec= " $(NSEC)
#	@echo $(FC) -I$(PWD) -c $<
	@$(FC) -I$(PWD) -c $<

clean : 
	@rm -f *.o
	@rm -f *.mod

