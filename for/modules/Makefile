FLAGS=
SRC := $(wildcard */*for)
MAIN = specyl.for cyl1.inc.for cyl2.inc.for files.inc.for eq_csha.blc.for ppcon.blc.for ppconcyl.blc.for pppoincare.blc.for settings.blc.for
SRCS = $(filter-out $(MAIN),$(SRC))
SCFUNAME = FOR0
PWD = `pwd`
NULL = "  "
NSEC = `pwd | awk '{sub(/for\//," "); print $$2}'`
PARDIR = $(shell dirname $(PWD)/../../)
BASENAME = "/ricercatori/ft/specyl/veranda/flow/"
SPDIR = $(shell readlink -f $(BASENAME))
BTCHFILE = $(BATCHFILE)_$(NSEC)
PROG_BIN = $(SIMDIR)bin/flwspc_
exec = $(wildcard $(PARDIR)/bin/*bin)
SC_OTHER_SECTIONS_MIN := $(foreach sc_minus, $(SC_OTHER_SECTIONS), $$(( $(sc_minus) - 1 )) )
.PHONY : append append2

include ../settings.inc.mk

define batch_templ
echo $(BATCHFILE)
echo "#!/bin/bash" >> $(BATCHFILE)
echo "#PBS -l select=1:mpiprocs=$(SC_NUM_THREADS):ncpus=$(SC_NUM_THREADS):mem=100mb:host=rat2" >> $(BATCHFILE)
#echo "#PBS -l select=1:mpiprocs=$(SC_NUM_THREADS):ncpus=$(SC_NUM_THREADS):mem=100mb:host=ra36" >> $(BATCHFILE)
echo "#PBS -l walltime=2000:00:00" >> $(BATCHFILE)
echo "#PBS -M $(USER)@igi.cnr.it" >> $(BATCHFILE)
echo "#PBS -j oe" >> $(BATCHFILE)
echo "#PBS -m ae" >> $(BATCHFILE)
echo "#PBS -r n" >> $(BATCHFILE)
echo "#PBS -q plong" >> $(BATCHFILE)
echo "#PBS -o $(ERRORFILE)" >> $(BATCHFILE)
echo export OMP_NUM_THREADS=$(SC_NUM_THREADS) >> $(BATCHFILE)
#echo time $(SIMDIR)bin/flwspc_$(1).bin '>>' $(LOGFILE) >> $(BATCHFILE)
echo time $(SIMDIR)bin/flwspc_$(1).bin $(1) '>>' $(LOGFILE) >> $(BATCHFILE)
echo "Sottometto il job per il run di specyl"
endef

define append_templ0

echo "       prefix="'"$(PREFIX)"'" " >> $(SIMFORDIR)/$(1)/files.inc.for
echo "       simdir="'"/$(SIM)"'" " >> $(SIMFORDIR)/$(1)/files.inc.for
echo "       un_all="'"specyl_$(1)_all.dat"'"" >> $(SIMFORDIR)/$(1)/files.inc.for
endef

define append_templ1

echo "       un_end="'"specyl_$(1)_end.dat"'"" >> $(SIMFORDIR)/$(1)/files.inc.for 
endef

append : 
	@$(foreach section,$(SC_SECTIONS),$(call append_templ0,$(section)))
	@$(foreach section,$(SC_SECTIONS),$(call append_templ1,$(section)))
	@if [ "$(SC_INPUT_1)" = "" ]; then \
	 echo "      un_prevend='dum_0_end.dat'" >> $(SIMFORDIR)/1/files.inc.for ; \
	 else echo "       un_prevend='specyl_0_end.dat'" >> $(SIMFORDIR)/1/files.inc.for ; \
	 fi;

append2 :
	@for x in $(SC_OTHER_SECTIONS) ; do \
	 echo "       un_prevend='specyl_$$(( $$x - 1 ))_end.dat'" >> $(SIMFORDIR)/$${x}/files.inc.for ; \
	done;

launch : 
	@$(foreach section,$(SC_SECTIONS),$(call batch_templ,$(section)))
	@qsub $(BATCHFILE)

