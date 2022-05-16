##;# Marco: definition of some input to the makefile

#;# SIM: relative path of the simulation folder
#;# MF: if MF==1 then the b,v,j field will also be written on disk using the modulus-phase notation, elseif MF==0 or not defined the fields will not be written on disk using the modulus-phase notation

#SHELL = /bin/sh
export
include make.inc

#;# Marco, deal with TOK
ifndef $(TOK)
	TOK = 0
endif

ifeq ($(TOK),1)
PREFIX = /0/rfx/ricercatori/ft/specyl/veranda/flow
SPECYLBIN=/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/specyl.xxx
POSTCODE = /0/rfx/ricercatori/ft/specyl/veranda/flow/for/postproc.F #!# code for postprocessing
POSTEXEC = /0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx #!# code for postprocessing
BINDIR = /0/rfx/ricercatori/ft/specyl/veranda/flow/bin
FORDIR = /0/rfx/ricercatori/ft/specyl/veranda/flow/for
SIMDIR = $(PREFIX)/$(SIM)
SETDIR = /0/rfx/ricercatori/ft/specyl/veranda/flow/for/settings
POSTFILE = $(PREFIX)/$(SIM)log/postproc.sh
endif
ifeq ($(TOK),0)
PREFIX = /0/rfx/ricercatori/ft/specyl/veranda/flow
#PREFIX = #/0/rfx/ricercatori/ft/specyl/veranda/flow
SPECYLBIN=/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/specyl.xxx
BINDIR = $(PREFIX)/bin
FORDIR = $(PREFIX)/for
POSTCODE = $(FORDIR)/postproc.F #!# code for postprocessing
POSTBIN=$(PREFIX)/bin/postproc.xxx
POSTEXEC = $(BINDIR)/postproc.xxx
SIMDIR = $(PREFIX)/$(SIM)
POSTFILE = $(PREFIX)/$(SIM)log/postproc.sh
SETDIR = $(PREFIX)/for/settings
endif


FC=ifort
OMP=-qopenmp
#!# Marco, gfortran
#FC=gfortran
#OMP=-fopenmp
#FC=pgf95
#OMP=-mp
#DEBUG=-g -traceback -warn interface
DEBUG=
SKYLAKE=-xcore-AVX512 -mtune=skylake
USER=veranda
#SPECYLBIN=$(PREFIX)/bin/specyl.xxx
POINCBIN=$(PREFIX)/bin/pppoincare.xxx
POSTBIN=$(PREFIX)/bin/postproc.xxx
ENDBIN=$(PREFIX)/bin/write_end.xxx
MODDIR = $(PREFIX)/for/modules
MODDIR2 = $(PREFIX)
PRODIR = $(PREFIX)/pro
SIMFORDIR = $(SIMDIR)for
SIMLOGDIR = $(SIMDIR)log
SIMDATDIR = $(SIMDIR)dat
SIMBINDIR = $(SIMDIR)bin
SC_BIN = $(wildcard $(SIMBINDIR)/*.bin )
SC_SAV = $(wildcard $(SIMDATDIR)/specyl_*_all.sav )
SC_DAT = $(wildcard $(SIMDATDIR)/*_end.dat )
SCFUNAME = FOR0
SIMPATH = $(strip $(subst /, ,$(SIM)))
DATE = $(shell date +%Y%m%d%H%M)
#SIMNAME = $(word $(words $(SIMPATH)), $(SIMPATH))
SIMNAMEDATE = flw_$(DATE)
POINCNAMEDATE = poinc_$(DATE)
LOGFILE = $(SIMLOGDIR)/$(SIMNAMEDATE).log
POINCLOGFILE = $(SIMLOGDIR)/$(POINCNAMEDATE).log
POSTERRORFILE = $(SIMLOGDIR)/postproc.err
ERRORFILE = $(SIMLOGDIR)/$(SIMNAMEDATE).err
POINCERRORFILE = $(SIMLOGDIR)/$(POINCNAMEDATE).err
BATCHFILE =  $(SIMLOGDIR)/$(SIMNAMEDATE).pbs
ENDFILE =  $(SIMLOGDIR)/write_end.sh
POINCFILE =  $(SIMLOGDIR)/ppoinc.fdryk
SIMDIR0 = $(PREFIX)/$(SIM0)
SIMDIR1 = $(PREFIX)/$(SIM1)
#vpath %.for $(MODDIR)

SC_NUM_THREADS = 6
SIMDIR = $(PREFIX)/$(SIM)

#SC_NUM_SECTIONS = $(words $(SC_SECTIONS))
#;# Marco, I want to change the definition of SC_NUM_SECTIONS to deal with the case when not all the sections written in settings.inc.mk have already been computed.
#NUMBER_OF_FILES = $(shell find $(SIMDATDIR) -name specyl_*_all.dat)
NUM_OF_FILES = $(wildcard $(SIMDATDIR)/specyl_*_all.dat)
NUMBER_OF_FILES = $(words $(NUM_OF_FILES))
#;# dovrei dire che se NUMBER_OF_FILES è minore di SC_SECTIONS allora bisogna usare NUMBER_OF_FILES


#SC_OUTPUTS_ALL = $(foreach section,$(SC_SECTIONS),$(SC_OUTPUT_$(section)_ALL))
#SC_OUTPUTS_END = $(foreach section,$(SC_SECTIONS),$(SC_OUTPUT_$(section)_END))
#ONE := 1
#SC_OTHER_SECTIONS := $(filter-out $(ONE),$(SC_SECTIONS))


#;# Marco, define newline
define nl


endef
define pv
;
endef
define andand
&&
endef

#;# Marco, deal with MF
ifndef $(MF)
	MF = 1	
endif
#;# Marco, deal with 3D
ifndef $(3D)
	3D = 1
endif
#;# Marco, deal with EQU
ifndef $(EQU)
	EQU = 1
endif
#;# Marco, deal with HEL
ifndef $(HEL)
	HEL = 10
endif
#;# Marco, deal with TMIN
ifndef $(TMIN)
	TMIN = 0
endif
#;# Marco, deal with TMAX
ifndef $(TMAX)
	TMAX = 1
endif
#;# Marco, deal with TT
ifndef $(TT)
	TT = 0
endif
#;# Marco, deal with MM
ifndef $(TMAX)
	MM = 1
endif
#;# Marco, deal with NN
ifndef $(TMAX)
	NN = -10
endif
#;# Marco, deal with NTH
ifndef $(NTH)
	NTH = 1
endif
#;# Marco, deal with NPHI
ifndef $(NPHI)
	TT = 1
endif
#;# Marco, deal with RDS
ifndef $(RDS)
	RDS = 69
endif
#;# Marco, deal with NSEC
ifndef $(NSEC)
	NSEC = 1
endif
#;# Marco, deal with ZED
ifndef $(ZED)
	ZED = 0
endif



rwildcard=$(wildcard imf_?profiles_?.sav)
TEMP_SAV = $(wildcard $(SIMDATDIR)/imf_?profiles_?.sav)
#!# Marco, molto potente, elimina tutti i file sav temporanei
ALL_DIR = $(wildcard $(CHOSEPATH)/*/dat/imf_?profiles_?.sav)

SC_LAST_SECTIONS = $(wordlist 2, $(SC_NUM_SECTIONS), $(SC_SECTIONS))

#!# Marco, parte in cui definisco i defaults-value per il numero di sezioni della simulazione
#ifndef $(NSEC)
#	NSEC = 1
#endif

SECSEC := $(shell seq 1 $(NSEC))
SC_SECTIONS = $(SECSEC)
SECSEC_LESS := $(shell echo ${NSEC}- 1 | bc)
SC_SECTIONS_LESS := $(shell seq 1 $(SECSEC_LESS))

#ifndef $(SC_SECTIONS)
#	SC_SECTIONS = $(SECSEC)
#endif

#!#! Marco, parte di esecuzione

define pbs_templ
echo $(BATCHFILE)
echo "#!/bin/bash" >> $(BATCHFILE)
echo "#PBS -l select=1:mpiprocs=$(SC_NUM_THREADS):ncpus=$(SC_NUM_THREADS):mem=100mb:host=rat2" >> $(BATCHFILE)
echo "#PBS -l walltime=700:00:00" >> $(BATCHFILE)
echo "#PBS -M $(USER)@igi.cnr.it" >> $(BATCHFILE)
echo "#PBS -j oe" >> $(BATCHFILE)
echo "#PBS -m ae" >> $(BATCHFILE)
echo "#PBS -r n" >> $(BATCHFILE)
echo "#PBS -q lm" >> $(BATCHFILE)
echo "#PBS -o $(ERRORFILE)" >> $(BATCHFILE)
echo export OMP_NUM_THREADS=$(SC_NUM_THREADS) >> $(BATCHFILE)
echo cd $(SIMDIR) >> $(BATCHFILE)
echo module load intel/pe-xe-2018 >> $(BATCHFILE)
echo time $(SPECYLBIN) $(1) '>>' $(LOGFILE) >> $(BATCHFILE)
echo "Sottometto il job per il run di specyl"
endef

define pbs_templ2
echo $(SPECYLBIN) $(1) '>>' $(LOGFILE) '&&' >> $(BATCHFILE)
endef
define pbs_templ3
echo $(SPECYLBIN) $(1) '>>' $(LOGFILE) >> $(BATCHFILE)
endef

define pbs_poinc
echo $(POINCFILE)
touch $(POINCFILE)
echo "ciao"
echo "#!/bin/bash" >> $(POINCFILE)
echo "#PBS -l select=1:mpiprocs=$(SC_NUM_THREADS):ncpus=$(SC_NUM_THREADS):mem=100mb:host=rat2" >> $(POINCFILE)
echo "#PBS -l walltime=70:00:00" >> $(POINCFILE)
echo "#PBS -M $(USER)@igi.cnr.it" >> $(POINCFILE)
echo "#PBS -j oe" >> $(POINCFILE)
echo "#PBS -m ae" >> $(POINCFILE)
echo "#PBS -r n" >> $(POINCFILE)
echo "#PBS -q plong" >> $(POINCFILE)
echo "#PBS -o $(POINCERRORFILE)" >> $(POINCFILE)
echo export OMP_NUM_THREADS=$(SC_NUM_THREADS) >> $(POINCFILE)
echo cd $(SIMDIR) >> $(POINCFILE)
echo time $(POINCBIN) '>>' $(POINCLOGFILE) >> $(POINCFILE)
echo "Sottometto il job per fare il poincare file"
endef

#SC_END = $(wildcard $(SIMDAT)/specyl_*_end.dat)
SC_END = $(foreach SECTION,$(SC_SECTIONS),$(SIMDAT)/specyl_$(SECTION)_end.dat)
SC_DAT = $(foreach SECTION,$(SC_SECTIONS),$(SIMDAT)/specyl_$(SECTION)_all.dat)

prova:
	@echo $(PREFIX)
	@echo $(POSTFILE)
	@echo $(POSTCODE)
	@echo $(SIMDIR)
	@echo $(MODDIR)
#	@echo 'Compiling' $(<F)
#eseguo : 
#	@echo $(SC_SECTIONS)
#	@echo $(SECSEC)
#	@$(foreach section,$(SC_SECTIONS),$(call pbs_templ,$(section)))

pppoinc : 
#	@$(foreach section,$(SC_SECTIONS),$(call pbs_templ,$(section)))
	@$(call pbs_poinc,$(section))

eseguo : build
	@echo $(SC_SECTIONS)
	@echo $(NSEC)
	@echo $(SECSEC_LESS)
	@echo $(SC_SECTIONS_LESS)
	@echo $(BATCHFILE)
	@echo "#!/bin/bash" >> $(BATCHFILE)
	@echo "#PBS -l select=1:mpiprocs=$(SC_NUM_THREADS):ncpus=$(SC_NUM_THREADS):mem=100mb:host=rat4" >> $(BATCHFILE)
	@echo "#PBS -l walltime=5:05:05" >> $(BATCHFILE)
	@echo "#PBS -M $(USER)@igi.cnr.it" >> $(BATCHFILE)
	@echo "#PBS -j oe" >> $(BATCHFILE)
	@echo "#PBS -m ae" >> $(BATCHFILE)
	@echo "#PBS -r n" >> $(BATCHFILE)
	@echo "#PBS -q lm" >> $(BATCHFILE)
	@echo "#PBS -o $(ERRORFILE)" >> $(BATCHFILE)
	@echo export OMP_NUM_THREADS=$(SC_NUM_THREADS) >> $(BATCHFILE)
	@echo module load intel/pe-xe-2018 >> $(BATCHFILE)
	@echo cd $(SIMDIR) >> $(BATCHFILE)
#ifeq ($(SECSEC),1)
#	@echo "ciao"
#else
#	@echo "nisba"
#endif
	@$(foreach section,$(SC_SECTIONS_LESS),$(call pbs_templ2,$(section))) $(call pbs_templ3,$(NSEC)) >> $(BATCHFILE)
	echo "Sottometto il job per il run di specyl"

#!#!#!# Marco, parte in cui compilo il codice 
OBJECTS = $(MODULES:.for=.o)
MODDS = $(MODULES:.for=.mod)
EXEC = $(BINDIR)/specyl.xxx
#POINCEXEC = $(BINDIR)/pppoincare_readold.xxx
POINCEXEC = $(BINDIR)/pppoincare.xxx
POINCEXP = $(BINDIR)/exp_poinc.xxx
ENDEXEC = $(BINDIR)/write_end.xxx
CODE = $(FORDIR)/specyl.for
POINCCODE = $(FORDIR)/pppoincare.F
EXPPOINCCODE = $(FORDIR)/exp_poinc.F
POINCPYTHONEXP = $(BINDIR)/exp_poinc_v2.xxx
EXPPOINCPYTHONCODE = $(FORDIR)/exp_poinc_v2.F
POSTCODE = $(FORDIR)/postproc.F #!# code for postprocessing
ENDCODE = $(FORDIR)/write_end.F #!# code for writing the specyl_end_file
#VPATH = $(MODDIR)/


compile : $(EXEC) 
	@echo $(CODE)
#

$(EXEC) : $(OBJECTS) $(CODE)
	@echo $(MODDIR)
	@echo $(OBJECTS)
#	$(FC) $(DEBUG) $(SKYLAKE) $(OMP) -r8 $(OBJECTS) -I$(MODDIR) -o $(BINDIR)/specyl.xxx $(FORDIR)/specyl.for
	$(FC) $(DEBUG) $(OMP) -r8 $(OBJECTS) -I$(MODDIR) -o $(BINDIR)/specyl.xxx $(FORDIR)/specyl.for
#!# Marco, per gfortran
#	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -I$(MODDIR) -o $(BINDIR)/specyl.xxx $(FORDIR)/specyl.for

poinc_compile : $(OBJECTS) $(POINCEXEC) 

$(POINCEXEC) : $(OBJECTS) $(POINCCODE)
	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -I$(MODDIR) -I$(MODDIR2) -o $(BINDIR)/pppoincare.xxx $(FORDIR)/pppoincare.F

exppoinc_compile : $(POINCEXP) 

$(POINCEXP) : $(OBJECTS) $(EXPPOINCCODE)
#	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -r8 -I$(MODDIR) -o $(BINDIR)/exp_poinc_v2.xxx $(FORDIR)/exp_poinc_v2.F
	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -r8 -I$(MODDIR) -o $(BINDIR)/exp_poinc.xxx $(FORDIR)/exp_poinc.F

exppython_compile : $(POINCPYTHONEXP) 

$(POINCPYTHONEXP) : $(OBJECTS) $(EXPPOINCPYTHONCODE)
	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -r8 -I$(MODDIR) -o $(POINCPYTHONEXP) $(EXPPOINCPYTHONCODE)

#;# lavoro per Italo
toksol_compile : toksol

toksol : $(OBJECTS) $(FORDIR)/exp_toksol.F
	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -r8 -I$(MODDIR) -o $(BINDIR)/exp_toksol.bin $(FORDIR)/exp_toksol.F
###HOWTO
### make -f Makefile toksol_compile
### ./bin/exp_toksol.bin

#;# Poincare di current density
jpoinc_compile : jpoinc

jpoinc : $(OBJECTS) $(FORDIR)/jpoinc.F
	$(FC) $(DEBUG) $(OMP) $(OBJECTS) -r8 -I$(MODDIR) -o $(BINDIR)/jpoinc.bin $(FORDIR)/jpoinc.F


post_compile : $(POSTEXEC) 

$(POSTEXEC) : $(POSTCODE)
#	$(FC) $(DEBUG) $(SKYLAKE) -o $(BINDIR)/postproc.xxx $(FORDIR)/postproc.F
#	#!# Marco, tolgo la keywork SKYLAKE perché non serve in questo frangente
	$(FC) $(DEBUG) -o $(BINDIR)/postproc.xxx $(FORDIR)/postproc.F

post_process : $(POSTEXEC)
	@echo $(POSTFILE)
	@echo $(SETDIR)
#	@echo $(SC_SECTIONS_LESS)
	@rm -f $(POSTFILE) #!# to remove the old script
	@echo '#!/bin/bash' >> $(POSTFILE)
	@echo cd $(SIMDIR) >> $(POSTFILE)
	@cat $(SETDIR)/postproc.sh >> $(POSTFILE)
	@chmod +x $(POSTFILE)
	bash $(POSTFILE)
#
##	@echo cd $(SIMDIR) >> $(POSTFILE)
##	@$(foreach section,$(SC_SECTIONS),$(call post_process_time,$(section))) >> $(POSTFILE)
##	@echo $(foreach section,$(SC_SECTIONS_LESS),$(POSTBIN) $(section) "&&")	$(POSTBIN) $(NSEC) >> $(POSTFILE)

#!# Marco, part to deal with the writing of specyl_?_end.dat files
end_compile : $(ENDEXEC) 

$(ENDEXEC) : $(ENDCODE)
	$(FC) $(DEBUG) -o $(ENDEXEC) $(FORDIR)/write_end.F

end_process : $(ENDEXEC)
	@echo $(SC_SECTIONS_LESS)
	@rm -f $(ENDFILE) #!# to remove the old script
	@echo "#!/bin/bash" >> $(ENDFILE)
	@echo cd $(SIMDIR) >> $(ENDFILE)
	@echo $(foreach section,$(SC_SECTIONS_LESS),$(ENDBIN) $(section) "&&")	$(ENDBIN) $(NSEC) >> $(ENDFILE)
	@chmod +x $(ENDFILE)
	@bash $(ENDFILE)

#!# example:
#!# make 


#$(MODULES) : $(OBJECTS)
#	@echo 'Compiling' $(<F)
#	$(FC) $(DEBUG) $(OMP) -c $< -o $(MODDIR)/$(@F) 

%.o : %.for 
	@echo 'Compiling' $(<F)
#The file-within-directory part of the file name of the target. If the value of ‘$@’ is dir/foo.o then ‘$(@F)’ is foo.o. ‘$(@F)’ is equivalent to ‘$(notdir $@)’.
	@echo $(@F)
#	$(FC) $(DEBUG) $(SKYLAKE) $(OMP) -r8 -c $< -o $(MODDIR)/$(@F) 

	$(FC) $(DEBUG) $(OMP) -r8 -c $< -o $(MODDIR)/$(@F) 
# Marco, per gfortran che non vuole -r8
#	$(FC) $(DEBUG) $(OMP) -c $< -o $(MODDIR)/$(@F) 
# -module $(MODDIR)

clean :
	rm -f $(OBJECTS)
	rm -f $(MODDS)

#!#!#!#!#! Marco, parte in cui costruisco la struttura della cartella della simulazione $(SIM) con una struttura standard con una sola sezione


provo :
	@echo $(NSEC)
	@echo $(SECSEC)

#build : tree spectrum spec2 spec3 cpspec 
#!# Marco, aprile 2018, scrivo una procedura automatica dentro specyl per avere lo spectrum.dat giusto nella cartella della simulazione
build : tree spec2 spec3 cpspec 


#-include $(SIMDIR)settings.inc.mk
tree : 
	@ mkdir -p $(SIMLOGDIR)
	@ mkdir -p $(SIMDATDIR)
	@ mkdir -p $(SIMBINDIR)
	@ mkdir -p $(SIMFORDIR)
	@echo "SC_SECTIONS = "$(SECSEC) >> $(SIMDIR)/settings.inc.mk
	@cp for/settings/settings.blc.for $(SIMFORDIR)
	@cp for/modules/nmt.in $(SIMFORDIR)
	@cp for/settings/nettings.in $(SIMFORDIR)
	@echo "Numero sezioni create: " $(SC_SECTIONS)
	@for x in $(SC_SECTIONS) ; do \
	 mkdir -p $(SIMFORDIR)/$${x} ; \
	done; \

#spec2 :
#	@if [ "$(3D)" = "0" ]; then \
#	 cp for/modules/spc2_2d.in $(SIMFORDIR)/spc2.in ; \
#	 else cp for/modules/spc2_3d.in $(SIMFORDIR)/spc2.in ;\
#	 fi; \

spec2 :
	@if [ "$(3D)" = "0" ]; then \
	  if [ "$(TOK)" = "1" ]; then \
	   cp for/modules/spc2_2dtok.in $(SIMFORDIR)/spc2.in ; \
	   else cp for/modules/spc2_2d.in $(SIMFORDIR)/spc2.in ;\
          fi; \
	 else cp for/modules/spc2_3d.in $(SIMFORDIR)/spc2.in ;\
	 fi; \

#spec3 :
#	@if [ "$(3D)" = "0" ]; then \
#	 cp for/modules/spc3_2d.in $(SIMFORDIR)/spc3.in ; \
#	 else cp for/modules/spc3_3d.in $(SIMFORDIR)/spc3.in ;\
#	 fi; \

spec3 :
	@if [ "$(3D)" = "0" ]; then \
	  if [ "$(TOK)" = "1" ]; then \
	   cp for/modules/spc3_2dtok.in $(SIMFORDIR)/spc3.in ; \
	   else cp for/modules/spc3_2d.in $(SIMFORDIR)/spc3.in ; \
	  fi; \
	 else cp for/modules/spc3_3d.in $(SIMFORDIR)/spc3.in ;\
	 fi; \

cpspec : spec2 spec3
	@for x in $(SC_SECTIONS) ; do \
	 cp $(SIMFORDIR)/spc2.in $(SIMFORDIR)/$${x}/spc2.in ; \
	 cp $(SIMFORDIR)/spc3.in $(SIMFORDIR)/$${x}/spc3.in ; \
	done; 
	@rm -f $(SIMFORDIR)/spc2.in
	@rm -f $(SIMFORDIR)/spc3.in

spectrum :
	@if [ "$(3D)" = "0" ]; then \
	  if [ "$(TOK)" = "1" ]; then \
	  cp for/modules/spectrumtok.sav $(SIMDATDIR)/spectrum.sav ; \
	  else cp for/modules/spectrum2d.sav $(SIMDATDIR)/spectrum.sav ; \
	  fi; \
	 else cp for/modules/spectrum3d.sav $(SIMDATDIR)/spectrum.sav ;\
	 fi; \

build_clean:
	@rm -fr $(SIMFORDIR)/*/*.in
	@rm -fr $(SIMLOGDIR)/*
	@rm -fr $(SIMBINDIR)/*

#!#!#! Marco, clean part
clean_simdir :
	@echo "pulisco la directory"
	@$(foreach AAA, $(MODNAMES), $(foreach SECTION, $(SC_SECTIONS), rm -f $(SIMFORDIR)/$(SECTION)/$(AAA).o ${nl}))
	@$(foreach AAA, $(MODNAMES), $(foreach SECTION, $(SC_SECTIONS), rm -f $(SIMFORDIR)/$(SECTION)/$(AAA).mod ${nl}))
	@rm -f $(SIMLOGDIR)/*
	@rm -f $(SIMDATDIR)/imf*sav
	@rm -f $(SIMDATDIR)/bf*sav
	@rm -f $(SIMDATDIR)/*energy.sav
	@rm -f $(SIMDATDIR)/specyl_?_all.sav
	@rm -f $(SIMDATDIR)/*eps
	@rm -f $(SIMDATDIR)/*png
	@rm -f $(SIMDATDIR)/*pdf
#

#!# Marco, IDL part

#;# Marco, deprecated, slow. Launch "end_process" instead
#i_end :
#	@echo "Create the specyl end file"
#	@cd $(PRODIR);				\
#	echo "i_writeend,'$(SIMDIR)',$(3D) & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup

i_change :
	@echo "Create the specyl end file with changed spectrum"
	@cd $(PRODIR);				\
	echo "i_changeend,'$(SIMDIR)',$(3D) & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

i_change_interp :
	@echo "Create the specyl end file with changed spectrum in a finer radial mesh to be specified in the source file"
	@cd $(PRODIR);				\
	echo "i_interpolateend,'$(SIMDIR)',$(3D) & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

post_begin22:
	@echo "ireadsc_fort22,'$(SIMDIR)','$(NUMBER_OF_FILES)',$(MF),$(3D) & exit"
	@cd $(PRODIR);				\
	echo "ireadsc_fort22,'$(SIMDIR)','$(NUMBER_OF_FILES)',$(MF),$(3D) & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

post_begin: post_process
	@echo "ireadsc_fort,'$(SIMDIR)','$(NUMBER_OF_FILES)',$(MF),$(3D) & exit"
	@cd $(PRODIR);				\
	echo "ireadsc_fort,'$(SIMDIR)','$(NUMBER_OF_FILES)',$(MF),$(3D) & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "   "
	@echo "Unisco i imf_bprofiles_?.sav in idl e produco imf_bprofiles.sav"
	@cd $(PRODIR);				\
	echo "imergebmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "   "
	@echo "Unisco i imf_vprofiles_?.sav in idl e produco imf_vprofiles.sav"
	@cd $(PRODIR);				\
	echo "imergevmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "Salvo il campo magnetico all'edge"
	@cd $(PRODIR);				\
	echo "i_saveedge,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

begin: 
	@echo "ireadsc2,'$(SIMDIR)','$(NUMBER_OF_FILES)',$(MF),$(3D) & exit"
	@cd $(PRODIR);				\
	echo "ireadsc2,'$(SIMDIR)','$(NUMBER_OF_FILES)',$(MF),$(3D) & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "   "
#	@echo "Unisco i imf_bprofiles_?.sav in idl e produco imf_bprofiles.sav"
#	@cd $(PRODIR);				\
#	echo "imergebmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup
#	@echo "   "
#	@echo "Unisco i imf_vprofiles_?.sav in idl e produco imf_vprofiles.sav"
#	@cd $(PRODIR);				\
#	echo "imergevmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup
#	@echo "Salvo il campo magnetico all'edge"
#	@cd $(PRODIR);				\
#	echo "i_saveedge,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup

#!# Marco, IDL plot part

#........................................   fields plot

i_plot: post_begin
	@echo "Plotting b,v in mf notation"
	@echo "Remind to set 3D={0,1}"
	@echo $(3D), $(RDS)
	@cd $(PRODIR);				\
	echo "plot_bmf,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@cd $(PRODIR);				\
	echo "plot_edge,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@cd $(PRODIR);				\
	echo "plot_vmf,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "Plotting magnetic field energy in time"
	@cd $(PRODIR);				\
	echo "plot_enb,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup


#........................................   equilibrium profiles in idl

i_eq :
	@echo "   "
	@echo "Unisco i bv_eq_?.sav in idl e produco bv_eq.sav"
	@cd $(PRODIR);				\
	echo "isaveeqprofiles,'$(SIMDIR)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
maj  :
	@echo "   "
	@echo "major radius"
	@cd $(PRODIR);				\
	echo "prova_mj,'$(SIMDIR)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

i_qprof :
	@echo "Save axisymmetric q profile in time"
	@cd $(PRODIR);				\
	echo "qprof,'$(SIMDIR)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   v profiles in idl

i_savea : 
	@echo "Computing the vector potential"
	@cd $(PRODIR);				\
	echo "idl_savea,'$(SIMDIR)','$(TT)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

i_savechi3 : 
	@echo "Computing the helical flux function"
	@cd $(PRODIR);				\
	echo "idl_savechi_v3,'$(SIMDIR)','$(3D)','$(TOK)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#;# EXAMPLE
#make -f Makefile SIM=archive/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/ 3D=1 TOK=0 HEL=-11 TMIN=1000 TMAX=1000 DELTA=5 NTH=128 NZZ=256 i_savechi3
#make -f Makefile SIM=archive/tokamak/3d/s1000m3e4_br1e-4m2n1_plus29/ 3D=1 TOK=1 HEL=-1 NTH=64 NZZ=64 TMIN=9000 TMAX=9200 DELTA=5 i_savechi3

i_mfv :
	@echo "   "
	@echo "Unisco i imf_vprofiles_?.sav in idl e produco imf_vprofiles.sav"
	@cd $(PRODIR);				\
	echo "imergevmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   b profiles in idl

i_mfb :
	@echo "   "
	@echo "Unisco i imf_bprofiles_?.sav in idl e produco imf_bprofiles.sav"
	@cd $(PRODIR);				\
	echo "imergebmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup


i_modeb :
	@echo "   "
	@echo "Unisco i imf_bprofiles_?.sav in idl e produco imf_bprofiles.sav"
	@cd $(PRODIR);				\
	echo "isavemodeprofiles,'$(SIMDIR)','$(MF)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup


i_edge:
	@echo "Salvo il campo magnetico all'edge"
	@cd $(PRODIR);				\
	echo "i_saveedge,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup


#........................................   average in time of b profiles 
i_averageb:
	@echo "Salvo la media dei profili delle autofunzioni di br"
	@cd $(PRODIR);				\
	echo "average_profiles,'$(SIMDIR)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   profiles in idl

i_mf :
	@echo "   "
	@echo "Unisco i imf_bprofiles_?.sav in idl e produco imf_bprofiles.sav"
	@cd $(PRODIR);				\
	echo "imergebmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "   "
	@echo "Unisco i imf_vprofiles_?.sav in idl e produco imf_vprofiles.sav"
	@cd $(PRODIR);				\
	echo "imergevmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   energy time dependece in idl

i_energy :
	@echo "   "
	@echo "Calcolo energia cinetica e magnetica nel tempo"
	@cd $(PRODIR);				\
	echo "idl_energy,'$(SIMDIR)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   current profiles in idl

i_jprofiles : 
	@echo "Computing ijprofiles.sav reading the magnetic field"
	@cd $(PRODIR);				\
	echo "isavejprofiles,'$(SIMDIR)','$(MF)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

i_vortprofiles : 
	@echo "Computing iomprofiles.sav (vorticity) reading the magnetic field"
	@cd $(PRODIR);				\
	echo "isaveomprofiles,'$(SIMDIR)','$(MF)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

i_writejfort : 
	@echo "Writing ijprofiles.sav in fortran format"
	@cd $(PRODIR);				\
	echo "iwritejprofiles,'$(SIMDIR)','$(MF)','$(3D)',$(TT) & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   current-related quantities in idl

i_jpar : 
	@echo "Computing jpar and j2 reading the magnetic field or current field"
	@cd $(PRODIR);				\
	echo "plot_etaj2,'$(SIMDIR)','$(3D)','$(TOK)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(HEL)','$(DELTA)','$(MYMIN)','$(MYMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE make -f Makefile SIM=archive/rfp/2d/theta16/h10/s3e3m3e3/ 3D=0 NTH=64 NZZ=256 TOK=0 TMIN=-1 TMAX=-1 DELTA=1 MYMIN=0 MYMAX=4 i_jpar
#!#make -f Makefile SIM=archive_flow/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/ 3D=1 TOK=0 NTH=128 NZZ=512 TMIN=2250 TMAX=2350 HEL=-7 DELTA=2 MYMIN=0 MYMAX=2 i_jpar

corsomhd : 
	@echo "Figure corso MHD"
	@cd $(PRODIR);				\
	echo "plot_corso,'$(SIMDIR)','$(3D)','$(TOK)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(HEL)','$(DELTA)','$(MYMIN)','$(MYMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE make -f Makefile SIM=archive/rfp/2d/theta16/h10/s3e3m3e3/ 3D=0 NTH=64 NZZ=256 TOK=0 TMIN=-1 TMAX=-1 DELTA=1 MYMIN=0 MYMAX=4 corso_mhd
#!#make -f Makefile SIM=archive_flow/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/ 3D=1 TOK=0 NTH=128 NZZ=512 TMIN=2250 TMAX=2350 HEL=-7 DELTA=2 MYMIN=0 MYMAX=2 corso_mhd

#........................................   electric field profiles in idl

i_el : 
	@echo "Computing el.sav (electric field)"
	@cd $(PRODIR);				\
	echo "i_el,'$(SIMDIR)','$(D3)','$(MF)','$(EQU)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   momentum in idl

i_momentum : 
	@echo "Computing momentum.sav reading the velocity field"
	@cd $(PRODIR);				\
	echo "momentum,'$(SIMDIR)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   save electromagnetic torque at a specific time
i_tau : 
	@echo "Compute torque "
	@cd $(PRODIR);				\
	echo "tau,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   save integrated electromagnetic torque z direction
i_integtau : 
	@echo "Compute integrated electromagnetic torque z direction"
#	@mkdir -p $(SIMDIR)'dat/itp/'$(ITP)'/'
	@cd $(PRODIR);				\
	echo "integrated_tau,'$(SIMDIR)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   angular momentum in idl

i_lprofiles.sav : 
	@echo "Computing ilprofiles.sav reading the velocity field"
	@echo "Remind to set TMIN, TMAX"
	@cd $(PRODIR);				\
	echo "isavelprofiles,'$(SIMDIR)','$(MF)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   magnetic helicity

i_helicity: 
	@echo "magnetic helicity in time"
	@cd $(PRODIR);				\
	echo "idl_helicity,'$(SIMDIR)','$(3D)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   F-Theta in time in idl

i_savef : 
	@echo "Computing F-Theta value in time reading the magnetic field"
	@cd $(PRODIR);				\
	echo "idl_savef,'$(SIMDIR)','$(3D)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   magnetic field in modulus&phase notation

i_saveeq :
	@echo "Read specyl_?_all.dat and save the equilibrium part of the fields"
	@cd $(PRODIR);				\
	echo "ireadeq,'$(SIMDIR)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   save fields to be read by SpeCyl

i_savespecyl : 
	@echo "Write on disk *changed* fields to be read by SpeCyl"
	@cd $(PRODIR);				\
	echo "idl_savespecyl,'$(SIMDIR)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup


i_antifft :
	@echo "Reconstruct fields in space"
	@echo "antifft_reconst_v3,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit";
	@echo "Remind to set the field, as in $(FLD)"
	@echo $(FLD)
	@cd $(PRODIR);				\
	echo "antifft_reconst_v3,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#;# EXAMPLE   make -f Makefile SIM=archive/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/ FLD=b 3D=1 HEL=-7 NTH=128 NZZ=1024 TMIN=1200 TMAX=1200 DELTA=1 i_antifft
i_antifftguscio :
	@echo "Reconstruct fields in space"
	@echo "antifft_guscio,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit";
	@echo "Remind to set the field, as in $(FLD)"
	@echo $(FLD)
	@cd $(PRODIR);				\
	echo "antifft_guscio,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#; example  make -f Makefile SIM=archive/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/ 3D=1 NTH=128 NZZ=512 TMIN=2230 TMAX=2232 DELTA=1 FLD=b i_antifftguscio
#;# make -f Makefile SIM=archive_flow/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/ 3D=1 NTH=128 NZZ=512 TMIN=1000 TMAX=2390 DELTA=1 FLD=b i_antifftguscio

#	echo "antifft_reconst_v2,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
#
tok_antifft :
	@echo "Reconstruct fields in space"
	@echo "antifft_reconst_v2,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit";
	@echo "Remind to set the field, as in $(FLD)"
	@echo $(FLD)
	@cd $(PRODIR);				\
	echo "antifft_reconst_tok,'$(SIMDIR)','$(FLD)','$(3D)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

flowxy :
	@echo "Plot flow"
	@cd $(PRODIR);				\
	echo "flow_xy,'$(SIMDIR)','$(3D)','$(TOK)','$(HEL)','$(TMIN)','$(TMAX)','$(DELTA)','$(MYMIN)','$(MYMAX)', '$(NTH)','$(NZZ)', '$(ZED)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE make -f Makefile SIM=archive/tokamak/2d/v2_s1000m3e4_br1e-4_plus29_Q4e5/ 3D=0 TOK=1 TMIN=-1 TMAX=-1 DELTA=1 HEL=0.5 MYMIN=1 MYMAX=4 flowxy
#!# EXAMPLE RFP make -f Makefile SIM=archive/rfp/2d/theta16/momsour/h35/s100m3e4/br1e-5Q3.33e5/ 3D=0 TOK=0 TMIN=-1 TMAX=-1 DELTA=1 HEL=-35 MYMIN=1 MYMAX=4 flowxy
#!# example make -f Makefile SIM=archive/rfp/3d/theta16/s100m100/v2_nmp7_mp6e-3_qinfty/ 3D=0 TOK=0 TMIN=2075 TMAX=2075 DELTA=1 HEL=-7 MYMIN=0 MYMAX=2 flowxy

i_timechihel :
	@echo "Reconstruct helical flux function in time"
	@echo $(HEL)
	@cd $(PRODIR);				\
	echo "time_chihel,'$(SIMDIR)','$(3D)','$(TOK)','$(HEL)','$(TMIN)','$(TMAX)','$(DELTA)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE: make -f Makefile SIM=archive/tokamak/2d/Q5e5s1000m3e4_vario_brext/br1e-2/ 3D=0 TOK=1 HEL=0.5 TMIN=980 TMAX=990 DELTA=5 i_timechihel
#!# EXAMPLE: make -f Makefile SIM=archive/rfp/2d/theta16/momsour/h35/s100m3e4/br1e-5Q8e6/ 3D=0 TOK=0 HEL=35 TMIN=-1 TMAX=-1 DELTA=1 i_timechihel

chihel3d :
	@echo "Reconstruct helical flux function in time"
	@echo $(HEL)
	@cd $(PRODIR);				\
	echo "chihel3d,'$(SIMDIR)','$(3D)','$(TOK)','$(HEL)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE: make -f Makefile SIM=archive/tokamak/2d/Q5e5s1000m3e4_vario_brext/br1e-2/ 3D=0 TOK=1 HEL=0.5 NTH=64 NZZ=64 TMIN=5000 TMAX=990 DELTA=5 chihel3d
#!# EXAMPLE: make -f Makefile SIM=archive/rfp/2d/theta16/momsour/h35/s100m3e4/br1e-5Q8e6/ 3D=0 TOK=0 HEL=35 TMIN=-1 TMAX=-1 DELTA=1 

show_timechihel :
	@echo "Reconstruct helical flux function in time"
	@cd $(PRODIR);				\
	echo "show_time_chihel,'$(SIMDIR)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE: make -f Makefile SIM=archive/tokamak/2d/Q5e5s1000m3e4_vario_brext/br1e-2/ show_timechihel

displacement :
	@echo "Reconstruct displacement"
	@echo "REMIND set displacement,'$(SIMDIR)','$(3D)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)'"
	@cd $(PRODIR);				\
	echo "displacement,'$(SIMDIR)','$(3D)','$(NTH)','$(NZZ)','$(TMIN)','$(TMAX)','$(DELTA)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE:  make -f Makefile SIM=archive/rfp/3d/theta16/s100m100/base_fsg/ 3D=1 NTH=128 NZZ=256 TMIN=300 TMAX=300 DELTA=1 displacement

#;# Marco, 15 maggio 2019
#;# Marco, 30 marzo 2022, unisco con vperp
flow_shear :
	@echo "Reconstruct flow-shear"
	@cd $(PRODIR);				\
	echo "flow_shear,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(DELTA)', '$(NTH)', '$(NZZ)', '$(MYMIN)', '$(MYMAX)' & exit" > idlstartup;\
	echo "vperp,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(DELTA)', '$(NTH)', '$(NZZ)', '$(MYMIN)', '$(MYMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE: 
#!# make -f Makefile SIM=archive/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/ 3D=1 TOK=0 FLD=v TMIN=1300 TMAX=1700 DELTA=10 flow_shear

#;# Marco, 31 maggio 2021
vperp :
	@echo "Reconstruct perpendicular and parallel components of velocity and electric field"
	@cd $(PRODIR);				\
	echo "vperp,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(DELTA)', '$(NTH)', '$(NZZ)', '$(MYMIN)', '$(MYMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#!# EXAMPLE: 
#!# make -f Makefile SIM=archive/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/ 3D=1 TOK=0 FLD=v TMIN=1300 TMAX=1700 DELTA=10 vperp
#;# Marco, plot part
#........................................   velocity field plot

i_plotv: 
	@echo "Plotting v in mf notation"
	@echo "Remind to set 3D={0,1}"
	@echo $(3D)
	@cd $(PRODIR);				\
	echo "plot_vmf,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   magnetic field plot

i_plotb: 
	@echo "Plotting b in mf notation"
	@echo "Remind to set 3D={0,1}"
	@echo $(SIMDIR)
	@echo $(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)','$(TOK)
	@cd $(PRODIR);				\
	echo "plot_bmf,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#	echo "plot_bmf2,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
#	@cd $(PRODIR);				\
#	echo "plot_edge2,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup

i_tok3d: post_begin
	@echo "Plotting b in mf notation"
	@echo "Remind to set 3D={0,1}"
	@echo $(3D)
	@cd $(PRODIR);				\
	echo "plot_tok3d,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@cd $(PRODIR);				\
	echo "plot_enbtok3d,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#........................................   fields plot

#........................................   fields plot2

i_plot2: ##begin 
	@echo "Plotting b,v in mf notation"
	@echo "Remind to set 3D={0,1}"
#	@echo $(3D)
#	@echo "Unisco i imf_vprofiles_?.sav in idl e produco imf_vprofiles.sav"
#	@cd $(PRODIR);				\
#	echo "imergevmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup
#	@cd $(PRODIR);				\
#	echo "imergebmf,'$(SIMDIR)','$(3D)','$(MF)' & exit" > idlstartup;\
#	/usr/local/exelis/idl85/bin/idl idlstartup
	@cd $(PRODIR);				\
	echo "plot_bmf,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@cd $(PRODIR);				\
	echo "plot_edge,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@cd $(PRODIR);				\
	echo "plot_vmf,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)','$(RDS)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
	@echo "Plotting magnetic field energy in time"
	@cd $(PRODIR);				\
	echo "plot_enb,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   edge magnetic field plot

i_edgeb: 
	@echo "Plotting b in mf notation"
	@echo "Remind to set 3D={0,1}"
	@echo $(3D)
	@cd $(PRODIR);				\
	echo "plot_edge,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   magnetic energy plot

i_enb: 
	@echo "Plotting magnetic field energy in time"
	@echo "Remind to set 3D={0,1}"
	@echo "3D= "$(3D)
	@cd $(PRODIR);				\
	echo "plot_enb,'$(SIMDIR)','$(3D)','$(TMIN)','$(TMAX)'  & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   magnetic field eigenfunctions at specific time

i_plottimeb: 
	@echo "Plotting b eigenfunctions at t="$(TT)
	@echo "Remind to set TT=?? and the wave number of the mode MM=? NN=?"
	@cd $(PRODIR);				\
	echo "idl_pltb,'$(SIMDIR)','$(TT)','$(MM)','$(NN)','$(3D)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
#;#   make -f Makefile SIM=archive/rfp/2d/theta16/momsour/h10/s10m3e3/v2_br0_Qinfty/ 3D=0 TOK=0 MM=1 NN=-10 TT=-1 i_plottimeb

#........................................   velocity field eigenfunctions at specific time

i_plottimev: 
	@echo "Plotting v eigenfunctions at t="$(TT)
	@echo "Remind to set TT=?? and the wave number of the mode MM=? NN=?"
	@cd $(PRODIR);				\
	echo "idl_pltv,'$(SIMDIR)','$(TT)','$(MM)','$(NN)','$(3D)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
# make -f Makefile SIM=archive/rfp/3d/theta16/s100m10/base_perturbazionerandomphases/ TT=-1 MM=0 NN=0 i_plottimev

#........................................   current field eigenfunctions at specific time

i_plottimej: 
	@echo "Plotting j eigenfunctions at t="$(TT)
	@echo "Remind to set TT=?? and the wave number of the mode MM=? NN=?"
	@cd $(PRODIR);				\
	echo "idl_pltj,'$(SIMDIR)','$(TT)','$(MM)','$(NN)','$(3D)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#........................................   safety factor at specific time

i_plotq: 
	@echo "Plotting axysimemtric safety factor at t="$(TT)
	@echo "Remind to set TT=??"
	@cd $(PRODIR);				\
	echo "idl_pltq,'$(SIMDIR)','$(TT)','$(3D)','$(TOK)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup
# make -f Makefile SIM=archive/rfp/3d/theta16/s100m10/base_perturbazionerandomphases/ TT=-1 i_plotq

#........................................   compare magnetic field in two simulations

SIMDIR0=$(PREFIX)"/"$(SIM0)
SIMDIR1=$(PREFIX)"/"$(SIM1)
i_compare: 
	@echo "Comparing sim1 and sim2 b-field"
	@echo "Remind to set 3D={0,1}"
	@echo $(3D)
	@echo $(SIM0)
	@echo $(SIMDIR0)
	@cd $(PRODIR);				\
	echo "compare,'$(SIMDIR0)','$(SIMDIR1)','$(3D)','$(TMIN)','$(TMAX)' & exit" > idlstartup;\
	/usr/local/exelis/idl85/bin/idl idlstartup

#make -f Makefile SIM0=archive/rfp/3d/theta16/s100m10/mp0q500/ SIM1=archive/rfp/3d/theta16/s100m10/base_perturbazionerandomphases/ i_compare
