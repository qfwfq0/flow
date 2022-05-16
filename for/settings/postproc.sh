#!# sec1
file_dat="dat/specyl_1_all.dat"
file_idl="dat/specyl_1_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_1_all.dat`
time_idl=`stat -c%Y dat/specyl_1_all.sav`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 1
else
echo "non rifare #1"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 1
fi;
else
echo "non esiste il file dat #1"
fi;

#!# sec2
file_dat="dat/specyl_2_all.dat"
file_idl="dat/specyl_2_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_2_all.dat`
time_idl=`stat -c%Y dat/specyl_2_all.sav`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 2
else
echo "non rifare #2"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 2
fi;
else
echo "non esiste il file dat #2"
fi;

#!# sec3
file_dat="dat/specyl_3_all.dat"
file_idl="dat/specyl_3_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_3_all.dat`
time_idl=`stat -c%Y dat/specyl_3_all.sav`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 3
else
echo "non rifare #3"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 3
fi;
else
echo "non esiste il file dat #3"
fi;

#!# sec4
file_dat="dat/specyl_4_all.dat"
file_idl="dat/specyl_4_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_4_all.dat`
time_idl=`stat -c%Y dat/specyl_4_all.sav`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 4
else
echo "non rifare #4"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 4
fi;
else
echo "non esiste il file dat #4"
fi;

#!# sec5
file_dat="dat/specyl_5_all.dat"
file_idl="dat/specyl_5_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_5_all.dat`
time_idl=`stat -c%Y dat/specyl_5_all.sav`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 5
else
echo "non rifare #5"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 5
fi;
else
echo "non esiste il file dat #5"
fi;

#!# sec6
file_dat="dat/specyl_6_all.dat"
file_idl="dat/specyl_6_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_6_all.dat`
time_idl=`stat -c%Y dat/specyl_6_all.sav`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 6
else
echo "non rifare #6"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 6
fi;
else
echo "non esiste il file dat #6"
fi;

#!# sec7
file_dat="dat/specyl_7_all.dat"
file_idl="dat/specyl_7_all.sav"
if [ -e ${file_dat} ];
then
if [ -e ${file_idl} ];
then
time_dat=`stat -c%Y dat/specyl_7_all.dat`
time_idl=`stat -c%Y dat/specyl_7_all.savt`
#echo ${time_dat}
#echo ${time_idl}
if [ ${time_dat} \> ${time_idl} ];
then
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 7
else
echo "non rifare #7"
fi;
else
/0/rfx/ricercatori/ft/specyl/veranda/flow/bin/postproc.xxx 7
fi;
else
echo "non esiste il file dat #7"
fi;
