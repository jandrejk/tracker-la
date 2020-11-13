#!/bin/csh -f 

if ( $#argv < 1 ) then
   echo "nothing to be done"
   exit
endif

set os = "6"
if ( $#argv == 2 ) then
 set os = $2
endif

set cfg = $1
set rootdir = `cat $cfg | grep -v '#' | grep "CalibTreeDirectory" | awk -F "=" '{print $2}'`
set datatype = `basename $rootdir`
set basedir = `dirname $rootdir`
set baseanadir = `echo $basedir | awk -F "calibTrees" '{print $2}'`

set version  = `cat $cfg | grep '#' | grep "VERSION" | awk -F "=" '{print $2}'`

set ptmin  = `cat $cfg | grep -v '#' |  grep "ptMin"  | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`
set ptmax  = `cat $cfg | grep -v '#' |  grep "ptMax"  | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`
set etamin = `cat $cfg | grep -v '#' |  grep "etaMin" | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`
set etamax = `cat $cfg | grep -v '#' |  grep "etaMax" | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`

set hitmin  = `cat $cfg | grep -v '#' |  grep "nHitsValidMin" | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`
set hitmax  = `cat $cfg | grep -v '#' |  grep "nHitsValidMax" | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`
set chi2min = `cat $cfg | grep -v '#' |  grep "chi2ndofMin"   | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`
set chi2max = `cat $cfg | grep -v '#' |  grep "chi2ndofMax"   | awk -F "=" '{print $2}' | sed -e "s/\./p/g"`

if ( $ptmin != "" ) then
   set ptmin = "_ptmin"$ptmin
endif
if ( $ptmax != "" ) then
   set ptmin = "_ptmax"$ptmax
endif
if ( $etamin != "" ) then
   set etamin = "_etamin"$etamin
endif
if ( $etamax != "" ) then
   set etamax = "_etamax"$etamax
endif
if ( $hitmin != "" ) then
   set hitmin = "_hitmin"$hitmin
endif
if ( $hitmax != "" ) then
   set hitmax = "_hitmax"$hitmax
endif
if ( $chi2max != "" ) then
   set chi2max = "_chi2max"$chi2max
endif
if ( $chi2min != "" ) then
   set chi2min = "_chi2min"$chi2min
endif

if ( $version != "" ) then
   set version = "v"$version"/"
endif

set ananame = $datatype$ptmin$ptmax$etamin$etamax$hitmin$hitmax$chi2min$chi2max

set anadir = /nfs/dust/cms/user/walsh/cms/tracker/dpg/la$baseanadir/$ananame

mkdir -p $anadir
#ln -s $anadir

set files = `/bin/ls -1 $rootdir/calibTree_*_*.root`

set current_run = "000000"
foreach f ( $files )
   set run = `basename $f .root | awk -F "_" '{print $2}'`
   if ( $run == $current_run ) then
      continue
   else
      set current_run = $run
   endif
   
   if ( $run == "318517" || $run == "318519" || $run == "318520" || $run == "318521"  ) then
      continue
   endif
   set command = "SiStripLAMonitor -c $cfg -r $run" 
   set condor_dir = "$anadir/Condor_SiStripLAMonitor/job_run"$run
   if ( -d $condor_dir ) then
      echo "Skipping run $run. Directory $condor_dir already exists ..." 
      continue
   endif
   mkdir -p $condor_dir
   cp -p $cfg $condor_dir
   cd $condor_dir
   echo "Submitting Lorentz Angle job for run $run ..."
   condor_macro_submit.sh $command
   
   cd -
   
#   if ( $run == "320442" ) then
#      break
#   endif
end

exit
