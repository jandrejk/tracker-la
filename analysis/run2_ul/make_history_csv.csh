#!/bin/csh -f

set subdet = $1
set trigger = L1MuOpen
set color = $2

set files = `find ./$trigger -name $subdet"_*.csv" | sort`

set outf = $subdet"_Run2_UL_"$trigger"_"$color".csv"
if ( -e $outf ) then
   rm -f $outf
endif

foreach f ( $files )
   set era = `echo $f | awk -F "/" '{print $4}'`
   set data = `cat $f | grep $color | awk -F ',' '{print $4","$5","$6","$7","$8","$9","$10","$11","$3}'`
   echo $era", "$data >> $outf
end

exit
