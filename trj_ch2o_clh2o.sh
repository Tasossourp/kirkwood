#!/bin/bash

export trj='/lustre/project/m2_trr146/asourpis/paper_ccn_h2o_0.25_0.75/75/NVTrun100ns'
#export nccn=1216
export nh2o=1220

#for ((n=1217;n<=$nh2o;n+=1))
#do

#start=1
#end=500
#incr=25

#for ((i=$start;i<=$end;i+=incr))

#do
var=$r
#var=`ls -l|awk "BEGIN {print $rad/100}"` #this is a way to use real numbers in bash shell scripts

m=0
for ((t=1;t<=25000;t+=10))
do

gmx trajectory -f ''$trj'/sim_100ns_0.0-V_nm.xtc' -s ''$trj/'sim_100ns_0.0-V_nm.tpr' -selrpos whole_mol_com -seltype atom -b $t -e $t  -n 'index_'$nh2o'h2o_h2o_r'$var'.ndx' -ox 'central_'$nh2o'h2o_h2o_r'$var'_t'$t'.xvg' <<EOF
0
EOF

gmx trajectory -f ''$trj'/sim_100ns_0.0-V_nm.xtc' -s ''$trj/'sim_100ns_0.0-V_nm.tpr' -selrpos whole_mol_com -seltype atom -b $t -e $t  -n 'index_'$nh2o'h2o_h2o_r'$var'.ndx' -ox 'collect_'$nh2o'h2o_h2o_r'$var'_t'$t'.xvg' <<EOF
$(($m+1))
EOF

done

#done

#done






























