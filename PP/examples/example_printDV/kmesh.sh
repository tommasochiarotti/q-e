q1=$1; q2=$2;q3=$3
nk_irr=`grep "number of k " SrVO3.scf.out | awk '{print $5}'`
echo "K_POINTS crystal"
echo `echo $nk_irr*2 | bc`
for i in `seq 1 $nk_irr`; do
kpoint=`grep -A $i "cryst. coord." SrVO3.scf.out | tail -1 | awk '{printf "%16.12f %16.12f %16.12f \n", $5, $6, $7}' | sed -e "s/),//g"`
qpoint=`grep -A $i "cryst. coord." SrVO3.scf.out | tail -1 | awk -v q1=$q1 -v q2=$q2 -v q3=$q3 '{printf "%16.12f %16.12f %16.12f \n", $5+q1, $6+q2, $7+q3}' | sed -e "s/),//g"`
weight=`grep -A $i "cryst. coord." SrVO3.scf.out | tail -1 | awk '{printf "%16.12f \n", $10}'`
printf "%16.12f %16.12f %16.12f %16.12f \n" $kpoint $weight
printf "%16.12f %16.12f %16.12f %16.12f \n" $qpoint 0.0000
done
