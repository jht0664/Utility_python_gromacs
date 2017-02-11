# THis is for box-length and enthalpy
# $1 : top folder
# $2 : p or c (pred or corr)
# Example: ./kofke-calc.sh Tr0.94 p

for PHASE in gas liq
do

cd $1/$PHASE-$2
echo 10 16 | gmx_mpi energy
tail -2500 energy.xvg > calc.xvg
echo "box-x"
awk '{ sum += $2 }END{ print sum / NR }' calc.xvg
echo "Enthalpy"
awk '{ sum += $3 }END{ print sum / NR }' calc.xvg

cd ../../

done
