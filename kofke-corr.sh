# THis is for correction part
# $1 : new temp.
# $2 : new pressure
# Example: ./kofke.sh Tr0.94 32.87335142
text_slurm="mpirun -np 48 /home/hjung52/gromacs5.1.2/bin/gmx_mpi mdrun -pin on -v >& log"
if [ ! -d "$1/gas-c" ]; then
  mkdir -p $1/{gas,liq}-c
fi

for PHASE in gas liq
do

cd $1/$PHASE-c
cp ../$PHASE-p/confout.gro conf.gro
cp ../$PHASE-p/grompp.mdp ./
head -6 ../$PHASE-p/$1-$PHASE.slurm > $1-$PHASE.slurm

workdir=$(pwd)

echo "cd ${workdir}" >> $1-$PHASE.slurm
echo "${text_slurm}" >> $1-$PHASE.slurm
echo "exit 0" >> $1-$PHASE.slurm

sed -i '/ref-p                    =/c\ref-p                    = '"$2"' ;' grompp.mdp

if [ "$PHASE" = "liq" ]
then
	gmx_mpi grompp -p ../../../FF.Argon/RNP/topol-4k.top
else 
	gmx_mpi grompp -p ../../../FF.Argon/White/topol-4k.top	
fi

qsub $1-$PHASE.slurm

cd ../..

done
