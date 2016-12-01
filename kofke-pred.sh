# THis is for prediction part
# $1 : new temp.
# $2 : old folder to get files
# $3 : new pressure
# $4 : new temperature
# Example: ./kofke.sh Tr0.94 Tr0.92 32.87335142 144.8256
text_slurm="mpirun -np 48 /home/hjung52/gromacs5.1.2/bin/gmx_mpi mdrun -pin on -v >& log"
if [ ! -d "$1" ]; then
  mkdir $1
  mkdir -p $1/{gas,liq}-p
fi

for PHASE in gas liq
do

cd $1/$PHASE-p
cp ../../$2/$PHASE-c/confout.gro conf.gro
cp ../../$2/$PHASE-c/grompp.mdp ./
head -6 ../../$2/$PHASE-c/$2-$PHASE.slurm > $1-$PHASE.slurm

workdir=$(pwd)

echo "cd ${workdir}" >> $1-$PHASE.slurm
echo "${text_slurm}" >> $1-$PHASE.slurm
echo "exit 0" >> $1-$PHASE.slurm

sed -i '/ref_t                    =/c\ref_t                    = '"$4"' ;' grompp.mdp
sed -i '/ref-p                    =/c\ref-p                    = '"$3"' ;' grompp.mdp

if [ "$PHASE" = "liq" ]
then
        gmx_mpi grompp -p ../../../FF.Argon/RNP/topol-4k.top
else
        gmx_mpi grompp -p ../../../FF.Argon/White/topol-4k.top
fi
qsub $1-$PHASE.slurm

cd ../..

done
