#/bin/bash
# $1 = # molecules of A
# $2 = # molecules of B
# $3 = density
# $4 = ratio of box-z
# $5 = foldername
# $6 = seperation or not (YES/NO)
# $7 = lattice insertion? (YES/NO)
# $8 = number fraction of A of one phase

if [ $6 = 'YES' ]
then
	name=$5-'sep'
else
	name=$5
fi

mkdir $name
cd $name

comm="python3 /home/htjung/NVT/initial.py -na $1 -nb $2 -d $3 -r $4"
space=" "
yes1="-sep YES -fr $8"
yes2='-mt 100000'
c=""
if [ $6 = 'YES' ]
then
	c=$comm$space$yes1
else
	c=$comm
fi
if [ $7 = 'NO' ]
then
	c=$c$space$yes2
else
	c=$c
fi

echo $c 

$c

mv init.ic composite.ic
cp ../../colloid.inp ./
ln -sf /home/htjung/NVT/colloid.x
cp ../../run.pbs ./
MYPWD=$(pwd)
pbsname=${MYPWD: -15}
sed -i "s#PWD#$MYPWD#g" run.pbs
sed -i "s#NAME#$pbsname#g" run.pbs
#qsub run.pbs
cd ..

