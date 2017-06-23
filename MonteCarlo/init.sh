#/bin/bash
# $1 = # molecules of A = # molecules of B
# $2 = density
# $3 = ratio of box-z
# $4 = foldername
# $5 = seperation or not (YES/NO)
# $6 = lattice insertion? (YES/NO)
# $7 = number fraction of A of one phase

if [ $5 = 'YES' ]
then
	name=$4-'sep'
else
	name=$4
fi

mkdir $name
cd $name

comm="python3 /home/htjung/NVT/initial.py -na $1 -nb $1 -d $2 -r $3"
space=" "
yes1="-sep YES -fr $7"
yes2='-mt 100000'
c=""
if [ $5 = 'YES' ]
then
	c=$comm$space$yes1
else
	c=$comm
fi
if [ $6 = 'NO' ]
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
qsub run.pbs
cd ..

