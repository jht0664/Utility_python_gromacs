#!/bin/bash
# traverse subdirectories in a current directory
# with command-line $1 (ex. ~/Utility/n-prof.sh) 
dirlist=$(find $(pwd) -mindepth 1 -maxdepth 1 -type d)

for dir in $dirlist
do 
# remaining suffolders
	pushd $dir
# show current directory
	echo $dir
# do a command-line
	$1 $2 $3
# show current pwd
	popd 
done >traverse-dir.log

