#!/bin/bash
# Create a new distribution version

# TODO before:
#   - Check version number in theo_header.py and doc/source/conf.py
#   - git tag set by this script

echo "Syntax: make_distrib.bash [Version]"
echo

if [ $# -eq 0 ]
then
   echo "Please enter version number!"
   exit 1
fi

SDIR=`pwd`

LTDIR="MATILDA_$1"
TDIR="$SDIR/../Versions/$LTDIR"
echo "Creating new version in $TDIR"

mkdir $TDIR || exit 2

echo git tag -a v$1

cp README.rst LICENSE setpaths.bash $TDIR

cp -r scripts $TDIR
cp -r matilda $TDIR
cp -r EXAMPLES $TDIR

echo "Removing binary pyc files"
find $TDIR -name '*pyc' -exec rm {} \;

# create tar with shorter relative paths
cd $SDIR/../Versions
tar -czf $LTDIR.tar.gz $LTDIR
