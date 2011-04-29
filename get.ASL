#!/bin/sh

echo " "
echo "Script for downloading ASL source code"
echo " "

echo "Check for wget"
wgetAvailable=`which wget 2>/dev/null | wc -w`
if test ! $wgetAvailable = 1; then
  echo "wget not found"
  exit -1;
fi


echo "Download ASL code from www.netlib.org"
rm -f solvers.tar
wget ftp://www.netlib.org/ampl/solvers.tar

echo "Unpack the files"
rm -rf solvers
tar xf solvers.tar

echo "Uncompress the files"
gunzip -fr solvers

echo "Delete the tar file"
rm solvers.tar

echo " "
echo "Script complete, ASL source retrieved successfully!"
echo " "

