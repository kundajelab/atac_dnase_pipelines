#!/bin/bash

SRC=$1
DEST=$2

if [ ! -d $SRC ]
then
  echo "Source directory ($SRC) doesn't exist!"
  exit 2
fi

if [ $# < 2 ]
then
  echo "Destination server not specified!"
  exit 3
fi

scp -pr $SRC "$DEST:$(dirname $SRC)/"
