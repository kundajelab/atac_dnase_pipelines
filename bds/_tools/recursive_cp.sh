#!/bin/bash

if [ "$#" -lt 3 ]
then
  echo "Usage: recursive_cp.sh [SOURCE_DIR] [DEST_DIR] [FILE_EXT]"
  exit 1
fi

SRC=$1
DEST=$2
EXT=$3

if [ ! -d $SRC ]
then
  echo "Source directory ($SRC) doesn't exist!"
  exit 2
fi

if [ ! -d $DEST ]
then
  echo "Destination directory ($DEST) doesn't exist!"
  exit 3
fi

cd $SRC

for f in $(find . -name "*.$EXT" )
do
  BASENAME=$(basename $f)
  DIRNAME=$(dirname $f)
  mkdir -p "$DEST/$DIRNAME"

  FULLPATH=$(readlink -f $f)
  TARGET="$DEST/$DIRNAME/$BASENAME"
  cp $FULLPATH $TARGET
done

