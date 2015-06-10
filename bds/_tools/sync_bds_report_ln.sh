#!/bin/bash

if [ "$#" -lt 2 ]
then
  echo "Usage: sync_bds_report_ln.sh [SOURCE_DIR] [DEST_DIR]"
  exit 1
fi

SRC=$1
DEST=$2

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

cd $DEST
find $DEST -type l -xtype l -delete
find $DEST -mindepth 1 -type d -empty -delete

cd $SRC

#for f in $(find . -type f ! -name "*.js" ! -path "*.bds.*" )
for f in $(find . -type f ! -name "*.sh" ! -name "*.cluster" ! -name "*.stderr" ! -name "*.stdout" ! -name "*bds.pid*" )
do
  BASENAME=$(basename $f)
  DIRNAME=$(dirname $f)
  mkdir -p "$DEST/$DIRNAME"

  FULLPATH=$(readlink -f $f)
  TARGET="$DEST/$DIRNAME/$BASENAME"
  if [ ! -f $TARGET ];
  then
    ln -s $FULLPATH $TARGET
  fi
done

