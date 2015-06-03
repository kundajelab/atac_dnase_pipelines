#!/bin/bash

if [ "$#" -lt 2 ]
then
  echo "Usage: sync_bds_report.sh [SOURCE_DIR] [DEST_DIR]"
  exit 1
fi

SRC=$1
DEST=$2

if [ ! -d $SRC ]
then
  echo "Source directory (bds working dir.) ($SRC) doesn't exist!"
  exit 2
fi

if [ ! -d $DEST ]
then
  echo "Destination directory (web dir.) ($DEST) doesn't exist!"
  exit 3
fi

cd $SRC

#for f in $(find . -name '*.html' -not -name '*_parallel_*' -or -name '*.js' -not -name '*_parallel_*' -or -name '*.pdf' -or -name '*.png' -or -name '*.log' -or -name '*.qc' )
for f in $(find . -name '*.html' -or -name '*.js' -or -name '*.pdf' -or -name '*.png' -or -name '*.log' -or -name '*.qc' -or -name '*.preseq.dat' )
do
  BASENAME=$(basename $f)
  DIRNAME=$(dirname $f)
  if [ ! -d "$DEST/$DIRNAME" ]; then
    mkdir -p "$DEST/$DIRNAME"
  fi

  FULLPATH=$(readlink -f $f)
  TARGET="$DEST/$DIRNAME/$BASENAME"

  if [ -f $TARGET ]; then
    SIZE1=$(stat -c%s $FULLPATH)
    SIZE2=$(stat -c%s $TARGET)
    if [[ "$SIZE1" -ne "$SIZE2" ]]; then
      cp $FULLPATH $TARGET
    fi
  else
    cp $FULLPATH $TARGET
  fi

done

