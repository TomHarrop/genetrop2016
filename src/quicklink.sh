#!/bin/bash

export TARGET="/home/tom/Projects/fiveAccessions/output"
export DESTINATION="data/fiveacc"

# make directories
find $TARGET -type d -print0 | xargs -0 bash -c 'for DIR in "$@"; 
do
  mkdir -p $DESTINATION${DIR#$TARGET}
  done' -

# make links
find $TARGET -type f -name "*.Rds" -print0 | xargs -0 bash -c 'for file in "$@"; 
do
  ln -s "$file" "$DESTINATION"${file#$TARGET}
  done' -