#!/bin/bash

# Usage: ./testdiff_all.sh <source_directory>
# Example: ./testdiff_all.sh ../modernFortran_Poseidon_stable181123/

#trap errors
set -eE -o functrace
failure() {
  local lineno=$1
  local msg=$2
  echo "Failed at $lineno: $msg"
}
trap 'failure ${LINENO} "$BASH_COMMAND"' ERR

# Check if source directory is provided and exists
if [ -z "$1" ]; then
  echo "Source directory not provided. Usage: $0 <source_directory>"
  exit 1
elif [ ! -d "$1" ]; then
  echo "Source directory $1 does not exist."
  exit 1
fi

srcdir=$1

# Clear or create diff_all.out at the start
> diff_all.out

#main script
cd ..
find . -type f -name '*' | while read -r file; do
    #echo "$file"
    #echo "$file" >>diff_all.out
    # Correctly form the source file path
    srcfile="${srcdir}/${file#./}"
    # Ensure srcfile path doesn't start with multiple slashes in case srcdir was given as a relative path
    srcfile=$(echo $srcfile | sed 's://*:/:g')
    # Check if the file exists in the source directory and is readable
    if [ -f "$srcfile" ]; then
        diff -q "$file" "$srcfile" >> scripts/diff_all.out || echo "Diff command failed for $file" >> scripts/diff_all.out
    else
        echo "Source file $srcfile not found." >> scripts/diff_all.out
    fi
done

