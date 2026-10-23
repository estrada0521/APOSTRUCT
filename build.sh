#!/bin/sh
# Build apo, the Source blob, the verification library and source_dump (what apo
# reads from the Source) into bin/.
# CC picks the compiler (a mingw-w64 one builds Windows executables, natively under MSYS2 or as a cross build);
# HOSTCC builds the blob compiler, which runs on this machine.
set -e
cd "$(dirname "$0")"
mkdir -p bin
CC=${CC:-cc}
HOSTCC=${HOSTCC:-cc}
SRC=$(find core -name '*.c' | LC_ALL=C sort)
INC=$(find core -type d | LC_ALL=C sort | sed 's/^/-I/')
FLAGS="-O2 -ffp-contract=off -fno-builtin-sin -fno-builtin-cos $INC"
EXE=
LIBS="-lm -lpthread"
case "$($CC -dumpmachine)" in
*mingw* | *windows*) EXE=.exe LIBS="-lm -lpthread -lws2_32 -static" ;;
esac
$HOSTCC -O2 $INC -o bin/source_blob_build source/blob_build.c -lm
./bin/source_blob_build source bin/source.blob
$CC $FLAGS -o bin/source_dump$EXE source/source_dump.c $SRC $LIBS
[ -z "$EXE" ] && $CC $FLAGS -fPIC -shared -o bin/_apostruct.so $SRC -lm
$CC $FLAGS -o bin/apo$EXE $SRC cli/apo_main.c gui/server.c $LIBS
cat $( (find core -name '*.[ch]'; echo source/blob_build.c) | LC_ALL=C sort) | { command -v shasum >/dev/null && shasum -a 256 || sha256sum; } | cut -d' ' -f1 > bin/_apostruct.source
