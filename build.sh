#!/bin/sh
# Build apo, the Source blob, the verification library and byu_export (the BYU-format tables) into bin/.
set -e
cd "$(dirname "$0")"
mkdir -p bin
FLAGS="-O2 -ffp-contract=off -fno-builtin-sin -fno-builtin-cos -Icore"
cc -O2 -Icore -o bin/source_blob_build source/blob_build.c -lm
./bin/source_blob_build source bin/source.blob
cc $FLAGS -o bin/byu_export source/byu_export.c core/*.c -lm -lpthread
cc $FLAGS -fPIC -shared -o bin/_apostruct.so core/*.c -lm
cc $FLAGS -o bin/apo core/*.c cli/apo_main.c gui/server.c -lm -lpthread
cat $(ls core/*.c core/*.h source/blob_build.c | LC_ALL=C sort) | { command -v shasum >/dev/null && shasum -a 256 || sha256sum; } | cut -d' ' -f1 > bin/_apostruct.source
