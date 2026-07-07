#!/usr/bin/env bash
# Fetch Eigen 3.3.9 headers into extern/eigen (gitignored).
# NMSS requires Eigen 3.3 -- 3.4 breaks the vec(double) indexing used here.
set -euo pipefail

VER=3.3.9
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$HERE/extern/eigen"

if [ -f "$DEST/Eigen/Dense" ]; then
    echo "Eigen already present at $DEST"
    exit 0
fi

mkdir -p "$HERE/extern"
TMP="$(mktemp -d)"
URL="https://gitlab.com/libeigen/eigen/-/archive/${VER}/eigen-${VER}.tar.gz"
echo "downloading $URL"
curl -sL --max-time 180 "$URL" -o "$TMP/eigen.tar.gz"
tar -xzf "$TMP/eigen.tar.gz" -C "$TMP"
rm -rf "$DEST"
mv "$TMP/eigen-${VER}" "$DEST"
rm -rf "$TMP"
echo "Eigen ${VER} installed at $DEST"
