#!/usr/bin/env bash

##########################################################################
# Pipeline:  ChloroFORGE
# Author:    Lucas Munnes
# 	     Institute for Biological Data Science, HHU
# GitHub:    https://github.com/usadellab/ChloroFORGE
##########################################################################

set -euo pipefail

# =========================================================
# BASE DIRECTORY (REPO ROOT)
# =========================================================
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="$BASE_DIR/dependencies"
BIN_DIR="$INSTALL_DIR/bin"

mkdir -p "$BIN_DIR"

echo "========================================"
echo "Installing tools into: $INSTALL_DIR"
echo "========================================"

export PATH="$BIN_DIR:$PATH"

# =========================================================
# HELPERS
# =========================================================
exists_local() {
    [[ -x "$BIN_DIR/$1" ]]
}

exists_global() {
    command -v "$1" >/dev/null 2>&1
}

link_global() {
    TOOL="$1"
    TOOL_PATH="$(command -v "$TOOL")"

    echo "→ Using system $TOOL at $TOOL_PATH"
    ln -sf "$TOOL_PATH" "$BIN_DIR/$TOOL"
}
chmod 555 chloroFORGE.sh
# =========================================================
# MINIMAP2 2.30
# =========================================================
if exists_local minimap2; then
    echo "[OK] minimap2 already in dependencies"
elif exists_global minimap2; then
    link_global minimap2
else
    echo "[INSTALL] minimap2"
    cd "$INSTALL_DIR"
    git clone https://github.com/lh3/minimap2.git
    cd minimap2
	git checkout v2.30
    make
    cp minimap2 "$BIN_DIR/"
fi

# =========================================================
# SEQKIT 2.13.0
# =========================================================
if exists_local seqkit; then
    echo "[OK] seqkit already in dependencies"
elif exists_global seqkit; then
    link_global seqkit
else
    echo "[INSTALL] seqkit"
    cd "$INSTALL_DIR"
    wget -q https://github.com/shenwei356/seqkit/releases/download/v2.13.0/seqkit_linux_amd64.tar.gz
    tar -xzf seqkit_linux_amd64.tar.gz
    mv seqkit "$BIN_DIR/"
    rm -f seqkit_linux_amd64.tar.gz
fi

# =========================================================
# BLAST+ 2.17.0
# =========================================================
if exists_local blastn; then
    echo "[OK] BLAST already in dependencies"
elif exists_global blastn; then
    link_global blastn
    link_global makeblastdb
else
    echo "[INSTALL] BLAST+"
    cd "$INSTALL_DIR"
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*-x64-linux.tar.gz
    tar -xzf ncbi-blast-*-x64-linux.tar.gz

    BLAST_DIR=$(find . -maxdepth 1 -type d -name "ncbi-blast-*")
    cp "$BLAST_DIR/bin/"* "$BIN_DIR/"
    rm -rf "$BLAST_DIR"
fi

# =========================================================
# FLYE 2.9.6
# =========================================================
INSTALL_FLYE=false

if exists_local flye; then
    echo "[OK] flye already in dependencies"
elif exists_global flye; then
    VERSION=$(flye --version 2>/dev/null || true)

    if [[ "$VERSION" == *"2.9.6"* ]]; then
        link_global flye
    else
        echo "[WARN] Flye found but wrong version: $VERSION"
        INSTALL_FLYE=true
    fi
else
    INSTALL_FLYE=true
fi

if [[ "$INSTALL_FLYE" == true ]]; then
    echo "[INSTALL] Flye 2.9.6"

    cd "$INSTALL_DIR"


    git clone https://github.com/fenderglass/Flye.git
    cd Flye

    git checkout 2.9.6

    make
    mv bin/flye "$BIN_DIR/"
fi

# =========================================================
# DONE
# =========================================================
echo "========================================"
echo "Setup complete"
echo ""
echo "Path to installation:"
echo "export PATH=$BIN_DIR:\$PATH"
echo "========================================"
