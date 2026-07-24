#!/bin/bash

##########################################################################
# Pipeline:  ChloroFORGE
# Author:    Lucas Munnes
# 			 Institute for Biological Data Science, HHU
# GitHub:    https://github.com/usadellab/ChloroFORGE
##########################################################################

set -e 

#####################################
# Usage function
#####################################
usage() {
    echo "   ________    __                  __________  ____  ____________"
    echo "  / ____/ /_  / /___  _________  / ____/ __ \/ __ \/ ____/ ____/"
    echo " / /   / __ \/ / __ \/ ___/ __ \/ /_  / / / / /_/ / / __/ __/   "
    echo " / /___/ / / / / /_/ / /  / /_/ / __/ / /_/ / _, _/ /_/ / /___   "
    echo " \____/_/ /_/_/\____/_/   \____/_/    \____/_/ |_|\____/_____/   "
    echo
    echo "ChloroFORGE: A Pipeline to reassemble Chloroplast contigs from precomputed assemblies"
    echo
    echo "Usage: $0 -g GENOME -c CHLOROPLAST -t THREADS -o OUTPUT [-s ESTIMATED_SIZE -l CHROM_LIST -x CP_COV -f MIN_OVERLAP]"
    echo
    echo "Options:"
    echo "  -o OUTPUT        Sample name / project identifier"
    echo "  -g GENOME        Path to genome assembly"
    echo "  -c CHLOROPLAST   Path to chloroplast reference genome"
		echo "  -s ESTIMATED_SIZE Estimated size of the chloroplast (default: 150k)"
    echo "  -l CHROM_LIST    Chromosome list for anchored contigs (optional)"
    echo "  -t THREADS       Number of threads"
    echo "  -x CP_COV        Target cp contig coverage for Flye (default: 50)"
    echo "  -f MIN_OVERLAP   Minimum overlap for Flye assembly (default: 5000)"
    echo "      --allow-lowcov   Allow Flye even if cp coverage < target"
    echo
    exit 1
}

#####################################
# Parse arguments
#####################################
CP_COV=50
MIN_OVERLAP=5000
ALLOW_LOWCOV=false
ESTIMATED_SIZE=150k
THREADS=1
while getopts "o:g:c:s:l:t:x:f:-:" opt; do
    case $opt in
        o) SAMPLE="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        c) CHLOROPLAST="$OPTARG" ;;
        s) ESTIMATED_SIZE="$OPTARG" ;;
				l) CHROM_LIST="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        x) CP_COV="$OPTARG" ;;
        f) MIN_OVERLAP="$OPTARG" ;;
        -)
            case "${OPTARG}" in
                allow-lowcov) ALLOW_LOWCOV=true ;;
                *) usage ;;
            esac ;;
        *) usage ;;
    esac
done

[[ -z "$SAMPLE" || -z "$GENOME" || -z "$CHLOROPLAST" ]] && usage

#####################################
# Setup Working Directory
#####################################
echo "--- Starting ChloroForge for $SAMPLE ---"

USE_CHROM_LIST=true
[[ -z "$CHROM_LIST" ]] && USE_CHROM_LIST=false && echo "No chromosome list provided. Using full genome."

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORKDIR="./${SAMPLE}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

DEP_DIR="$SCRIPT_DIR/dependencies"
BIN_DIR="$DEP_DIR/bin"
FLYE_BIN="$DEP_DIR/Flye/bin"

export PATH="$PATH:$FLYE_BIN:$BIN_DIR:"

check_tool() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "[ERROR] Tool not found: $1"
        exit 1
    }
}

check_tool minimap2
check_tool flye
check_tool seqkit
check_tool blastn

echo "[INFO] All required tools found"

#####################################
# Step 1: Generate unanchored contigs
#####################################
UNANCHORED="unanchored_${SAMPLE}.fasta"

if [[ "$USE_CHROM_LIST" == false ]]; then
    ln -sf "$(realpath "$GENOME")" "$UNANCHORED"
    touch .unanchored.done
else
    if [[ ! -f ".unanchored.done" ]]; then
        echo "Step 1: Extracting unanchored contigs..."
        seqkit grep -v -f "$CHROM_LIST" "$GENOME" -o "$UNANCHORED"
        touch .unanchored.done
    fi
fi

#####################################
# Step 2: Make BLAST database
#####################################
if [[ ! -f ".blast_db.done" ]]; then
    echo "Step 2: Creating BLAST database..."
    makeblastdb -in "$CHLOROPLAST" -dbtype nucl -out "$SAMPLE"
    touch .blast_db.done
fi

#####################################
# Step 3: Run BLAST
#####################################
if [[ ! -f ".blast_hits.done" ]]; then
    echo "Step 3: Running BLAST..."
    blastn -query "$UNANCHORED" -db "$SAMPLE" -task blastn -evalue 0.001 \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
        -num_threads "$THREADS" > blast_hits.tsv
    touch .blast_hits.done
fi

#####################################
# Step 4: Filter high-identity contigs
#####################################
if [[ ! -f ".filter_cp.done" ]]; then
    echo "Step 4: Filtering chloroplast contigs..."
    awk -F'\t' '
    {
        q=$1; qlen=$5; s=$7; e=$8;
        if(s>e){t=s;s=e;e=t}
        cov[q]+=e-s; len[q]=qlen
    }
    END{
        for(i in cov){
            perc=(cov[i]/len[i])*100
            if(perc>=95) print i
        }
    }' blast_hits.tsv > cp_contigs.txt
    touch .filter_cp.done
fi

#####################################
# Step 5: Extract & Downsample
#####################################
CP_CONTIGS="cp_${SAMPLE}_contigs.fasta"
CP_DOWNSAMPLED="cp_${SAMPLE}_contigs_${CP_COV}x.fasta"
DONE_EXTRACT=".extract_cp_${CP_COV}x.done"

if [[ ! -f "$DONE_EXTRACT" ]]; then
    echo "Step 5: New Target Coverage ${CP_COV}x detected or first run. Processing..."
    rm -f "$CP_DOWNSAMPLED"

    seqkit grep -f cp_contigs.txt "unanchored_${SAMPLE}.fasta" -o "$CP_CONTIGS"

    CP_REF_LEN=$(seqkit stats -T "$CHLOROPLAST" | awk 'NR==2{print $5}')
    CP_TOTAL_LEN=$(seqkit stats -T "$CP_CONTIGS" | awk 'NR==2{print $5}')
    
    REAL_COV=$(awk "BEGIN{printf \"%.2f\", $CP_TOTAL_LEN/$CP_REF_LEN}")
    FRAC=$(awk "BEGIN{printf \"%.6f\", ($CP_REF_LEN*$CP_COV)/$CP_TOTAL_LEN}")

    echo "Real coverage: ${REAL_COV}x | Target: ${CP_COV}x"

    if awk "BEGIN{exit !($REAL_COV < $CP_COV)}"; then
        if [[ "$ALLOW_LOWCOV" != "true" ]]; then
            echo "[ERROR] Coverage ${REAL_COV}x < ${CP_COV}x. Use --allow-lowcov to proceed."
            exit 1
        fi
    fi

    if awk "BEGIN{exit !($FRAC < 1)}"; then
        seqkit sample -p "$FRAC" -s 100 "$CP_CONTIGS" > "$CP_DOWNSAMPLED"
    else
        cp "$CP_CONTIGS" "$CP_DOWNSAMPLED"
    fi
    touch "$DONE_EXTRACT"
fi


#####################################
# Step 6: Flye assembly
#####################################
if [[ ! -f ".flye.done" ]]; then
    echo "Step 6: Running Flye assembly..."
    flye \
        --subassemblies "$CP_DOWNSAMPLED" \
        --out-dir flye_cp_out \
        --genome-size "$ESTIMATED_SIZE" \
        --threads "$THREADS" \
        --plasmids \
        --min-overlap "$MIN_OVERLAP" \
        --asm-coverage "$CP_COV"
    touch .flye.done
fi

#####################################
# Step 7: Process assembly graph
#####################################
if [[ ! -f ".process_graph.done" ]]; then
    echo "Step 7: Processing assembly graph..."
    mkdir -p results_cp
    awk '/^S/{print ">"$2"\n"$3}' flye_cp_out/assembly_graph.gfa | fold > results_cp/edges.fa
    awk -F ":" '/^S/{print substr($1,3,7)"\t"$3}' flye_cp_out/assembly_graph.gfa | sort -k2 -n -r > results_cp/edges_depth
    touch .process_graph.done
fi

#####################################
# Step 8: Orientation/Flipping
#####################################
if [[ ! -f ".flip.done" ]]; then
    echo "Step 8: Running flipping and orientation..."
    SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
    
    PY_SCRIPT="${SCRIPT_DIR}/edges.py"
    if [[ -f "results_cp/edges.fa" ]]; then
        python3 "$PY_SCRIPT" "results_cp/edges.fa" "$CHLOROPLAST" "$SAMPLE"
        touch .flip.done
    fi
fi

#####################################
# Step 9: Validation
#####################################
CHLORO="results_cp/chloroplast_final_assembly/${SAMPLE}_chloroplast.fasta"

if [[ ! -f "$CHLORO" ]]; then
    echo "ERROR: Final assembly file not found at $CHLORO"
    exit 1
fi

echo "Step 9: Validating assembly length..."

REF_LEN=$(seqkit stats -T "$CHLOROPLAST" | awk 'NR==2 {print $5}')
ASM_LEN=$(seqkit stats -T "$CHLORO" | awk 'NR==2 {print $5}')

echo "   Reference: $REF_LEN bp | Assembly: $ASM_LEN bp"

TOL=0.10
MIN_LEN=$(awk -v len="$REF_LEN" -v tol="$TOL" 'BEGIN {printf "%.0f", len*(1-tol)}')
MAX_LEN=$(awk -v len="$REF_LEN" -v tol="$TOL" 'BEGIN {printf "%.0f", len*(1+tol)}')

if (( ASM_LEN < MIN_LEN || ASM_LEN > MAX_LEN )); then
    echo "ERROR: Chloroplast assembly length outside acceptable range."
    exit 1
else
    echo "Length in range — concatenating final genome."
fi

#####################################
# Step 10: Final Genome Construction
#####################################
echo "Step 10: Constructing final genome..."
FINAL="final_genome_${SAMPLE}.fasta"

seqkit grep -v -f "cp_contigs.txt" "$GENOME" -o "${SAMPLE}_genomic_contigs.fasta"
cat "${SAMPLE}_genomic_contigs.fasta" "$CHLORO" > "$FINAL"

echo "Pipeline completed successfully! Output: $FINAL"
