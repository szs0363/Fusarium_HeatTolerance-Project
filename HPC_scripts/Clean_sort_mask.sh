#!/bin/bash
#Load funannotate
module load funannotate/1.8.13

# Define directories
INPUT_DIR="/scratch/aubsxs002/Fusarium_HTproject/FSSCgenomes"
CLEAN_DIR="/scratch/aubsxs002/Fusarium_HTproject/genomes_clean"
SORT_DIR="/scratch/aubsxs002/Fusarium_HTproject/genomes_sorted"
MASK_DIR="/scratch/aubsxs002/Fusarium_HTproject/genomes_masked"

# Create directories
mkdir -p "$CLEAN_DIR" "$SORT_DIR" "$MASK_DIR"

echo "=== FUNANNOTATE GENOME PREPROCESSING ==="

for genome in "$INPUT_DIR"/*.fasta; do
    basename=$(basename "$genome" .fasta)
    echo "Processing: $basename"
    
    # STEP 1: Clean genome (remove short contigs)
    echo "  Step 1: Cleaning genome..."
    funannotate clean -i "$genome" \
                     -o "$CLEAN_DIR/${basename}_clean.fasta" \
                     --minlen 1000
    
    # STEP 2: Sort contigs by size
    echo "  Step 2: Sorting contigs..."
    funannotate sort -i "$CLEAN_DIR/${basename}_clean.fasta" \
                    -o "$SORT_DIR/${basename}_sorted.fasta"
    
    # STEP 3: Mask repetitive elements
    echo "  Step 3: Masking repeats..."
    funannotate mask -i "$SORT_DIR/${basename}_sorted.fasta" \
                    -o "$MASK_DIR/${basename}_masked.fasta" \
                    --cpus 16
    
    echo "  ✓ Completed: $basename"
    echo ""
done

echo "=== PREPROCESSING COMPLETE ==="
echo "Clean genomes: $CLEAN_DIR"
echo "Sorted genomes: $SORT_DIR"
echo "Masked genomes: $MASK_DIR"
echo ""
echo "Next step: Run gene prediction on masked genomes"