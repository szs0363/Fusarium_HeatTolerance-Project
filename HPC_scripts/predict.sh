#!/bin/bash


module load funannotate/1.8.13

# Define directories
MASK_DIR="/scratch/aubsxs002/Fusarium_HTproject/genomes_masked"
PREDICT_DIR="/scratch/aubsxs002/Fusarium_HTproject/predictions_all"
PROTEIN_DIR="/scratch/aubsxs002/Fusarium_HTproject/proteins_all"

# Create directories
mkdir -p "$PREDICT_DIR" "$PROTEIN_DIR"

echo "=== FUNANNOTATE GENE PREDICTION ==="

for genome in "$MASK_DIR"/*_masked.fasta; do
    basename=$(basename "$genome" _masked.fasta)
    echo "Predicting genes for: $basename"
    
    funannotate predict -i "$genome" \
                       -o "$PREDICT_DIR/${basename}" \
                       -s "$basename" \
                       --species "Fusarium oxysporum" \
                       --busco_db /home/aubsxs002/funannotate_db/pezizomycotina \
                       --cpus 32 \
                       --optimize_augustus
    
    # Extract proteins
    if [ -f "$PREDICT_DIR/${basename}/predict_results/${basename}.proteins.fa" ]; then
        cp "$PREDICT_DIR/${basename}/predict_results/${basename}.proteins.fa" "$PROTEIN_DIR/${basename}.faa"
        echo "  ✓ Proteins extracted: ${basename}.faa"
    fi
    
    echo "  ✓ Completed: $basename"
    echo ""
done

echo "Gene prediction complete!"

