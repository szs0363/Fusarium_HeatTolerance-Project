module load busco

GENOME_DIR="/scratch/aubsxs002/Fusarium_HTproject/FSSCgenomes"
OUTPUT_DIR="/scratch/aubsxs002/Fusarium_HTproject/results/busco_results"

mkdir -p "$OUTPUT_DIR"

for genome in "$GENOME_DIR"/*.fasta; do
    basename=$(basename "$genome" .fasta)
    echo "Running BUSCO for: $basename"
    
    busco -i "$genome" \
          -o "${basename}_busco" \
          -l fungi_odb10 \
          -m genome \
          --cpu 16 \
          --out_path "$OUTPUT_DIR"
done

echo "BUSCO analysis complete!"