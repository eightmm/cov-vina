#!/bin/bash
# Benchmark suite for 6lu7 with various warheads

set -e

PROTEIN="examples/6lu7/6lu7_pocket.pdb"
RESIDUE="CYS145"
OUTPUT_BASE="examples/6lu7/results"

# Create output directory
mkdir -p "$OUTPUT_BASE"

echo "=== CovVina Benchmark Suite ==="
echo "Protein: $PROTEIN"
echo "Target residue: $RESIDUE"
echo ""

# Test ligands with different warheads
declare -A LIGANDS=(
    ["vinyl_alanine"]="C=CC(=O)N[C@@H](C)C(=O)O"
    ["acrylamide_phenyl"]="C=CC(=O)Nc1ccccc1"
    ["chloroacetamide"]="ClCC(=O)Nc1ccccc1"
    ["vinyl_sulfonamide"]="C=CS(=O)(=O)Nc1ccccc1"
    ["bromoacetamide"]="BrCC(=O)Nc1ccccc1"
    ["acrylamide_methoxy"]="C=CC(=O)Nc1ccc(OC)cc1"
)

# Run each ligand
for name in "${!LIGANDS[@]}"; do
    smiles="${LIGANDS[$name]}"
    output_dir="$OUTPUT_BASE/$name"

    echo "[$name] SMILES: $smiles"
    echo "  Running docking..."

    ~/.local/bin/uv run python scripts/run_covalent_pipeline.py \
        -p "$PROTEIN" \
        -q "$smiles" \
        -r "$RESIDUE" \
        -o "$output_dir" \
        --optimize \
        -n 500 \
        --opt_steps 200 \
        --save_all

    echo "  ✓ Complete: $output_dir"
    echo ""
done

echo "=== Benchmark Complete ==="
echo "Results saved to: $OUTPUT_BASE"
