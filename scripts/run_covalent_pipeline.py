import argparse

from cov_vina import run_covalent_pipeline


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CovVina: Covalent Docking Pipeline with Warhead Anchor & Vina Scoring"
    )

    # Required
    parser.add_argument("-p", "--protein", required=True,
                        help="Path to the protein pocket PDB file")
    parser.add_argument("-q", "--query_ligand", required=True,
                        help="SMILES string or path to an SDF file of the query ligand")

    # Covalent-specific
    parser.add_argument("-r", "--reactive_residue", default=None,
                        help="Reactive residue specifier, e.g. 'CYS145' or 'CYS145:A'. "
                             "If omitted, auto-detect the first supported reactive residue.")
    parser.add_argument("--pocket_cutoff", type=float, default=12.0,
                        help="Distance cutoff (Å) for pocket extraction around residue. "
                             "Default 12Å covers Vina's effective range. Use 10Å for small molecules, "
                             "15-20Å for large ligands/peptides.")

    # Output
    parser.add_argument("-o", "--out_dir", default="output_predictions",
                        help="Directory to save the resulting SDF file (default: output_predictions)")

    # Conformer generation
    parser.add_argument("-n", "--num_confs", type=int, default=1000,
                        help="Number of conformers to generate (default: 1000)")
    parser.add_argument("--rmsd_threshold", type=float, default=1.0,
                        help="RMSD threshold (Å) for clustering (default: 1.0)")

    # Relaxation
    parser.add_argument("--no_mmff", action="store_true",
                        help="Disable MMFF94 force field relaxation after anchor placement")

    # Optimization
    parser.add_argument("--optimize", action="store_true",
                        help="Enable gradient-based torsion optimization")
    parser.add_argument("--opt_steps", type=int, default=100,
                        help="Number of optimization steps (default: 100)")
    parser.add_argument("--opt_lr", type=float, default=0.05,
                        help="Learning rate for optimization (default: 0.05)")
    parser.add_argument("--opt_batch_size", type=int, default=128,
                        help="Batch size for optimization (default: 128)")
    parser.add_argument("--optimizer", type=str, choices=["adam", "adamw", "lbfgs"],
                        default="adam",
                        help="Optimizer for torsion optimization (default: adam)")

    # Scoring
    parser.add_argument("--weight_preset", type=str,
                        choices=["vina", "vina_lp", "vinardo"], default="vina",
                        help="Vina scoring weight preset (default: vina)")
    parser.set_defaults(torsion_penalty=True)
    torsion_group = parser.add_mutually_exclusive_group()
    torsion_group.add_argument("--torsion_penalty", dest="torsion_penalty",
                               action="store_true",
                               help="Include torsional entropy penalty (default)")
    torsion_group.add_argument("--no_torsion_penalty", dest="torsion_penalty",
                               action="store_false",
                               help="Disable torsional entropy penalty")

    # Output control
    parser.add_argument("--save_all", action="store_true",
                        help="Save all poses instead of top-k")
    parser.add_argument("--top_k", type=int, default=None,
                        help="Number of top poses to save (default: 3 without --optimize, all with --optimize)")

    args = parser.parse_args()

    results = run_covalent_pipeline(
        protein_pdb=args.protein,
        query_ligand=args.query_ligand,
        reactive_residue=args.reactive_residue,
        output_dir=args.out_dir,
        pocket_cutoff=args.pocket_cutoff,
        num_confs=args.num_confs,
        rmsd_threshold=args.rmsd_threshold,
        mmff_optimize=not args.no_mmff,
        optimize=args.optimize,
        optimizer=args.optimizer,
        opt_steps=args.opt_steps,
        opt_lr=args.opt_lr,
        opt_batch_size=args.opt_batch_size,
        weight_preset=args.weight_preset,
        torsion_penalty=args.torsion_penalty,
        save_all_poses=args.save_all if args.save_all else None,
        top_k=args.top_k,
    )

    print(f"\nBest score: {results['best_score']:.3f} kcal/mol")
    print(f"Warhead: {results['warhead_type']}")
    print(f"Anchor: {results['anchor_residue']} ({results['anchor_atom']})")
