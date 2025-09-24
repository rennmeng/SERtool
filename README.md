# SERtool
This tool is used to search for proteins with repetitive sequence amino acids and can also be applied to the search of nucleic acid sequences.


# Usage:

python SERtool.py [TARGET] --input FILE --start INT --point INT INT [INT INT]

High-performance sequence scanning tool powered by Rust. Scans protein sequences for amino acid patterns with dynamic hit thresholds.

Prerequisites:
  This tool is written in Rust and requires:
    1. Rust toolchain: cd SERtool
    2. maturin (Python): pip install maturin
  Build the extension first:
    maturin develop    # for development
    or
    maturin build --release && pip install target/wheels/*.whl

Positional Arguments:
  [TARGET]              Amino acid pattern(s) to scan for.
                        Formats:
                          X,Y    : independent tasks (e.g., D,E)
                          X-Y    : combined task (e.g., D-E → polyDE)
                          .X     : wildcard (e.g., .E → DE, AE, RE...), excludes XX

Required Arguments:
  --input FILE          Input FASTA file (e.g., clinvar_mutant.fasta)
  --start INT           Starting window size (e.g., 20)
  --point x1 y1 [x2 y2] Points defining hit threshold vs window size.
                        Example: 18 20 30 50 → linear model from (18,20) to (30,50)
                        Use two points for dynamic threshold, or one point (e.g. 50 30) for fixed.

Optional Arguments:
  --mode fix            Use fixed-mode filtering (e.g., exactly 30 hits in 50aa)
  -h, --help            Show this help message and exit

Examples:
  Independent scan for 20 amino acids:
  python SERtool.py A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y --input clinvar_mutant.fasta --start 20 --point 18 20 30 50

  Combined D-E analysis:
  python SERtool.py D-E --input clinvar_mutant.fasta --start 20 --point 18 20 30 50

  Find ≥10 consecutive matches:
  python SERtool.py D,E --input clinvar_mutant.fasta --start 10 --point 20 20 30 30

  Wildcard: all X-S motifs (e.g., DS, RS)
  python SERtool .S --input clinvar_mutant.fasta --start 15 --point 5 10 10 30

  Fixed mode: find regions with at least 30 hits in 50aa
  python SERtool.py A --input clinvar_mutant.fasta --mode fix --point 30 50
