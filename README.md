# SERtool
This tool is used to search for proteins with repetitive sequence amino acids and can also be applied to the search of nucleic acid sequences.

📚 SERtool - High-Performance Sequence Scanning and Hit Calculation
A fast, Rust-accelerated tool for scanning protein sequences and calculating hit rates based on sliding windows and amino acid patterns. Built with PyO3 for seamless Python integration.

🔧 Command Syntax
python SERtool.py [TARGET] --input <fasta_file> --start <int> --point <int>... [OPTIONS]

🎯 1. Target Amino Acid Patterns ([TARGET])
Define the amino acid pattern(s) to scan for. Supports three syntax types:

Pattern      Example   Meaning

Combined     D-E       Treat D and E as a single unit (e.g., polyDE). Combined analysis.
Independent  D,E       Analyze D and E as separate tasks (e.g., polyD and polyE independently).
Wildcard     .E        Match all dipeptides ending in E (e.g., DE, AE, RE, PE), excluding EE. The '.' acts as a wildcard. Useful for motif screening (e.g., .S, .R).

✅ Examples:

A,C,D,E → Run independent analyses for polyA, polyC, polyD, polyE.
D-E → Combine D and E into a single polyDE-like region analysis.
.S → Scan for all X-S motifs (e.g., DS, RS, GS), useful for phosphorylation site context studies.
📁 Input File
--input <filename.fasta>
Specify the input FASTA file containing protein or gene sequences to analyze (e.g., clinvar_mutant.fasta).

🪟 Sliding Window Configuration
--start <int>
Set the starting window size for scanning.
The tool begins calculation from this window length and increases incrementally.

Example: --start 20 → start scanning with a window of 20 amino acids.

📈 Hit Calculation Model (--point)
--point x1 y1 x2 y2
Define two reference points (x1, y1) and (x2, y2) to establish a linear relationship between window size and hit count threshold.

The line formed by these points determines how many hits are required at any given window size.
Used to dynamically adjust the hit threshold as window size changes.
Common Use Cases:
Goal                                      Command Snippet

Fixed threshold within 50aa if 30 hits    --point 30 30 50 30
Require proportional hits                 --point 18 20 30 50
Detect runs with ≥10 consecutive matches  --start 10 --point 20 20 30 30

🔍 Example:
--point 18 20 30 50 → calculates slope between (18,20) and (30,50), then maps window size to required hit count using linear interpolation.

⚙️ Optional Arguments
--mode fix
Use fixed-mode analysis. Only evaluate results at exact point(s) defined by --point, typically used for strict filtering (e.g., find regions where at least 30 out of 50 residues match).

Tip: Combine with --point 1 100 to find sequences containing a specific segment of interest.

✅ Example Commands
Scan for multiple homopolymeric tracts independently:
python SERtool.py A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y --input clinvar_mutant.fasta --start 20 --point 18 20 30 50

Analyze polyDE combined regions:
python SERtool.py D-E --input clinvar_mutant.fasta --start 20 --point 18 20 30 50

Find regions with ≥10 consecutive matching residues:
python SERtool.py D,E --input clinvar_mutant.fasta --start 10 --point 20 20 30 30

Screen for serine-rich motifs (e.g., DS, RS, GS):
python SERtool.py .S --input clinvar_mutant.fasta --start 15 --point 5 15 10 30

Extract sequences containing a specific subsequence (fixed mode):
python SERtool.py A --input clinvar_mutant.fasta --mode fix --point 1 100

💡 Tips
Use smaller --start values for detecting short motifs.
Adjust --point values to tune sensitivity: steeper slope = stricter matching.
Combine wildcards (.X) with biological context (e.g., phosphorylation, disorder) for targeted screening.
🚀 Performance Note
This tool is powered by Rust via PyO3 and uses parallel processing (rayon) for high-performance scanning of large sequence datasets.
Ensure you've built the extension with maturin develop or installed the wheel before running.
