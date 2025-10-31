# SERtool

This tool is used to search for proteins with repetitive sequence amino acids and can also be applied to the search of nucleic acid sequences.

## Usage
```bash
python SERtool.py [TARGET] --input FILE --start INT --point INT INT [INT INT]
```

High-performance sequence scanning tool powered by Rust. Scans protein sequences for amino acid patterns with dynamic hit thresholds.

### Prerequisites
This tool is written in Rust and requires:
  1. Install requirements: 
     ```bash
     git clone https://github.com/rennmeng/SERtool
     cd SERtool
     pip install -r requirements.txt
     ```
  2. Install the Rust core: 
     ```bash
     sudo apt update
     sudo apt install build-essential
     curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
     ```
  3. Build maturin (Install SERtool): 
     ```bash
     . "$HOME/.cargo/env"
     maturin build --release && pip install target/wheels/*.whl
     ```

### Positional Arguments
PATTERN formats:<br>
&nbsp;&nbsp;<code>X,Y</code>&nbsp;&nbsp;&nbsp;: Independent scan (e.g., <code>D,E</code> → poly-D and poly-E)<br>
&nbsp;&nbsp;<code>XY</code>&nbsp;&nbsp;&nbsp;&nbsp;: Tandem motif (e.g., <code>DE</code> → DEDEDE)<br>
&nbsp;&nbsp;<code>X-Y</code>&nbsp;&nbsp;&nbsp;: Mixed class (e.g., <code>D-E</code> → DDDEEDE)<br>
&nbsp;&nbsp;<code>.X</code>&nbsp;&nbsp;&nbsp;&nbsp;: Wildcard (e.g., <code>.E</code> → any-E)

### Required Arguments
`--input FILE`  
&nbsp;&nbsp;Input FASTA file (e.g., test.fasta)  
`--start INT`  
&nbsp;&nbsp;Start window size (e.g., 20) 
`--end INT`  
&nbsp;&nbsp;End window size (e.g., 100)  
`--mode STR`  
&nbsp;&nbsp;mode (Options: formula (default) or fix)  
`--point x1 y1 [x2 y2]`  
&nbsp;&nbsp;Points defining hit threshold vs window size.  
&nbsp;&nbsp;Example: 18 20 30 50 → linear model from (18,20) to (30,50)  
&nbsp;&nbsp;Use two points for dynamic threshold, or one point (e.g. 50 30) for fixed.

### Examples:

Independent scan for 20 amino acids:
```bash
python SERtool.py A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y --input test.fasta --start 20 --point 18 20 30 50
```

Combined D-E analysis:
```bash
python SERtool.py D-E --input test.fasta --start 20 --point 18 20 30 50
```

Wildcard: all X-S motifs (e.g., DS, RS)
```bash
python SERtool.py .S --input test.fasta --start 15 --point 5 10 10 30
```

Fixed mode: find regions with at least 30 hits in 50aa
```bash
python SERtool.py E --input test.fasta --mode fix --point 30 50
```

Consecutive Mode: Screen consecutive target sequences
Format: --point hit1 hit1 hit2 hit2
Example: Find ≥10 consecutive matches: ***start=10***
```bash
python SERtool.py E --input test.fasta --start 10 --point 20 20 30 30
```

Score Mode: Use weighted scoring (e.g., each hit counts as 1, but some positions contribute extra score).
Format: --point hit1 [2*hit1 + score] hit2 [2*hit2 + score]
Example: Interpreted as: window = 2×hit + 10. Enables non-linear sensitivity: ***score=10***
```bash
python SERtool.py E --input test.fasta --start 15 --point 10 30 20 50
```

Rate Mode: Use proportional threshold (e.g., k×hit) for density-based filtering.
Format: --point hit1 [k*hit1] hit2 [k*hit2]
Example: Interpreted as: window= 2×hit (k=2). Ensures high-density regions (e.g., ≥50% occupancy): ***k=2***
```bash
python SERtool.py E --input test.fasta --start 15 --point 10 20 20 40
```

## Dataset
We generated a dataset of HGVS mutation sequences and developed a web-based tool for analysis(manuscript in preparation), both of which are available at https://www.sertool.net.</br>
The dataset was processed, filtered, and standardized to construct a comprehensive collection of mutation sequences for downstream analysis of repetitive motifs and functional impact assessment, VCF files were obtained by August 2, 2025.

## About This Tool
https://www.sertool.net</br>https://bio.tools/sertool</br>
If you encounter any problems with the Rust installation, you can use pySERtool.py instead, which is a pure Python implementation. Although it is relatively slower, approximately 30–40 times slower than the Rust version based on our previous benchmarks, it still meets the requirements for routine sequence searches.</br>
Notably, pySERtool supports configurable sequence step (step) and window sliding step (wstep), which were omitted from the Rust implementation due to limited performance gains and increased computational overhead. Dynamic stepping in the Rust version showed negligible improvement over simple iteration.</br>
Example usage: python pySERtool.py D-E --input test.fasta --start 20 --point 18 20 30 50</br>
Contact Us: rennmeng@smail.nju.edu.cn
