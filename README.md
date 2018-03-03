# Sequence Alignment Algorithms
Scripts to align sequences together with penalties. Includes global and local alignment algorithms

TODO: Fix and test affine gap alignment

## Instructions
infile.txt is a text file with two strings of amino acids; see test files for formatting.
Global alignment:
```bash
./global_aligment.py <infile.txt> <('blosum62', 'pam250', 'basic')> <sigma value>
```
Local alignment:
```bash
./local_alignment.py <infile.txt> <('blosum62', 'pam250', 'basic')> <sigma value>
```
