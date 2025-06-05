# MPI-Based Distance Matrix Computation for Multiple Sequence Alignment (MSA)

This project implements a parallel algorithm using MPI to compute a distance matrix for Multiple Sequence Alignment (MSA). It draws inspiration from progressive alignment strategies used in tools like MAFFT, with a focus on correctness, parallel efficiency, and scalability.

## Objective

- Perform pairwise alignments between biological sequences (DNA/protein)
- Calculate alignment scores and convert them into a distance matrix
- Use MPI to parallelize the computation across multiple processes
- Optionally generate a guide tree to assist with progressive alignment

## File Structure


## Algorithm Overview

### Pairwise Sequence Alignment

- Each MPI rank receives a subset of sequence pairs.
- Global alignment is performed using the Needleman–Wunsch algorithm.
- Alignment score is calculated for each pair and converted into a distance.

### Distance Matrix Construction

- Pairwise distances are gathered by the root process (rank 0).
- A complete symmetric matrix is printed or optionally saved to a file.



