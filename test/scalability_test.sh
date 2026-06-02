#!/bin/bash

# Nome eseguibile
EXEC=./main

# File output risultati
OUTFILE=results.txt

# Pulisce file precedente
echo "MPI_PROCS TIME(s)" > $OUTFILE

echo "Test running on a 100 element grid, with tolerance 1e-06 and max iterations set to 1e06"

# Lista dei processi MPI da testare
PROCS_LIST="1 2 4 8 16"

echo "Starting scalability test"

for N in $PROCS_LIST
do
    echo "Test using $N MPI processes"

    mpirun -np $N $EXEC > /dev/null

    echo "$N $avg" >> $OUTFILE

done

echo "Test completato. Risultati in $OUTFILE"