#!/bin/bash

# Nome eseguibile
EXEC=./test/scalability_test

# File output risultati
OUTFILE=test/results.txt

# Pulisce file precedente e scrive l'intestazione allineata a sinistra
printf "%-12s %-15s %-15s %-12s\n" "MPI_PROCS" "TIME(s)" "ERROR" "ITER" > $OUTFILE

echo "Test running on a 256 element grid, with tolerance 1e-06 and max iterations set to 1e06"

# Lista dei processi MPI da testare
PROCS_LIST="1 2 4 8 16"

RUNS=3                 # ← quante volte ripetere ogni test

SRC=test/scalability_test.cpp         # ← il tuo file sorgente

echo "Compiling $SRC..."
mpicxx -O2 -fopenmp -Iinclude -o $EXEC $SRC src/solver.cpp

if [ $? -ne 0 ]; then
    echo "Errore di compilazione. Esco."
    exit 1
fi

echo "Compilazione completata."
echo "Starting scalability test"

for N in $PROCS_LIST
do
    echo "----------------------------------------"
    echo "Running test with $N MPI processes..."

    # 1. Esegue il programma MPI e cattura l'output
    OUTPUT=$(mpirun -np $N $EXEC)

    # STAMPA DI DEBUG: Vediamo a schermo cosa sputa fuori il C++
    echo "Program output:"
    echo "$OUTPUT"

    # 2. Isola la riga contenente i dati
    DATA_LINE=$(echo "$OUTPUT" | grep "TIME=")

    # 3. Estrae i singoli valori numerici rimuovendo le etichette di testo
    TIME_VAL=$(echo "$DATA_LINE" | awk '{print $1}' | sed 's/TIME=//')
    ERR_VAL=$(echo "$DATA_LINE" | awk '{print $2}' | sed 's/ERROR=//')
    ITER_VAL=$(echo "$DATA_LINE" | awk '{print $3}' | sed 's/ITER=//')

    # Se per qualche motivo le variabili sono vuote, mettiamo un flag di errore per non far rompere printf
    if [ -z "$TIME_VAL" ]; then TIME_VAL="N/A"; fi
    if [ -z "$ERR_VAL" ]; then ERR_VAL="N/A"; fi
    if [ -z "$ITER_VAL" ]; then ITER_VAL="N/A"; fi

    # 4. Scrive i dati formattati nel file dei risultati usando spaziature fisse
    printf "%-12s %-15s %-15s %-12s\n" "$N" "$TIME_VAL" "$ERR_VAL" "$ITER_VAL" >> $OUTFILE
done

echo "----------------------------------------"
echo "Test completato. Risultati in $OUTFILE"