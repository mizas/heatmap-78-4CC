#!/bin/bash

# Directorio donde se encuentran los archivos
directory="03-proteinortho-graph"

# Iterar sobre cada archivo que comienza con "file"
for file in "$directory"/file*; do
    # Nombre del archivo de salida
    output_file="${file%.txt}_cut.txt"
    
    # Cortar la primera y segunda columna y guardarlas en el archivo de salida
    cut -f1,2 "$file" > "$output_file"
done
