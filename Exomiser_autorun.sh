#!/bin/bash

set -eou pipefail

# Directorios
exomiser_input="$HOME/Documentos/Fertility/Exomas/exomiser_input"
exomiser_output="$HOME/Documentos/Fertility/Exomas/exomiser_output"

mkdir -p "${exomiser_input}" "${exomiser_output}"

# Files
template="$HOME/Documentos/Fertility/scripts/analysis-exome-modified.yml"

# Moverte a la carpeta donde esta el exomiser
cd "$HOME/<carpeta>"

for i in $(seq -w 1 4); do 
  sample="EX26${i}"
  results="$HOME/Documentos/Fertility/Exomas/${sample}/results"
 
  # copiar todos los archivos a un solo contenedor
  cp "${results}/merge_${sample}.vcf" ${exomiser_input}
  
  # Crear yml para la corrida de exomiser
  sed "s/sample/${sample}/g" ${template} > "${exomiser_input}/analysis-exome-${sample}.yml"
  
  java -jar exomiser-cli-15.1.0.jar analyse --analysis "${exomiser_input}/analysis-exome-${sample}.yml"

done
