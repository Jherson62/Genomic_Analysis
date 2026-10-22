#!/bin/bash
set -eou pipefail

# Entorno NGStools
source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate NGStools

# Directorios
exomiser_input="$HOME/Documents/Exomas/exomiser_input"
exomiser_output="$HOME/Documents/Exomas/exomiser_output"

mkdir -p "${exomiser_input}" "${exomiser_output}"

# Files
template="$HOME/Documents/Exomas/templates/analysis-exome-modified.yml"

# Entramos a la carpeta
cd "$HOME/exomiser-cli-15.1.0"

for i in $(seq -w 1); do 
  sample="EX26650"
  
  echo ${sample}
  
  results="$HOME/Documents/Exomas/${sample}/results"
 
  # Mover todos los archivos a un solo contenedor
  cp "${results}/merge_${sample}.vcf" ${exomiser_input}
  
  # Crear yml para la corrida de exomiser
  sed "s/sample/${sample}/g" ${template} > "${exomiser_input}/analysis-exome-${sample}.yml"
  
  java -jar exomiser-cli-15.1.0.jar analyse --analysis "${exomiser_input}/analysis-exome-${sample}.yml"; done
