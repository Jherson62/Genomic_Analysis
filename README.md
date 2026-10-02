# PROYECTO LONGEVIDAD (LIFESPAN)

La longevidad es un rasgo complejo, muchos genes interactúan para brindarte la información fenotípica
por esta razón, se utilizará un Poligenic Score. Se toma una lista de variantes a los que se les asigna
un "peso" o tamaño de efecto y se suma para analizar su predisposición genética.

* Directorio: **LifeSpan**
  Se encuentra los datos descargados de la pagina [PGS Catalog](https://www.pgscatalog.org/)
*  Script: **life_span.R**

## Primer intento 
Polygenic Score ID: [PGS000906](https://www.pgscatalog.org/score/PGS000906/)

- Rasgo para la determinación de la esperanza de vida
- Publicación: "Polygenic Risk Score of Longevity Predicts Longer Survival Across an Age Continuum"
- DOI (2021): [10.1093/gerona/glaa289](10.1093/gerona/glaa289)
- Creado con GRCh37 (hg19) pero contiene archivos armonizado para ambos hg19 y hg38
- 330 variantes
- Ancestria: Europea
- Survival analysis performed in the control samples used to develop the PRS

**Archivos**

- PGS000906.txt.gz: Archivo original. Contiene variantes y su peso.
- PGS000906_hmPOS_GRCh37/GRCh38.txt: Archivo con posiciones estadarizadas según genoma de ref.

**Metadata**

- Modelo: Aditivo lineal simple -> los pesos son coeficientes de regresión. Solo se suman.
- *Sumatoria snp x peso.*
- No hay interacción entre SNPs. No hay correlación entre SNPs.
- ***Problema***: Es una data europea, las interaciones los ligamientos entre SNPs cambian cuando 
hablamos de población Peruana por ejemplo.



# PROYECTO TRISOMIA EN EXOMAS
objetivo: Analizar si con exomas puedo tener información sobre trisomia o monosomia del crh22.

# DORADO BASECALLER

## METILACIÓN
Un archivo POD5 almacena la señal eléctrica cruda de la secuenciación. Para obtener los datos de
modificaciones, dorado analiza la señal eléctrica utilizando un modelo de redes neuronales diferentes
al convencional.
Escribe las etiquetas MM (Modified base) y ML (Modification Likelihood/Probability) en el archivo BAM.

Para ver la lista total de modelos `dorado download --list`

Dentro del información de la secuenciación encontramos los siguientes datos:

**protocol=sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-NBD114-24:400**

*   MIN114: tecnologia MinION con química V14
*   e8_2: Versión de la enzima motor (E8.2 protein)
*   400K: velocidad de traslocación del poro (400bps)
*   barcode kit: SQK-NBD114-24 (Native Barcoding Kit 24 V14)
*   FLO-MIN114: celda de flujo R10.4.1
*   dna_r10.4.1: muestra de ADN corrido en esa celda.

**Duplex calling**

Estrategia bioinformatica diseñada para elevar la precisión de la lectura (ayuda a superar Q30 o 99.9% accuracy). 
En celdas con químicas actuales (kit V14 y celdas R10.4.1) una porcion de ADN se secuencia a doble hebra, primero pasa una e inmediatamente pasa a secuenciarse la otra cadena. Este artificio permite obtener información de ambas hebras para el basecalling. No aplicable en Rapid Sequencing kit.

*Nota: No puedes utilizar la demultiplexación (separar por barcodes) al mismo tiempo que el duplex calling. Se haría primero el duplex y luego la demultiplexación*

`dorado demux --kit-name SQK-NBD114-24 --output-dir ./demux_duplex duplex_calls.bam`

El estándar moderno de Oxford Nanopore para almacenar lecturas demultiplexadas es el formato uBAM (BAM sin alinear).

***Native Barcoding***: La ligación ocurre en ambos extremos de la molécula. 
***Rapid Barcoding***: Solo se etiqueta un extremos con el barcode.


