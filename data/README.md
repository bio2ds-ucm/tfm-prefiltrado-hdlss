# Datos

Los datos **no se incluyen** en este repositorio por su tamaño (~340 MB) y porque
proceden de una fuente pública.

## Origen

Cohorte **SCAN-B** (Dalal et al., 2022), expresión génica (RNA-seq) de tumores
primarios de mama. Disponible públicamente en el repositorio **GEO** con el
identificador **GSE202203**.

## Cómo preparar los datos

1. Descargar los datos de GSE202203 desde GEO.
2. Generar (o colocar) el objeto procesado utilizado en el TFM con el nombre exacto:

   ```
   data/GSE202203_rnaseq_data_orig_obj_tumorSize.rds
   ```

   Es un `data.frame` donde las filas son muestras y las columnas son genes, más una
   columna respuesta `obj` (1 = subtipo Luminal A, 0 = resto).

Los scripts de `scripts/ejecuciones/` esperan encontrar el fichero en esta ruta
relativa (`data/...`).
