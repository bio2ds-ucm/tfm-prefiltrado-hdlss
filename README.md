# Estudio de técnicas de pre-filtrado en modelos de machine learning para manejar datos de alta-dimensionalidad

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/status-in%20preparation-orange.svg)

Trabajo Fin de Máster — Máster en Bioestadística, Facultad de Estudios Estadísticos, Universidad Complutense de Madrid.

**Autoras:** Cristina Ordaz 
**Tutoras:** Aida Calviño, Pedro Contró y Silvia Pineda

## Descripción

Este repositorio contiene el código y los resultados de un estudio sobre técnicas de prefiltrado de variables (p-valor, MRMR, CMIM) y aumento de datos (SMOTE, Bootstrap, Ruido) aplicadas a datos transcriptómicos de expresión génica de tumores primarios de mama (cohorte SCAN-B), en un escenario HDLSS (*High Dimension, Low Sample Size*) donde p ≫ n. La variable respuesta indica si el tumor es de subtipo Luminal A (1) o no (0).

El modelo base es una regresión logística con regularización Lasso (L1).

## Procedimientos evaluados

Siguiendo el diagrama de repeticiones del estudio:

| Procedimiento | Descripción |
|---|---|
| A | Baseline: train completo, sin prefiltrado |
| B | Train completo + prefiltrado de variables |
| C | Submuestra (n = 100), sin prefiltrado |
| D | Submuestra (n = 100) + prefiltrado |
| F | Submuestra + aumento con SMOTE |
| G | Submuestra + aumento con Bootstrap |
| H | Submuestra + aumento con Ruido |

- Procedimientos **A y B**: 5 repeticiones (una por semilla, train completo).
- Procedimientos **C, D, F, G, H**: 100 repeticiones (5 semillas × 20 submuestras de n = 100).

## Estructura del repositorio

```
.
├── data/                  Instrucciones para obtener los datos (los .rds NO se versionan)
├── scripts/
│   ├── ejecuciones/       5 scripts, uno por semilla
│   └── graficos/          Scripts que generan las figuras del póster
├── resultados/            Se rellena al ejecutar los scripts (los .rds NO se versionan)
├── figuras/               Figuras finales (PNG)
└── doc/                   Póster en PDF
```

## Datos

Los datos provienen de la cohorte SCAN-B (Dalal et al., 2022), disponibles públicamente en el repositorio GEO (identificador **GSE202203**). Por su tamaño (~340 MB), el archivo `.rds` procesado **no se incluye** en este repositorio. En `data/README.md` se explica cómo obtenerlo y con qué nombre colocarlo para que los scripts lo encuentren.

## Cómo reproducir los resultados

### Requisitos

- R (versión utilizada: **completar con la versión de R con la que se ejecutó**).
- Paquetes principales: `glmnet`, `praznik` (filtros MRMR y CMIM), `smotefamily` o `themis` (SMOTE), `dplyr`, `ggplot2`, `tidyr`.


### Pasos

1. Colocar el archivo de datos en `data/` siguiendo las instrucciones de `data/README.md`.
2. Ejecutar los 5 scripts de `scripts/ejecuciones/`. Cada uno crea una subcarpeta `resultados/semilla_X/` con sus salidas (`elastic_net_results_FINAL.rds`, `coefs_results_FINAL.rds`, `vars_input_results_FINAL.rds`).
3. Ejecutar los scripts de `scripts/graficos/` para generar las figuras (en `figuras/`) a partir de esos `.rds`.

### Nota sobre los scripts de ejecución

Los 5 scripts de `scripts/ejecuciones/` son prácticamente idénticos y solo difieren en la semilla utilizada. Se han mantenido como archivos separados para facilitar la trazabilidad de cada ejecución. Si se modifica algo en uno, conviene replicar el cambio en los otros cuatro.

## Resultados principales

Resumen de lo observado en el estudio (ver póster en `doc/` para detalle):

- En **train completo**, el prefiltrado (Procedimiento B) alcanza un rendimiento equiparable al baseline (A) cuando se retienen ≥ 1000 variables. CMIM es el filtro más eficaz; el p-valor el más débil con pocas variables.
- El **solapamiento entre filtros es muy bajo**: con 50 variables ninguna es compartida por los tres métodos; con 1000, solo 89 lo son.
- Con **muestra reducida** (n = 100), ninguna configuración alcanza el rendimiento del baseline. El prefiltrado no compensa la pérdida de tamaño muestral.
- Las técnicas de **aumento de datos** evaluadas (SMOTE, Bootstrap, Ruido) producen AUCs en torno a 0.60–0.61, sin mejorar de forma clara el escenario de submuestra.

## Contacto

cordaz@ucm.es — Facultad de Estudios Estadísticos, UCM.
pcontro@ucm.es - Facultad de Estudios Estadísticos, UCM.
aida.calvino@ucm.es - Facultad de Estudios Estadísticos, UCM.
sipineda@ucm.es - Facultad de Estudios Estadísticos, UCM.
