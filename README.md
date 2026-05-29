# Estudio de técnicas de pre-filtrado en modelos de machine learning para manejar datos de alta-dimensionalidad

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20443008.svg)](https://doi.org/10.5281/zenodo.20443008)

Trabajo Fin de Máster — Máster en Bioestadística  
Facultad de Estudios Estadísticos · Universidad Complutense de Madrid

|  |  |
|---|---|
| 🎓 **Autora** | Cristina Ordaz |
| 🧭 **Tutores** | Pedro Contró, Aida Calviño, Silvia Pineda |

## 📌 Descripción

Este repositorio contiene el código y los resultados de un estudio sobre técnicas de prefiltrado de variables (p-valor, MRMR, CMIM) y aumento de datos (SMOTE, Bootstrap, Ruido) aplicadas a datos transcriptómicos de expresión génica de tumores primarios de mama (cohorte SCAN-B), en un escenario HDLSS (*High Dimension, Low Sample Size*) donde p ≫ n. La variable respuesta indica si el tumor es de subtipo Luminal A (1) o no (0).

El modelo base es una regresión logística con regularización Lasso (L1).

## ⚙️ Procedimientos evaluados

Siguiendo el diagrama de repeticiones del estudio:

| Procedimiento | Descripción |
|---|---|
| A | Baseline: train completo, sin prefiltrado |
| B | Train completo + prefiltrado de variables |
| C | Submuestra (n = 100), sin prefiltrado |
| D | Submuestra (n = 100) + prefiltrado |
| E | Submuestra + aumento con datos reales |
| F | Submuestra + aumento con SMOTE |
| G | Submuestra + aumento con Bootstrap |
| H | Submuestra + aumento con Ruido |

- Procedimientos **A y B**: 5 repeticiones (una por semilla, train completo).
- Procedimientos **C, D, E, F, G, H**: 100 repeticiones (5 semillas × 20 submuestras de n = 100).

## 📁 Estructura del repositorio

```
.
├── data/                  Instrucciones para descargar el dataset
├── scripts/
│   ├── ejecuciones/       5 scripts, uno por semilla
│   └── graficos/          Scripts que generan las figuras
├── resultados/            Se rellena al ejecutar los scripts (los .rds NO se versionan)
├── figuras/               Figuras finales (PNG)
└── doc/                   Póster en PDF y TFM en PDF
```

## 🧬 Datos

Los datos provienen de la cohorte SCAN-B (Dalal et al., 2022), disponibles públicamente en el repositorio GEO (identificador **GSE202203**). Debido al tamaño del archivo `.rds` (~340 MB) procesado se incluye en este repositorio las instrucciones para acceder al dataset en la carpeta `data/`. 

## ▶️ Cómo reproducir los resultados

### Requisitos

- R (versión utilizada: 4.5.2.
- Paquetes principales: `glmnet`, `praznik`, `smotefamily`, `dplyr`, `ggplot2`, `tidyr`, `caret`, `pROC`, `doParallel`, `foreach`, `mRMRe`.


### Pasos

1. Colocar el archivo de datos en `data/` siguiendo las instrucciones de `data/README.md`.
2. Ejecutar los 5 scripts de `scripts/ejecuciones/`. Cada uno crea una subcarpeta `resultados/semilla_X/` con sus salidas (`elastic_net_results_FINAL.rds`, `coefs_results_FINAL.rds`, `vars_input_results_FINAL.rds`).
3. Ejecutar los scripts de `scripts/graficos/` para generar las figuras (en `figuras/`) a partir de esos `.rds`.

### Nota sobre los scripts de ejecución

Los 5 scripts de `scripts/ejecuciones/` son prácticamente idénticos y solo difieren en la semilla utilizada. Se han mantenido como archivos separados para facilitar la trazabilidad de cada ejecución. Si se modifica algo en uno, conviene replicar el cambio en los otros cuatro.

## 📊 Resultados principales

Resumen de lo observado en el estudio:

- **p-valor** y **MRMR** son los filtros más sólidos. El p-valor rinde de forma consistente y se adapta mejor a la escasez de observaciones; MRMR es la alternativa multivariante más estable. La elección depende del objetivo.
- **CMIM** es el mejor con muestra grande pero frágil con muestra reducida, y su combinación con SMOTE se descarta. La **selección aleatoria** queda descartada por rendir siempre peor.
- La **selección de variables** es muy inestable y los filtros son complementarios. Pocas variables consistentes con n=100 y cada filtro favorece variables distintas.
- Las **observaciones reales** (E) dan las mayores ganancias; cuando no son viables, **SMOTE** es la mejor opción sintética, sobre todo con p-valor o MRMR. Bootstrap y ruido fueron descartados.
- **Como recomendación**, prefiltrar con p-valor o MRMR junto con datos reales; si no es posible, SMOTE con cualquiera de ellos. Evitar CMIM × SMOTE.

## 📑 Citación

```bibtex
@software{ordaz2026code,
  author    = {Ordaz, Cristina},
  title     = {Code for: Estudio de técnicas de pre-filtrado en modelos de machine learning para manejar datos de alta-dimensionalidad},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {[Zenodo DOI]},
  url       = {https://github.com/bio2ds-ucm/tfm-prefiltrado-hdlss}
}
```

See [`CITATION.cff`](CITATION.cff).

## 📄 License

This code is released under the [MIT License](LICENSE).

## 📬 Contacto

cordaz@ucm.es — Facultad de Estudios Estadísticos, UCM.
pcontro@ucm.es - Facultad de Estudios Estadísticos, UCM.
aida.calvino@ucm.es - Facultad de Estudios Estadísticos, UCM.
sipineda@ucm.es - Facultad de Estudios Estadísticos, UCM.
