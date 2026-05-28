# Estudio de técnicas de pre-filtrado en modelos de machine learning para manejar datos de alta-dimensionalidad

Trabajo Fin de Máster — Máster en Bioestadística, Universidad Complutense de Madrid.

**Autores:** Cristina Ordaz, Pedro Contró
**Tutoras:** Aida Calviño, Silvia Pineda

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

## Bibliografía

1. Fernández-Casal R, Costa J, Oviedo M. *Métodos predictivos de aprendizaje estadístico*. Servizo de Publicacións, Universidade da Coruña; 2024. <https://rubenfcasal.github.io/aprendizaje_estadistico/>
2. Bühlmann P, Van de Geer S. *Statistics for high-dimensional data: methods, theory and applications*. Springer Series in Statistics. Berlin, Heidelberg: Springer Berlin Heidelberg; 2011.
3. Giraud C. *Introduction to high-dimensional statistics*. 2nd ed. Chapman & Hall/CRC Monographs on Statistics and Applied Probability. Chapman and Hall/CRC; 2021.
4. Fernández-Tortolero Á, Reigosa-Yániz A. Subtipos del carcinoma luminal de mama según el consenso de Saint Gallen en un grupo de pacientes venezolanas. *Biomédica*. 2021;41(3):531-40.
5. Gaona-Romero C, Domínguez-Recio ME, Comino-Méndez I, Ortega-Jiménez MV, Lavado-Valenzuela R, Alba E. Luminal and basal subtypes across carcinomas: molecular programs beyond tissue of origin. *Cancers*. 2025;17(16):2720.
6. Chawla NV, Bowyer KW, Hall LO, Kegelmeyer WP. SMOTE: synthetic minority over-sampling technique. *Journal of Artificial Intelligence Research*. 2002;16:321-57.
7. Calviño Martínez A, Alonso Revenga JM. *Introducción a la ciencia de datos con R*. García Maroto Editores; 2025.
8. Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. A simulation study of the number of events per variable in logistic regression analysis. *Journal of Clinical Epidemiology*. 1996;49(12):1373-9.
9. Sur P, Candès EJ. A modern maximum-likelihood theory for high-dimensional logistic regression. *Proceedings of the National Academy of Sciences of the United States of America*. 2019;116(29):14516-25.
10. Albert A, Anderson JA. On the existence of maximum likelihood estimates in logistic regression models. *Biometrika*. 1984;71(1):1-10.
11. Lacan A. *Transcriptomics data generation with deep generative models* [Thèse de doctorat]. Université Paris-Saclay; 2025. <https://theses.hal.science/tel-04996930>
12. García-Vicente C, Chushig-Muzo D, Mora-Jiménez I, Fabelo H, Gram IT, Løchen ML, et al. Evaluation of synthetic categorical data generation techniques for predicting cardiovascular diseases and post-hoc interpretability of the risk factors. *Applied Sciences*. 2023;13(7):4119.
13. Gómez Deraves Á, Gómez Marquina K. *Muestreo estadístico para docentes y estudiantes*. 1st ed. Tecana American University; 2019. <https://tauniversity.org/libros/muestreo-estadistico-para-docentes-y-estudiantes>
14. Belío Miranda J. *Métodos Bootstrap y sus aplicaciones* [Trabajo Fin de Grado]. Universidad de Zaragoza; 2020. <https://zaguan.unizar.es/record/98153>
15. Tsai TI, Li DC. Utilize bootstrap in small data set learning for pilot run modeling of manufacturing systems. *Expert Systems with Applications*. 2008;35(3):1293-300.
16. Bolón-Canedo V, Sánchez-Maroño N, Alonso-Betanzos A. *Feature Selection for High-Dimensional Data*. Artificial Intelligence: Foundations, Theory, and Algorithms. Cham: Springer International Publishing; 2015.
17. Bommert A, Welchowski T, Schmid M, Rahnenführer J. Benchmark of filter methods for feature selection in high-dimensional gene expression survival data. *Briefings in Bioinformatics*. 2022;23(1):bbab354.
18. Peng H, Long F, Ding C. Feature selection based on mutual information: criteria of max-dependency, max-relevance, and min-redundancy. *IEEE Transactions on Pattern Analysis and Machine Intelligence*. 2005;27(8):1226-38.
19. Méndez López ÁJ. *Selección de características en entornos de dimensionalidad masiva y muestra reducida: aplicación al estudio de metilación del ADN en la enfermedad de Alzheimer* [Tesis doctoral]. Universidad Autónoma de Madrid; 2024. <https://hdl.handle.net/10486/715205>
20. EL-Manzalawy Y, Hsieh TY, Shivakumar M, Kim D, Honavar V. Min-redundancy and max-relevance multi-view feature selection for predicting ovarian cancer survival using multi-omics data. *BMC Medical Genomics*. 2018;11(S3):71.
21. Fleuret F. Fast binary feature selection with conditional mutual information. *Journal of Machine Learning Research*. 2004;5:1531-55.
22. Dalal H, Dahlgren M, Gladchuk S, Brueffer C, Gruvberger-Saal SK, Saal LH. Clinical associations of ESR2 (estrogen receptor beta) expression across thousands of primary breast tumors. *Scientific Reports*. 2022;12(1):4696.
23. Hastie T, Tibshirani R, Friedman J. *The elements of statistical learning: data mining, inference, and prediction*. 2nd ed. Springer Series in Statistics. New York: Springer; 2009.
24. Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. *Radiology*. 1982;143(1):29-36.
25. Dong X, Yang Y, Xu G, Tian Z, Yang Q, Gong Y, et al. The initial expression alterations occurring to transcription factors during the formation of breast cancer: evidence from bioinformatics. *Cancer Medicine*. 2022;11(5):1371-95.
26. Binato R, Corrêa S, Panis C, Ferreira G, Petrone I, da Costa IR, et al. NRIP1 is activated by C-JUN/C-FOS and activates the expression of PGR, ESR1 and CCND1 in luminal A breast cancer. *Scientific Reports*. 2021;11(1):21159.

## Contacto

cordaz@ucm.es — Facultad de Estudios Estadísticos, UCM.
