# Resultados

Esta carpeta **no se versiona** (los `.rds` se excluyen en `.gitignore`) porque los
resultados se regeneran ejecutando los scripts.

Al ejecutar cada script de `scripts/ejecuciones/` se crea automáticamente una
subcarpeta por semilla con sus salidas:

```
resultados/
├── semilla_12345/
│   ├── elastic_net_results_FINAL.rds
│   ├── coefs_results_FINAL.rds
│   └── vars_input_results_FINAL.rds
├── semilla_17689/
├── semilla_48271/
├── semilla_60532/
└── semilla_93054/
```

Los scripts de `scripts/graficos/` leen estos `.rds` desde aquí y escriben las
figuras en `figuras/`.
