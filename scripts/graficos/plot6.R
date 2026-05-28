library(dplyr)
library(ggplot2)
library(tidyr)


RESULTS_ROOT <- "resultados"
OUT_DIR      <- "figuras"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEMILLAS_SPLIT <- c(12345, 48271, 93054, 17689, 60532)


UMBRAL_FREC <- 0.5


NIVEL_AUMENTO <- 130


cargar_datos <- function(semilla, archivo) {
  path <- file.path(RESULTS_ROOT, paste0("semilla_", semilla), archivo)
  df <- readRDS(path)
  df$semilla_split <- semilla
  df
}


summary_all <- bind_rows(lapply(SEMILLAS_SPLIT, function(sem)
  cargar_datos(sem, "elastic_net_results_FINAL.rds")))


coefs_all <- bind_rows(lapply(SEMILLAS_SPLIT, function(sem)
  cargar_datos(sem, "coefs_results_FINAL.rds")))

cat("Summary filas:", nrow(summary_all), "\n")
cat("Coefs filas:  ", nrow(coefs_all), "\n")


identificar_p_optimo <- function(proc_letra, col_nivel) {
  summary_all %>%
    filter(procedimiento == proc_letra,
           filtro %in% c("p_valor", "mrmr", "cmim"),
           !!sym(col_nivel) == NIVEL_AUMENTO) %>%
    group_by(filtro, n_vars) %>%
    summarise(median_auc = median(auc_test, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(filtro) %>%
    slice_max(median_auc, n = 1) %>%
    ungroup() %>%
    select(filtro, n_vars_opt = n_vars, median_auc)
}

p_opt_E <- identificar_p_optimo("E", "n_obs")
p_opt_F <- identificar_p_optimo("F", "n_s")

cat("\n=== p óptimo para E (n=130) ===\n"); print(p_opt_E)
cat("\n=== p óptimo para F (n=130) ===\n"); print(p_opt_F)


obtener_vars_lasso <- function(proc_letra, filtro_v, n_vars_v, col_nivel = NULL,
                               nivel_v = NULL) {
  s <- summary_all %>%
    filter(procedimiento == proc_letra) %>%
    filter(filtro == filtro_v | is.na(filtro)) %>%
    filter(n_vars == n_vars_v | is.na(n_vars))
  
  if (!is.null(col_nivel) && !is.null(nivel_v)) {
    s <- s %>% filter(!!sym(col_nivel) == nivel_v)
  }
  
  ids_relevantes <- unique(s$execution_id)
  
  c <- coefs_all %>%
    filter(execution_id %in% ids_relevantes,
           s == "lambda.min",
           variable != "(Intercept)",
           coef != 0)
  

  N_ejec <- length(ids_relevantes)
  
  c %>%
    group_by(variable) %>%
    summarise(n_apariciones = n(),
              .groups = "drop") %>%
    mutate(frecuencia = n_apariciones / N_ejec,
           N_total    = N_ejec)
}

vars_C <- obtener_vars_lasso("C", filtro_v = NA, n_vars_v = NA)
cat(sprintf("\nVariables en C (N_total=%d): %d distintas\n",
            unique(vars_C$N_total), nrow(vars_C)))

vars_E_pvalor <- obtener_vars_lasso("E", "p_valor",
                                    p_opt_E %>% filter(filtro == "p_valor") %>% pull(n_vars_opt),
                                    "n_obs", NIVEL_AUMENTO)
vars_E_mrmr   <- obtener_vars_lasso("E", "mrmr",
                                    p_opt_E %>% filter(filtro == "mrmr") %>% pull(n_vars_opt),
                                    "n_obs", NIVEL_AUMENTO)
vars_E_cmim   <- obtener_vars_lasso("E", "cmim",
                                    p_opt_E %>% filter(filtro == "cmim") %>% pull(n_vars_opt),
                                    "n_obs", NIVEL_AUMENTO)

vars_F_pvalor <- obtener_vars_lasso("F", "p_valor",
                                    p_opt_F %>% filter(filtro == "p_valor") %>% pull(n_vars_opt),
                                    "n_s", NIVEL_AUMENTO)
vars_F_mrmr   <- obtener_vars_lasso("F", "mrmr",
                                    p_opt_F %>% filter(filtro == "mrmr") %>% pull(n_vars_opt),
                                    "n_s", NIVEL_AUMENTO)
vars_F_cmim   <- obtener_vars_lasso("F", "cmim",
                                    p_opt_F %>% filter(filtro == "cmim") %>% pull(n_vars_opt),
                                    "n_s", NIVEL_AUMENTO)



aplicar_umbral <- function(df, umbral = UMBRAL_FREC) {
  df %>% filter(frecuencia >= umbral) %>% pull(variable)
}

set_C        <- aplicar_umbral(vars_C)
set_E_pvalor <- aplicar_umbral(vars_E_pvalor)
set_E_mrmr   <- aplicar_umbral(vars_E_mrmr)
set_E_cmim   <- aplicar_umbral(vars_E_cmim)
set_F_pvalor <- aplicar_umbral(vars_F_pvalor)
set_F_mrmr   <- aplicar_umbral(vars_F_mrmr)
set_F_cmim   <- aplicar_umbral(vars_F_cmim)

cat(sprintf("\n=== Variables seleccionadas (umbral %.0f%%) ===\n", UMBRAL_FREC * 100))
cat(sprintf("  C:           %d variables\n", length(set_C)))
cat(sprintf("  E p-valor:   %d variables\n", length(set_E_pvalor)))
cat(sprintf("  E MRMR:      %d variables\n", length(set_E_mrmr)))
cat(sprintf("  E CMIM:      %d variables\n", length(set_E_cmim)))
cat(sprintf("  F p-valor:   %d variables\n", length(set_F_pvalor)))
cat(sprintf("  F MRMR:      %d variables\n", length(set_F_mrmr)))
cat(sprintf("  F CMIM:      %d variables\n", length(set_F_cmim)))

tabla_solapamiento <- function(sets, nombres) {
  n <- length(sets)
  m <- matrix(NA, n, n)
  rownames(m) <- nombres
  colnames(m) <- nombres
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) m[i, j] <- length(sets[[i]])
      else m[i, j] <- length(intersect(sets[[i]], sets[[j]]))
    }
  }
  m
}

cat("\n=== TABLA SOLAPAMIENTO C vs E ===\n")
tabla_E <- tabla_solapamiento(
  list(set_C, set_E_pvalor, set_E_mrmr, set_E_cmim),
  c("C", "E_p_valor", "E_MRMR", "E_CMIM")
)
print(tabla_E)
write.csv(tabla_E, file.path(OUT_DIR, "tabla_solapamiento_C_vs_E.csv"))

cat("\n=== TABLA SOLAPAMIENTO C vs F ===\n")
tabla_F <- tabla_solapamiento(
  list(set_C, set_F_pvalor, set_F_mrmr, set_F_cmim),
  c("C", "F_p_valor", "F_MRMR", "F_CMIM")
)
print(tabla_F)
write.csv(tabla_F, file.path(OUT_DIR, "tabla_solapamiento_C_vs_F.csv"))


construir_heatmap_data <- function(vars_list, nombres) {
  df_long <- bind_rows(lapply(seq_along(vars_list), function(i) {
    vars_list[[i]] %>% mutate(procedimiento = nombres[i])
  }))
  
  vars_top <- df_long %>%
    filter(frecuencia >= UMBRAL_FREC) %>%
    pull(variable) %>%
    unique()
  
  df_long %>%
    filter(variable %in% vars_top) %>%
    complete(variable = vars_top, procedimiento = nombres,
             fill = list(frecuencia = 0))
}

plot_heatmap_vars <- function(df_long, titulo, file) {
  
  orden_vars <- df_long %>%
    group_by(variable) %>%
    summarise(mean_frec = mean(frecuencia, na.rm = TRUE)) %>%
    arrange(desc(mean_frec)) %>%
    pull(variable)
  
  df_long$variable <- factor(df_long$variable, levels = rev(orden_vars))
  
  p <- ggplot(df_long, aes(x = procedimiento, y = variable, fill = frecuencia)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", frecuencia)),
              size = 3, color = "grey20") +
    scale_fill_gradient(low = "#FFFFFF", high = "#EC407A",
                        limits = c(0, 1),
                        name = "Frecuencia") +
    labs(
      title    = titulo,
      subtitle = sprintf("Variables con frecuencia ≥ %.0f%% en al menos un procedimiento",
                         UMBRAL_FREC * 100),
      x        = NULL,
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      panel.grid       = element_blank(),
      legend.position  = "right"
    )
  
  altura <- max(4, 0.25 * length(unique(df_long$variable)) + 2)
  ggsave(file, p, width = 8, height = altura, dpi = 150)
  cat("  →", file, "\n")
}

cat("\n[Generando heatmaps]\n")

df_E <- construir_heatmap_data(
  list(vars_C, vars_E_pvalor, vars_E_mrmr, vars_E_cmim),
  c("C", "E_p_valor", "E_MRMR", "E_CMIM")
)
plot_heatmap_vars(df_E, "Variables seleccionadas — C vs E",
                  file.path(OUT_DIR, "P6_heatmap_vars_C_vs_E.png"))

df_F <- construir_heatmap_data(
  list(vars_C, vars_F_pvalor, vars_F_mrmr, vars_F_cmim),
  c("C", "F_p_valor", "F_MRMR", "F_CMIM")
)
plot_heatmap_vars(df_F, "Variables seleccionadas — C vs F",
                  file.path(OUT_DIR, "P6_heatmap_vars_C_vs_F.png"))

cat("\n===== ANÁLISIS COMPLETADO =====\n")
cat("Directorio:", OUT_DIR, "\n")
