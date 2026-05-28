library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

RESULTS_ROOT <- "resultados"
OUT_DIR      <- "figuras"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEMILLAS <- c(12345, 48271, 93054, 17689, 60532)

PREFILTER_LABELS <- c(
  "p_valor" = "p-valor",
  "mrmr"    = "MRMR",
  "cmim"    = "CMIM"
)

PREFILTER_ORDER <- c("p_valor", "mrmr", "cmim")


PREFILTER_GRADIENT <- list(
  "p_valor" = c(light = "#FCE4EC", dark = "#880E4F"),  
  "mrmr"    = c(light = "#EDE7F6", dark = "#311B92"),  
  "cmim"    = c(light = "#E3F2FD", dark = "#0D47A1")   
)


PARES <- list(
  "E" = list(
    col_nivel = "n_obs",
    tec_label = "E (observaciones reales)",
    x_label   = "Número de observaciones",
    sufijo    = "E"
  ),
  "F" = list(
    col_nivel = "n_s",
    tec_label = "F (SMOTE)",
    x_label   = "Tamaño total tras SMOTE",
    sufijo    = "F"
  ),
  "G" = list(
    col_nivel = "n_b",
    tec_label = "G (bootstrap)",
    x_label   = "Tamaño total tras bootstrap",
    sufijo    = "G"
  ),
  "H" = list(
    col_nivel = "n_n",
    tec_label = "H (ruido)",
    x_label   = "Tamaño total tras ruido",
    sufijo    = "H"
  )
)


cargar_semilla <- function(semilla) {
  path <- file.path(RESULTS_ROOT, paste0("semilla_", semilla),
                    "elastic_net_results_FINAL.rds")
  df <- readRDS(path)
  df$semilla_split <- semilla
  df
}

datos <- bind_rows(lapply(SEMILLAS, cargar_semilla)) %>%
  rename(
    procedure    = procedimiento,
    preselection = filtro,
    n_variables  = n_vars
  )

cat("Filas totales cargadas:", nrow(datos), "\n")


df_D_seed <- datos %>%
  filter(procedure == "D",
         preselection %in% PREFILTER_ORDER) %>%
  group_by(semilla_split, preselection, n_variables) %>%
  summarise(auc_D = median(auc_test, na.rm = TRUE),
            .groups = "drop")


calcular_ganancia <- function(proc, col_nivel) {
  
  df_proc_seed <- datos %>%
    filter(procedure == proc,
           preselection %in% PREFILTER_ORDER) %>%
    rename(nivel = !!sym(col_nivel)) %>%
    group_by(semilla_split, preselection, n_variables, nivel) %>%
    summarise(auc_proc = median(auc_test, na.rm = TRUE),
              .groups = "drop")
  
  # Diferencia emparejada por semilla
  df_diff_seed <- df_proc_seed %>%
    left_join(df_D_seed,
              by = c("semilla_split", "preselection", "n_variables")) %>%
    mutate(diff = auc_proc - auc_D)
  
  # Resumen entre semillas: mediana + Q1-Q3 de las diferencias
  df_diff_seed %>%
    group_by(preselection, n_variables, nivel) %>%
    summarise(
      ganancia = median(diff, na.rm = TRUE),
      q25      = quantile(diff, 0.25, na.rm = TRUE),
      q75      = quantile(diff, 0.75, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(p_factor = factor(n_variables, levels = sort(unique(n_variables))))
}

ganancias <- lapply(names(PARES), function(proc) {
  calcular_ganancia(proc, PARES[[proc]]$col_nivel)
})
names(ganancias) <- names(PARES)


y_max <- max(unlist(lapply(ganancias, function(d) {
  max(c(abs(d$q25), abs(d$q75)), na.rm = TRUE)
})))
y_lim <- c(-y_max, y_max) * 1.05
cat(sprintf("\nRango Y común (los 12 PNG): [%.4f, %.4f]\n", y_lim[1], y_lim[2]))


plot_filtro <- function(df, filtro_clave, par_cfg) {
  
  df_f <- df %>% filter(preselection == filtro_clave)
  filtro_label <- PREFILTER_LABELS[[filtro_clave]]
  
  gradient_colors <- PREFILTER_GRADIENT[[filtro_clave]]
  n_p <- length(levels(df_f$p_factor))
  p_colors <- colorRampPalette(c(gradient_colors["light"],
                                 gradient_colors["dark"]))(n_p)
  names(p_colors) <- levels(df_f$p_factor)
  
  ggplot(df_f,
         aes(x = nivel, y = ganancia,
             color = p_factor, fill = p_factor, group = p_factor)) +
    geom_hline(yintercept = 0, color = "grey40",
               linetype = "dashed", linewidth = 0.7) +
    geom_ribbon(aes(ymin = q25, ymax = q75),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.6) +
    scale_color_manual(values = p_colors,
                       name   = "Nº variables (p)") +
    scale_fill_manual(values = p_colors, guide = "none") +
    scale_x_continuous(breaks = sort(unique(df_f$nivel))) +
    coord_cartesian(ylim = y_lim) +
    labs(
      title    = paste0("Impacto del aumento de observaciones — ",
                        filtro_label),
      subtitle = paste0("Diferencia AUC del procedimiento ",
                        par_cfg$tec_label,
                        " respecto a D"),
      x        = par_cfg$x_label,
      y        = sprintf("AUC(%s) − AUC(D)", par_cfg$sufijo)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 10, color = "grey30")
    )
}


guardar_png <- function(p, file, w = 8, h = 5.5) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = 150)
  cat("  →", file, "\n")
}

cat("\n[Generando los 12 gráficos]\n")

for (proc in names(PARES)) {
  par_cfg <- PARES[[proc]]
  df_g    <- ganancias[[proc]]
  
  for (filtro in PREFILTER_ORDER) {
    p <- plot_filtro(df_g, filtro, par_cfg)
    guardar_png(
      p,
      file.path(OUT_DIR,
                sprintf("P3_ganancia_%s_vs_D_%s.png",
                        par_cfg$sufijo, filtro))
    )
  }
}

cat("\n===== FIGURAS GENERADAS =====\n")
cat("Directorio:", OUT_DIR, "\n")
cat("Archivos: P3_ganancia_{E,F,G,H}_vs_D_{p_valor,mrmr,cmim}.png\n")
