library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)



RESULTS_ROOT <- "resultados"
OUT_DIR      <- "figuras"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEMILLAS <- c(12345, 48271, 93054, 17689, 60532)

N_VARS_B <- c(50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000)

PREFILTER_COLORS <- c(
  "mrmr"      = "#B39DDB",
  "cmim"      = "#90CAF9",
  "p_valor"   = "#F4A7B9",
  "aleatorio" = "#A5D6A7"
)

PREFILTER_LABELS <- c(
  "mrmr"      = "MRMR",
  "cmim"      = "CMIM",
  "p_valor"   = "p-valor",
  "aleatorio" = "Aleatorio"
)

PREFILTER_ORDER <- c("aleatorio", "p_valor", "cmim", "mrmr")



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

# n_variables que efectivamente usa D (puede diferir de B)
N_VARS_D <- datos %>%
  filter(procedure == "D") %>%
  pull(n_variables) %>%
  unique() %>%
  sort()

cat("n_variables en B:", paste(N_VARS_B, collapse = ", "), "\n")
cat("n_variables en D:", paste(N_VARS_D, collapse = ", "), "\n")



guardar_png <- function(p, file, w = 10, h = 5) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = 150)
  cat("  →", file, "\n")
}

resumen_lineas <- function(df, proc, n_vars_levels) {
  df_seed <- df %>%
    filter(procedure == proc,
           preselection %in% PREFILTER_ORDER) %>%
    group_by(semilla_split, preselection, n_variables) %>%
    summarise(median_auc = median(auc_test, na.rm = TRUE),
              .groups    = "drop")

  df_seed %>%
    group_by(preselection, n_variables) %>%
    summarise(
      median_val = median(median_auc, na.rm = TRUE),
      q25        = quantile(median_auc, 0.25, na.rm = TRUE),
      q75        = quantile(median_auc, 0.75, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      n_variables  = factor(n_variables, levels = n_vars_levels),
      preselection = factor(preselection, levels = PREFILTER_ORDER)
    )
}

resumen_barras <- function(df, proc) {
  df_seed <- df %>%
    filter(procedure == proc,
           preselection %in% PREFILTER_ORDER) %>%
    group_by(semilla_split, preselection) %>%
    summarise(median_auc = median(auc_test, na.rm = TRUE),
              .groups    = "drop")

  df_seed %>%
    group_by(preselection) %>%
    summarise(
      median_val = median(median_auc, na.rm = TRUE),
      q25        = quantile(median_auc, 0.25, na.rm = TRUE),
      q75        = quantile(median_auc, 0.75, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    arrange(median_val) %>%
    mutate(preselection = factor(preselection, levels = preselection))
}

stats_baseline <- function(df, proc) {
  df_seed <- df %>%
    filter(procedure == proc) %>%
    group_by(semilla_split) %>%
    summarise(median_auc = median(auc_test, na.rm = TRUE),
              .groups    = "drop")

  list(
    median = median(df_seed$median_auc, na.rm = TRUE),
    q25    = quantile(df_seed$median_auc, 0.25, na.rm = TRUE),
    q75    = quantile(df_seed$median_auc, 0.75, na.rm = TRUE)
  )
}

plot_lineas <- function(df_summary, baseline, ylim_extra = 0.01,
                        title, subtitle, caption,
                        x_lab = "Número de variables") {
  ylim_y <- c(
    min(c(df_summary$q25, baseline$q25), na.rm = TRUE) - ylim_extra,
    max(c(df_summary$q75, baseline$q75), na.rm = TRUE) + ylim_extra
  )
  pd <- position_dodge(width = 0.4)

  ggplot(df_summary,
         aes(x = n_variables, y = median_val,
             color = preselection, group = preselection)) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = baseline$q25, ymax = baseline$q75,
             fill = "grey70", alpha = 0.3) +
    geom_hline(yintercept = baseline$median,
               color = "grey40", linetype = "dashed", linewidth = 0.7) +
    geom_errorbar(aes(ymin = q25, ymax = q75),
                  width = 0.25, linewidth = 0.7, position = pd) +
    geom_line(linewidth = 1.1, position = pd) +
    geom_point(size = 2.8, position = pd) +
    scale_color_manual(
      values = PREFILTER_COLORS[PREFILTER_ORDER],
      labels = PREFILTER_LABELS[PREFILTER_ORDER]
    ) +
    coord_cartesian(ylim = ylim_y) +
    labs(
      title    = title,
      subtitle = subtitle,
      caption  = caption,
      x        = x_lab,
      y        = "AUC test",
      color    = "Filtro"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 10, color = "grey30")
    )
}

plot_barras <- function(df_summary, baseline,
                        title, subtitle, caption) {
  xlim_x <- c(
    min(c(df_summary$q25, baseline$q25), na.rm = TRUE) - 0.01,
    max(c(df_summary$q75, baseline$q75), na.rm = TRUE) + 0.025
  )

  ggplot(df_summary,
         aes(x = median_val, y = preselection, fill = preselection)) +
    annotate("rect",
             xmin = baseline$q25, xmax = baseline$q75,
             ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.3) +
    geom_vline(xintercept = baseline$median,
               color = "grey40", linetype = "dashed", linewidth = 0.7) +
    geom_col(width = 0.65, color = NA) +
    geom_errorbarh(aes(xmin = q25, xmax = q75),
                   height = 0.25, linewidth = 0.7, color = "grey30") +
    geom_text(aes(label = sprintf("%.3f", median_val),
                  x = q75 + 0.003),
              hjust = 0, size = 4, color = "grey20") +
    scale_fill_manual(values = PREFILTER_COLORS,
                      labels = PREFILTER_LABELS,
                      guide  = "none") +
    scale_y_discrete(labels = PREFILTER_LABELS) +
    coord_cartesian(xlim = xlim_x) +
    labs(
      title    = title,
      subtitle = subtitle,
      caption  = caption,
      x        = "AUC test",
      y        = "Filtro"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 12),
      plot.caption       = element_text(size = 10, color = "grey30")
    )
}



A <- stats_baseline(datos, "A")
C <- stats_baseline(datos, "C")

cat(sprintf("\nProcedimiento A: mediana = %.4f, Q1 = %.4f, Q3 = %.4f\n",
            A$median, A$q25, A$q75))
cat(sprintf("Procedimiento C: mediana = %.4f, Q1 = %.4f, Q3 = %.4f\n",
            C$median, C$q25, C$q75))



cat("\n[Gráfico 1] A vs B — líneas\n")

df_lineas_B <- resumen_lineas(datos, "B", N_VARS_B)

p_B_lineas <- plot_lineas(
  df_summary = df_lineas_B,
  baseline   = A,
  title      = "Comparación entre métodos de prefiltrado",
  subtitle   = "Procedimiento B — mediana entre semillas y barras Q1–Q3",
  caption    = "Banda gris: Q1–Q3 del procedimiento A. Línea discontinua: mediana de A"
)

guardar_png(p_B_lineas, file.path(OUT_DIR, "P1_B_lineas.png"), w = 10, h = 6)



cat("\n[Gráfico 2] C vs D — líneas\n")

df_lineas_D <- resumen_lineas(datos, "D", N_VARS_D)

p_D_lineas <- plot_lineas(
  df_summary = df_lineas_D,
  baseline   = C,
  title      = "Comparación entre métodos de prefiltrado",
  subtitle   = "Procedimiento D — mediana entre semillas y barras Q1–Q3",
  caption    = "Banda gris: Q1–Q3 del procedimiento C (submuestra sin filtrado). Línea discontinua: mediana de C"
)

guardar_png(p_D_lineas, file.path(OUT_DIR, "P1_D_lineas.png"), w = 10, h = 6)




cat("\n===== FIGURAS GENERADAS =====\n")
cat("Directorio:", OUT_DIR, "\n")
cat("  - P1_B_lineas.png  (A vs B, líneas)\n")
cat("  - P1_D_lineas.png  (C vs D, líneas)\n")
