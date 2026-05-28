library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

RESULTS_ROOT <- "resultados"
OUT_DIR      <- "figuras"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEMILLAS <- c(12345, 48271, 93054, 17689, 60532)

PREFILTER_LABELS <- c(
  "p_valor"   = "p-valor",
  "mrmr"      = "MRMR",
  "cmim"      = "CMIM",
  "aleatorio" = "Aleatorio"
)

PREFILTER_ORDER <- c("p_valor", "mrmr", "cmim", "aleatorio")


TEC_COLORS <- c(
  "F" = "#FFE082",  
  "G" = "#FFAB91",  
  "H" = "#80CBC4"   
)

TEC_LABELS <- c(
  "F" = "SMOTE (F)",
  "G" = "Bootstrap (G)",
  "H" = "Ruido (H)"
)

NIVELES <- c(110, 120, 130)


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


datos_FGH <- datos %>%
  filter(procedure %in% c("F", "G", "H"),
         preselection %in% PREFILTER_ORDER) %>%
  mutate(
    nivel = case_when(
      procedure == "F" ~ n_s,
      procedure == "G" ~ n_b,
      procedure == "H" ~ n_n
    )
  ) %>%
  filter(nivel %in% NIVELES)

cat("\nFilas F/G/H tras filtrar nivel ∈ {110,120,130}:", nrow(datos_FGH), "\n")
cat("Distribución por procedimiento × nivel:\n")
print(table(datos_FGH$procedure, datos_FGH$nivel))


datos_D <- datos %>%
  filter(procedure == "D",
         preselection %in% PREFILTER_ORDER)

df_seed_D <- datos_D %>%
  group_by(semilla_split, procedure, preselection, n_variables) %>%
  summarise(median_auc = median(auc_test, na.rm = TRUE),
            .groups = "drop")

ref_D <- df_seed_D %>%
  summarise(
    median_val = median(median_auc, na.rm = TRUE),
    q25        = quantile(median_auc, 0.25, na.rm = TRUE),
    q75        = quantile(median_auc, 0.75, na.rm = TRUE)
  )

cat("\nReferencia Procedimiento D:\n")
print(ref_D)


df_seed <- datos_FGH %>%
  group_by(semilla_split, procedure, preselection, n_variables, nivel) %>%
  summarise(median_auc = median(auc_test, na.rm = TRUE),
            .groups = "drop")

df_summary <- df_seed %>%
  group_by(procedure, preselection, n_variables, nivel) %>%
  summarise(
    median_val = median(median_auc, na.rm = TRUE),
    q25        = quantile(median_auc, 0.25, na.rm = TRUE),
    q75        = quantile(median_auc, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )


guardar_png <- function(p, file, w = 9, h = 4.5) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = 150)
  cat("  →", file, "\n")
}

df_bars_global <- df_seed %>%
  group_by(procedure, nivel) %>%
  summarise(
    median_val = median(median_auc, na.rm = TRUE),
    q25        = quantile(median_auc, 0.25, na.rm = TRUE),
    q75        = quantile(median_auc, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

BAR_COLORS <- c(
  "F_110" = "#FFECB3", "F_120" = "#FFD54F", "F_130" = "#FFA000",
  "G_110" = "#FFCCBC", "G_120" = "#FF8A65", "G_130" = "#E64A19",
  "H_110" = "#B2DFDB", "H_120" = "#4DB6AC", "H_130" = "#00796B"
)

plot_barras_FGH_unico <- function(df_bars, ref_D) {
  
  df_p <- df_bars %>%
    mutate(
      proc_letra = procedure,
      procedure  = factor(procedure,
                          levels = c("F", "G", "H"),
                          labels = TEC_LABELS[c("F", "G", "H")]),
      nivel      = factor(nivel, levels = rev(NIVELES)),
      grupo      = paste(proc_letra, as.character(nivel), sep = "_")
    )
  
  xmin_data <- min(df_p$q25, ref_D$q25, na.rm = TRUE)
  xmax_data <- max(df_p$q75, ref_D$q75, na.rm = TRUE)
  xlim_bars <- c(xmin_data - 0.005, xmax_data + 0.012)
  
  ggplot(df_p,
         aes(x = median_val, y = nivel, fill = grupo)) +
    geom_rect(data = data.frame(xmin = ref_D$q25, xmax = ref_D$q75),
              aes(xmin = xmin, xmax = xmax),
              ymin = -Inf, ymax = Inf,
              inherit.aes = FALSE,
              fill = "grey70", alpha = 0.25) +
    geom_vline(xintercept = ref_D$median_val,
               linetype = "dashed", color = "grey30",
               linewidth = 0.5) +
    geom_col(width = 0.75, color = NA) +
    geom_errorbarh(aes(xmin = q25, xmax = q75),
                   height = 0.25, linewidth = 0.6, color = "grey30") +
    geom_text(aes(label = sprintf("%.3f", median_val),
                  x = q75 + 0.003),
              hjust = 0, size = 3.6, color = "grey20") +
    scale_fill_manual(values = BAR_COLORS, guide = "none") +
    scale_y_discrete(labels = function(x) paste0("n=", x)) +
    coord_cartesian(xlim = xlim_bars) +
    facet_grid(rows = vars(procedure), scales = "free_y", switch = "y") +
    labs(
      title    = "Comparación entre técnicas de aumento de datos",
      x        = "AUC test",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 11),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, face = "bold",
                                        size = 12, hjust = 1,
                                        margin = ggplot2::margin(r = 8)),
      panel.spacing.y    = unit(0, "lines"),
      panel.background   = element_rect(fill = "transparent", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      strip.background   = element_rect(fill = "transparent", color = NA),
      plot.caption       = element_text(size = 10, color = "grey30",
                                        hjust = 1)
    )
}

cat("\n[Parte 1] Barras agrupadas F/G/H × nivel (gráfico único)\n")
p_unico <- plot_barras_FGH_unico(df_bars_global, ref_D)
guardar_png(p_unico,
            file.path(OUT_DIR, "P4_barras_FGH_agrupadas.png"),
            w = 9, h = 6)

cat("\n===== FIGURAS GENERADAS =====\n")
cat("Directorio:", OUT_DIR, "\n")
