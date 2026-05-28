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


COLOR_LOW   <- "#C62828"   
COLOR_MID   <- "#FFFFFF"  
COLOR_HIGH  <- "#2E7D32"   


POWER_EXP <- 0.55



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


resumir_heatmap <- function(df) {
  df %>%
    group_by(semilla_split, preselection, n_variables) %>%
    summarise(median_auc = median(auc_test, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(preselection, n_variables) %>%
    summarise(median_val = median(median_auc, na.rm = TRUE),
              .groups = "drop")
}

df_B <- datos %>%
  filter(procedure == "B",
         preselection %in% PREFILTER_ORDER) %>%
  resumir_heatmap()

df_D <- datos %>%
  filter(procedure == "D",
         preselection %in% PREFILTER_ORDER) %>%
  resumir_heatmap()




aplicar_transformacion <- function(df) {
  ancla <- median(df$median_val, na.rm = TRUE)
  vmin  <- min(df$median_val, na.rm = TRUE)
  vmax  <- max(df$median_val, na.rm = TRUE)

  df %>%
    mutate(
      ancla_block = ancla,
      vmin_block  = vmin,
      vmax_block  = vmax,
      raw_t = if_else(
        median_val <= ancla,
        -((ancla - median_val) / (ancla - vmin)),
        ((median_val - ancla) / (vmax - ancla))
      ),
      fill_val = sign(raw_t) * abs(raw_t) ^ POWER_EXP
    )
}

df_B <- aplicar_transformacion(df_B)
df_D <- aplicar_transformacion(df_D)

cat("\nValores en bloque B:\n"); print(df_B %>% select(preselection, n_variables, median_val, fill_val))
cat("\nValores en bloque D:\n"); print(df_D %>% select(preselection, n_variables, median_val, fill_val))



plot_heatmap <- function(df, titulo, subtitulo) {

  ancla <- unique(df$ancla_block)
  vmin  <- unique(df$vmin_block)
  vmax  <- unique(df$vmax_block)

  cat(sprintf("  → ancla: %.4f | min: %.4f | max: %.4f\n",
              ancla, vmin, vmax))

  df <- df %>%
    mutate(
      preselection = factor(preselection, levels = rev(PREFILTER_ORDER)),
      n_variables  = factor(n_variables, levels = sort(unique(n_variables))),
      label        = sprintf("%.3f", median_val)
    )

  legend_aucs <- c(vmin, ancla, vmax)
  legend_fills <- sapply(legend_aucs, function(v) {
    raw_t <- if (v <= ancla) -((ancla - v) / (ancla - vmin)) else ((v - ancla) / (vmax - ancla))
    sign(raw_t) * abs(raw_t) ^ POWER_EXP
  })

  ggplot(df, aes(x = n_variables, y = preselection, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 3.6,
              fontface = "bold", color = "grey15") +
    scale_y_discrete(labels = PREFILTER_LABELS) +
    scale_fill_gradient2(
      low      = COLOR_LOW,
      mid      = COLOR_MID,
      high     = COLOR_HIGH,
      midpoint = 0,
      limits   = c(-1, 1),
      breaks   = legend_fills,
      labels   = sprintf("%.3f", legend_aucs),
      name     = "AUC test",
      guide    = guide_colorbar(
        barheight       = unit(4, "cm"),
        ticks           = TRUE,
        ticks.colour    = "grey40",
        ticks.linewidth = 0.5,
        frame.colour    = NA,
        draw.ulim       = FALSE,
        draw.llim       = FALSE
      )
    ) +
    labs(
      title    = titulo,
      subtitle = subtitulo,
      caption  = sprintf("Escala divergente anclada en la mediana del bloque (%.3f): verde = por encima, rosa = por debajo",
                         ancla),
      x        = "Número de variables retenidas (p)",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      panel.grid      = element_blank(),
      axis.text.y     = element_text(size = 12),
      axis.text.x     = element_text(angle = 0),
      legend.position = "right",
      legend.title    = element_text(margin = ggplot2::margin(b = 10)),
      plot.caption    = element_text(size = 10, color = "grey30")
    )
}


guardar_png <- function(p, file, w = 11, h = 4.5) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = 150)
  cat("  →", file, "\n")
}

cat("\n[Heatmap B]\n")
p_B <- plot_heatmap(
  df_B,
  titulo    = "Comportamiento de cada filtro en función de p",
  subtitulo = "Procedimiento B — train completo (n grande)"
)
guardar_png(p_B, file.path(OUT_DIR, "P2_heatmap_B.png"),
            w = 12, h = 4.5)

cat("\n[Heatmap D]\n")
p_D <- plot_heatmap(
  df_D,
  titulo    = "Comportamiento de cada filtro en función de p",
  subtitulo = "Procedimiento D — submuestra"
)
guardar_png(p_D, file.path(OUT_DIR, "P2_heatmap_D.png"),
            w = 9, h = 4.5)

cat("\n===== FIGURAS GENERADAS =====\n")
cat("Directorio:", OUT_DIR, "\n")
cat("  - P2_heatmap_B.png  (Procedimiento B)\n")
cat("  - P2_heatmap_D.png  (Procedimiento D)\n")
