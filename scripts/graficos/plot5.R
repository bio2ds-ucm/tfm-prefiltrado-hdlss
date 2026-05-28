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
  "H" = "#80CBC4",  
  "D" = "#c4a2fc"   
)

TEC_LABELS <- c(
  "F" = "SMOTE (F)",
  "G" = "Bootstrap (G)",
  "H" = "Ruido (H)",
  "D" = "Sin aumento (D)"
)

NIVELES     <- c(110, 120, 130)
EMPATE_TOL  <- 0.001  


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



datos_D <- datos %>%
  filter(procedure == "D",
         preselection %in% PREFILTER_ORDER)


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

cat("\nFilas D:",   nrow(datos_D), "\n")
cat("Filas F/G/H:", nrow(datos_FGH), "\n")


df_seed_D <- datos_D %>%
  group_by(semilla_split, preselection, n_variables) %>%
  summarise(median_auc = median(auc_test, na.rm = TRUE),
            .groups = "drop")

df_summary_D <- df_seed_D %>%
  group_by(preselection, n_variables) %>%
  summarise(auc_D = median(median_auc, na.rm = TRUE),
            .groups = "drop")


df_seed_FGH <- datos_FGH %>%
  group_by(semilla_split, procedure, preselection, n_variables, nivel) %>%
  summarise(median_auc = median(auc_test, na.rm = TRUE),
            .groups = "drop")

df_summary_FGH <- df_seed_FGH %>%
  group_by(procedure, preselection, n_variables, nivel) %>%
  summarise(auc_tec = median(median_auc, na.rm = TRUE),
            .groups = "drop")


df_delta <- df_summary_FGH %>%
  left_join(df_summary_D, by = c("preselection", "n_variables")) %>%
  mutate(delta = auc_tec - auc_D)

cat("\nFilas df_delta:", nrow(df_delta), "\n")



calcular_ganador_vs_D <- function(df_delta_nivel) {
  df_delta_nivel %>%
    group_by(preselection, n_variables) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    summarise(
      auc_D    = auc_D[1],
      delta_max = delta[1],
      ganadores = list({
        if (delta[1] <= EMPATE_TOL) {
          "D"
        } else {
          sort(procedure[delta >= delta_max - EMPATE_TOL])
        }
      }),
      delta_lbl_val = if_else(delta[1] <= EMPATE_TOL,
                              -delta_max, delta_max),
      .groups = "drop"
    ) %>%
    mutate(
      n_ganadores = lengths(ganadores),
      empate      = n_ganadores > 1,
      gana_D      = vapply(ganadores, function(x) identical(x, "D"),
                           logical(1)),
      etiqueta    = vapply(ganadores, paste, character(1), collapse = "="),
      delta_lbl   = sprintf("+%.3f", delta_lbl_val)
    )
}


construir_poligonos <- function(df_ganador) {
  x_levels <- levels(df_ganador$n_variables)
  y_levels <- levels(df_ganador$preselection)

  half <- 0.5

  poly_list <- lapply(seq_len(nrow(df_ganador)), function(i) {
    fila    <- df_ganador[i, ]
    gans    <- fila$ganadores[[1]]
    ng      <- length(gans)
    cx      <- match(as.character(fila$n_variables),  x_levels)
    cy      <- match(as.character(fila$preselection), y_levels)
    cell_id <- paste0(cx, "_", cy)

    if (ng == 1) {
      data.frame(
        cell_id = cell_id, poly_id = paste0(cell_id, "_1"),
        ganador = gans,
        x = cx + c(-half,  half, half, -half),
        y = cy + c(-half, -half, half,  half)
      )
    } else if (ng == 2) {
      do.call(rbind, list(
        data.frame(
          cell_id = cell_id, poly_id = paste0(cell_id, "_1"),
          ganador = gans[1],
          x = cx + c(-half,  half,  half),
          y = cy + c(-half, -half,  half)
        ),
        data.frame(
          cell_id = cell_id, poly_id = paste0(cell_id, "_2"),
          ganador = gans[2],
          x = cx + c(-half,  half, -half),
          y = cy + c(-half,  half,  half)
        )
      ))
    } else {
      TL <- c(cx - half, cy + half)
      TR <- c(cx + half, cy + half)
      BR <- c(cx + half, cy - half)
      BL <- c(cx - half, cy - half)
      C  <- c(cx, cy)
      do.call(rbind, list(
        data.frame(
          cell_id = cell_id, poly_id = paste0(cell_id, "_1"),
          ganador = gans[1],
          x = c(C[1], TL[1], TR[1]),
          y = c(C[2], TL[2], TR[2])
        ),
        data.frame(
          cell_id = cell_id, poly_id = paste0(cell_id, "_2"),
          ganador = gans[2],
          x = c(C[1], TR[1], BR[1]),
          y = c(C[2], TR[2], BR[2])
        ),
        data.frame(
          cell_id = cell_id, poly_id = paste0(cell_id, "_3"),
          ganador = gans[3],
          x = c(C[1], BR[1], BL[1], TL[1]),
          y = c(C[2], BR[2], BL[2], TL[2])
        )
      ))
    }
  })

  do.call(rbind, poly_list)
}


guardar_png <- function(p, file, w = 9, h = 4.5) {
  ggsave(filename = file, plot = p, width = w, height = h, dpi = 150)
  cat("  →", file, "\n")
}

plot_heatmap_vs_D <- function(df_ganador, nivel_actual) {

  df_ganador <- df_ganador %>%
    mutate(
      preselection = factor(preselection,
                            levels = rev(PREFILTER_ORDER),
                            labels = PREFILTER_LABELS[rev(PREFILTER_ORDER)]),
      n_variables  = factor(n_variables, levels = sort(unique(n_variables)))
    )

  df_poly <- construir_poligonos(df_ganador)

  x_levels <- levels(df_ganador$n_variables)
  y_levels <- levels(df_ganador$preselection)

  df_lbl <- df_ganador %>%
    mutate(
      cx = match(as.character(n_variables),  x_levels),
      cy = match(as.character(preselection), y_levels)
    )


  df_empate <- df_lbl %>% filter(empate)

  if (nrow(df_empate) > 0) {
    df_overlay_poly <- df_empate %>%
      mutate(cell_id = paste0(cx, "_", cy)) %>%
      rowwise() %>%
      do({
        r <- .
        data.frame(
          cell_id = r$cell_id,
          x = r$cx + c(-0.5,  0.5, 0.5, -0.5),
          y = r$cy + c(-0.5, -0.5, 0.5,  0.5)
        )
      }) %>%
      ungroup()
  } else {
    df_overlay_poly <- data.frame(cell_id = character(0),
                                  x = numeric(0),
                                  y = numeric(0))
  }

  ggplot() +
    geom_polygon(data = df_poly,
                 aes(x = x, y = y, group = poly_id, fill = ganador),
                 color = "white", linewidth = 0.8) +
    geom_polygon(data = df_overlay_poly,
                 aes(x = x, y = y, group = cell_id),
                 fill = "grey60", alpha = 0.45,
                 color = "white", linewidth = 0.8) +
    geom_text(data = df_lbl,
              aes(x = cx, y = cy, label = etiqueta),
              size = 4.4, fontface = "bold", color = "grey15",
              vjust = -0.4) +
    geom_text(data = df_lbl,
              aes(x = cx, y = cy, label = delta_lbl),
              size = 3, color = "grey25", vjust = 1.6) +
    scale_fill_manual(values = TEC_COLORS,
                      labels = TEC_LABELS,
                      breaks = c("F", "G", "H", "D"),
                      name   = "Ganador") +
    scale_x_continuous(breaks = seq_along(x_levels), labels = x_levels,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(y_levels), labels = y_levels,
                       expand = c(0, 0)) +
    coord_equal(clip = "off") +
    labs(
      title    = "Mejora respecto al modelo sin aumento (D)",
      x        = "Número de variables retenidas (p)",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      panel.grid       = element_blank(),
      axis.text.y      = element_text(size = 12),
      legend.position  = "right",
      legend.title     = element_text(margin = ggplot2::margin(b = 10)),
      plot.caption     = element_text(size = 10, color = "grey30")
    )
}

cat("\n[Heatmaps] Mejora vs D por nivel\n")
for (n in NIVELES) {
  df_n <- df_delta %>% filter(nivel == n)
  df_g <- calcular_ganador_vs_D(df_n)

  cat(sprintf("\n  Nivel %d:\n", n))
  print(df_g %>%
          count(etiqueta) %>%
          arrange(desc(n)))

  p <- plot_heatmap_vs_D(df_g, n)
  guardar_png(p, file.path(OUT_DIR,
                           sprintf("P5_heatmap_vs_D_n%d.png", n)),
              w = 10, h = 4.5)
}

cat("\n===== FIGURAS GENERADAS =====\n")
cat("Directorio:", OUT_DIR, "\n")
