# =============================================================
# SCRIPT 1 (DATASET B) - Semilla split: 12345
# Dataset: GSE202203_rnaseq_data_orig_obj_tumorSize
# =============================================================
# .libPaths("RUTA/A/TUS/Rlibs")  # opcional: solo si usas una libreria de R personalizada
# ---- Librerías ----
library(caret)
library(pROC)
library(dplyr)
library(glmnet)
library(doParallel)
library(foreach)
library(mRMRe)
library(tidyr)
library(praznik)
library(smotefamily)

SEMILLA_SCRIPT <- 12345
RESULTS_DIR    <- "resultados/semilla_12345"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

n_cores <- as.integer(Sys.getenv("N_CORES", detectCores() - 1))
cluster <- parallel::makePSOCKcluster(n_cores)
doParallel::registerDoParallel(cluster)
# Forzar inicializacion de .doSnowGlobals en workers
parallel::clusterCall(cluster, function() {
  # .libPaths("RUTA/A/TUS/Rlibs")  # opcional: solo si usas una libreria de R personalizada
  library(doParallel)
  assign(".doSnowGlobals", new.env(), envir = globalenv())
})

# =============================================================
# FUNCIONES
# =============================================================

combinar_resultados <- function(lista, procedimiento) {
  lista <- Filter(function(x) is.list(x) && !is.null(x$summary), lista)
  if (length(lista) == 0) {
    return(list(summary = NULL, coefs = NULL, vars_input = NULL))
  }
  sum_df  <- do.call(dplyr::bind_rows, lapply(lista, `[[`, "summary"))
  coef_df <- do.call(dplyr::bind_rows, lapply(lista, `[[`, "coefs"))
  var_df  <- do.call(dplyr::bind_rows, lapply(lista, `[[`, "vars_input"))
  if (!is.null(sum_df))  sum_df$procedimiento  <- procedimiento
  if (!is.null(coef_df)) coef_df$procedimiento <- procedimiento
  if (!is.null(var_df))  var_df$procedimiento  <- procedimiento
  list(summary = sum_df, coefs = coef_df, vars_input = var_df)
}

seleccionar_variables <- function(data, semilla = 12345) {
  variable_fija <- "obj"
  if (!variable_fija %in% colnames(data)) stop("La variable fija no existe en el dataset.")
  if (!is.null(semilla)) set.seed(semilla)
  total_vars <- ncol(data)
  num_seleccion_total <- floor(0.80 * total_vars)
  variables_restantes <- setdiff(colnames(data), variable_fija)
  num_aleatorias <- num_seleccion_total - 1
  if (num_aleatorias < 0) stop("El dataset tiene muy pocas variables.")
  variables_aleatorias <- sample(variables_restantes, num_aleatorias)
  nuevo_df <- data[, c(variable_fija, variables_aleatorias)]
  return(nuevo_df)
}

generar_rankings_prefiltros <- function(data, target, top_max = 5000, min_obs = 5,
                                        n_bins_cmim = 3,
                                        threads_cmim = parallel::detectCores() - 1) {
  if(!is.data.frame(data)) stop("`data` debe ser un data.frame")
  if(!target %in% names(data)) stop("Target no encontrada")
  
  df <- data
  df$sample_id <- NULL
  df[[target]] <- as.factor(df[[target]])
  predictoras <- setdiff(names(df), target)
  
  p_values <- sapply(predictoras, function(var){
    x <- df[[var]]; y <- df[[target]]
    if(!is.numeric(x)) return(NA_real_)
    tmp <- na.omit(data.frame(x = x, y = y))
    if(nrow(tmp) < min_obs || length(unique(tmp$y)) < 2 || length(unique(tmp$x)) < 2) return(NA_real_)
    test <- try(t.test(tmp$x ~ tmp$y), silent = TRUE)
    if(inherits(test, "try-error")) return(NA_real_)
    test$p.value
  })
  rank_p <- data.frame(variable = names(p_values), p_value = p_values)
  rank_p <- rank_p[!is.na(rank_p$p_value),]
  rank_p <- rank_p[order(rank_p$p_value),]
  rownames(rank_p) <- NULL
  if(nrow(rank_p) > top_max) rank_p <- rank_p[seq_len(top_max),]
  
  df_m <- df
  pred_num <- predictoras[sapply(df_m[predictoras], is.numeric)]
  df_m <- df_m[, c(pred_num, target), drop = FALSE]
  df_m[[target]] <- as.numeric(as.factor(df_m[[target]]))
  rank_m <- data.frame(variable = character(0), mrmr_rank = integer(0))
  if(length(pred_num) >= 1) {
    data_mrmr <- mRMRe::mRMR.data(data = df_m)
    target_index <- which(colnames(df_m) == target)
    modelo <- mRMRe::mRMR.classic(data = data_mrmr, target_indices = target_index,
                                  feature_count = min(top_max, length(pred_num)))
    idx <- mRMRe::solutions(modelo)[[1]]
    vars_m <- colnames(df_m)[idx]
    vars_m <- vars_m[vars_m != target]
    rank_m <- data.frame(variable = vars_m, mrmr_rank = seq_along(vars_m))
    rownames(rank_m) <- NULL
    if(nrow(rank_m) > top_max) rank_m <- rank_m[seq_len(top_max),]
  }
  
  pred_num_cmim <- predictoras[sapply(df[predictoras], is.numeric)]
  rank_c <- data.frame(variable = character(0), cmim_rank = integer(0))
  if(length(pred_num_cmim) >= 1) {
    x_cmim <- df[, pred_num_cmim, drop = FALSE]
    y_cmim <- df[[target]]
    completos <- stats::complete.cases(x_cmim, y_cmim)
    x_cmim <- x_cmim[completos,]; y_cmim <- y_cmim[completos]
    vars_validas <- vapply(x_cmim, function(z) length(unique(z)) > 1, logical(1))
    x_cmim <- x_cmim[, vars_validas, drop = FALSE]
    if(ncol(x_cmim) >= 1 && length(unique(y_cmim)) >= 2) {
      x_train_disc <- as.data.frame(apply(x_cmim, 2, function(col)
        infotheo::discretize(col, nbins = n_bins_cmim)[, 1]))
      y_cmim <- droplevels(as.factor(y_cmim))
      k_cmim <- min(top_max, ncol(x_train_disc))
      sel_cmim <- praznik::CMIM(x_train_disc, y_cmim, k = k_cmim,
                                threads = threads_cmim)$selection
      vars_c <- colnames(x_train_disc)[sel_cmim]
      rank_c <- data.frame(variable = vars_c, cmim_rank = seq_along(vars_c))
      rownames(rank_c) <- NULL
    }
  }
  
  return(list(p_valor = rank_p, mrmr = rank_m, cmim = rank_c,
              top_max = top_max, target = target))
}

filtro_pvalor <- function(data, target, rankings, n_vars) {
  r <- rankings$p_valor
  vars_ok <- intersect(r$variable, setdiff(names(data), target))
  n <- min(as.integer(n_vars), length(vars_ok))
  data[, c(vars_ok[seq_len(n)], target), drop = FALSE]
}

filtro_mrmr <- function(data, target, rankings, n_vars) {
  r <- rankings$mrmr
  vars_ok <- intersect(r$variable, setdiff(names(data), target))
  n <- min(as.integer(n_vars), length(vars_ok))
  data[, c(vars_ok[seq_len(n)], target), drop = FALSE]
}

filtro_cmim <- function(data, target, rankings, n_vars) {
  r <- rankings$cmim
  vars_ok <- intersect(r$variable, setdiff(names(data), target))
  n <- min(as.integer(n_vars), length(vars_ok))
  data[, c(vars_ok[seq_len(n)], target), drop = FALSE]
}

filtro_aleatorio <- function(data, target, n_vars, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  vars_disponibles <- setdiff(names(data), target)
  n <- min(as.integer(n_vars), length(vars_disponibles))
  vars_sel <- sample(vars_disponibles, size = n, replace = FALSE)
  data[, c(vars_sel, target), drop = FALSE]
}

muestreo <- function(data, group_var, n = 100, seed = 12345) {
  set.seed(seed)
  y_all <- data[[group_var]]
  tab <- table(y_all); prop <- tab / sum(tab)
  n_str <- floor(n * prop)
  faltan <- n - sum(n_str)
  if (faltan > 0) {
    restos <- (n * prop) - n_str
    add_levels <- names(sort(restos, decreasing = TRUE))[seq_len(faltan)]
    n_str[add_levels] <- n_str[add_levels] + 1
  }
  idx_sample <- unlist(lapply(names(n_str), function(level) {
    rows <- which(y_all == level)
    sample(rows, size = n_str[[level]], replace = FALSE)
  }))
  list(muestra = data[idx_sample,, drop = FALSE],
       n_por_clase_muestra = table(data[idx_sample, group_var]))
}

elastic_net_single <- function(train_data, test_data, target, execution_id,
                               alpha = 1, nfolds = 5, seed = 1234) {
  set.seed(seed)
  y_train <- as.factor(train_data[[target]])
  x_train <- model.matrix(~ . - 1,
                          train_data[, setdiff(colnames(train_data), target), drop = FALSE])
  cv_model <- glmnet::cv.glmnet(x = x_train, y = y_train, alpha = alpha,
                                family = "binomial", nfolds = nfolds,
                                type.measure = "auc", standardize = TRUE)
  x_test_tmp <- model.matrix(~ . - 1,
                             test_data[, setdiff(colnames(test_data), target), drop = FALSE])
  x_test <- x_test_tmp[, colnames(x_train), drop = FALSE]
  pred_test <- predict(cv_model, newx = x_test, s = "lambda.min", type = "response")
  roc_obj   <- pROC::roc(test_data[[target]], as.numeric(pred_test))
  auc_test  <- as.numeric(roc_obj$auc)
  
  coef_min <- as.matrix(coef(cv_model, s = "lambda.min"))
  coef_1se <- as.matrix(coef(cv_model, s = "lambda.1se"))
  nnz_min  <- sum(coef_min[-which(rownames(coef_min) == "(Intercept)"),] != 0)
  nnz_1se  <- sum(coef_1se[-which(rownames(coef_1se) == "(Intercept)"),] != 0)
  
  # --- Tabla de coeficientes no cero (incluye intercepto) ---
  coefs_df <- rbind(
    data.frame(execution_id = execution_id,
               s            = "lambda.min",
               variable     = rownames(coef_min),
               coef         = as.numeric(coef_min),
               stringsAsFactors = FALSE),
    data.frame(execution_id = execution_id,
               s            = "lambda.1se",
               variable     = rownames(coef_1se),
               coef         = as.numeric(coef_1se),
               stringsAsFactors = FALSE)
  )
  coefs_df <- coefs_df[coefs_df$coef != 0, , drop = FALSE]
  rownames(coefs_df) <- NULL
  
  # --- Variables de entrada al modelo ---
  vars_input_df <- data.frame(
    execution_id = execution_id,
    variable     = setdiff(colnames(train_data), target),
    stringsAsFactors = FALSE
  )
  
  summary_df <- data.frame(
    execution_id  = execution_id,
    alpha         = alpha,
    lambda_min    = cv_model$lambda.min,
    lambda_1se    = cv_model$lambda.1se,
    auc_test      = auc_test,
    n_nonzero_min = nnz_min,
    n_nonzero_1se = nnz_1se,
    stringsAsFactors = FALSE
  )
  
  list(summary = summary_df, coefs = coefs_df, vars_input = vars_input_df)
}

bootstrap_sample_stratified <- function(data, group_var, n, seed = NULL) {
  n  <- as.integer(n)
  n0 <- nrow(data)
  if (!is.null(seed)) set.seed(seed)
  
  # Sin aumento: devolver la muestra tal cual (equivalente a D)
  if (n <= n0) return(data)
  
  n_extra <- n - n0
  g <- as.factor(data[[group_var]])
  props <- table(g) / length(g)
  
  n_per_class <- floor(props * n_extra)
  remainder <- n_extra - sum(n_per_class)
  if (remainder > 0) {
    frac <- (props * n_extra) - n_per_class
    add_classes <- names(sort(frac, decreasing = TRUE))[seq_len(remainder)]
    n_per_class[add_classes] <- n_per_class[add_classes] + 1
  }
  
  idx_extra <- integer(0)
  for (cls in names(n_per_class)) {
    cls_idx <- which(g == cls)
    k <- as.integer(n_per_class[[cls]])
    if (k > 0) idx_extra <- c(idx_extra, sample(cls_idx, size = k, replace = TRUE))
  }
  
  out <- rbind(data, data[idx_extra,, drop = FALSE])
  rownames(out) <- NULL
  out
}

add_noise_sample_stratified <- function(data, group_var, n, noise_factor = 0.01, seed = NULL) {
  n <- as.integer(n)
  if (!is.null(seed)) set.seed(seed)
  g <- as.factor(data[[group_var]])
  n_originales <- nrow(data)
  n_nuevas <- n - n_originales
  if (n_nuevas <= 0) return(data)
  props <- table(g) / length(g)
  n_per_class <- floor(props * n_nuevas)
  remainder <- n_nuevas - sum(n_per_class)
  if (remainder > 0) {
    frac <- (props * n_nuevas) - n_per_class
    add_classes <- names(sort(frac, decreasing = TRUE))[seq_len(remainder)]
    n_per_class[add_classes] <- n_per_class[add_classes] + 1
  }
  num_cols <- setdiff(names(data), group_var)
  num_cols <- num_cols[vapply(data[num_cols], is.numeric, logical(1))]
  sd_cols <- vapply(data[num_cols], function(x) { s <- sd(x, na.rm=TRUE); if(is.na(s)||s==0) 0 else s }, numeric(1))
  sinteticas_list <- lapply(names(n_per_class), function(cls) {
    k <- as.integer(n_per_class[[cls]])
    if (k == 0) return(NULL)
    cls_idx <- which(g == cls)
    base_rows <- data[sample(cls_idx, size = k, replace = TRUE),, drop = FALSE]
    for (col in num_cols) {
      sd_ruido <- sd_cols[col] * noise_factor
      if (sd_ruido > 0) base_rows[[col]] <- base_rows[[col]] + rnorm(k, 0, sd_ruido)
    }
    base_rows
  })
  sinteticas <- do.call(rbind, Filter(Negate(is.null), sinteticas_list))
  rownames(sinteticas) <- NULL
  out <- rbind(data, sinteticas)
  rownames(out) <- NULL
  out
}

smote <- function(data, group_var, n, seed = 12345, drop_cols = c("sample_id")) {
  set.seed(seed)
  df <- as.data.frame(data)

  # Sin aumento: devolver la muestra tal cual (equivalente a D)
  if (n == nrow(df)) return(df)

  df[intersect(drop_cols, names(df))] <- NULL
  df[[group_var]] <- as.factor(df[[group_var]])
  X <- df[, setdiff(names(df), group_var), drop = FALSE]
  X <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
  dat <- cbind(X, df[group_var])
  names(dat)[ncol(dat)] <- group_var
  n0 <- nrow(dat)
  
  # Caso submuestreo: igual que antes
  if (n <= n0) {
    out <- dat[sample.int(n0, n), , drop = FALSE]
    out$.smote_source <- "Original"
    rownames(out) <- NULL
    return(out)
  }
  
  # Reparto de sintéticas por clase respetando proporciones (idéntico a la versión UBL)
  props      <- table(dat[[group_var]]) / n0
  n_nuevas   <- n - n0
  n_sint_cls <- floor(props * n_nuevas)
  remainder  <- n_nuevas - sum(n_sint_cls)
  if (remainder > 0) {
    frac <- (props * n_nuevas) - n_sint_cls
    add_classes <- names(sort(frac, decreasing = TRUE))[seq_len(remainder)]
    n_sint_cls[add_classes] <- n_sint_cls[add_classes] + 1
  }
  
  tab       <- table(dat[[group_var]])
  feat_cols <- setdiff(names(dat), group_var)
  
  # Bucle por clase: smotefamily::SMOTE es binario, así que para cada clase
  # construimos y_bin = "yes"/"no" (yes = clase objetivo) y pedimos el dup_size
  # necesario para generar n_target sintéticas de esa clase.
  synthetic_list <- list()
  
  for (cls in names(tab)) {
    n_target <- as.integer(n_sint_cls[[cls]])
    if (is.na(n_target) || n_target <= 0) next
    
    n_pos <- as.integer(tab[[cls]])       # muestras de la clase objetivo
    n_cls <- n0 - n_pos                    # muestras del resto (clase mayoritaria relativa)
    
    # SMOTE necesita al menos K+1 vecinos de la clase positiva
    k <- min(5L, max(1L, n_pos - 1L))
    if (n_pos < 2L) next  # no se puede interpolar con <2 puntos
    
    # dup_size: cuántas veces queremos duplicar la clase positiva.
    # Queremos generar ~ n_target sintéticas => dup_est = ceiling(n_target / n_pos).
    dup_est <- max(1L, as.integer(ceiling(n_target / n_pos)))
    
    y_bin <- ifelse(as.character(dat[[group_var]]) == cls, "yes", "no")
    
    sm <- smotefamily::SMOTE(
      X        = dat[, feat_cols, drop = FALSE],
      target   = y_bin,
      K        = k,
      dup_size = dup_est
    )
    
    # smotefamily devuelve $syn_data con las filas sintéticas y la columna "class"
    syn <- sm$syn_data
    if (is.null(syn) || nrow(syn) == 0) next
    
    # Quitar la columna "class" que añade smotefamily y poner la clase real
    syn[["class"]] <- NULL
    syn[[group_var]] <- factor(cls, levels = levels(dat[[group_var]]))
    
    # Reordenar columnas para que coincidan con dat
    syn <- syn[, names(dat), drop = FALSE]
    
    # Asegurar tipos numéricos en las features (smotefamily a veces devuelve character)
    for (cc in feat_cols) syn[[cc]] <- as.numeric(syn[[cc]])
    
    # Ajustar al número exacto pedido para esta clase
    if (nrow(syn) > n_target) {
      syn <- syn[sample.int(nrow(syn), n_target), , drop = FALSE]
    } else if (nrow(syn) < n_target) {
      syn <- syn[sample.int(nrow(syn), n_target, replace = TRUE), , drop = FALSE]
    }
    
    synthetic_list[[cls]] <- syn
  }
  
  synthetic <- if (length(synthetic_list)) do.call(rbind, synthetic_list) else dat[0, , drop = FALSE]
  
  # Fallback idéntico al original por si no se pudo sintetizar nada
  needed <- n - n0
  if (nrow(synthetic) == 0) {
    synthetic <- dat[sample.int(n0, needed, replace = TRUE), , drop = FALSE]
  } else if (nrow(synthetic) < needed) {
    extra <- dat[sample.int(n0, needed - nrow(synthetic), replace = TRUE), , drop = FALSE]
    synthetic <- rbind(synthetic, extra)
  } else if (nrow(synthetic) > needed) {
    synthetic <- synthetic[sample.int(nrow(synthetic), needed), , drop = FALSE]
  }
  
  dat$.smote_source       <- "Original"
  synthetic$.smote_source <- "SMOTE"
  out <- rbind(dat, synthetic)
  rownames(out) <- NULL
  out
}

# =============================================================
# DATOS
# =============================================================

datos <- readRDS("data/GSE202203_rnaseq_data_orig_obj_tumorSize.rds")
datos$obj <- as.factor(make.names(datos$obj))
datos_nuevos <- seleccionar_variables(datos, semilla = SEMILLA_SCRIPT)

set.seed(SEMILLA_SCRIPT)
trainIndex <- caret::createDataPartition(datos_nuevos$obj, p = 0.90, list = FALSE)
data_train <- datos_nuevos[trainIndex,]
data_test  <- datos_nuevos[-trainIndex,]

SEMILLAS <- c(12345, 48271, 93054, 17689, 60532,
              21983, 75410, 34867, 82159, 46073,
              59314, 28746, 91025, 63481, 37952,
              14608, 70293, 55871, 42136, 86749)

elastic_net_results <- NULL
coefs_results       <- NULL
vars_input_results  <- NULL

# =============================================================
# PROCEDIMIENTO A
# =============================================================

cat("\n===== PROCEDIMIENTO A =====\n")
res_a <- elastic_net_single(
  train_data   = data_train,
  test_data    = data_test,
  target       = "obj",
  execution_id = "Procedimiento A",
  alpha        = 1,
  nfolds       = 5
)

res_a$summary$procedimiento    <- "A"
res_a$coefs$procedimiento      <- "A"
res_a$vars_input$procedimiento <- "A"

elastic_net_results <- res_a$summary
coefs_results       <- res_a$coefs
vars_input_results  <- res_a$vars_input

cat("Procedimiento A completado. AUC:", res_a$summary$auc_test, "\n")

# =============================================================
# PROCEDIMIENTO B
# =============================================================

cat("\n===== PROCEDIMIENTO B =====\n")

rankings_b <- generar_rankings_prefiltros(data_train, target = "obj",
                                          top_max = 5000, n_bins_cmim = 3,
                                          threads_cmim = 1)
saveRDS(rankings_b, file.path(RESULTS_DIR, "rankings_B.rds"))

n_vars_b  <- c(50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000)
filtros_b <- c("p_valor", "mrmr", "cmim", "aleatorio")

combinaciones_b <- expand.grid(
  ejec   = seq_along(SEMILLAS),
  n_vars = n_vars_b,
  filtro = filtros_b,
  stringsAsFactors = FALSE
)

res_b_list <- foreach(
  row       = seq_len(nrow(combinaciones_b)),
  .packages = c("glmnet", "pROC", "mRMRe", "praznik", "infotheo"),
  .export   = c("combinaciones_b", "SEMILLAS", "data_train", "data_test",
                "rankings_b", "filtro_pvalor", "filtro_mrmr", "filtro_cmim",
                "filtro_aleatorio", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  i      <- combinaciones_b$ejec[row]
  n_vars <- combinaciones_b$n_vars[row]
  filtro <- combinaciones_b$filtro[row]
  seed_i <- SEMILLAS[i]
  
  datos_filtrados <- tryCatch({
    if (filtro == "aleatorio")    filtro_aleatorio(data_train, "obj", n_vars, seed = seed_i)
    else if (filtro == "p_valor") filtro_pvalor  (data_train, "obj", rankings_b, n_vars)
    else if (filtro == "mrmr")    filtro_mrmr    (data_train, "obj", rankings_b, n_vars)
    else                          filtro_cmim    (data_train, "obj", rankings_b, n_vars)
  }, error = function(e) NULL)
  if (is.null(datos_filtrados)) return(NULL)
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = datos_filtrados,
      test_data    = data_test,
      target       = "obj",
      execution_id = sprintf("Procedimiento B %d variables %s Ejecución %d",
                             n_vars, filtro, i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$filtro    <- filtro
  out$summary$n_vars    <- n_vars
  out$summary$ejec      <- i
  out$coefs$filtro      <- filtro
  out$coefs$n_vars      <- n_vars
  out$coefs$ejec        <- i
  out$vars_input$filtro <- filtro
  out$vars_input$n_vars <- n_vars
  out$vars_input$ejec   <- i
  out
}

res_b <- combinar_resultados(res_b_list, "B")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_b$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_b$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_b$vars_input)
cat("Procedimiento B completado.\n")

# =============================================================
# PROCEDIMIENTO C
# =============================================================

cat("\n===== PROCEDIMIENTO C =====\n")

res_c_list <- foreach(
  i         = seq_along(SEMILLAS),
  .packages = c("glmnet", "pROC"),
  .export   = c("SEMILLAS", "data_train", "data_test",
                "muestreo", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  seed_i  <- SEMILLAS[i]
  muestra <- muestreo(data_train, "obj", n = 100, seed = seed_i)$muestra
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = muestra,
      test_data    = data_test,
      target       = "obj",
      execution_id = paste("Procedimiento C Ejecución", i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$ejec    <- i
  out$coefs$ejec      <- i
  out$vars_input$ejec <- i
  out
}

muestras_c <- lapply(seq_along(SEMILLAS), function(i)
  muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra)
names(muestras_c) <- paste0("semilla_", SEMILLAS)
saveRDS(muestras_c, file.path(RESULTS_DIR, "muestras_C.rds"))

res_c <- combinar_resultados(res_c_list, "C")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_c$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_c$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_c$vars_input)
cat("Procedimiento C completado.\n")

# =============================================================
# PROCEDIMIENTO D
# =============================================================

cat("\n===== PROCEDIMIENTO D =====\n")

n_vars_d  <- c(50, 100, 200, 300, 500, 1000)
filtros_d <- c("p_valor", "mrmr", "cmim", "aleatorio")

# Rankings por ejecución (secuencial, no se puede paralelizar)
rankings_d <- lapply(seq_along(SEMILLAS), function(i) {
  cat(sprintf("  Rankings D ejecución %d / %d\n", i, length(SEMILLAS)))
  muestra <- muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra
  generar_rankings_prefiltros(muestra, target = "obj", top_max = 1000,
                              n_bins_cmim = 3, threads_cmim = 1)
})
muestras_d <- lapply(seq_along(SEMILLAS), function(i)
  muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra)
names(rankings_d) <- paste0("semilla_", SEMILLAS)
names(muestras_d) <- paste0("semilla_", SEMILLAS)
saveRDS(rankings_d, file.path(RESULTS_DIR, "rankings_D.rds"))
saveRDS(muestras_d, file.path(RESULTS_DIR, "muestras_D.rds"))

combinaciones_d <- expand.grid(
  ejec   = seq_along(SEMILLAS),
  n_vars = n_vars_d,
  filtro = filtros_d,
  stringsAsFactors = FALSE
)

res_d_list <- foreach(
  row       = seq_len(nrow(combinaciones_d)),
  .packages = c("glmnet", "pROC"),
  .export   = c("combinaciones_d", "SEMILLAS", "data_test", "muestras_d",
                "rankings_d", "filtro_pvalor", "filtro_mrmr", "filtro_cmim",
                "filtro_aleatorio", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  i      <- combinaciones_d$ejec[row]
  n_vars <- combinaciones_d$n_vars[row]
  filtro <- combinaciones_d$filtro[row]
  seed_i <- SEMILLAS[i]
  muestra  <- muestras_d[[i]]
  rankings <- rankings_d[[i]]
  
  datos_filtrados <- tryCatch({
    if (filtro == "aleatorio")    filtro_aleatorio(muestra, "obj", n_vars, seed = seed_i)
    else if (filtro == "p_valor") filtro_pvalor  (muestra, "obj", rankings, n_vars)
    else if (filtro == "mrmr")    filtro_mrmr    (muestra, "obj", rankings, n_vars)
    else                          filtro_cmim    (muestra, "obj", rankings, n_vars)
  }, error = function(e) NULL)
  if (is.null(datos_filtrados)) return(NULL)
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = datos_filtrados,
      test_data    = data_test,
      target       = "obj",
      execution_id = sprintf("Procedimiento D %d variables %s Ejecución %d",
                             n_vars, filtro, i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$filtro    <- filtro
  out$summary$n_vars    <- n_vars
  out$summary$ejec      <- i
  out$coefs$filtro      <- filtro
  out$coefs$n_vars      <- n_vars
  out$coefs$ejec        <- i
  out$vars_input$filtro <- filtro
  out$vars_input$n_vars <- n_vars
  out$vars_input$ejec   <- i
  out
}

res_d <- combinar_resultados(res_d_list, "D")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_d$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_d$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_d$vars_input)
cat("Procedimiento D completado.\n")

# =============================================================
# PROCEDIMIENTO E
# =============================================================

cat("\n===== PROCEDIMIENTO E =====\n")

n_vars_e   <- c(50, 100, 200, 300, 500, 1000)
aumentos_e <- c(0, 0.10, 0.20, 0.30)
n_obs_e    <- round(100 * (1 + aumentos_e))
names(n_obs_e) <- paste0(aumentos_e * 100, "%")

rankings_e <- lapply(seq_along(SEMILLAS), function(i) {
  cat(sprintf("  Rankings E ejecución %d / %d\n", i, length(SEMILLAS)))
  muestra <- muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra
  generar_rankings_prefiltros(muestra, target = "obj", top_max = 1000,
                              n_bins_cmim = 3, threads_cmim = 1)
})
muestras_e <- lapply(seq_along(SEMILLAS), function(i)
  muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra)
names(rankings_e) <- paste0("semilla_", SEMILLAS)
names(muestras_e) <- paste0("semilla_", SEMILLAS)
saveRDS(rankings_e, file.path(RESULTS_DIR, "rankings_E.rds"))
saveRDS(muestras_e, file.path(RESULTS_DIR, "muestras_E.rds"))

combinaciones_e <- expand.grid(
  ejec   = seq_along(SEMILLAS),
  n_vars = n_vars_e,
  n_obs  = n_obs_e,
  filtro = c("p_valor", "mrmr", "cmim", "aleatorio"),
  stringsAsFactors = FALSE
)
combinaciones_e <- combinaciones_e  # sin restriccion n_obs > n_vars

res_e_list <- foreach(
  row       = seq_len(nrow(combinaciones_e)),
  .packages = c("glmnet", "pROC"),
  .export   = c("combinaciones_e", "SEMILLAS", "data_train", "data_test",
                "rankings_e", "muestreo", "filtro_pvalor", "filtro_mrmr",
                "filtro_cmim", "filtro_aleatorio", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  i      <- combinaciones_e$ejec[row]
  n_vars <- combinaciones_e$n_vars[row]
  n_obs  <- combinaciones_e$n_obs[row]
  filtro <- combinaciones_e$filtro[row]
  seed_i <- SEMILLAS[i]
  rankings <- rankings_e[[i]]
  
  muestra_aug <- tryCatch(
    muestreo(data_train, "obj", n = n_obs, seed = seed_i)$muestra,
    error = function(e) NULL)
  if (is.null(muestra_aug)) return(NULL)
  
  datos_filtrados <- tryCatch({
    if (filtro == "aleatorio")    filtro_aleatorio(muestra_aug, "obj", n_vars, seed = seed_i)
    else if (filtro == "p_valor") filtro_pvalor  (muestra_aug, "obj", rankings, n_vars)
    else if (filtro == "mrmr")    filtro_mrmr    (muestra_aug, "obj", rankings, n_vars)
    else                          filtro_cmim    (muestra_aug, "obj", rankings, n_vars)
  }, error = function(e) NULL)
  if (is.null(datos_filtrados)) return(NULL)
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = datos_filtrados,
      test_data    = data_test,
      target       = "obj",
      execution_id = sprintf("Procedimiento E %d variables %s %d obs Ejecución %d",
                             n_vars, filtro, n_obs, i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$filtro    <- filtro
  out$summary$n_vars    <- n_vars
  out$summary$n_obs     <- n_obs
  out$summary$ejec      <- i
  out$coefs$filtro      <- filtro
  out$coefs$n_vars      <- n_vars
  out$coefs$n_obs       <- n_obs
  out$coefs$ejec        <- i
  out$vars_input$filtro <- filtro
  out$vars_input$n_vars <- n_vars
  out$vars_input$n_obs  <- n_obs
  out$vars_input$ejec   <- i
  out
}

res_e <- combinar_resultados(res_e_list, "E")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_e$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_e$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_e$vars_input)
cat("Procedimiento E completado.\n")

# =============================================================
# PROCEDIMIENTO F
# =============================================================

cat("\n===== PROCEDIMIENTO F =====\n")

aumentos_f <- c(0, 0.10, 0.20, 0.30)
n_smote_f  <- round(100 * (1 + aumentos_f))
names(n_smote_f) <- paste0(aumentos_f * 100, "%")

rankings_f <- lapply(seq_along(SEMILLAS), function(i) {
  cat(sprintf("  Rankings F ejecución %d / %d\n", i, length(SEMILLAS)))
  muestra <- muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra
  generar_rankings_prefiltros(muestra, target = "obj", top_max = 1000,
                              n_bins_cmim = 3, threads_cmim = 1)
})
muestras_f <- lapply(seq_along(SEMILLAS), function(i)
  muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra)
names(rankings_f) <- paste0("semilla_", SEMILLAS)
names(muestras_f) <- paste0("semilla_", SEMILLAS)
saveRDS(rankings_f, file.path(RESULTS_DIR, "rankings_F.rds"))
saveRDS(muestras_f, file.path(RESULTS_DIR, "muestras_F.rds"))

combinaciones_f <- expand.grid(
  ejec   = seq_along(SEMILLAS),
  n_vars = c(50, 100, 200, 300, 500, 1000),
  n_s    = n_smote_f,
  filtro = c("p_valor", "mrmr", "cmim", "aleatorio"),
  stringsAsFactors = FALSE
)

res_f_list <- foreach(
  row       = seq_len(nrow(combinaciones_f)),
  .packages = c("glmnet", "pROC", "smotefamily"),
  .export   = c("combinaciones_f", "SEMILLAS", "data_test", "muestras_f",
                "rankings_f", "smote", "filtro_pvalor", "filtro_mrmr",
                "filtro_cmim", "filtro_aleatorio", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  i      <- combinaciones_f$ejec[row]
  n_vars <- combinaciones_f$n_vars[row]
  n_s    <- combinaciones_f$n_s[row]
  filtro <- combinaciones_f$filtro[row]
  seed_i <- SEMILLAS[i]
  muestra  <- muestras_f[[i]]
  rankings <- rankings_f[[i]]
  
  smote_dat <- tryCatch(
    smote(muestra, "obj", n = n_s, seed = seed_i),
    error = function(e) NULL)
  if (is.null(smote_dat)) return(NULL)
  
  datos_filtrados <- tryCatch({
    if (filtro == "aleatorio")    filtro_aleatorio(smote_dat, "obj", n_vars, seed = seed_i)
    else if (filtro == "p_valor") filtro_pvalor  (smote_dat, "obj", rankings, n_vars)
    else if (filtro == "mrmr")    filtro_mrmr    (smote_dat, "obj", rankings, n_vars)
    else                          filtro_cmim    (smote_dat, "obj", rankings, n_vars)
  }, error = function(e) NULL)
  if (is.null(datos_filtrados)) return(NULL)
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = datos_filtrados,
      test_data    = data_test,
      target       = "obj",
      execution_id = sprintf("Procedimiento F %d variables %s %d obs Ejecución %d",
                             n_vars, filtro, n_s, i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$filtro    <- filtro
  out$summary$n_vars    <- n_vars
  out$summary$n_s       <- n_s
  out$summary$ejec      <- i
  out$coefs$filtro      <- filtro
  out$coefs$n_vars      <- n_vars
  out$coefs$n_s         <- n_s
  out$coefs$ejec        <- i
  out$vars_input$filtro <- filtro
  out$vars_input$n_vars <- n_vars
  out$vars_input$n_s    <- n_s
  out$vars_input$ejec   <- i
  out
}

res_f <- combinar_resultados(res_f_list, "F")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_f$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_f$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_f$vars_input)
cat("Procedimiento F completado.\n")

# =============================================================
# PROCEDIMIENTO G
# =============================================================

cat("\n===== PROCEDIMIENTO G =====\n")

aumentos_g <- c(0, 0.10, 0.20, 0.30)
n_boot_g   <- round(100 * (1 + aumentos_g))
names(n_boot_g) <- paste0(aumentos_g * 100, "%")

rankings_g <- lapply(seq_along(SEMILLAS), function(i) {
  cat(sprintf("  Rankings G ejecución %d / %d\n", i, length(SEMILLAS)))
  muestra <- muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra
  generar_rankings_prefiltros(muestra, target = "obj", top_max = 1000,
                              n_bins_cmim = 3, threads_cmim = 1)
})
muestras_g <- lapply(seq_along(SEMILLAS), function(i)
  muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra)
names(rankings_g) <- paste0("semilla_", SEMILLAS)
names(muestras_g) <- paste0("semilla_", SEMILLAS)
saveRDS(rankings_g, file.path(RESULTS_DIR, "rankings_G.rds"))
saveRDS(muestras_g, file.path(RESULTS_DIR, "muestras_G.rds"))

combinaciones_g <- expand.grid(
  ejec   = seq_along(SEMILLAS),
  n_vars = c(50, 100, 200, 300, 500, 1000),
  n_b    = n_boot_g,
  filtro = c("p_valor", "mrmr", "cmim", "aleatorio"),
  stringsAsFactors = FALSE
)

res_g_list <- foreach(
  row       = seq_len(nrow(combinaciones_g)),
  .packages = c("glmnet", "pROC"),
  .export   = c("combinaciones_g", "SEMILLAS", "data_test", "muestras_g",
                "rankings_g", "bootstrap_sample_stratified", "filtro_pvalor",
                "filtro_mrmr", "filtro_cmim", "filtro_aleatorio", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  i      <- combinaciones_g$ejec[row]
  n_vars <- combinaciones_g$n_vars[row]
  n_b    <- combinaciones_g$n_b[row]
  filtro <- combinaciones_g$filtro[row]
  seed_i <- SEMILLAS[i]
  muestra  <- muestras_g[[i]]
  rankings <- rankings_g[[i]]
  
  boot_dat <- tryCatch(
    bootstrap_sample_stratified(muestra, "obj", n = n_b, seed = seed_i),
    error = function(e) NULL)
  if (is.null(boot_dat)) return(NULL)
  
  datos_filtrados <- tryCatch({
    if (filtro == "aleatorio")    filtro_aleatorio(boot_dat, "obj", n_vars, seed = seed_i)
    else if (filtro == "p_valor") filtro_pvalor  (boot_dat, "obj", rankings, n_vars)
    else if (filtro == "mrmr")    filtro_mrmr    (boot_dat, "obj", rankings, n_vars)
    else                          filtro_cmim    (boot_dat, "obj", rankings, n_vars)
  }, error = function(e) NULL)
  if (is.null(datos_filtrados)) return(NULL)
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = datos_filtrados,
      test_data    = data_test,
      target       = "obj",
      execution_id = sprintf("Procedimiento G %d variables %s %d obs Ejecución %d",
                             n_vars, filtro, n_b, i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$filtro    <- filtro
  out$summary$n_vars    <- n_vars
  out$summary$n_b       <- n_b
  out$summary$ejec      <- i
  out$coefs$filtro      <- filtro
  out$coefs$n_vars      <- n_vars
  out$coefs$n_b         <- n_b
  out$coefs$ejec        <- i
  out$vars_input$filtro <- filtro
  out$vars_input$n_vars <- n_vars
  out$vars_input$n_b    <- n_b
  out$vars_input$ejec   <- i
  out
}

res_g <- combinar_resultados(res_g_list, "G")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_g$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_g$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_g$vars_input)
cat("Procedimiento G completado.\n")

# =============================================================
# PROCEDIMIENTO H
# =============================================================

cat("\n===== PROCEDIMIENTO H =====\n")

aumentos_h   <- c(0, 0.10, 0.20, 0.30)
n_noise_h    <- round(100 * (1 + aumentos_h))
names(n_noise_h) <- paste0(aumentos_h * 100, "%")
noise_factor <- 0.01

rankings_h <- lapply(seq_along(SEMILLAS), function(i) {
  cat(sprintf("  Rankings H ejecución %d / %d\n", i, length(SEMILLAS)))
  muestra <- muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra
  generar_rankings_prefiltros(muestra, target = "obj", top_max = 1000,
                              n_bins_cmim = 3, threads_cmim = 1)
})
muestras_h <- lapply(seq_along(SEMILLAS), function(i)
  muestreo(data_train, "obj", n = 100, seed = SEMILLAS[i])$muestra)
names(rankings_h) <- paste0("semilla_", SEMILLAS)
names(muestras_h) <- paste0("semilla_", SEMILLAS)
saveRDS(rankings_h, file.path(RESULTS_DIR, "rankings_H.rds"))
saveRDS(muestras_h, file.path(RESULTS_DIR, "muestras_H.rds"))

combinaciones_h <- expand.grid(
  ejec   = seq_along(SEMILLAS),
  n_vars = c(50, 100, 200, 300, 500, 1000),
  n_n    = n_noise_h,
  filtro = c("p_valor", "mrmr", "cmim", "aleatorio"),
  stringsAsFactors = FALSE
)

res_h_list <- foreach(
  row       = seq_len(nrow(combinaciones_h)),
  .packages = c("glmnet", "pROC"),
  .export   = c("combinaciones_h", "SEMILLAS", "data_test", "muestras_h",
                "rankings_h", "noise_factor", "add_noise_sample_stratified",
                "filtro_pvalor", "filtro_mrmr", "filtro_cmim",
                "filtro_aleatorio", "elastic_net_single"),
  .errorhandling = "pass"
) %dopar% {
  i      <- combinaciones_h$ejec[row]
  n_vars <- combinaciones_h$n_vars[row]
  n_n    <- combinaciones_h$n_n[row]
  filtro <- combinaciones_h$filtro[row]
  seed_i <- SEMILLAS[i]
  muestra  <- muestras_h[[i]]
  rankings <- rankings_h[[i]]
  
  noise_dat <- tryCatch(
    add_noise_sample_stratified(muestra, "obj", n = n_n,
                                noise_factor = noise_factor, seed = seed_i),
    error = function(e) NULL)
  if (is.null(noise_dat)) return(NULL)
  
  datos_filtrados <- tryCatch({
    if (filtro == "aleatorio")    filtro_aleatorio(noise_dat, "obj", n_vars, seed = seed_i)
    else if (filtro == "p_valor") filtro_pvalor  (noise_dat, "obj", rankings, n_vars)
    else if (filtro == "mrmr")    filtro_mrmr    (noise_dat, "obj", rankings, n_vars)
    else                          filtro_cmim    (noise_dat, "obj", rankings, n_vars)
  }, error = function(e) NULL)
  if (is.null(datos_filtrados)) return(NULL)
  
  out <- tryCatch(
    elastic_net_single(
      train_data   = datos_filtrados,
      test_data    = data_test,
      target       = "obj",
      execution_id = sprintf("Procedimiento H %d variables %s %d obs Ejecución %d",
                             n_vars, filtro, n_n, i),
      alpha = 1, nfolds = 5, seed = seed_i
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  
  out$summary$filtro    <- filtro
  out$summary$n_vars    <- n_vars
  out$summary$n_n       <- n_n
  out$summary$ejec      <- i
  out$coefs$filtro      <- filtro
  out$coefs$n_vars      <- n_vars
  out$coefs$n_n         <- n_n
  out$coefs$ejec        <- i
  out$vars_input$filtro <- filtro
  out$vars_input$n_vars <- n_vars
  out$vars_input$n_n    <- n_n
  out$vars_input$ejec   <- i
  out
}

res_h <- combinar_resultados(res_h_list, "H")
elastic_net_results <- dplyr::bind_rows(elastic_net_results, res_h$summary)
coefs_results       <- dplyr::bind_rows(coefs_results,       res_h$coefs)
vars_input_results  <- dplyr::bind_rows(vars_input_results,  res_h$vars_input)
cat("Procedimiento H completado.\n")

# =============================================================
# GUARDADO FINAL
# =============================================================

saveRDS(elastic_net_results, file.path(RESULTS_DIR, "elastic_net_results_FINAL.rds"))
saveRDS(coefs_results,       file.path(RESULTS_DIR, "coefs_results_FINAL.rds"))
saveRDS(vars_input_results,  file.path(RESULTS_DIR, "vars_input_results_FINAL.rds"))

cat("\n===== SCRIPT 1 COMPLETADO =====\n")
cat("Resultados guardados en:", RESULTS_DIR, "\n")

stopCluster(cluster)