Análisis Kappa de Cohen para hotspots P70
================

- [Descripción del análisis](#descripción-del-análisis)
- [Script original](#script-original)

# Descripción del análisis

Este documento reproduce el análisis del **índice Kappa de Cohen**
aplicado a hotspots definidos mediante el percentil 70 (P70). El
objetivo es evaluar el grado de coincidencia entre las celdas
identificadas como hotspots de riqueza y la presencia de instrumentos de
conservación marina o áreas protegidas.

El índice Kappa permite cuantificar si el acuerdo observado es mayor al
esperado por azar.

Interpretación general:

- $Kappa < 0$: acuerdo menor al esperado por azar
- $Kappa pprox 0: acuerdo similar al azar
- $0.01 - 0.20: acuerdo leve
- $0.21 - 0.40: acuerdo aceptable
- $0.41 - 0.60: acuerdo moderado
- $0.61 - 0.80: acuerdo sustancial
- $0.81 - 1.00: acuerdo casi perfecto

# Script original

El siguiente bloque contiene el script original convertido a un
documento reproducible para GitHub.

``` r
library(sf)
library(dplyr)

# --- AJUSTA AQUÍ ---
gpkg_path  <- "hex60km_todo.gpkg"
layer_name <- "hex60_con_metricas___hex_metricas"
campo_riqueza <- "rich_raw"
campo_use     <- "use"
campo_binario <- "amp_bin"   # tu campo 0/1 de ANP

hex <- st_read(gpkg_path, layer = layer_name, quiet = TRUE) %>% st_drop_geometry()
stopifnot(all(c(campo_riqueza, campo_use, campo_binario) %in% names(hex)))

df <- hex %>% 
  filter(.data[[campo_use]] == 1) %>%
  transmute(rich = .data[[campo_riqueza]],
            binB = as.integer(.data[[campo_binario]])) %>%
  filter(!is.na(rich), !is.na(binB))


library(psych)  # cohen.kappa
library(irr)    # kappa2
library(dplyr)
library(readr)  # write_csv (opcional)

# --- función que calcula TODO para un percentil p ---
eval_concordancia <- function(d, p = 0.70, Btimes = 2000, seed = 123){
  stopifnot(all(c("rich","binB") %in% names(d)))
  d <- d %>% filter(!is.na(rich), !is.na(binB))
  if(nrow(d) < 5) stop("Muy pocas filas tras el filtro.")
  
  thr <- as.numeric(quantile(d$rich, probs = p, na.rm = TRUE, type = 7))
  A <- as.integer(d$rich >= thr)   # hotspot binario
  B <- as.integer(d$binB)          # anp_bin
  
  # Matriz de confusión con niveles fijos para no fallar cuando falta una celda
  conf <- table(factor(A, levels = c(0,1)), factor(B, levels = c(0,1)))
  TN <- conf["0","0"]; FP <- conf["1","0"]; FN <- conf["0","1"]; TP <- conf["1","1"]
  N  <- sum(conf)
  
  # Métricas clásicas
  acc  <- (TP + TN) / N
  prec <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)      # PPV
  rec  <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)      # Sensibilidad
  f1   <- ifelse(is.na(prec) | is.na(rec) | (prec+rec)==0, NA_real_, 2*prec*rec/(prec+rec))
  jacc <- ifelse((TP + FP + FN) > 0, TP / (TP + FP + FN), NA_real_)
  prevA <- mean(A == 1); prevB <- mean(B == 1)
  
  # Kappa (psych) con IC
  ck <- psych::cohen.kappa(cbind(A,B))
  kap  <- ck$kappa
  ci_l <- as.numeric(ck$confid[1]); ci_h <- as.numeric(ck$confid[2])
  
  # Bootstrap para IC de Kappa (percentil)
  set.seed(seed)
  vals <- numeric(Btimes)
  n <- length(A)
  for(i in seq_len(Btimes)){
    idx <- sample.int(n, n, replace = TRUE)
    vals[i] <- irr::kappa2(cbind(A[idx], B[idx]), weight = "unweighted")$value
  }
  boot_est <- mean(vals, na.rm = TRUE)
  boot_l   <- as.numeric(quantile(vals, 0.025, na.rm = TRUE))
  boot_h   <- as.numeric(quantile(vals, 0.975, na.rm = TRUE))
  
  tibble(
    percentil = p,
    umbral_rich = thr,
    n_total = n,
    n_pos = sum(A==1),
    n_neg = sum(A==0),
    TP = as.integer(TP), FP = as.integer(FP), FN = as.integer(FN), TN = as.integer(TN),
    accuracy = acc,
    precision = prec,
    recall = rec,
    F1 = f1,
    jaccard = jacc,
    prev_A = prevA, prev_B = prevB,
    kappa = kap, kappa_ci_low = ci_l, kappa_ci_high = ci_h,
    kappa_boot = boot_est, kappa_boot_low = boot_l, kappa_boot_high = boot_h
  )
}

# --- correr para P70, P80, P90 ---
percentiles <- c(0.70, 0.80, 0.90)
res <- bind_rows(lapply(percentiles, function(p) eval_concordancia(df, p = p, Btimes = 2000)))

# ver resultados
res
```

<div class="kable-table">

| percentil | umbral_rich | n_total | n_pos | n_neg | TP | FP | FN | TN | accuracy | precision | recall | F1 | jaccard | prev_A | prev_B | kappa | kappa_ci_low | kappa_ci_high | kappa_boot | kappa_boot_low | kappa_boot_high |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.7 | 14.0 | 105 | 36 | 69 | 20 | 16 | 42 | 27 | 0.4476190 | 0.5555556 | 0.3225806 | 0.4081633 | 0.2564103 | 0.3428571 | 0.5904762 | -0.0453141 | -0.2152403 | -0.2152403 | -0.0485854 | -0.2218020 | 0.1226630 |
| 0.8 | 19.4 | 105 | 21 | 84 | 11 | 10 | 51 | 33 | 0.4190476 | 0.5238095 | 0.1774194 | 0.2650602 | 0.1527778 | 0.2000000 | 0.5904762 | -0.0481100 | -0.1864082 | -0.1864082 | -0.0489479 | -0.1918750 | 0.0842294 |
| 0.9 | 29.6 | 105 | 11 | 94 | 6 | 5 | 56 | 38 | 0.4190476 | 0.5454545 | 0.0967742 | 0.1643836 | 0.0895522 | 0.1047619 | 0.5904762 | -0.0165053 | -0.1188178 | -0.1188178 | -0.0169375 | -0.1184456 | 0.0815442 |

</div>

``` r
# (opcional) exportar a CSV
write_csv(res, "kappa_jaccard_f1_bootstrap_P70_P80_P90.csv")


#============TABLA COMPACTA==================
library(dplyr)
tabla_compacta <- res %>%
  mutate(
    pct_hotspots_en_ANP = TP / (TP + FP),
    pct_ANP_que_son_hot = TP / (TP + FN)
  ) %>%
  transmute(
    percentil,
    umbral_rich = round(umbral_rich,1),
    n_pos, n_neg,
    kappa = round(kappa,3),
    kappa_boot_CI = sprintf("[%.3f, %.3f]", kappa_boot_low, kappa_boot_high),
    F1 = round(F1,3),
    Jaccard = round(jaccard,3),
    hotspots_en_ANP = scales::percent(pct_hotspots_en_ANP),
    ANP_que_son_hot = scales::percent(pct_ANP_que_son_hot)
  )
tabla_compacta
```

<div class="kable-table">

| percentil | umbral_rich | n_pos | n_neg | kappa | kappa_boot_CI | F1 | Jaccard | hotspots_en_ANP | ANP_que_son_hot |
|---:|---:|---:|---:|---:|:---|---:|---:|:---|:---|
| 0.7 | 14.0 | 36 | 69 | -0.045 | \[-0.222, 0.123\] | 0.408 | 0.256 | 55.6% | 32.3% |
| 0.8 | 19.4 | 21 | 84 | -0.048 | \[-0.192, 0.084\] | 0.265 | 0.153 | 52.4% | 17.7% |
| 0.9 | 29.6 | 11 | 94 | -0.017 | \[-0.118, 0.082\] | 0.164 | 0.090 | 54.5% | 9.7% |

</div>
