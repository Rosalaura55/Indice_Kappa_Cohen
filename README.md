Índice Kappa de Cohen para evaluar coincidencia entre hotspots y áreas
protegidas
================

- [Descripción del análisis](#descripción-del-análisis)
- [Script original](#script-original)
- [Notas para reproducibilidad](#notas-para-reproducibilidad)
- [Paquetes sugeridos](#paquetes-sugeridos)

# Descripción del análisis

Este documento reproduce el análisis del **índice Kappa de Cohen**,
utilizado para evaluar el grado de coincidencia entre dos
clasificaciones espaciales. En este caso, el análisis puede usarse para
comparar la presencia de **hotspots de riqueza** con la presencia de
**instrumentos de conservación marina** o áreas protegidas.

El índice Kappa permite estimar si el acuerdo observado entre ambas
clasificaciones es mayor al esperado por azar. Valores cercanos a 1
indican alto acuerdo, valores cercanos a 0 indican acuerdo similar al
azar y valores negativos indican menor acuerdo que el esperado por azar.

# Script original

El siguiente bloque contiene el script original convertido a un
documento reproducible para GitHub.

``` r
# ==== KAPPA entre dos capas binarias por hex ====
# Ejemplo: hotspot (riqueza) vs anp_bin (cobertura ANP)

library(sf)
library(dplyr)
library(janitor)   # para limpiar nombres y tablas  (opcional)
library(psych)     # cohen.kappa con IC
library(irr)       # kappa2 (alternativa)

# -------- 1) Cargar capa de hex con métricas --------
gpkg_path <- "hex_final.gpkg"   # <-- ajusta tu ruta
st_layers(gpkg_path)
```

<div class="kable-table">

| name | geomtype | driver | features | fields | crs |
|:---|:---|:---|---:|---:|:---|
| hex_con_metricas\_\_hex_metricas | Multi Polygon | GPKG | 10002 | 22 | WGS 84 / NSIDC EASE-Grid 2.0 Global , PROJCRS\[“WGS 84 / NSIDC EASE-Grid 2.0 Global”, |

    BASEGEOGCRS["WGS 84",
        ENSEMBLE["World Geodetic System 1984 ensemble",
            MEMBER["World Geodetic System 1984 (Transit)"],
            MEMBER["World Geodetic System 1984 (G730)"],
            MEMBER["World Geodetic System 1984 (G873)"],
            MEMBER["World Geodetic System 1984 (G1150)"],
            MEMBER["World Geodetic System 1984 (G1674)"],
            MEMBER["World Geodetic System 1984 (G1762)"],
            MEMBER["World Geodetic System 1984 (G2139)"],
            MEMBER["World Geodetic System 1984 (G2296)"],
            ELLIPSOID["WGS 84",6378137,298.257223563,
                LENGTHUNIT["metre",1]],
            ENSEMBLEACCURACY[2.0]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4326]],
    CONVERSION["US NSIDC EASE-Grid 2.0 Global",
        METHOD["Lambert Cylindrical Equal Area",
            ID["EPSG",9835]],
        PARAMETER["Latitude of 1st standard parallel",30,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8823]],
        PARAMETER["Longitude of natural origin",0,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8802]],
        PARAMETER["False easting",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8806]],
        PARAMETER["False northing",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8807]]],
    CS[Cartesian,2],
        AXIS["easting (X)",east,
            ORDER[1],
            LENGTHUNIT["metre",1]],
        AXIS["northing (Y)",north,
            ORDER[2],
            LENGTHUNIT["metre",1]],
    USAGE[
        SCOPE["Environmental science - used as basis for EASE grid."],
        AREA["World between 86°S and 86°N."],
        BBOX[-86,-180,86,180]],
    ID["EPSG",6933]] |

</div>

``` r
layer_name <- "hex_con_metricas__hex_metricas"         

hex <- st_read(gpkg_path, layer = layer_name, quiet = TRUE) %>%
  janitor::clean_names()

# Checa qué campos tienes
names(hex)
```

    ##  [1] "id"           "left"         "top"          "right"        "bottom"      
    ##  [6] "row_index"    "col_index"    "area"         "id2"          "n_records"   
    ## [11] "rich_raw"     "use"          "hotspot"      "id2_str"      "area_anp_m2" 
    ## [16] "area_hex_km2" "area_anp_km2" "prop_anp"     "pct_anp"      "cover_amp"   
    ## [21] "anp_bin"      "amp_bin_01"   "geom"

``` r
# -------- 2) Seleccionar columnas binarias y filtrar --------
# Ajusta los nombres si difieren en tu base:
campo_a <- "hotspot"   # binario 0/1 
campo_b <- "anp_bin"   # binario 0/1 (p.ej. >= 30% ANP)

usar_filtro_use <- TRUE   # pon FALSE si no quieres filtrar por 'use'
umbral_use <- 1

df <- hex %>%
  st_drop_geometry() %>%
  { if (usar_filtro_use && "use" %in% names(.)) filter(., .data[["use"]] >= umbral_use) else . } %>%
  transmute(
    a = .data[[campo_a]],
    b = .data[[campo_b]]
  ) %>%
  mutate(
    a = as.integer(a),
    b = as.integer(b)
  ) %>%
  filter(!is.na(a), !is.na(b))

# Comprobación rápida
table(df$a, df$b)
```

    ##    
    ##      0  1
    ##   0 40 50
    ##   1  7  5

``` r
# -------- 3) Matriz de confusión y métricas básicas --------
conf <- table(A = df$a, B = df$b)      # A = hotspot, B = anp_bin
conf
```

    ##    B
    ## A    0  1
    ##   0 40 50
    ##   1  7  5

``` r
# Métricas descriptivas (tratando A=1 como "positivo"):
TP <- ifelse("1" %in% rownames(conf) && "1" %in% colnames(conf), conf["1","1"], 0)
FP <- ifelse("1" %in% rownames(conf) && "0" %in% colnames(conf), conf["1","0"], 0)
FN <- ifelse("0" %in% rownames(conf) && "1" %in% colnames(conf), conf["0","1"], 0)
TN <- ifelse("0" %in% rownames(conf) && "0" %in% colnames(conf), conf["0","0"], 0)

acc  <- (TP + TN) / sum(conf)                      # exactitud
sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)  # sensibilidad (recall de A=1 respecto a B)
spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)  # especificidad
ppv  <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)  # precisión (PPV)
npv  <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)  # NPV
prev_a <- mean(df$a == 1)                          # prevalencia de A=1
prev_b <- mean(df$b == 1)                          # prevalencia de B=1

list(
  confusion_matrix = conf,
  accuracy = acc,
  sensitivity = sens,
  specificity = spec,
  ppv = ppv,
  npv = npv,
  prevalence_a = prev_a,
  prevalence_b = prev_b
)
```

    ## $confusion_matrix
    ##    B
    ## A    0  1
    ##   0 40 50
    ##   1  7  5
    ## 
    ## $accuracy
    ## [1] 0.4411765
    ## 
    ## $sensitivity
    ## [1] 0.09090909
    ## 
    ## $specificity
    ## [1] 0.8510638
    ## 
    ## $ppv
    ## [1] 0.4166667
    ## 
    ## $npv
    ## [1] 0.4444444
    ## 
    ## $prevalence_a
    ## [1] 0.1176471
    ## 
    ## $prevalence_b
    ## [1] 0.5392157

``` r
# -------- 4) Cohen's Kappa (con IC) --------
# Opción 1: psych::cohen.kappa (acepta matriz de 2 columnas o tabla)
ck <- psych::cohen.kappa(cbind(df$a, df$b))
ck
```

    ## Call: cohen.kappa1(x = x, w = w, n.obs = n.obs, alpha = alpha, levels = levels, 
    ##     w.exp = w.exp)
    ## 
    ## Cohen Kappa and Weighted Kappa correlation coefficients and confidence boundaries 
    ##                  lower estimate upper
    ## unweighted kappa -0.17   -0.054 0.065
    ## weighted kappa   -0.17   -0.054 0.065
    ## 
    ##  Number of subjects = 102


