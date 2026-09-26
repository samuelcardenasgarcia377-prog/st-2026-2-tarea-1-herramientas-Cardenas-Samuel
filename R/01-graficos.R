# R/01-graficos.R
# Gráficos, ACF manual y correlograma
# ============================================================


# ============================================================
# graficar_serie()
# ============================================================

graficar_serie <- function(
    datos,
    titulo
) {
  
  stopifnot(
    is.data.frame(datos),
    all(
      c("t", "fecha", "y") %in%
        names(datos)
    ),
    is.character(titulo),
    length(titulo) == 1,
    !is.na(titulo)
  )
  
  if (nrow(datos) == 0) {
    
    stop(
      "No hay observaciones para graficar."
    )
  }
  
  if (!inherits(
    datos$fecha,
    "Date"
  )) {
    
    stop(
      "La columna fecha debe ser de clase Date."
    )
  }
  
  if (!is.numeric(datos$y)) {
    
    stop(
      "La columna y debe ser numérica."
    )
  }
  
  if (anyNA(datos$y) ||
      any(!is.finite(datos$y))) {
    
    stop(
      "La columna y contiene valores no finitos."
    )
  }
  
  
  fuente <- attr(
    datos,
    "fuente"
  )
  
  unidad <- attr(
    datos,
    "unidad"
  )
  
  
  if (is.null(fuente) ||
      length(fuente) == 0) {
    
    fuente <- "Fuente no especificada"
  }
  
  
  if (is.null(unidad) ||
      length(unidad) == 0) {
    
    unidad <- "unidades"
  }
  
  
  n <- nrow(datos)
  
  
  # ----------------------------------------------------------
  # Número de etiquetas de fecha
  # ----------------------------------------------------------
  
  numero_cortes <- min(
    8,
    n
  )
  
  numero_cortes <- max(
    2,
    numero_cortes
  )
  
  posiciones <- unique(
    round(
      seq(
        1,
        n,
        length.out = numero_cortes
      )
    )
  )
  
  fechas_corte <- datos$fecha[
    posiciones
  ]
  
  
  # ----------------------------------------------------------
  # Gráfico
  # ----------------------------------------------------------
  
  ggplot(
    datos,
    aes(
      x = fecha,
      y = y
    )
  ) +
    
    geom_line(
      linewidth = 0.7
    ) +
    
    labs(
      title = titulo,
      x = "Fecha",
      y = unidad,
      caption = paste0(
        "Fuente: ",
        fuente,
        " | n = ",
        n
      )
    ) +
    
    scale_x_date(
      breaks = fechas_corte,
      date_labels = "%Y-%m"
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


# ============================================================
# acf_manual()
#
# ACF muestral calculada manualmente:
#
#                 sum_{t=h+1}^T
# (Y_t-Ybar)(Y_{t-h}-Ybar)
#
# r_h = ---------------------------------
#                 sum_{t=1}^T
#              (Y_t-Ybar)^2
#
# Se utiliza UN ÚNICO divisor T en la normalización.
# ============================================================

acf_manual <- function(
    y,
    m = max(
      1,
      min(
        floor(length(y) / 4),
        24
      )
    )
) {
  
  stopifnot(
    is.numeric(y)
  )
  
  if (length(y) < 2) {
    
    stop(
      "La serie debe tener al menos dos observaciones."
    )
  }
  
  if (anyNA(y) ||
      any(!is.finite(y))) {
    
    stop(
      "y contiene valores NA, NaN o infinitos."
    )
  }
  
  T <- length(y)
  
  
  # ----------------------------------------------------------
  # Validación de m
  # ----------------------------------------------------------
  
  if (
    length(m) != 1 ||
    !is.numeric(m) ||
    !is.finite(m) ||
    m != floor(m)
  ) {
    
    stop(
      "m debe ser un único entero finito."
    )
  }
  
  m <- as.integer(m)
  
  if (m < 1) {
    
    stop(
      "m debe ser al menos 1."
    )
  }
  
  if (m >= T) {
    
    m <- T - 1
  }
  
  
  # ----------------------------------------------------------
  # Media
  # ----------------------------------------------------------
  
  media <- mean(y)
  
  
  # ----------------------------------------------------------
  # Denominador
  # ----------------------------------------------------------
  
  denominador <- sum(
    (y - media)^2
  )
  
  if (denominador <= 0) {
    
    stop(
      "La serie es constante; la ACF no está definida."
    )
  }
  
  
  # ----------------------------------------------------------
  # ACF manual
  # ----------------------------------------------------------
  
  r <- numeric(m)
  
  for (h in seq_len(m)) {
    
    numerador <- sum(
      (y[(h + 1):T] - media) *
        (y[1:(T - h)] - media)
    )
    
    r[h] <-
      numerador /
      denominador
  }
  
  
  tibble::tibble(
    lag = seq_len(m),
    acf = r
  )
}


# ============================================================
# correlograma()
#
# Devuelve un gráfico ACF/PACF vertical.
#
# ACF:
#   calculada manualmente.
#
# PACF:
#   stats::pacf(plot = FALSE)
#
# Banda:
#   qnorm((1 + .95)/2) / sqrt(T)
#
# Lag 0:
#   excluido del gráfico.
# ============================================================

correlograma <- function(
    datos,
    m = max(
      1,
      min(
        floor(nrow(datos) / 4),
        24
      )
    )
) {
  
  stopifnot(
    is.data.frame(datos),
    all(
      c("t", "fecha", "y") %in%
        names(datos)
    )
  )
  
  y <- datos$y
  
  
  if (!is.numeric(y)) {
    
    stop(
      "La columna y debe ser numérica."
    )
  }
  
  if (anyNA(y) ||
      any(!is.finite(y))) {
    
    stop(
      "La serie contiene valores no finitos."
    )
  }
  
  T <- length(y)
  
  
  if (T < 2) {
    
    stop(
      "La serie debe tener al menos dos observaciones."
    )
  }
  
  
  # ----------------------------------------------------------
  # Validación de m
  # ----------------------------------------------------------
  
  if (
    length(m) != 1 ||
    !is.numeric(m) ||
    !is.finite(m) ||
    m != floor(m)
  ) {
    
    stop(
      "m debe ser un único entero finito."
    )
  }
  
  m <- as.integer(m)
  
  if (m < 1) {
    
    stop(
      "m debe ser al menos 1."
    )
  }
  
  m <- min(
    m,
    T - 1
  )
  
  
  # ----------------------------------------------------------
  # ACF manual
  # ----------------------------------------------------------
  
  acf_datos <- acf_manual(
    y,
    m = m
  )
  
  
  # ----------------------------------------------------------
  # PACF
  # ----------------------------------------------------------
  
  pacf_resultado <- stats::pacf(
    y,
    lag.max = m,
    plot = FALSE
  )
  
  pacf_datos <- tibble::tibble(
    lag = seq_len(m),
    pacf = as.numeric(
      pacf_resultado$acf[
        1:m
      ]
    )
  )
  
  
  # ----------------------------------------------------------
  # Banda de ruido blanco al 95 %
  #
  # Equivalente a la utilizada por plot.acf()
  # bajo la aproximación normal.
  # ----------------------------------------------------------
  
  banda <- stats::qnorm(
    (1 + 0.95) / 2
  ) / sqrt(T)
  
  
  # ----------------------------------------------------------
  # Gráfico ACF
  # ----------------------------------------------------------
  
  grafico_acf <- ggplot(
    acf_datos,
    aes(
      x = lag,
      y = acf
    )
  ) +
    
    geom_hline(
      yintercept = 0
    ) +
    
    geom_hline(
      yintercept = c(
        -banda,
        banda
      ),
      linetype = "dashed"
    ) +
    
    geom_segment(
      aes(
        xend = lag,
        yend = 0
      ),
      linewidth = 0.6
    ) +
    
    labs(
      title = "Función de autocorrelación (ACF)",
      x = "Rezago",
      y = "ACF"
    ) +
    
    theme_minimal()
  
  
  # ----------------------------------------------------------
  # Gráfico PACF
  # ----------------------------------------------------------
  
  grafico_pacf <- ggplot(
    pacf_datos,
    aes(
      x = lag,
      y = pacf
    )
  ) +
    
    geom_hline(
      yintercept = 0
    ) +
    
    geom_hline(
      yintercept = c(
        -banda,
        banda
      ),
      linetype = "dashed"
    ) +
    
    geom_segment(
      aes(
        xend = lag,
        yend = 0
      ),
      linewidth = 0.6
    ) +
    
    labs(
      title = "Función de autocorrelación parcial (PACF)",
      x = "Rezago",
      y = "PACF"
    ) +
    
    theme_minimal()
  
  
  # ----------------------------------------------------------
  # Panel final
  # ----------------------------------------------------------
  
  grafico_acf /
    grafico_pacf
}

# ============================================================
