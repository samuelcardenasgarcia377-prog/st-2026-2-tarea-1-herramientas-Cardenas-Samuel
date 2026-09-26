# R/03-evaluacion.R
# ============================================================

# ============================================================
# 03-evaluacion.R
# Métricas, pruebas y validación
# ============================================================


# ============================================================
# VALIDACIÓN DE VECTORES
# ============================================================

.validar_vectores <- function(
    y,
    yhat
) {
  
  stopifnot(
    is.numeric(y),
    is.numeric(yhat),
    length(y) == length(yhat)
  )
  
  if (
    anyNA(y) ||
    anyNA(yhat) ||
    any(!is.finite(y)) ||
    any(!is.finite(yhat))
  ) {
    
    stop(
      "y y yhat deben contener únicamente valores finitos."
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# MÉTRICAS
# ============================================================

medidas <- function(
    y,
    yhat,
    escala_mase = NULL
) {
  
  .validar_vectores(
    y,
    yhat
  )
  
  errores <-
    y - yhat
  
  mse <-
    mean(
      errores^2
    )
  
  mad <-
    mean(
      abs(errores)
    )
  
  
  # ----------------------------------------------------------
  # MAPE
  # ----------------------------------------------------------
  
  indices_mape <-
    abs(y) > 0
  
  if (any(indices_mape)) {
    
    mape <-
      mean(
        abs(
          errores[
            indices_mape
          ] /
            y[
              indices_mape
            ]
        )
      ) * 100
    
  } else {
    
    mape <- NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Escala MASE
  # ----------------------------------------------------------
  
  if (is.null(escala_mase)) {
    
    if (length(y) < 2) {
      
      escala_mase <- NA_real_
      
    } else {
      
      escala_mase <-
        mean(
          abs(
            diff(y)
          )
        )
    }
  }
  
  
  if (
    is.na(escala_mase) ||
    !is.finite(escala_mase) ||
    escala_mase <= 0
  ) {
    
    mase <- NA_real_
    
  } else {
    
    mase <-
      mad /
      escala_mase
  }
  
  
  tibble::tibble(
    MSE = mse,
    MAD = mad,
    MAPE = mape,
    MASE = mase
  )
}


# ============================================================
# INGENUO
# ============================================================

pronostico_ingenuo <- function(
    y,
    h
) {
  
  .validar_y(y)
  
  h <-
    .validar_h(h)
  
  rep(
    tail(y, 1),
    h
  )
}


pronostico_ingenuo_estacional <- function(
    y,
    frecuencia,
    h
) {
  
  .validar_y(y)
  
  stopifnot(
    is.numeric(frecuencia),
    length(frecuencia) == 1,
    is.finite(frecuencia),
    frecuencia == as.integer(frecuencia),
    frecuencia >= 1
  )
  
  frecuencia <-
    as.integer(
      frecuencia
    )
  
  h <-
    .validar_h(h)
  
  if (
    length(y) <
    frecuencia
  ) {
    
    stop(
      "No hay suficientes observaciones para el ingenuo estacional."
    )
  }
  
  rep(
    tail(
      y,
      frecuencia
    ),
    length.out = h
  )
}


escala_ingenuo <- function(y) {
  
  .validar_y(y)
  
  mean(
    abs(
      diff(y)
    )
  )
}


escala_ingenuo_estacional <- function(
    y,
    frecuencia
) {
  
  .validar_y(y)
  
  stopifnot(
    is.numeric(frecuencia),
    length(frecuencia) == 1,
    is.finite(frecuencia),
    frecuencia == as.integer(frecuencia),
    frecuencia >= 1
  )
  
  frecuencia <-
    as.integer(
      frecuencia
    )
  
  if (
    frecuencia >=
    length(y)
  ) {
    
    stop(
      "frecuencia debe ser menor que la longitud de y."
    )
  }
  
  mean(
    abs(
      y[
        (frecuencia + 1):
          length(y)
      ] -
        y[
          1:
            (length(y) - frecuencia)
        ]
    )
  )
}


# ============================================================
# LJUNG-BOX
# ============================================================

ljung_box <- function(
    r,
    T,
    m,
    p = 0
) {
  
  stopifnot(
    is.numeric(r),
    is.numeric(T),
    length(T) == 1,
    is.finite(T),
    is.numeric(m),
    length(m) == 1,
    is.finite(m),
    is.numeric(p),
    length(p) == 1,
    is.finite(p)
  )
  
  if (
    T <= m ||
    m < 1 ||
    m != as.integer(m) ||
    p < 0 ||
    p != as.integer(p) ||
    m - p <= 0
  ) {
    
    stop(
      "Se requiere T > m >= 1 y m-p > 0."
    )
  }
  
  m <-
    as.integer(m)
  
  p <-
    as.integer(p)
  
  if (
    length(r) <
    m
  ) {
    
    stop(
      "r no contiene suficientes autocorrelaciones."
    )
  }
  
  r <-
    r[
      seq_len(m)
    ]
  
  if (
    anyNA(r) ||
    any(!is.finite(r))
  ) {
    
    stop(
      "Las autocorrelaciones deben ser finitas."
    )
  }
  
  h <-
    seq_len(m)
  
  
  # ----------------------------------------------------------
  # Estadístico:
  #
  # Q_m =
  # T(T+2) sum_{h=1}^m r_h^2/(T-h)
  # ----------------------------------------------------------
  
  Q <-
    T *
    (T + 2) *
    sum(
      r^2 /
        (T - h)
    )
  
  gl <-
    m - p
  
  valor_critico <-
    stats::qchisq(
      0.95,
      df = gl
    )
  
  p_valor <-
    stats::pchisq(
      Q,
      df = gl,
      lower.tail = FALSE
    )
  
  list(
    estadistico = Q,
    gl = gl,
    valor_critico_5 =
      valor_critico,
    p_valor = p_valor,
    rechazar_5 =
      Q > valor_critico
  )
}


# ============================================================
# JARQUE-BERA
# ============================================================

jarque_bera <- function(e) {
  
  stopifnot(
    is.numeric(e),
    length(e) >= 3
  )
  
  if (
    anyNA(e) ||
    any(!is.finite(e))
  ) {
    
    stop(
      "Los residuos deben ser finitos."
    )
  }
  
  N <-
    length(e)
  
  media <-
    mean(e)
  
  m2 <-
    mean(
      (e - media)^2
    )
  
  m3 <-
    mean(
      (e - media)^3
    )
  
  m4 <-
    mean(
      (e - media)^4
    )
  
  if (m2 <= 0) {
    
    stop(
      "La varianza de los residuos es cero."
    )
  }
  
  A <-
    m3 /
    m2^(3 / 2)
  
  K <-
    m4 /
    m2^2
  
  JB <-
    N / 6 *
    (
      A^2 +
        (K - 3)^2 / 4
    )
  
  valor_critico <-
    stats::qchisq(
      0.95,
      df = 2
    )
  
  p_valor <-
    stats::pchisq(
      JB,
      df = 2,
      lower.tail = FALSE
    )
  
  list(
    estadistico = JB,
    gl = 2,
    valor_critico_5 =
      valor_critico,
    p_valor = p_valor,
    asimetria = A,
    curtosis = K,
    rechazar_5 =
      JB > valor_critico
  )
}


# ============================================================
# DURBIN-WATSON
# ============================================================

durbin_watson <- function(e) {
  
  stopifnot(
    is.numeric(e),
    length(e) >= 2
  )
  
  if (
    anyNA(e) ||
    any(!is.finite(e))
  ) {
    
    stop(
      "Los residuos deben ser finitos."
    )
  }
  
  denominador <-
    sum(
      e^2
    )
  
  if (
    denominador <= 0
  ) {
    
    stop(
      "El denominador del estadístico Durbin-Watson es cero."
    )
  }
  
  sum(
    diff(e)^2
  ) /
    denominador
}


# ============================================================
# PRUEBA DE MEDIA CERO
# ============================================================

prueba_media_cero <- function(e) {
  
  stopifnot(
    is.numeric(e),
    length(e) >= 2
  )
  
  if (
    anyNA(e) ||
    any(!is.finite(e))
  ) {
    
    stop(
      "Los residuos deben ser finitos."
    )
  }
  
  n <-
    length(e)
  
  media <-
    mean(e)
  
  s <-
    stats::sd(e)
  
  if (
    s <= 0 ||
    !is.finite(s)
  ) {
    
    stop(
      "La desviación estándar es cero."
    )
  }
  
  estadistico <-
    media /
    (
      s /
        sqrt(n)
    )
  
  gl <-
    n - 1
  
  valor_critico <-
    stats::qt(
      0.975,
      df = gl
    )
  
  p_valor <-
    2 *
    stats::pt(
      -abs(estadistico),
      df = gl
    )
  
  list(
    estadistico = estadistico,
    gl = gl,
    valor_critico_5 =
      valor_critico,
    p_valor = p_valor,
    rechazar_5 =
      abs(estadistico) >
      valor_critico
  )
}


# ============================================================
# ACF INDIVIDUAL
# ============================================================

autocorrelaciones_individuales <- function(
    x,
    m = NULL,
    ci = 0.95
) {
  
  stopifnot(
    is.numeric(x),
    length(x) >= 2,
    is.numeric(ci),
    length(ci) == 1,
    is.finite(ci),
    ci > 0,
    ci < 1
  )
  
  if (
    anyNA(x) ||
    any(!is.finite(x))
  ) {
    
    stop(
      "x debe contener valores finitos."
    )
  }
  
  T <-
    length(x)
  
  if (is.null(m)) {
    
    m <-
      max(
        1,
        min(
          floor(T / 4),
          24
        )
      )
  }
  
  stopifnot(
    is.numeric(m),
    length(m) == 1,
    is.finite(m),
    m == as.integer(m)
  )
  
  m <-
    as.integer(m)
  
  if (
    m < 1 ||
    m >= T
  ) {
    
    stop(
      "Debe cumplirse 1 <= m < length(x)."
    )
  }
  
  media <-
    mean(x)
  
  denominador <-
    sum(
      (x - media)^2
    )
  
  if (
    denominador <= 0
  ) {
    
    stop(
      "La serie es constante."
    )
  }
  
  r <-
    numeric(m)
  
  for (h in seq_len(m)) {
    
    r[h] <-
      sum(
        (x[(h + 1):T] - media) *
          (x[1:(T - h)] - media)
      ) /
      denominador
  }
  
  banda <-
    stats::qnorm(
      (1 + ci) / 2
    ) /
    sqrt(T)
  
  tibble::tibble(
    lag = seq_len(m),
    acf = r,
    limite_inferior = -banda,
    limite_superior = banda,
    significativo =
      abs(r) > banda
  )
}


# ============================================================
# VALIDACIÓN INTEGRAL DE ERRORES
# ============================================================

validar_errores <- function(
    errores,
    p = 0,
    m = NULL,
    titulo = "Validación de errores",
    ajustados = NULL
) {
  
  stopifnot(
    is.numeric(errores),
    length(errores) >= 3,
    is.numeric(p),
    length(p) == 1,
    p >= 0,
    p == as.integer(p)
  )
  
  if (
    anyNA(errores) ||
    any(!is.finite(errores))
  ) {
    
    stop(
      "Los errores deben ser finitos."
    )
  }
  
  T <-
    length(errores)
  
  p <-
    as.integer(p)
  
  
  if (is.null(m)) {
    
    m <-
      max(
        1,
        min(
          floor(T / 4),
          24
        )
      )
  }
  
  stopifnot(
    is.numeric(m),
    length(m) == 1,
    is.finite(m),
    m == as.integer(m)
  )
  
  m <-
    as.integer(m)
  
  m <-
    min(
      m,
      T - 1
    )
  
  if (
    m < 1 ||
    m - p <= 0
  ) {
    
    stop(
      "Debe cumplirse m-p > 0."
    )
  }
  
  
  acf_errores <-
    autocorrelaciones_individuales(
      errores,
      m = m
    )
  
  lb <-
    ljung_box(
      acf_errores$acf,
      T = T,
      m = m,
      p = p
    )
  
  media <-
    prueba_media_cero(
      errores
    )
  
  jb <-
    jarque_bera(
      errores
    )
  
  dw <-
    durbin_watson(
      errores
    )
  
  
  resultados <- list(
    T = T,
    m = m,
    p = p,
    acf = acf_errores,
    ljung_box = lb,
    media_cero = media,
    jarque_bera = jb,
    durbin_watson = dw
  )
  
  
  # ----------------------------------------------------------
  # ACF de residuos
  # ----------------------------------------------------------
  
  grafico_acf <-
    ggplot2::ggplot(
      acf_errores,
      ggplot2::aes(
        x = lag,
        y = acf
      )
    ) +
    
    ggplot2::geom_hline(
      yintercept = 0
    ) +
    
    ggplot2::geom_hline(
      ggplot2::aes(
        yintercept =
          limite_inferior
      ),
      linetype = "dashed"
    ) +
    
    ggplot2::geom_hline(
      ggplot2::aes(
        yintercept =
          limite_superior
      ),
      linetype = "dashed"
    ) +
    
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = lag,
        yend = 0
      )
    ) +
    
    ggplot2::labs(
      title =
        paste(
          titulo,
          "- ACF de errores"
        ),
      x = "Rezago",
      y = "ACF"
    ) +
    
    ggplot2::theme_minimal()
  
  resultados$grafico_acf <-
    grafico_acf
  
  
  # ----------------------------------------------------------
  # Errores sobre el tiempo
  # ----------------------------------------------------------
  
  datos_errores_tiempo <-
    tibble::tibble(
      t = seq_along(errores),
      errores = errores
    )
  
  grafico_errores_tiempo <-
    ggplot2::ggplot(
      datos_errores_tiempo,
      ggplot2::aes(
        x = t,
        y = errores
      )
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste(titulo, "- Errores sobre el tiempo"),
      x = "Tiempo",
      y = "Error"
    ) +
    ggplot2::theme_minimal()
  
  resultados$grafico_errores_tiempo <-
    grafico_errores_tiempo
  
  
  # ----------------------------------------------------------
  # Residuos vs ajustados
  # ----------------------------------------------------------
  
  if (!is.null(ajustados)) {
    
    if (
      length(ajustados) !=
      length(errores)
    ) {
      
      stop(
        "ajustados y errores deben tener la misma longitud."
      )
    }
    
    if (
      anyNA(ajustados) ||
      any(!is.finite(ajustados))
    ) {
      
      stop(
        "ajustados debe contener únicamente valores finitos."
      )
    }
    
    datos_residuos <-
      tibble::tibble(
        ajustados = ajustados,
        residuos = errores
      )
    
    grafico_residuos <-
      ggplot2::ggplot(
        datos_residuos,
        ggplot2::aes(
          x = ajustados,
          y = residuos
        )
      ) +
      
      ggplot2::geom_point() +
      
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      
      ggplot2::labs(
        title =
          paste(
            titulo,
            "- Residuos vs ajustados"
          ),
        x = "Valores ajustados",
        y = "Residuos"
      ) +
      
      ggplot2::theme_minimal()
    
    resultados$grafico_residuos <-
      grafico_residuos
  }
  
  
  resultados
}


# ============================================================
# VERIFICACIÓN DE LJUNG-BOX
# ============================================================

verificar_ljung_box <- function(
    x,
    m,
    p = 0
) {
  
  .validar_y(x)
  
  m <-
    as.integer(m)
  
  r <-
    acf_manual(
      x,
      m = m
    )$acf
  
  resultado_manual <-
    ljung_box(
      r = r,
      T = length(x),
      m = m,
      p = p
    )
  
  resultado_R <-
    stats::Box.test(
      x,
      lag = m,
      type = "Ljung-Box",
      fitdf = p
    )
  
  diferencia_Q <-
    abs(
      resultado_manual$estadistico -
        unname(
          resultado_R$statistic
        )
    )
  
  diferencia_p <-
    abs(
      resultado_manual$p_valor -
        resultado_R$p.value
    )
  
  list(
    manual = resultado_manual,
    R = resultado_R,
    diferencia_Q = diferencia_Q,
    diferencia_p = diferencia_p,
    coincide =
      diferencia_Q < 1e-10 &&
      diferencia_p < 1e-10
  )
}


# ============================================================
# VERIFICACIÓN SES
# ============================================================

verificar_ses <- function(
    y,
    alpha
) {
  
  .validar_y(y)
  
  ajuste <-
    ajustar_ses(
      y,
      alpha
    )
  
  n <-
    length(y)
  
  nivel_ce <-
    numeric(n)
  
  nivel_ce[1] <-
    y[1]
  
  for (t in 2:n) {
    
    error_t <-
      y[t] -
      nivel_ce[t - 1]
    
    nivel_ce[t] <-
      nivel_ce[t - 1] +
      alpha *
      error_t
  }
  
  diferencia <-
    max(
      abs(
        ajuste$parametros$nivel -
          nivel_ce
      )
    )
  
  list(
    diferencia_maxima =
      diferencia,
    coincide =
      diferencia < 1e-12
  )
}


# ============================================================
# VERIFICACIÓN HOLT
# ============================================================

verificar_holt <- function(
    y,
    alpha,
    beta
) {
  
  .validar_y(y)
  
  ajuste <-
    ajustar_holt(
      y,
      alpha,
      beta
    )
  
  n <-
    length(y)
  
  nivel_ce <-
    numeric(n)
  
  tendencia_ce <-
    numeric(n)
  
  nivel_ce[1] <-
    y[1]
  
  tendencia_ce[1] <-
    0
  
  for (t in 2:n) {
    
    error_t <-
      y[t] -
      (
        nivel_ce[t - 1] +
          tendencia_ce[t - 1]
      )
    
    nivel_ce[t] <-
      nivel_ce[t - 1] +
      tendencia_ce[t - 1] +
      alpha *
      error_t
    
    tendencia_ce[t] <-
      tendencia_ce[t - 1] +
      alpha *
      beta *
      error_t
  }
  
  diferencia_nivel <-
    max(
      abs(
        ajuste$parametros$nivel -
          nivel_ce
      )
    )
  
  diferencia_tendencia <-
    max(
      abs(
        ajuste$parametros$tendencia -
          tendencia_ce
      )
    )
  
  list(
    diferencia_nivel =
      diferencia_nivel,
    diferencia_tendencia =
      diferencia_tendencia,
    coincide =
      diferencia_nivel < 1e-12 &&
      diferencia_tendencia < 1e-12
  )
}


# ============================================================
# MÉTRICAS DE ESTIMACIÓN
# ============================================================

medidas_estimacion <- function(
    y,
    ajuste,
    escala_mase = NULL
) {
  
  validos <-
    is.finite(
      ajuste$yhat
    )
  
  if (!any(validos)) {
    
    stop(
      "No existen pronósticos de un paso válidos."
    )
  }
  
  medidas(
    y[validos],
    ajuste$yhat[validos],
    escala_mase
  )
}


# ============================================================
# COMPARACIÓN CON INGENUO
# ============================================================

comparar_con_ingenuo <- function(
    y_estimacion,
    y_validacion,
    pronostico_metodo,
    frecuencia = NULL
) {
  
  .validar_y(
    y_estimacion
  )
  
  .validar_y(
    y_validacion
  )
  
  h <-
    length(
      y_validacion
    )
  
  
  if (
    !is.null(frecuencia) &&
    is.finite(frecuencia) &&
    frecuencia > 1
  ) {
    
    frecuencia <-
      as.integer(
        frecuencia
      )
    
    pronostico_referencia <-
      pronostico_ingenuo_estacional(
        y_estimacion,
        frecuencia,
        h
      )
    
    escala <-
      escala_ingenuo_estacional(
        y_estimacion,
        frecuencia
      )
    
    referencia <-
      "ingenuo estacional"
    
  } else {
    
    pronostico_referencia <-
      pronostico_ingenuo(
        y_estimacion,
        h
      )
    
    escala <-
      escala_ingenuo(
        y_estimacion
      )
    
    referencia <-
      "ingenuo"
  }
  
  
  metodo <-
    medidas(
      y_validacion,
      pronostico_metodo,
      escala
    )
  
  ingenuo <-
    medidas(
      y_validacion,
      pronostico_referencia,
      escala
    )
  
  
  list(
    referencia = referencia,
    metodo = metodo,
    ingenuo = ingenuo,
    escala_mase = escala,
    pronostico_referencia =
      pronostico_referencia
  )
}


# ============================================================
# FIN DE 03-EVALUACION.R
# ============================================================



# ============================================================
