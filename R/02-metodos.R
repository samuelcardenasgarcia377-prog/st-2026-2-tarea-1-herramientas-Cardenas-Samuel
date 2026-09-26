# R/02-metodos.R
# ============================================================

# ============================================================
# 02-metodos.R
# Ocho métodos de pronóstico implementados manualmente
# ============================================================


# ============================================================
# VALIDACIONES AUXILIARES
# ============================================================

.validar_y <- function(y) {
  
  stopifnot(is.numeric(y))
  
  if (!is.numeric(y)) {
    stop("y debe ser un vector numérico.")
  }
  
  if (length(y) < 2) {
    stop("y debe contener al menos dos observaciones.")
  }
  
  if (anyNA(y) || any(!is.finite(y))) {
    stop("y debe contener únicamente valores finitos.")
  }
  
  invisible(TRUE)
}


.validar_k <- function(k, T, minimo = 2) {
  
  stopifnot(
    is.numeric(k),
    length(k) == 1,
    is.finite(k),
    k == as.integer(k)
  )
  
  k <- as.integer(k)
  
  if (k < minimo) {
    stop(
      paste0(
        "k debe ser un entero mayor o igual que ",
        minimo,
        "."
      )
    )
  }
  
  if (k >= T) {
    stop(
      "k debe ser menor que la longitud de la serie."
    )
  }
  
  k
}


.validar_constante <- function(x, nombre) {
  
  stopifnot(
    is.numeric(x),
    length(x) == 1,
    is.finite(x)
  )
  
  if (x <= 0 || x >= 1) {
    stop(
      paste0(
        nombre,
        " debe estar estrictamente entre 0 y 1."
      )
    )
  }
  
  invisible(TRUE)
}


.validar_h <- function(h) {
  
  stopifnot(
    is.numeric(h),
    length(h) == 1,
    is.finite(h),
    h >= 1,
    h == as.integer(h)
  )
  
  as.integer(h)
}


# ============================================================
# 1. MEDIA HISTÓRICA
# ============================================================

ajustar_media <- function(y) {
  
  .validar_y(y)
  
  T <- length(y)
  
  yhat <- rep(
    NA_real_,
    T
  )
  
  # Yhat_2 = Y_1
  yhat[2] <- y[1]
  
  # Pronósticos históricos mediante media acumulada
  if (T >= 3) {
    
    for (t in 3:T) {
      
      yhat[t] <-
        mean(
          y[1:(t - 1)]
        )
    }
  }
  
  nivel_final <- mean(y)
  
  pronosticar <- function(h) {
    
    h <- .validar_h(h)
    
    rep(
      nivel_final,
      h
    )
  }
  
  parametros <- list(
    metodo = "media",
    nivel_final = nivel_final,
    
    # La media histórica no estima un parámetro
    # de ajuste libre mediante rejilla.
    p_modelo = 1L
  )
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = parametros
  )
}


# ============================================================
# 2. MEDIA MÓVIL
# ============================================================

ajustar_mm <- function(y, k) {
  
  .validar_y(y)
  
  T <- length(y)
  
  k <- .validar_k(
    k,
    T
  )
  
  yhat <- rep(
    NA_real_,
    T
  )
  
  # Primer pronóstico disponible:
  # media de las k observaciones anteriores.
  for (t in (k + 1):T) {
    
    yhat[t] <-
      mean(
        y[(t - k):(t - 1)]
      )
  }
  
  nivel_final <-
    mean(
      tail(y, k)
    )
  
  pronosticar <- function(h) {
    
    h <- .validar_h(h)
    
    rep(
      nivel_final,
      h
    )
  }
  
  parametros <- list(
    metodo = "mm",
    k = k,
    nivel_final = nivel_final,
    
    # k es seleccionado mediante la rejilla de optimización.
    p_modelo = 1L
  )
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = parametros
  )
}


# ============================================================
# 3. SUAVIZAMIENTO EXPONENCIAL SIMPLE
# ============================================================

ajustar_ses <- function(
    y,
    alpha
) {
  
  .validar_y(y)
  
  .validar_constante(
    alpha,
    "alpha"
  )
  
  T <- length(y)
  
  yhat <- rep(
    NA_real_,
    T
  )
  
  nivel <- rep(
    NA_real_,
    T
  )
  
  # Inicialización:
  # L_1 = Y_1
  nivel[1] <- y[1]
  
  # Yhat_2 = Y_1
  yhat[2] <- nivel[1]
  
  for (t in 2:T) {
    
    error_t <-
      y[t] -
      nivel[t - 1]
    
    # Forma de corrección del error:
    #
    # L_t = L_{t-1} + alpha e_t
    nivel[t] <-
      nivel[t - 1] +
      alpha * error_t
    
    if (t < T) {
      
      yhat[t + 1] <-
        nivel[t]
    }
  }
  
  nivel_final <- nivel[T]
  
  pronosticar <- function(h) {
    
    h <- .validar_h(h)
    
    rep(
      nivel_final,
      h
    )
  }
  
  parametros <- list(
    metodo = "ses",
    alpha = alpha,
    nivel = nivel,
    nivel_final = nivel_final,
    p_modelo = 1L
  )
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = parametros
  )
}


# ============================================================
# 4. DOBLE MEDIA MÓVIL
# ============================================================

ajustar_dmm <- function(
    y,
    k
) {
  
  .validar_y(y)
  
  T <- length(y)
  
  k <- .validar_k(
    k,
    T
  )
  
  # Para calcular la DMM final se necesitan 2k - 1
  # observaciones. Para obtener además al menos un
  # pronóstico de un paso dentro de la muestra se requiere
  # T >= 2k.
  if (T < 2 * k) {
    
    stop(
      paste0(
        "Para DMM se requieren al menos 2k observaciones ",
        "para obtener al menos un error de un paso."
      )
    )
  }
  
  MM <- rep(
    NA_real_,
    T
  )
  
  DMM <- rep(
    NA_real_,
    T
  )
  
  Ehat <- rep(
    NA_real_,
    T
  )
  
  beta1 <- rep(
    NA_real_,
    T
  )
  
  yhat <- rep(
    NA_real_,
    T
  )
  
  
  # ----------------------------------------------------------
  # Primera media móvil
  # ----------------------------------------------------------
  
  for (t in k:T) {
    
    MM[t] <-
      mean(
        y[(t - k + 1):t]
      )
  }
  
  
  # ----------------------------------------------------------
  # Segunda media móvil
  # ----------------------------------------------------------
  
  inicio_dmm <- 2 * k - 1
  
  for (t in inicio_dmm:T) {
    
    DMM[t] <-
      mean(
        MM[(t - k + 1):t]
      )
    
    Ehat[t] <-
      2 * MM[t] -
      DMM[t]
    
    beta1[t] <-
      2 / (k - 1) *
      (
        MM[t] -
          DMM[t]
      )
    
    # Pronóstico de un paso:
    #
    # Yhat_{t+1} = Ehat_t + beta1_t
    if (t < T) {
      
      yhat[t + 1] <-
        Ehat[t] +
        beta1[t]
    }
  }
  
  E_final <- Ehat[T]
  
  beta_final <- beta1[T]
  
  
  pronosticar <- function(h) {
    
    h <- .validar_h(h)
    
    E_final +
      beta_final *
      seq_len(h)
  }
  
  
  parametros <- list(
    metodo = "dmm",
    k = k,
    media_movil = MM,
    doble_media_movil = DMM,
    nivel = Ehat,
    tendencia = beta1,
    nivel_final = E_final,
    tendencia_final = beta_final,
    
    # El nivel y la pendiente se estiman en la DMM.
    p_modelo = 2L
  )
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = parametros
  )
}


# ============================================================
# MATRIZ DE DISEÑO
# ============================================================

.construir_X <- function(
    t,
    tipo
) {
  
  tipo <- match.arg(
    tipo,
    c(
      "lineal",
      "cuadratica",
      "exponencial"
    )
  )
  
  if (tipo == "lineal") {
    
    X <- cbind(
      intercepto = 1,
      t = t
    )
    
  } else if (tipo == "cuadratica") {
    
    X <- cbind(
      intercepto = 1,
      t = t,
      t2 = t^2
    )
    
  } else {
    
    # Para la tendencia exponencial:
    #
    # ln(Y_t) = a + theta t + e_t
    #
    # La matriz de diseño es lineal en los parámetros
    # transformados.
    X <- cbind(
      intercepto = 1,
      t = t
    )
  }
  
  X
}


# ============================================================
# HAC / BARTLETT
#
# Var(beta_hat) =
# (X'X)^(-1) S (X'X)^(-1)
#
# donde S es la matriz HAC de Bartlett.
# ============================================================

.hac_bartlett <- function(
    X,
    e
) {
  
  T <- nrow(X)
  
  q <- floor(
    4 * (T / 100)^(2 / 9)
  )
  
  q <- min(
    q,
    T - 1
  )
  
  XtX <- crossprod(X)
  
  XtX_inv <- solve(
    XtX
  )
  
  S <- matrix(
    0,
    nrow = ncol(X),
    ncol = ncol(X)
  )
  
  
  # ----------------------------------------------------------
  # Gamma_0
  # ----------------------------------------------------------
  
  for (t in seq_len(T)) {
    
    xt <- X[t, , drop = FALSE]
    
    S <- S +
      as.numeric(
        e[t]^2
      ) *
      crossprod(
        xt
      )
  }
  
  
  # ----------------------------------------------------------
  # Gamma_h + Gamma_h'
  # ----------------------------------------------------------
  
  if (q >= 1) {
    
    for (h in seq_len(q)) {
      
      peso <-
        1 -
        h / (q + 1)
      
      for (t in (h + 1):T) {
        
        xt <-
          X[t, , drop = FALSE]
        
        xlag <-
          X[t - h, , drop = FALSE]
        
        producto <-
          e[t] *
          e[t - h]
        
        S <- S +
          peso *
          producto *
          (
            crossprod(
              xt,
              xlag
            ) +
              crossprod(
                xlag,
                xt
              )
          )
      }
    }
  }
  
  
  V_HAC <-
    XtX_inv %*%
    S %*%
    XtX_inv
  
  se_HAC <- sqrt(
    pmax(
      diag(V_HAC),
      0
    )
  )
  
  list(
    V_HAC = V_HAC,
    se_HAC = se_HAC,
    lag_HAC = q
  )
}


# ============================================================
# TABLA DE COEFICIENTES
# ============================================================

.tabla_tendencia <- function(
    X,
    beta_hat,
    e
) {
  
  T <- nrow(X)
  
  p <- ncol(X)
  
  gl <- T - p
  
  if (gl <= 0) {
    stop(
      "No hay grados de libertad suficientes para estimar la varianza."
    )
  }
  
  XtX_inv <- solve(
    crossprod(X)
  )
  
  sigma2 <-
    sum(e^2) /
    gl
  
  se_OLS <- sqrt(
    pmax(
      diag(
        sigma2 * XtX_inv
      ),
      0
    )
  )
  
  hac <- .hac_bartlett(
    X,
    e
  )
  
  se_HAC <- hac$se_HAC
  
  t_HAC <-
    beta_hat /
    se_HAC
  
  p_HAC <-
    2 *
    stats::pt(
      -abs(t_HAC),
      df = gl
    )
  
  tabla <- tibble::tibble(
    termino = colnames(X),
    estimacion = as.numeric(beta_hat),
    SE_OLS = se_OLS,
    SE_HAC_Bartlett = se_HAC,
    t_robusto = t_HAC,
    p_robusto = p_HAC
  )
  
  
  # ----------------------------------------------------------
  # R2 sobre la variable realmente modelada.
  #
  # En lineal/cuadrática:
  #     variable modelada = y
  #
  # En exponencial:
  #     variable modelada = log(y)
  # ----------------------------------------------------------
  
  y_modelo <- as.numeric(
    X %*%
      beta_hat +
      e
  )
  
  SSE <- sum(
    e^2
  )
  
  SST <- sum(
    (
      y_modelo -
        mean(y_modelo)
    )^2
  )
  
  R2 <-
    if (SST > 0) {
      
      1 -
        SSE / SST
      
    } else {
      
      NA_real_
    }
  
  
  # ----------------------------------------------------------
  # Durbin-Watson
  # ----------------------------------------------------------
  
  DW <-
    sum(
      diff(e)^2
    ) /
    sum(
      e^2
    )
  
  list(
    tabla = tabla,
    sigma2 = sigma2,
    R2 = R2,
    DW = DW,
    lag_HAC = hac$lag_HAC
  )
}


# ============================================================
# 5, 6 Y 7. TENDENCIAS
# ============================================================

ajustar_tendencia <- function(
    y,
    tipo = c(
      "lineal",
      "cuadratica",
      "exponencial"
    ),
    corregir_sesgo = FALSE
) {
  
  .validar_y(y)
  
  tipo <- match.arg(
    tipo
  )
  
  stopifnot(
    is.logical(corregir_sesgo),
    length(corregir_sesgo) == 1
  )
  
  T <- length(y)
  
  t <- seq_len(T)
  
  
  # ----------------------------------------------------------
  # Transformación para tendencia exponencial
  # ----------------------------------------------------------
  
  if (tipo == "exponencial") {
    
    if (any(y <= 0)) {
      
      stop(
        "La tendencia exponencial requiere y > 0."
      )
    }
    
    y_modelo <- log(y)
    
  } else {
    
    y_modelo <- y
  }
  
  
  # ----------------------------------------------------------
  # Matriz de diseño
  # ----------------------------------------------------------
  
  X <- .construir_X(
    t,
    tipo
  )
  
  
  # ----------------------------------------------------------
  # Ecuaciones normales
  #
  # beta_hat = (X'X)^(-1) X'y
  # ----------------------------------------------------------
  
  beta_hat <- solve(
    crossprod(X),
    crossprod(
      X,
      y_modelo
    )
  )
  
  
  # ----------------------------------------------------------
  # Ajuste en escala del modelo
  # ----------------------------------------------------------
  
  ajuste_modelo <- as.numeric(
    X %*%
      beta_hat
  )
  
  e <- y_modelo -
    ajuste_modelo
  
  
  # ----------------------------------------------------------
  # Tabla estadística
  # ----------------------------------------------------------
  
  info <- .tabla_tendencia(
    X,
    beta_hat,
    e
  )
  
  
  # ----------------------------------------------------------
  # Tendencia lineal
  # ----------------------------------------------------------
  
  if (tipo == "lineal") {
    
    p <- 2L
    
    yhat <- ajuste_modelo
    
    pronosticar <- function(h) {
      
      h <- .validar_h(h)
      
      Xf <- .construir_X(
        T + seq_len(h),
        "lineal"
      )
      
      as.numeric(
        Xf %*%
          beta_hat
      )
    }
    
    coeficientes <- beta_hat
  }
  
  
  # ----------------------------------------------------------
  # Tendencia cuadrática
  # ----------------------------------------------------------
  
  else if (tipo == "cuadratica") {
    
    p <- 3L
    
    yhat <- ajuste_modelo
    
    pronosticar <- function(h) {
      
      h <- .validar_h(h)
      
      Xf <- .construir_X(
        T + seq_len(h),
        "cuadratica"
      )
      
      as.numeric(
        Xf %*%
          beta_hat
      )
    }
    
    coeficientes <- beta_hat
  }
  
  
  # ----------------------------------------------------------
  # Tendencia exponencial
  # ----------------------------------------------------------
  
  else {
    
    p <- 2L
    
    a_hat <- beta_hat[1]
    
    theta_hat <- beta_hat[2]
    
    beta0 <- exp(
      a_hat
    )
    
    beta1 <- exp(
      theta_hat
    )
    
    sigma2_log <-
      info$sigma2
    
    # --------------------------------------------------------
    # Bajo errores normales en escala logarítmica:
    #
    # E(Y|t) =
    # exp(a + theta t + sigma^2/2)
    #
    # Sin corrección, exp(a + theta t) representa
    # la mediana condicional.
    # --------------------------------------------------------
    
    factor_sesgo <-
      if (corregir_sesgo) {
        
        exp(
          sigma2_log / 2
        )
        
      } else {
        
        1
      }
    
    
    yhat <-
      exp(
        ajuste_modelo
      ) *
      factor_sesgo
    
    
    pronosticar <- function(h) {
      
      h <- .validar_h(h)
      
      tf <-
        T +
        seq_len(h)
      
      factor_sesgo *
        beta0 *
        beta1^tf
    }
    
    
    # Se conservan los parámetros en escala original
    # para facilitar la interpretación.
    coeficientes <- c(
      beta0 = beta0,
      beta1 = beta1
    )
  }
  
  
  # ----------------------------------------------------------
  # Parámetros
  # ----------------------------------------------------------
  
  parametros <- list(
    metodo = tipo,
    tipo = tipo,
    coeficientes = coeficientes,
    tabla_coeficientes =
      info$tabla,
    sigma2 = info$sigma2,
    R2 = info$R2,
    DW = info$DW,
    lag_HAC = info$lag_HAC,
    p_modelo = p,
    corregir_sesgo =
      corregir_sesgo
  )
  
  
  if (tipo == "exponencial") {
    
    parametros$escala_modelo <-
      "logaritmica"
    
    parametros$interpretacion <-
      paste(
        "La regresión linealizada estima la tendencia",
        "en escala logarítmica. Sin corrección, la",
        "transformación inversa corresponde a la mediana",
        "condicional. La corrección exp(sigma2/2) aproxima",
        "la media condicional bajo normalidad logarítmica."
      )
    
    parametros$beta0 <-
      beta0
    
    parametros$beta1 <-
      beta1
    
    parametros$sigma2_log <-
      sigma2_log
    
    parametros$factor_sesgo <-
      factor_sesgo
  }
  
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = parametros
  )
}


# ============================================================
# 8. HOLT
# ============================================================

ajustar_holt <- function(
    y,
    alpha,
    beta
) {
  
  .validar_y(y)
  
  .validar_constante(
    alpha,
    "alpha"
  )
  
  .validar_constante(
    beta,
    "beta"
  )
  
  T <- length(y)
  
  yhat <- rep(
    NA_real_,
    T
  )
  
  nivel <- rep(
    NA_real_,
    T
  )
  
  tendencia <- rep(
    NA_real_,
    T
  )
  
  
  # ----------------------------------------------------------
  # Inicialización
  # ----------------------------------------------------------
  
  nivel[1] <- y[1]
  
  tendencia[1] <- 0
  
  
  # ----------------------------------------------------------
  # Recursión
  #
  # Yhat_t = L_{t-1} + T_{t-1}
  #
  # L_t =
  #   L_{t-1} + T_{t-1} + alpha e_t
  #
  # T_t =
  #   T_{t-1} + alpha beta e_t
  #
  # Esta forma es algebraicamente equivalente a:
  #
  # L_t =
  # alpha Y_t +
  # (1-alpha)Yhat_t
  #
  # T_t =
  # beta(L_t-L_{t-1}) +
  # (1-beta)T_{t-1}
  # ----------------------------------------------------------
  
  for (t in 2:T) {
    
    yhat[t] <-
      nivel[t - 1] +
      tendencia[t - 1]
    
    error_t <-
      y[t] -
      yhat[t]
    
    nivel[t] <-
      nivel[t - 1] +
      tendencia[t - 1] +
      alpha * error_t
    
    tendencia[t] <-
      tendencia[t - 1] +
      alpha * beta * error_t
  }
  
  
  nivel_final <- nivel[T]
  
  tendencia_final <-
    tendencia[T]
  
  
  pronosticar <- function(h) {
    
    h <- .validar_h(h)
    
    nivel_final +
      tendencia_final *
      seq_len(h)
  }
  
  
  parametros <- list(
    metodo = "holt",
    alpha = alpha,
    beta = beta,
    nivel = nivel,
    tendencia = tendencia,
    nivel_final = nivel_final,
    tendencia_final =
      tendencia_final,
    p_modelo = 2L
  )
  
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = parametros
  )
}


# ============================================================
# OPTIMIZACIÓN
# ============================================================

optimizar <- function(
    y,
    metodo,
    rejilla = NULL
) {
  
  .validar_y(y)
  
  metodo <- match.arg(
    metodo,
    c(
      "mm",
      "dmm",
      "ses",
      "holt"
    )
  )
  
  T <- length(y)
  
  
  # ----------------------------------------------------------
  # Rejillas exigidas
  # ----------------------------------------------------------
  
  if (is.null(rejilla)) {
    
    if (metodo == "mm") {
      
      rejilla <- 2:12
      
    } else if (metodo == "dmm") {
      
      rejilla <- 2:12
      
    } else if (metodo == "ses") {
      
      rejilla <- seq(
        0.02,
        0.98,
        by = 0.02
      )
      
    } else {
      
      valores <- seq(
        0.05,
        0.95,
        by = 0.05
      )
      
      rejilla <- expand.grid(
        alpha = valores,
        beta = valores
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Validación MM / DMM
  # ----------------------------------------------------------
  
  if (
    metodo %in%
    c(
      "mm",
      "dmm"
    )
  ) {
    
    rejilla <- as.integer(
      rejilla
    )
    
    if (
      any(
        rejilla < 2
      ) ||
      any(
        rejilla > 12
      )
    ) {
      
      stop(
        "Para MM y DMM la rejilla debe estar entre 2 y 12."
      )
    }
    
    if (
      metodo == "mm" &&
      any(
        rejilla >= T
      )
    ) {
      
      stop(
        "Algún k de MM no es válido para la longitud de y."
      )
    }
    
    if (
      metodo == "dmm" &&
      any(
        2 * rejilla > T
      )
    ) {
      
      stop(
        "Algún k de DMM no permite obtener errores de un paso."
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Validación SES
  # ----------------------------------------------------------
  
  if (metodo == "ses") {
    
    if (
      any(
        rejilla <= 0 |
        rejilla >= 1
      )
    ) {
      
      stop(
        "La rejilla de SES debe estar estrictamente entre 0 y 1."
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Validación Holt
  # ----------------------------------------------------------
  
  if (metodo == "holt") {
    
    if (
      !all(
        c(
          "alpha",
          "beta"
        ) %in%
        names(rejilla)
      )
    ) {
      
      stop(
        "La rejilla de Holt debe contener alpha y beta."
      )
    }
    
    if (
      any(
        rejilla$alpha <= 0 |
        rejilla$alpha >= 1 |
        rejilla$beta <= 0 |
        rejilla$beta >= 1
      )
    ) {
      
      stop(
        "Los valores alpha y beta deben estar entre 0 y 1."
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Función interna de evaluación
  #
  # El MSE se calcula únicamente con los pronósticos
  # de un paso que existen dentro del segmento de
  # estimación.
  # ----------------------------------------------------------
  
  evaluar <- function(parametros) {
    
    if (metodo == "mm") {
      
      ajuste <-
        ajustar_mm(
          y,
          k = parametros
        )
      
    } else if (metodo == "dmm") {
      
      ajuste <-
        ajustar_dmm(
          y,
          k = parametros
        )
      
    } else if (metodo == "ses") {
      
      ajuste <-
        ajustar_ses(
          y,
          alpha = parametros
        )
      
    } else {
      
      ajuste <-
        ajustar_holt(
          y,
          alpha = parametros[1],
          beta = parametros[2]
        )
    }
    
    validos <-
      is.finite(
        ajuste$yhat
      )
    
    if (!any(validos)) {
      
      return(
        NA_real_
      )
    }
    
    mean(
      (
        y[validos] -
          ajuste$yhat[validos]
      )^2
    )
  }
  
  
  # ----------------------------------------------------------
  # MM / DMM / SES
  # ----------------------------------------------------------
  
  if (
    metodo %in%
    c(
      "mm",
      "dmm",
      "ses"
    )
  ) {
    
    MSE <- vapply(
      rejilla,
      evaluar,
      numeric(1)
    )
    
    resultados <-
      tibble::tibble(
        parametro = rejilla,
        MSE = MSE
      )
    
    if (
      all(
        !is.finite(MSE)
      )
    ) {
      
      stop(
        "No fue posible calcular MSE para ningún valor de la rejilla."
      )
    }
    
    indice <-
      which.min(
        MSE
      )
    
    optimo <-
      resultados[
        indice,
        ,
        drop = FALSE
      ]
  }
  
  
  # ----------------------------------------------------------
  # Holt
  # ----------------------------------------------------------
  
  else {
    
    MSE <- vapply(
      seq_len(
        nrow(rejilla)
      ),
      function(i) {
        
        evaluar(
          c(
            rejilla$alpha[i],
            rejilla$beta[i]
          )
        )
      },
      numeric(1)
    )
    
    resultados <-
      tibble::tibble(
        alpha = rejilla$alpha,
        beta = rejilla$beta,
        MSE = MSE
      )
    
    if (
      all(
        !is.finite(MSE)
      )
    ) {
      
      stop(
        "No fue posible calcular MSE para ningún punto de la rejilla."
      )
    }
    
    indice <-
      which.min(
        MSE
      )
    
    optimo <-
      resultados[
        indice,
        ,
        drop = FALSE
      ]
  }
  
  
  en_borde <-
    if (metodo == "holt") {
      optimo$alpha == min(rejilla$alpha) ||
        optimo$alpha == max(rejilla$alpha) ||
        optimo$beta == min(rejilla$beta) ||
        optimo$beta == max(rejilla$beta)
    } else {
      optimo$parametro == min(rejilla) ||
        optimo$parametro == max(rejilla)
    }
  
  list(
    metodo = metodo,
    resultados = resultados,
    optimo = optimo,
    indice_optimo = indice,
    en_borde = en_borde
  )
}


# ============================================================
# MÉTODOS DISPONIBLES
# ============================================================

metodos_disponibles <- function() {
  
  c(
    "media",
    "mm",
    "ses",
    "dmm",
    "lineal",
    "cuadratica",
    "exponencial",
    "holt"
  )
}


# ============================================================
# FIN DE 02-METODOS.R
# ============================================================



# ============================================================
