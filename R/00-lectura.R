# R/00-lectura.R
# Herramientas de lectura y estandarización de series
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)


# ============================================================
# Funciones auxiliares internas
# ============================================================

.validar_fechas <- function(fecha) {
  
  if (!inherits(fecha, "Date")) {
    stop("fecha debe ser un vector de clase Date.")
  }
  
  if (length(fecha) == 0) {
    stop("No hay fechas para validar.")
  }
  
  if (anyNA(fecha)) {
    stop("La columna de fechas contiene valores NA.")
  }
  
  if (length(fecha) >= 2) {
    
    if (any(diff(fecha) <= 0)) {
      stop(
        "Las fechas deben estar ordenadas estrictamente de forma creciente."
      )
    }
  }
  
  invisible(TRUE)
}


.inferir_frecuencia_fecha <- function(fecha) {
  
  .validar_fechas(fecha)
  
  if (length(fecha) < 2) {
    stop(
      "Se requieren al menos dos fechas para verificar la frecuencia."
    )
  }
  
  # ----------------------------------------------------------
  # Se trabaja con periodos calendario, no con diferencias
  # en días. Esto evita problemas como:
  #
  # enero  -> febrero = 31 días
  # febrero -> marzo = 28/29 días
  # marzo  -> abril = 31 días
  #
  # Todos ellos siguen siendo meses consecutivos.
  # ----------------------------------------------------------
  
  anio <- as.integer(
    format(fecha, "%Y")
  )
  
  mes <- as.integer(
    format(fecha, "%m")
  )
  
  indice_mensual <- 12 * anio + mes
  
  diferencias_mensuales <- diff(indice_mensual)
  
  if (all(diferencias_mensuales == 1)) {
    return(12)
  }
  
  if (all(diferencias_mensuales == 3)) {
    return(4)
  }
  
  if (all(diferencias_mensuales == 12)) {
    return(1)
  }
  
  stop(
    "Las fechas no están igualmente espaciadas con frecuencia anual, trimestral o mensual."
  )
}


.validar_frecuencia_ts <- function(frecuencia) {
  
  if (!is.numeric(frecuencia) ||
      length(frecuencia) != 1 ||
      !is.finite(frecuencia) ||
      frecuencia <= 0) {
    
    stop(
      "La frecuencia de la serie ts debe ser un número positivo."
    )
  }
  
  frecuencia <- as.integer(frecuencia)
  
  if (!frecuencia %in% c(1L, 4L, 12L)) {
    
    stop(
      paste0(
        "Frecuencia ts no implementada automáticamente: ",
        frecuencia,
        ". Se esperan frecuencias 1, 4 o 12."
      )
    )
  }
  
  frecuencia
}


# ============================================================
# leer_serie()
#
# Acepta:
#   1. un objeto ts
#   2. la ruta a un CSV con columnas:
#        fecha
#        valor
#
# Devuelve:
#   t      -> índice entero 1,...,n
#   fecha  -> clase Date
#   y      -> valores numéricos
#
# Atributos:
#   frecuencia
#   fuente
#   unidad
# ============================================================

leer_serie <- function(
    x,
    fuente,
    unidad
) {
  
  stopifnot(
    is.character(fuente),
    length(fuente) == 1,
    !is.na(fuente),
    is.character(unidad),
    length(unidad) == 1,
    !is.na(unidad)
  )
  
  
  # ==========================================================
  # CASO 1: objeto ts
  # ==========================================================
  
  if (inherits(x, "ts")) {
    
    y <- as.numeric(x)
    
    n <- length(y)
    
    if (n == 0) {
      stop("La serie no contiene observaciones.")
    }
    
    if (anyNA(y) ||
        any(!is.finite(y))) {
      
      stop(
        "La serie contiene valores NA, NaN o infinitos."
      )
    }
    
    frecuencia <- stats::frequency(x)
    
    frecuencia <- .validar_frecuencia_ts(
      frecuencia
    )
    
    inicio <- stats::start(x)
    
    # --------------------------------------------------------
    # Frecuencia anual
    # --------------------------------------------------------
    
    if (frecuencia == 1) {
      
      anio <- inicio[1] +
        seq_len(n) - 1
      
      fecha <- as.Date(
        sprintf(
          "%d-01-01",
          anio
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Frecuencia trimestral
    # --------------------------------------------------------
    
    else if (frecuencia == 4) {
      
      indice <- seq_len(n) - 1
      
      anio <- inicio[1] +
        floor(
          (inicio[2] - 1 + indice) / 4
        )
      
      trimestre <- (
        (inicio[2] - 1 + indice) %% 4
      ) + 1
      
      mes <- 1 +
        3 * (trimestre - 1)
      
      fecha <- as.Date(
        sprintf(
          "%d-%02d-01",
          anio,
          mes
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Frecuencia mensual
    # --------------------------------------------------------
    
    else if (frecuencia == 12) {
      
      indice <- seq_len(n) - 1
      
      anio <- inicio[1] +
        floor(
          (inicio[2] - 1 + indice) / 12
        )
      
      mes <- (
        (inicio[2] - 1 + indice) %% 12
      ) + 1
      
      fecha <- as.Date(
        sprintf(
          "%d-%02d-01",
          anio,
          mes
        )
      )
    }
    
    
    else {
      
      stop(
        paste0(
          "Frecuencia ts no implementada automáticamente: ",
          frecuencia
        )
      )
    }
  }
  
  
  # ==========================================================
  # CASO 2: archivo CSV
  # ==========================================================
  
  else if (
    is.character(x) &&
    length(x) == 1
  ) {
    
    datos <- utils::read.csv(
      x,
      stringsAsFactors = FALSE
    )
    
    if (!all(
      c("fecha", "valor") %in% names(datos)
    )) {
      
      stop(
        "El CSV debe contener las columnas 'fecha' y 'valor'."
      )
    }
    
    fecha <- as.Date(
      datos$fecha
    )
    
    y <- suppressWarnings(
      as.numeric(datos$valor)
    )
    
    if (anyNA(fecha)) {
      
      stop(
        "La columna 'fecha' contiene fechas no válidas."
      )
    }
    
    if (anyNA(y) ||
        any(!is.finite(y))) {
      
      stop(
        paste0(
          "La columna 'valor' contiene valores NA, ",
          "NaN, infinitos o no numéricos."
        )
      )
    }
    
    if (length(fecha) == 0) {
      
      stop(
        "El archivo no contiene observaciones."
      )
    }
    
    frecuencia <- .inferir_frecuencia_fecha(
      fecha
    )
  }
  
  
  # ==========================================================
  # CASO NO VÁLIDO
  # ==========================================================
  
  else {
    
    stop(
      "x debe ser un objeto ts o la ruta a un archivo CSV."
    )
  }
  
  
  # ==========================================================
  # Validaciones comunes
  # ==========================================================
  
  if (length(fecha) != length(y)) {
    
    stop(
      "La cantidad de fechas y observaciones no coincide."
    )
  }
  
  .validar_fechas(fecha)
  
  if (anyNA(y) ||
      any(!is.finite(y))) {
    
    stop(
      "La serie contiene valores no finitos."
    )
  }
  
  
  # ==========================================================
  # Construcción del objeto estandarizado
  # ==========================================================
  
  datos_salida <- tibble::tibble(
    
    t = seq_along(y),
    
    fecha = fecha,
    
    y = y
  )
  
  
  # ==========================================================
  # Atributos exigidos por la tarea
  # ==========================================================
  
  attr(
    datos_salida,
    "frecuencia"
  ) <- as.numeric(frecuencia)
  
  attr(
    datos_salida,
    "fuente"
  ) <- fuente
  
  attr(
    datos_salida,
    "unidad"
  ) <- unidad
  
  
  datos_salida
}


# ============================================================
