# ejemplos/ejemplos.R
# ============================================================

# ============================================================
# ejemplos.R
# Ejecución reproducible de los ocho métodos
# ============================================================


library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)


# ============================================================
# CARGA DE FUNCIONES
# ============================================================

source(
  "R/00-lectura.R"
)

source(
  "R/01-graficos.R"
)

source(
  "R/02-metodos.R"
)

source(
  "R/03-evaluacion.R"
)


# ============================================================
# CARPETA DE SALIDA
# ============================================================

if (!dir.exists("figs")) {
  
  dir.create(
    "figs",
    recursive = TRUE
  )
}


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

calcular_h <- function(
    datos,
    estacional = FALSE
) {
  
  T <-
    nrow(datos)
  
  if (T < 3) {
    
    stop(
      "Se requieren al menos tres observaciones."
    )
  }
  
  h <-
    min(
      12,
      floor(0.2 * T)
    )
  
  h <-
    max(
      h,
      1
    )
  
  frecuencia <-
    attr(
      datos,
      "frecuencia"
    )
  
  # Para series estacionales se exige al menos
  # un ciclo completo de validación.
  if (
    estacional &&
    is.finite(frecuencia) &&
    frecuencia > 1
  ) {
    
    h <-
      max(
        h,
        as.integer(frecuencia)
      )
  }
  
  max_h <- T - 2
  
  if (estacional &&
      is.finite(frecuencia) &&
      frecuencia > 1 &&
      frecuencia > max_h) {
    stop(
      "La serie no tiene suficientes observaciones para reservar un ciclo completo de validación."
    )
  }
  
  h <- min(h, max_h)
  
  if (h < 1) {
    stop(
      "La serie no tiene suficientes observaciones para separar la validación."
    )
  }
  
  as.integer(h)
}


guardar_figura <- function(
    grafico,
    nombre,
    ancho = 8,
    alto = 5
) {
  
  ggplot2::ggsave(
    filename =
      file.path(
        "figs",
        nombre
      ),
    plot = grafico,
    width = ancho,
    height = alto,
    units = "in"
  )
}


tabla_autocorrelaciones <- function(
    datos
) {
  
  m <-
    max(
      1,
      min(
        floor(
          nrow(datos) / 4
        ),
        24
      )
    )
  
  autocorrelaciones_individuales(
    datos$y,
    m = m
  )
}


# ============================================================
# GRÁFICO DE OPTIMIZACIÓN
# ============================================================

graficar_optimizacion <- function(
    optimizacion,
    nombre
) {
  
  tabla <-
    optimizacion$resultados
  
  metodo <-
    optimizacion$metodo
  
  
  if (metodo == "holt") {
    
    grafico <-
      ggplot2::ggplot(
        tabla,
        ggplot2::aes(
          x = alpha,
          y = beta,
          fill = MSE
        )
      ) +
      
      ggplot2::geom_tile() +
      
      ggplot2::geom_point(
        data =
          optimizacion$optimo,
        ggplot2::aes(
          x = alpha,
          y = beta
        ),
        inherit.aes = FALSE,
        size = 3
      ) +
      
      ggplot2::labs(
        title =
          "Optimización de Holt mediante MSE",
        x =
          expression(alpha),
        y =
          expression(beta),
        fill = "MSE"
      ) +
      
      ggplot2::theme_minimal()
    
  } else {
    
    etiqueta_x <-
      if (
        metodo %in%
        c(
          "mm",
          "dmm"
        )
      ) {
        
        "k"
        
      } else {
        
        "Parámetro"
      }
    
    
    grafico <-
      ggplot2::ggplot(
        tabla,
        ggplot2::aes(
          x = parametro,
          y = MSE
        )
      ) +
      
      ggplot2::geom_line() +
      
      ggplot2::geom_point() +
      
      ggplot2::geom_point(
        data =
          optimizacion$optimo,
        ggplot2::aes(
          x = parametro,
          y = MSE
        ),
        inherit.aes = FALSE,
        size = 3
      ) +
      
      ggplot2::labs(
        title =
          paste(
            "Optimización:",
            metodo
          ),
        x = etiqueta_x,
        y = "MSE"
      ) +
      
      ggplot2::theme_minimal()
  }
  
  
  guardar_figura(
    grafico,
    nombre
  )
  
  grafico
}


# ============================================================
# AJUSTAR CON EL ÓPTIMO
# ============================================================

ajustar_con_optimo <- function(
    y,
    metodo,
    optimizacion
) {
  
  optimo <-
    optimizacion$optimo
  
  
  if (metodo == "mm") {
    
    return(
      ajustar_mm(
        y,
        k =
          as.integer(
            optimo$parametro
          )
      )
    )
  }
  
  
  if (metodo == "dmm") {
    
    return(
      ajustar_dmm(
        y,
        k =
          as.integer(
            optimo$parametro
          )
      )
    )
  }
  
  
  if (metodo == "ses") {
    
    return(
      ajustar_ses(
        y,
        alpha =
          optimo$parametro
      )
    )
  }
  
  
  if (metodo == "holt") {
    
    return(
      ajustar_holt(
        y,
        alpha =
          optimo$alpha,
        beta =
          optimo$beta
      )
    )
  }
  
  
  stop(
    "Método no reconocido para optimización."
  )
}


# ============================================================
# EJECUCIÓN DE UN EJEMPLO
# ============================================================

ejecutar_ejemplo <- function(
    nombre,
    datos,
    metodo,
    justificacion,
    estacional = FALSE,
    corregir_sesgo = FALSE,
    archivo = NULL
) {
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "Ejemplo:",
    nombre,
    "\n"
  )
  
  cat(
    "Método:",
    metodo,
    "\n"
  )
  
  cat(
    "Justificación:",
    justificacion,
    "\n"
  )
  
  cat(
    "====================================================\n"
  )
  
  
  y <-
    datos$y
  
  T <-
    length(y)
  
  frecuencia <-
    attr(
      datos,
      "frecuencia"
    )
  
  
  # ==========================================================
  # HORIZONTE Y PARTICIÓN
  # ==========================================================
  
  h <-
    calcular_h(
      datos,
      estacional
    )
  
  T_est <-
    T - h
  
  y_estimacion <-
    y[
      seq_len(T_est)
    ]
  
  y_validacion <-
    y[
      (T_est + 1):T
    ]
  
  
  # ==========================================================
  # DESCRIPCIÓN
  # ==========================================================
  
  descripcion <- list(
    
    serie = nombre,
    
    fuente =
      attr(
        datos,
        "fuente"
      ),
    
    unidad =
      attr(
        datos,
        "unidad"
      ),
    
    n = T,
    
    frecuencia =
      frecuencia,
    
    fecha_inicio =
      min(
        datos$fecha
      ),
    
    fecha_final =
      max(
        datos$fecha
      ),
    
    h_validacion = h,
    
    observaciones_estimacion =
      T_est,
    
    justificacion =
      justificacion
  )
  
  
  # ==========================================================
  # GRÁFICO DE LA SERIE
  # ==========================================================
  
  grafico_serie <-
    graficar_serie(
      datos,
      paste(
        nombre,
        "- Serie"
      )
    )
  
  if (!is.null(archivo)) {
    
    guardar_figura(
      grafico_serie,
      paste0(
        archivo,
        "_serie.png"
      )
    )
  }
  
  
  # ==========================================================
  # CORRELOGRAMA
  # ==========================================================
  
  grafico_corr <-
    correlograma(
      datos
    )
  
  if (!is.null(archivo)) {
    
    guardar_figura(
      grafico_corr,
      paste0(
        archivo,
        "_correlograma.png"
      ),
      ancho = 8,
      alto = 7
    )
  }
  
  
  # ==========================================================
  # ACF INDIVIDUAL Y LJUNG-BOX DE LA SERIE
  # ==========================================================
  
  acf_serie <-
    tabla_autocorrelaciones(
      datos
    )
  
  m_serie <-
    nrow(
      acf_serie
    )
  
  lb_serie <-
    ljung_box(
      acf_serie$acf,
      T = T,
      m = m_serie,
      p = 0
    )
  
  
  # ==========================================================
  # AJUSTE
  # ==========================================================
  
  optimizacion <- NULL
  
  
  if (
    metodo %in%
    c(
      "mm",
      "dmm",
      "ses",
      "holt"
    )
  ) {
    
    optimizacion <-
      optimizar(
        y_estimacion,
        metodo
      )
    
    ajuste <-
      ajustar_con_optimo(
        y_estimacion,
        metodo,
        optimizacion
      )
    
    
    if (!is.null(archivo)) {
      
      graficar_optimizacion(
        optimizacion,
        paste0(
          archivo,
          "_optimizacion.png"
        )
      )
    }
    
  }
  
  
  else if (
    metodo == "media"
  ) {
    
    ajuste <-
      ajustar_media(
        y_estimacion
      )
    
  }
  
  
  else if (
    metodo %in%
    c(
      "lineal",
      "cuadratica",
      "exponencial"
    )
  ) {
    
    ajuste <-
      ajustar_tendencia(
        y_estimacion,
        tipo = metodo,
        corregir_sesgo =
          corregir_sesgo
      )
    
  }
  
  
  else {
    
    stop(
      paste(
        "Método no reconocido:",
        metodo
      )
    )
  }
  
  
  # ==========================================================
  # ESCALA MASE
  #
  # Si la serie se trata como estacional, se utiliza como
  # referencia el ingenuo estacional.
  # ==========================================================
  
  if (
    estacional &&
    is.finite(frecuencia) &&
    frecuencia > 1
  ) {
    
    escala <-
      escala_ingenuo_estacional(
        y_estimacion,
        frecuencia
      )
    
  } else {
    
    escala <-
      escala_ingenuo(
        y_estimacion
      )
  }
  
  
  # ==========================================================
  # MÉTRICAS DE ESTIMACIÓN
  # ==========================================================
  
  indicadores_estimacion <-
    medidas_estimacion(
      y_estimacion,
      ajuste,
      escala_mase =
        escala
    )
  
  
  # ==========================================================
  # PRONÓSTICO DE VALIDACIÓN
  # ==========================================================
  
  pronostico <-
    ajuste$pronosticar(
      h
    )
  
  if (
    length(pronostico) != h ||
    anyNA(pronostico) ||
    any(!is.finite(pronostico))
  ) {
    
    stop(
      paste(
        "El método",
        metodo,
        "no produjo un pronóstico válido."
      )
    )
  }
  
  
  indicadores_validacion <-
    medidas(
      y_validacion,
      pronostico,
      escala_mase =
        escala
    )
  
  
  # ==========================================================
  # COMPARACIÓN CON INGENUO
  # ==========================================================
  
  comparacion <-
    comparar_con_ingenuo(
      y_estimacion,
      y_validacion,
      pronostico,
      frecuencia =
        if (
          estacional &&
          is.finite(frecuencia) &&
          frecuencia > 1
        ) {
          
          frecuencia
          
        } else {
          
          NULL
        }
    )
  
  
  # ==========================================================
  # RESIDUOS DE ESTIMACIÓN
  # ==========================================================
  
  indices_validos <-
    is.finite(
      ajuste$yhat
    )
  
  if (!any(indices_validos)) {
    
    stop(
      "No hay errores de un paso disponibles."
    )
  }
  
  residuos <-
    y_estimacion[
      indices_validos
    ] -
    ajuste$yhat[
      indices_validos
    ]
  
  ajustados <-
    ajuste$yhat[
      indices_validos
    ]
  
  
  # ==========================================================
  # NÚMERO DE PARÁMETROS
  # ==========================================================
  
  p_modelo <-
    as.integer(
      ajuste$parametros$p_modelo
    )
  
  
  # ==========================================================
  # VALIDACIÓN DE RESIDUOS
  # ==========================================================
  
  m_residuos <-
    max(
      1,
      min(
        floor(
          length(residuos) / 4
        ),
        24
      )
    )
  
  # Garantizar grados de libertad positivos en Ljung-Box.
  m_residuos <-
    max(
      m_residuos,
      p_modelo + 1
    )
  
  m_residuos <-
    min(
      m_residuos,
      length(residuos) - 1
    )
  
  if (
    m_residuos -
    p_modelo <= 0
  ) {
    
    stop(
      "No hay suficientes rezagos para aplicar Ljung-Box con el p_modelo declarado."
    )
  }
  
  
  validacion_residuos <-
    validar_errores(
      residuos,
      p = p_modelo,
      m = m_residuos,
      titulo = nombre,
      ajustados = ajustados
    )
  
  
  # ==========================================================
  # GRÁFICOS DE RESIDUOS
  # ==========================================================
  
  if (!is.null(archivo)) {
    
    guardar_figura(
      validacion_residuos$grafico_acf,
      paste0(
        archivo,
        "_acf_residuos.png"
      )
    )
    
    guardar_figura(
      validacion_residuos$grafico_errores_tiempo,
      paste0(
        archivo,
        "_errores_tiempo.png"
      )
    )
    
    if (
      !is.null(
        validacion_residuos$grafico_residuos
      )
    ) {
      
      guardar_figura(
        validacion_residuos$grafico_residuos,
        paste0(
          archivo,
          "_residuos_ajustados.png"
        )
      )
    }
  }
  
  
  # ==========================================================
  # SUPERPOSICIÓN FINAL
  # ==========================================================
  
  datos_grafico <-
    tibble::tibble(
      t = seq_len(T),
      fecha = datos$fecha,
      observado = y
    )
  
  datos_pronostico <-
    tibble::tibble(
      fecha = datos$fecha[(T_est + 1):T],
      pronostico = pronostico
    )
  
  grafico_final <-
    ggplot(
      datos_grafico,
      aes(
        x = fecha
      )
    ) +
    
    geom_line(
      aes(
        y = observado
      ),
      linewidth = 0.7
    ) +
    
    geom_line(
      data = datos_pronostico,
      aes(
        x = fecha,
        y = pronostico
      ),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    
    geom_vline(
      xintercept =
        datos$fecha[T_est],
      linetype = "dotted"
    ) +
    
    labs(
      title =
        paste(
          nombre,
          "- Observado vs pronóstico"
        ),
      x = "Fecha",
      y =
        attr(
          datos,
          "unidad"
        ),
      caption =
        "Línea vertical: inicio del período de validación"
    ) +
    
    theme_minimal()
  
  if (!is.null(archivo)) {
    
    guardar_figura(
      grafico_final,
      paste0(
        archivo,
        "_final.png"
      )
    )
  }
  
  
  # ==========================================================
  # TEXTO DEL PARÁMETRO
  # ==========================================================
  
  if (
    is.null(
      optimizacion
    )
  ) {
    
    if (
      metodo == "media"
    ) {
      
      parametro_texto <-
        "ninguno"
      
    } else {
      
      parametro_texto <-
        metodo
    }
    
  }
  
  else if (
    metodo == "holt"
  ) {
    
    parametro_texto <-
      paste0(
        "alpha=",
        format(
          optimizacion$optimo$alpha,
          digits = 4
        ),
        ", beta=",
        format(
          optimizacion$optimo$beta,
          digits = 4
        )
      )
    
  }
  
  else {
    
    parametro_texto <-
      paste0(
        "parametro=",
        format(
          optimizacion$optimo$parametro,
          digits = 4
        )
      )
  }
  
  
  # ==========================================================
  # RESULTADO RESUMIDO
  # ==========================================================
  
  resultado <-
    tibble::tibble(
      
      serie = nombre,
      
      metodo = metodo,
      
      parametro =
        parametro_texto,
      
      MSE_estimacion =
        indicadores_estimacion$MSE,
      
      MAD_estimacion =
        indicadores_estimacion$MAD,
      
      MAPE_estimacion =
        indicadores_estimacion$MAPE,
      
      MASE_estimacion =
        indicadores_estimacion$MASE,
      
      MSE_validacion =
        indicadores_validacion$MSE,
      
      MAD_validacion =
        indicadores_validacion$MAD,
      
      MAPE_validacion =
        indicadores_validacion$MAPE,
      
      MASE_validacion =
        indicadores_validacion$MASE,
      
      MASE_ingenuo =
        comparacion$ingenuo$MASE,
      
      optimo_en_borde =
        if (is.null(optimizacion)) NA else optimizacion$en_borde,
      
      LB_p =
        validacion_residuos$ljung_box$p_valor,
      
      JB_p =
        validacion_residuos$jarque_bera$p_valor,
      
      DW =
        validacion_residuos$durbin_watson,
      
      LB_serie_p =
        lb_serie$p_valor
    )
  
  
  # ==========================================================
  # OBJETO COMPLETO
  # ==========================================================
  
  list(
    
    nombre = nombre,
    
    datos = datos,
    
    descripcion = descripcion,
    
    ajuste = ajuste,
    
    optimizacion =
      optimizacion,
    
    acf_serie =
      acf_serie,
    
    ljung_box_serie =
      lb_serie,
    
    medidas_estimacion =
      indicadores_estimacion,
    
    pronostico =
      pronostico,
    
    medidas_validacion =
      indicadores_validacion,
    
    comparacion =
      comparacion,
    
    residuos =
      residuos,
    
    validacion_residuos =
      validacion_residuos,
    
    grafico_final =
      grafico_final,
    
    resumen =
      resultado
  )
}


# ============================================================
# SERIES DE LOS EJEMPLOS
# ============================================================

data(
  Nile,
  package = "datasets"
)

data(
  nottem,
  package = "datasets"
)

data(
  LakeHuron,
  package = "datasets"
)

data(
  WWWusage,
  package = "datasets"
)

data(
  airmiles,
  package = "datasets"
)

data(
  UKDriverDeaths,
  package = "datasets"
)

data(
  JohnsonJohnson,
  package = "datasets"
)

data(
  austres,
  package = "datasets"
)

data(
  AirPassengers,
  package = "datasets"
)


# ============================================================
# LECTURA ESTANDARIZADA
# ============================================================

datos_Nile <-
  leer_serie(
    Nile,
    fuente =
      "R: datasets::Nile (help('Nile'))",
    unidad =
      "flujo del río"
  )


datos_nottem <-
  leer_serie(
    nottem,
    fuente =
      "R: datasets::nottem (help('nottem'))",
    unidad =
      "temperatura"
  )


datos_LakeHuron <-
  leer_serie(
    LakeHuron,
    fuente =
      "R: datasets::LakeHuron (help('LakeHuron'))",
    unidad =
      "nivel del lago"
  )


datos_WWWusage <-
  leer_serie(
    WWWusage,
    fuente =
      "R: datasets::WWWusage (help('WWWusage'))",
    unidad =
      "usuarios"
  )


datos_airmiles <-
  leer_serie(
    airmiles,
    fuente =
      "R: datasets::airmiles (help('airmiles'))",
    unidad =
      "millas aéreas"
  )


datos_UKDriverDeaths <-
  leer_serie(
    UKDriverDeaths,
    fuente =
      "R: datasets::UKDriverDeaths (help('UKDriverDeaths'))",
    unidad =
      "muertes de conductores"
  )


datos_JohnsonJohnson <-
  leer_serie(
    JohnsonJohnson,
    fuente =
      "R: datasets::JohnsonJohnson (help('JohnsonJohnson'))",
    unidad =
      "ganancias por acción"
  )


datos_austres <-
  leer_serie(
    austres,
    fuente =
      "R: datasets::austres (help('austres'))",
    unidad =
      "población"
  )


datos_AirPassengers <-
  leer_serie(
    AirPassengers,
    fuente =
      "R: datasets::AirPassengers (help('AirPassengers'))",
    unidad =
      "miles de pasajeros"
  )


# ============================================================
# OCHO EJEMPLOS
# ============================================================

resultados <-
  list()


# ------------------------------------------------------------
# 1. MEDIA HISTÓRICA
# ------------------------------------------------------------

resultados[[1]] <-
  ejecutar_ejemplo(
    
    nombre =
      "Nile",
    
    datos =
      datos_Nile,
    
    metodo =
      "media",
    
    justificacion =
      paste(
        "Serie anual sin estacionalidad explícita.",
        "La media histórica proporciona un pronóstico",
        "constante basado en toda la información disponible."
      ),
    
    estacional =
      FALSE,
    
    archivo =
      "01_Nile_media"
  )


# ------------------------------------------------------------
# 2. MEDIA MÓVIL
# ------------------------------------------------------------

resultados[[2]] <-
  ejecutar_ejemplo(
    
    nombre =
      "nottem",
    
    datos =
      datos_nottem,
    
    metodo =
      "mm",
    
    justificacion =
      paste(
        "Serie mensual. Se utiliza una media móvil",
        "para evaluar el comportamiento de corto plazo",
        "mediante un parámetro k seleccionado por MSE."
      ),
    
    estacional =
      TRUE,
    
    archivo =
      "02_nottem_mm"
  )


# ------------------------------------------------------------
# 3. SES
# ------------------------------------------------------------

resultados[[3]] <-
  ejecutar_ejemplo(
    
    nombre =
      "LakeHuron",
    
    datos =
      datos_LakeHuron,
    
    metodo =
      "ses",
    
    justificacion =
      paste(
        "Serie anual. El suavizamiento exponencial simple",
        "permite actualizar el nivel con un peso alpha",
        "seleccionado mediante MSE."
      ),
    
    estacional =
      FALSE,
    
    archivo =
      "03_LakeHuron_ses"
  )


# ------------------------------------------------------------
# 4. DOBLE MEDIA MÓVIL
# ------------------------------------------------------------

resultados[[4]] <-
  ejecutar_ejemplo(
    
    nombre =
      "WWWusage",
    
    datos =
      datos_WWWusage,
    
    metodo =
      "dmm",
    
    justificacion =
      paste(
        "Serie con evolución temporal y tendencia local.",
        "La doble media móvil permite incorporar nivel",
        "y tendencia mediante k."
      ),
    
    estacional =
      FALSE,
    
    archivo =
      "04_WWWusage_dmm"
  )


# ------------------------------------------------------------
# 5. TENDENCIA LINEAL
# ------------------------------------------------------------

resultados[[5]] <-
  ejecutar_ejemplo(
    
    nombre =
      "airmiles",
    
    datos =
      datos_airmiles,
    
    metodo =
      "lineal",
    
    justificacion =
      paste(
        "Serie anual con crecimiento temporal.",
        "Se utiliza una tendencia lineal estimada",
        "mediante las ecuaciones normales."
      ),
    
    estacional =
      FALSE,
    
    archivo =
      "05_airmiles_lineal"
  )


# ------------------------------------------------------------
# 6. TENDENCIA CUADRÁTICA
# ------------------------------------------------------------

resultados[[6]] <-
  ejecutar_ejemplo(
    
    nombre =
      "UKDriverDeaths",
    
    datos =
      datos_UKDriverDeaths,
    
    metodo =
      "cuadratica",
    
    justificacion =
      paste(
        "Serie mensual. Se utiliza una tendencia",
        "cuadrática para representar una posible curvatura",
        "en la evolución temporal."
      ),
    
    estacional =
      TRUE,
    
    archivo =
      "06_UKDriverDeaths_cuadratica"
  )


# ------------------------------------------------------------
# 7. TENDENCIA EXPONENCIAL
# ------------------------------------------------------------

resultados[[7]] <-
  ejecutar_ejemplo(
    
    nombre =
      "JohnsonJohnson",
    
    datos =
      datos_JohnsonJohnson,
    
    metodo =
      "exponencial",
    
    justificacion =
      paste(
        "Serie trimestral con crecimiento aproximadamente",
        "multiplicativo. Se utiliza la transformación",
        "logarítmica para estimar una tendencia exponencial."
      ),
    
    estacional =
      TRUE,
    
    corregir_sesgo =
      FALSE,
    
    archivo =
      "07_JohnsonJohnson_exponencial"
  )


# ------------------------------------------------------------
# 8. HOLT
# ------------------------------------------------------------

resultados[[8]] <-
  ejecutar_ejemplo(
    
    nombre =
      "austres",
    
    datos =
      datos_austres,
    
    metodo =
      "holt",
    
    justificacion =
      paste(
        "Serie trimestral con nivel y tendencia.",
        "Holt actualiza ambos componentes mediante",
        "alpha y beta seleccionados por MSE."
      ),
    
    estacional =
      TRUE,
    
    archivo =
      "08_austres_holt"
  )


# ============================================================
# CONTRAEJEMPLO
# ============================================================

contraejemplo <-
  ejecutar_ejemplo(
    
    nombre =
      "AirPassengers",
    
    datos =
      datos_AirPassengers,
    
    metodo =
      "media",
    
    justificacion =
      paste(
        "Contraejemplo deliberado. La serie presenta",
        "estructura estacional y crecimiento, mientras",
        "que la media histórica no representa explícitamente",
        "esos componentes."
      ),
    
    estacional =
      TRUE,
    
    archivo =
      "09_AirPassengers_contraejemplo"
  )


# ============================================================
# TABLA FINAL
# ============================================================

tabla_resultados <-
  purrr::map_dfr(
    resultados,
    function(x) {
      x$resumen
    }
  )

print(
  tabla_resultados
)

cat(
  "\nOptimización: un TRUE en optimo_en_borde indica que el mínimo de la rejilla está en uno de sus bordes; en el informe debe interpretarse antes de ampliar o cambiar la rejilla.\n"
)


# ============================================================
# TABLA DEL CONTRAEJEMPLO
# ============================================================

print(
  contraejemplo$resumen
)


# ============================================================
# VERIFICACIÓN ACF MANUAL VS acf()
# ============================================================

m_verificacion <-
  max(
    1,
    min(
      floor(
        length(
          datos_Nile$y
        ) / 4
      ),
      24
    )
  )

acf_manual_Nile <-
  acf_manual(
    datos_Nile$y,
    m =
      m_verificacion
  )

acf_R_Nile <-
  stats::acf(
    datos_Nile$y,
    lag.max =
      m_verificacion,
    plot = FALSE
  )

diferencia_acf <-
  max(
    abs(
      acf_manual_Nile$acf -
        as.numeric(
          acf_R_Nile$acf[-1]
        )
    )
  )

cat(
  "\nDiferencia máxima ACF manual vs acf():",
  diferencia_acf,
  "\n"
)

stopifnot(
  diferencia_acf < 1e-12
)


# ============================================================
# VERIFICACIÓN DE TENDENCIAS MANUALES VS lm()
# ============================================================
# lm() aparece únicamente aquí, como exige la tarea.

verificacion_lineal <-
  ajustar_tendencia(
    datos_airmiles$y,
    tipo = "lineal"
  )

modelo_lineal_R <-
  stats::lm(
    y ~ t,
    data = datos_airmiles
  )

stopifnot(
  max(
    abs(
      verificacion_lineal$parametros$tabla_coeficientes$estimacion -
        stats::coef(modelo_lineal_R)
    )
  ) < 1e-10
)

verificacion_cuadratica <-
  ajustar_tendencia(
    datos_UKDriverDeaths$y,
    tipo = "cuadratica"
  )

modelo_cuadratico_R <-
  stats::lm(
    y ~ t + I(t^2),
    data = datos_UKDriverDeaths
  )

stopifnot(
  max(
    abs(
      verificacion_cuadratica$parametros$tabla_coeficientes$estimacion -
        stats::coef(modelo_cuadratico_R)
    )
  ) < 1e-10
)

verificacion_exponencial <-
  ajustar_tendencia(
    datos_JohnsonJohnson$y,
    tipo = "exponencial"
  )

modelo_exponencial_R <-
  stats::lm(
    log(y) ~ t,
    data = datos_JohnsonJohnson
  )

coef_exponencial_R <-
  stats::coef(modelo_exponencial_R)

stopifnot(
  max(
    abs(
      verificacion_exponencial$parametros$tabla_coeficientes$estimacion -
        coef_exponencial_R
    )
  ) < 1e-10
)

cat(
  "\nVerificación de tendencias contra lm(): OK\n"
)


# ============================================================
# VERIFICACIÓN LJUNG-BOX MANUAL VS Box.test()
# ============================================================

verificacion_ljung_box_resultado <-
  verificar_ljung_box(
    datos_Nile$y,
    m =
      m_verificacion,
    p = 0
  )

print(
  verificacion_ljung_box_resultado
)

stopifnot(
  verificacion_ljung_box_resultado$coincide
)


# ============================================================
# VERIFICACIÓN SES
# ============================================================

verificacion_ses_resultado <-
  verificar_ses(
    LakeHuron,
    alpha =
      0.30
  )

print(
  verificacion_ses_resultado
)

stopifnot(
  verificacion_ses_resultado$coincide
)


# ============================================================
# VERIFICACIÓN HOLT
# ============================================================

verificacion_holt_resultado <-
  verificar_holt(
    austres,
    alpha =
      0.30,
    beta =
      0.20
  )

print(
  verificacion_holt_resultado
)

stopifnot(
  verificacion_holt_resultado$coincide
)


# ============================================================
# GUARDADO DE RESULTADOS
# ============================================================

saveRDS(
  resultados,
  file =
    file.path(
      "figs",
      "resultados_ejemplos.rds"
    )
)

saveRDS(
  contraejemplo,
  file =
    file.path(
      "figs",
      "contraejemplo.rds"
    )
)

utils::write.csv(
  tabla_resultados,
  file =
    file.path(
      "figs",
      "tabla_resultados.csv"
    ),
  row.names = FALSE
)


# ============================================================
# MENSAJE FINAL
# ============================================================

cat(
  "\n====================================================\n"
)

cat(
  "Ejecución finalizada correctamente.\n"
)

cat(
  "Número de ejemplos:",
  length(resultados),
  "\n"
)

cat(
  "Contraejemplo:",
  contraejemplo$nombre,
  "\n"
)

cat(
  "Diferencia máxima ACF:",
  diferencia_acf,
  "\n"
)

cat(
  "====================================================\n"
)

# ============================================================
# FIN DE ejemplos/ejemplos.R
# ============================================================
