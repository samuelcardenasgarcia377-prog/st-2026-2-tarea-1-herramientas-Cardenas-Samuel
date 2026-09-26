[PEGAR AQUÍ EL ENLACE AL REPOSITORIO DE GITHUB]

# Tarea 1 — Caja de herramientas de pronóstico

## Contenido del repositorio

- `R/00-lectura.R`: lectura y validación de series.
- `R/01-graficos.R`: gráficos de series y correlogramas.
- `R/02-metodos.R`: implementación manual de los ocho métodos.
- `R/03-evaluacion.R`: métricas y validación de residuos.
- `ejemplos/ejemplos.R`: ejecución de los ocho ejemplos y generación de figuras.
- `informe/informe.qmd`: fuente del informe.
- `informe/informe.html`: informe renderizado para entrega.
- `figs/`: figuras generadas durante el análisis.
- `sesion-info.txt`: información de la sesión de R.

## Ejecución

Desde la raíz del repositorio, con R 4.3 o superior:

```r
source("ejemplos/ejemplos.R")
```

El script ejecuta los ejemplos y genera las figuras correspondientes.

Para generar el informe, abrir `informe/informe.qmd` en RStudio y seleccionar **Render**. El resultado es `informe/informe.html`.

## Convenciones

- Las series se trabajan como vectores numéricos ordenados y sin faltantes.
- Los métodos reciben una serie y devuelven al menos `yhat`, `pronosticar(h)` y los parámetros estimados.
- Los parámetros de los métodos que requieren optimización se seleccionan mediante las rejillas definidas en la tarea.
- Las métricas utilizadas son MSE, MAD, MAPE y MASE.
- La validación incluye análisis de los residuos mediante media cero, Ljung–Box, Jarque–Bera y Durbin–Watson.
- No se utilizan funciones de paquetes de pronóstico prohibidos por la tarea.

## Resultados

Los resultados completos de los ocho ejemplos, incluyendo parámetros, medidas de error, diagnósticos y gráficos, se encuentran en `informe/informe.html`.

| Ejemplo | Serie | Método | Parámetros | MASE método | MASE referente |
|---|---|---|---|---:|---:|
| 1 | Nile | Media histórica | — | Ver informe | Ver informe |
| 2 | nottem | Media móvil | k optimizado | Ver informe | Ver informe |
| 3 | LakeHuron | SES | α optimizado | Ver informe | Ver informe |
| 4 | WWWusage | Doble media móvil | k optimizado | Ver informe | Ver informe |
| 5 | airmiles | Tendencia lineal | — | Ver informe | Ver informe |
| 6 | UKDriverDeaths | Tendencia cuadrática | — | Ver informe | Ver informe |
| 7 | JohnsonJohnson | Tendencia exponencial | — | Ver informe | Ver informe |
| 8 | austres | Holt | α, β optimizados | Ver informe | Ver informe |

## Uso de inteligencia artificial

Se utilizó ChatGPT como apoyo para revisar la estructura del repositorio, detectar y corregir errores de código, organizar el informe y facilitar la documentación. La implementación, ejecución y verificación de los procedimientos se realizaron en R y fueron revisadas antes de la entrega.
