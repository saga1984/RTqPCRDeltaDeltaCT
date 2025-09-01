#' Pipeline para graficar el ΔΔCt (método de Livak) para datos provenientes del equipo CFX96 (BioRad)' 
#'
#' Lee un CSV con Cq de RT-qPCR (CFX96/CFX Manager), filtra NTC/NA, selecciona las
#' mejores réplicas por grupo, solicita de forma interactiva el gen de interés,
#' genes de referencia y condiciones (control y tratamiento), calcula ΔCt/ΔΔCt y
#' Fold Change (2^-ΔΔCt), genera un gráfico de barras con error y guarda salidas.
#'
#' @param ruta_rtqpcr `character(1)` Ruta al CSV con columnas mínimas:
#'   `Target`, `Content`, `Sample`, `Cq`.
#' @param numero_replicas `integer(1)` Réplicas a conservar por grupo (se eligen
#'   las más cercanas a la media). Por defecto `3`.
#' @param resolucion `numeric(1)` Resolución (ppi) para dispositivos gráficos al
#'   exportar figuras. Por defecto `600`.
#' @param formatos `character()` Vector de formatos base a exportar (p. ej.
#'   `"jpeg"`, `"tiff"`, etc.). Deben existir funciones homónimas (`jpeg()`,
#'   `tiff()`, …). Por defecto `"jpeg"`.
#'
#' @details
#' La función es **interactiva** (usa `readline()`): solicita `gen_interes`,
#' `gen_referencia` (uno o varios separados por coma o espacio; se promedian),
#' `calibrador` (muestra control) y `tratamiento` (muestra tratada).
#'
#' Escribe en `dirname(ruta_rtqpcr)`:
#' - `Tabla_Fold_Change_Data_<archivo>_<tratamiento>_<GOI>_<REFS>.csv` (por réplica)
#' - `Fold_Change_Data_<archivo>_<tratamiento>_<GOI>_<REFS>.csv` (promedio y sd)
#' - `Tabla_filtrada_<archivo>_...csv` (datos usados)
#' - Figura `Boxplot_<archivo>.<ext>` y `fold_Livak_<archivo>_<...>.<ext>`
#'
#' **Supuestos de entrada**: el CSV contiene todos los targets y muestras
#' necesarios. Los NTC con Cq ≥ 31 se marcan como sospechosos; filas con `Cq`
#' faltante se excluyen. Para varias referencias se calcula la media de Cq por réplica.
#'
#' @return (invisible) una lista con:
#'   - `datos_filtrados`: data frame usado para ΔΔCt
#'   - `Fold_Change_df`: resultados por réplica (ΔCt, ΔΔCt, FC)
#'   - `Fold_Change_final`: promedios y desviaciones por grupo
#'   - `ggfold_Livak`: objeto `ggplot` del barplot final
#'
#' @section Efectos colaterales:
#' Abre dispositivos gráficos base, imprime/guarda figuras y escribe múltiples CSV.
#'
#' @seealso \link[ggplot2]{ggplot}, \link[psych]{describe}, \link[tidyr]{pivot_wider}
#'
#' @examples
#' \dontrun{
#' RTqPCR_data_analysis("2025/articulo_Rosy/target_genes.csv")
#' }
#'
#' @export
#' @import psych dplyr ggplot2 tidyr

# definir funcion
RTqPCR_data_analysis <- function(ruta_rtqpcr,
                                 numero_replicas = 3,
                                 resolucion = 600,
                                 formatos = "jpeg"){
  
  # librerias
  library(psych)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  # crear objeto de R
  datos <- read.csv(file = ruta_rtqpcr,
                    header = TRUE,
                    na.strings = "")
  
  # filtrar columnas
  datos <- datos[,c("Target",	"Content", "Sample",	"Cq")]
  
  # ordenar 
  datos <- datos[order(datos$Target, 
                       datos$Content, 
                       datos$Sample),]
  
  ############################ evaluar NAs vs genes ##############################
  
  # ver NAs
  NAs <- datos[datos$Content == "NTC",]
  
  # NAs que no amplificaron
  NAs_buenos <- NAs[is.na(NAs$Cq),]
  
  # NAs sospechosos
  NAs_malos <- NAs[NAs$Cq >= 31,]
  
  # datos sospechosos
  datos_malos <- datos[is.na(datos$Cq),]
  
  # eliminar datos sospechosos
  datos_buenos <- datos[!is.na(datos$Cq),]
  
  # datos_buenos_limpios
  datos_buenos_limpios <- datos_buenos[!datos_buenos$Content == "NTC",]
  
  # graficar boxplot
  boxplot <- ggplot(data = datos_buenos_limpios, 
                    aes(x = Target, 
                        y = Cq, 
                        fill = Sample)) +
    geom_boxplot() +
    facet_grid(Target~.) + # grid en filas
    theme_minimal() 
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    # crear y guardar los heatmaps
    match.fun(i)(paste(dirname(ruta_rtqpcr), 
                       "/Boxplot_", 
                       gsub("\\..+$", "",basename(ruta_rtqpcr)),
                       ".",i, 
                       sep = ""),
                 res = resolucion,
                 width = 2000,
                 height = 3000)
    
    # crear y guardar el heatmpat euclidean
    print(
      boxplot
    )
    dev.off()
  }
  
  # obtener promedios y desviaciones estandar (solo necesario cuando hay varias replicas)
  datos_agrupados <- datos_buenos_limpios %>% 
    group_by(Target, Sample) %>% mutate(
      Mean = mean(Cq),
      Dif_Abs = abs(Cq - Mean)) %>%
    arrange(Target, Sample, Dif_Abs) %>%
    slice_head(n = numero_replicas)
  
  ############################ pedir datos a ususarios ###########################
  
  # imprimir mensaje de genes
  cat("Los genes disponibles en el spreadsheet son:\n ")
  
  # imprimir en lista genes
  for(gene in unique(datos_agrupados$Target)){
    cat("\n", gene)
  }
  
  # pedir a usuario gen problema
  gen_interes <- readline(prompt ="Ingrese gen de interes: ")
  
  # pedir a usuario gen control
  gen_referencia <- readline(prompt = "Ingerese gen de referencia (como valores separados por espacio o coma): ")
  
  # imprimir mensaje de tratamientos
  cat("Los tratamientos disponible en el spreadsheet son:\n ")
  
  # imprimir en lista tratamientos/condiciones
  for(tratamientos in unique(datos_agrupados$Sample)){
    cat("\n", tratamientos)
  }
  
  # pedir a usuario ingresar tratmiento/condicion control
  calibrador <- readline(prompt = "Ingrese condicion control/calibrador: ")
  
  # pedir a usuario ingresar tratmiento/condicion
  tratamiento <- readline(prompt = "Ingrese tratamiento/condicion: ")
  
  ########################## comenzar a filtrar ##################################
  
  # filtrar en base a las 4 variables recie creadas
  datos_filtrados <- datos_agrupados[grepl(paste(gsub("[ ,]","|", gen_referencia), gen_interes, sep = "|"), 
                                           datos_agrupados$Target) &
                                       grepl(paste(tratamiento, calibrador, sep = "|"), 
                                             datos_agrupados$Sample), ] 
  
  # remover columnas usadas para obtener mejores replicas
  datos_filtrados[, c("Mean", "Dif_Abs")] <- NULL
  
  print("Evaluando mas de un gen de referencia")
  
  # # crear nuevo data frame para guardar Deltas-CTs y Fold Change
  # fold_change <- data.frame(
  #   "Sample_Type" = c("Control 1",
  #                "Control 2",
  #                "Control 3",
  #                "Treated 1",
  #                "Treated 2",
  #                "Treated 3"),
  #   "Sample" = 1:6,
  #   "DCt" = 1:6,
  #   "DDCt" = 1:6,
  #   "Fold_Change" =1:6)
  # 
  # llenar columnas 1 y 2
  # fold_change[grepl("Treated",fold_change$Sample_Type), "Sample"] <- paste(tratamiento, 1:3, sep = "_")
  # fold_change[grepl("Control",fold_change$Sample_Type), "Sample"] <- paste(gsub("[ ,]","|", gen_referencia), 1:3, sep = "_")
  
  #################### obtener columnas de Fold Change ###########################
  
  # convertir genes de referencia en un vector
  ref_vec <- strsplit(gen_referencia, ",|[ ]")[[1]]
  
  # filtrar genes de referencia calibrador
  cal_ref <- datos_filtrados[grepl(gsub(",|[ ]", "|", gen_referencia), datos_filtrados$Target) & 
                               datos_filtrados$Sample == calibrador,]
  
  # emparejar genes para obtener medias con pivot_wider
  cal_ref <- cal_ref %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    unnest(cols = ref_vec) %>% 
    ungroup() 
  
  # obtener nueva columna de medias
  cal_ref[, "Cq_Means"] <- rowMeans(cal_ref[,ref_vec], na.rm = T) 
  
  # agregar columna de nombre de genes de referencia
  cal_ref[, "Target"] <- gsub(",|[ ]", "|", gen_referencia) 
  
  # filtrar columnas de interes
  cal_ref <- cal_ref[, c("Target", "Sample", "Cq_Means")]
  
  # filtrar genes de referencia tratamiento
  trat_ref <- datos_filtrados[grepl(gsub(",|[ ]", "|", gen_referencia), datos_filtrados$Target) & 
                                datos_filtrados$Sample == tratamiento,]
  
  # emparejar genes para obtener medias con pivot_wider
  trat_ref <- trat_ref %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    unnest(ref_vec) %>% ungroup()
  
  # obtener nueva columna de medias
  trat_ref[, "Cq_Means"] <- rowMeans(trat_ref[,ref_vec], na.rm = T)
  
  # agregar columna de nombre de genes de referencia
  trat_ref[, "Target"] <- gsub(",|[ ]", "|", gen_referencia) 
  
  # filtrar columnas de interes
  trat_ref <- trat_ref[, c("Target", "Sample", "Cq_Means")]
  
  # obtener control goi
  cal_goi <- datos_filtrados[datos_filtrados$Target == gen_interes &
                               datos_filtrados$Sample == calibrador,]
  
  # eliminar columna innecesaria
  cal_goi[,"Content"] <- NULL
  
  # obtener tratmiento goi
  trat_goi <- datos_filtrados[datos_filtrados$Target == gen_interes &
                                datos_filtrados$Sample == tratamiento,]
  
  # eliminar columna innecesaria
  trat_goi[, "Content"] <- NULL
  
  # ordenar data frames
  cal_goi <- cal_goi[order(cal_goi$Cq),]
  trat_goi <- trat_goi[order(trat_goi$Cq),]
  cal_ref <- cal_ref[order(cal_ref$Cq_Means),]
  trat_ref <- trat_ref[order(trat_ref$Cq_Means),]
  
  # Delta CT Calibrador
  DCT_C13 <- cal_goi$Cq - cal_ref$Cq_Means
  DCT_C13 <- data.frame("Delta_CT" = DCT_C13)
  
  # EXTRA: CALIBRATOR
  calibrator <- colMeans(DCT_C13)
  
  # Delta CT Tratamiento
  DCT_T13 <- trat_goi$Cq - trat_ref$Cq_Means
  DCT_T13 <- data.frame("Delta_CT" = DCT_T13)
  
  # unir Delta CTs
  DCT <- rbind(DCT_C13, DCT_T13)
  
  # Delta Delta CT Calibrador
  DDCT_C13 <- DCT_C13 - calibrator
  DDCT_C13 <- data.frame("DeltaDelta_CT" = DDCT_C13)
  
  # Delta Delta CT Tratamiento
  DDCT_T13 <- DCT_T13 - calibrator
  DDCT_T13 <- data.frame("DeltaDelta_CT" = DDCT_T13)
  
  # unir Delta CTs
  DDCT <- rbind(DDCT_C13, DDCT_T13)
  
  # Fold change
  Fold_Change_Calibrador <- 2^-(DDCT_C13)
  Fold_Change_Tratamiento <- 2^-(DDCT_T13)
  
  # agregar columna de Sample
  Fold_Change_df <- data.frame("Sample" = c(rep("Control", 3), 
                                            rep("Tratamiento", 3))
  )
  
  # crear columna de fold change
  fold_change <- rbind(Fold_Change_Calibrador,
                       Fold_Change_Tratamiento)
  
  # unir columnas de deltas CTs y Fold Change
  Fold_Change_df[, "Sample_Real"] <- c(rep(gsub(",|[ ]","|", gen_referencia), 3),
                                       rep(tratamiento, 3))
  Fold_Change_df[,"DeltaCT"] <- DCT
  Fold_Change_df[,"DeltaDeltaCT"] <- DDCT
  Fold_Change_df[,"log2FC"] <- fold_change
  
  # agregar Desvest y Average
  Fold_Change_final <- Fold_Change_df %>% 
    group_by(Sample, Sample_Real) %>%
    summarise(Average = mean(log2FC),
              Desvest = sd(log2FC))
  
  ################################################################################
  
  # guardar data table de resultados de fold change
  write.csv(x = Fold_Change_df,
            file = paste(dirname(ruta_rtqpcr), 
                         "/Tabla_Fold_Change_Data_", 
                         gsub("\\..+$", "",basename(ruta_rtqpcr)), "_",
                         paste(tratamiento, "_",
                               gen_interes, "_",
                               gsub("[ ,]", "_", gen_referencia),
                               sep = ""),
                         ".csv", 
                         sep = ""),
            row.names = FALSE)
  
  
  # guardar data table de resultados finales de fold change
  write.csv(x = Fold_Change_final,
            file = paste(dirname(ruta_rtqpcr),
                         "/Fold_Change_Data_", 
                         gsub("\\..+$", "",basename(ruta_rtqpcr)), "_",
                         paste(tratamiento, "_",
                               gen_interes, "_",
                               gsub("[ ,]", "_", gen_referencia),
                               sep = ""),
                         ".csv", 
                         sep = ""),
            row.names = FALSE)
  
  # guardar archivo filtrado
  write.csv(x = datos_filtrados,
            file = paste(dirname(ruta_rtqpcr), 
                         "/Tabla_filtrada_", 
                         gsub("\\..+$", "",basename(ruta_rtqpcr)), "_",
                         paste(tratamiento, "_",
                               gen_interes, "_",
                               gsub("[ ,]", "_", gen_referencia),
                               sep = ""),
                         ".csv", 
                         sep = ""),
            row.names = FALSE)
  
  ################################################################################
  
  # graficar
  ggfold_Livak <- ggplot(data = Fold_Change_final, 
                         mapping = aes(x = Sample_Real,
                                       y = Average, 
                                       fill = Sample)) +
    geom_bar(stat = "identity") + 
    geom_errorbar(mapping = aes(ymin = Average - Desvest, 
                                ymax = Average + Desvest),
                  width=.2,
                  position=position_dodge(.9)) +
    theme_minimal() +
    labs(title = paste(gen_interes, 
                       " (", gsub("[ ,]"," & ", gen_referencia),") ",
                       sep = ""),
         x = "",
         y = "Fold Change")
  
  # mostrar en R
  print(ggfold_Livak)
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    # crear y guardar los heatmaps
    match.fun(i)(paste(dirname(ruta_rtqpcr), 
                       "/fold_Livak_", 
                       gsub("\\..+$", "",basename(ruta_rtqpcr)), "_",
                       paste(tratamiento, "_",
                             gen_interes, "_",
                             gsub("[ ,]","_", gen_referencia),
                             sep = ""),
                       ".",i, 
                       sep = ""),
                 res = resolucion,
                 width = 3000,
                 height = 2000)
    
    # crear y guardar el heatmpat euclidean
    print(
      ggfold_Livak
    )
    dev.off()
  }
  
}