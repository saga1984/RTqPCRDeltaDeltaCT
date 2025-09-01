#' Unir y graficar resultados de ΔΔCt ya procesados
#'
#' Reúne automáticamente los CSV generados por \code{RTqPCR_data_analysis()} para un
#' gen de interés y una condición, consolida en data frames “grandes” y produce
#' barplots y boxplots (con errores estándar) para \emph{Fold Change} y \emph{log2(FC)}.
#' Además, persiste objetos útiles en \code{.GlobalEnv}.
#'
#' @param gen_interes `character(1)` Gen de interés (p. ej. `"ERF03"`).
#' @param condicion `character(1)` Nombre corto de la condición/tratamiento
#'   (se usa en títulos y nombres de archivos).
#' @param gen_referencia `character()` Uno o varios genes de referencia. Si son
#'   varios, se asume que los CSV usan el nombre combinado con guiones bajos
#'   (p. ej. `Ubiq_Rps1`).
#' @param ruta_tablas `character(1)` Carpeta donde residen los CSV producidos por
#'   \code{RTqPCR_data_analysis()}.
#' @param resolucion `numeric(1)` ppp para exportar figuras. Default: 600.
#' @param formatos `character()` Formatos a exportar (funciones base homónimas:
#'   `tiff()`, `jpeg()`, …). Default: `c("tiff","jpeg")`.
#' @param fig_ancho,fig_alto `numeric(1)` Tamaño en píxeles de las figuras.
#'   Defaults: 5000 × 4000.
#'
#' @details
#' Busca archivos que comienzan con:
#' \itemize{
#'   \item \code{Fold_Change_Data_*.csv} (promedios/SE) para barplots.
#'   \item \code{Tabla_Fold_Change_Data_*.csv} (réplicas) para boxplots.
#' }
#' El patrón final incluye \code{<gen_interes>_<gen_referencia>.csv}. Se espera que
#' los CSV contengan, entre otras, las columnas:
#' \code{Sample}, \code{Treatment}, \code{Reference_Gene},
#' \code{FC_Average}, \code{SE_FC}, \code{log2FC_Average}, \code{SE_log2FC},
#' \code{FoldChange}, \code{log2FC}.
#'
#' Objetos creados en \code{.GlobalEnv}:
#' \itemize{
#'   \item \code{big_df_bar_<condicion>_<gen>_<refs>}
#'   \item \code{big_df_box_<condicion>_<gen>_<refs>}
#'   \item \code{FC_barplot_agrupado_*}, \code{log2FC_barplot_agrupado_*} (ggplot)
#' }
#'
#' Archivos escritos en \code{ruta_tablas}:
#' \itemize{
#'   \item \code{Fold_Change_agrupado_<cond>_<gen>_<refs>.csv}
#'   \item \code{Tabla_Fold_Change_agrupado_<cond>_<gen>_<refs>.csv}
#'   \item Figuras: \code{FC_barplot_agrupado_*}, \code{log2FC_barplot_agrupado_*},
#'         \code{log2FC/FC_boxplot_agrupado_*} en los formatos solicitados.
#' }
#'
#' @return (invisible) una lista con:
#' \itemize{
#'   \item \code{bar_df}, \code{box_df}: data frames consolidados
#'   \item \code{p_fc_bar}, \code{p_log2fc_bar}, \code{p_fc_box}, \code{p_log2fc_box}: objetos \code{ggplot}
#' }
#'
#' @section Notas:
#' - Los límites de \code{ylim()} se calculan a partir de medias ± SE.
#' - Los colores de muestra están fijados a \code{Control} y \code{Treatment}.
#'
#' @examples
#' \dontrun{
#' RTqPCR_data_analysis_unir_datos(
#'   gen_interes = "ERF03",
#'   condicion = "Dehydration",
#'   gen_referencia = c("Ubiq","Rps1"),
#'   ruta_tablas = "2025/articulo_Rosy/"
#' )
#' }
#'
#' @export
#' @import dplyr ggplot2 ggrepel


# funcion para unir datos por tratamiento/condicion y por gen de interes
RTqPCR_data_analysis_unir_datos <- function(gen_interes,
                                            condicion,
                                            gen_referencia,
                                            ruta_tablas,
                                            resolucion = 600,
                                            formatos = c("tiff", "jpeg"),
                                            fig_ancho = 5000,
                                            fig_alto = 4000){

  
  # librerias
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  # listar documentos en directorio de windows en base a metodo
  archivos_filtrados <- list.files(ruta_tablas, pattern = paste("^Fold_Change_Data_", sep = ""), full.names = TRUE)
  
  # determinar longitud de vector de genes de referencia, para distinguir 1 solo gen versus combinaciones
  n_ref_genes <- length(gen_referencia)
  
  if (n_ref_genes == 1) {
    # Solo un gen de referencia: buscar archivos que terminen exactamente con ese nombre
    patron <- paste0(gen_interes, "_", gen_referencia, "\\.csv$")
  } else {
    # Varios genes de referencia: buscar archivos que contengan todos los nombres combinados con "_"
    # definir variable
    gen_referencia <- paste(gen_referencia, collapse = "_")
    patron <- paste(gen_interes, gen_referencia, sep = "_")
  }
  
  
  # filtrar documentos de interes por tratamiento y gen de interes
  archivos_interes <- archivos_filtrados[grepl(patron, 
                                               archivos_filtrados)]
  
  # crear objetos de R: varios data frames unidos en formato de lista
  data_frames_list <- lapply(archivos_interes, read.csv)
  
  # convertir a data frame
  big_data_frame <- as.data.frame(do.call(rbind, data_frames_list))
  
  assign(paste("big_df_bar_", condicion, "_", gen_interes, "_", gen_referencia, sep = ""),
         big_data_frame, 
         envir = .GlobalEnv)
  
  # guardar data frame
  write.csv(x = big_data_frame,
            file = paste(paste(ruta_tablas, 
                               paste("/Fold_Change_agrupado_", 
                                     condicion, "_", gen_interes, "_", gen_referencia, sep = ""), 
                               sep = ""),
                         ".csv", 
                         sep = ""),
            row.names = FALSE)
  
  ############################# graficar #########################################
  
  # Vector de colores correcto (fuera del dataframe)
  colores_fijos <- c("Control" = "#1f77b4", 
                     "Treatment" = "#ff7f0e")
  
  
  # volver factor para ordenar eje X
  big_data_frame$Treatment <- factor(big_data_frame$Treatment,
                                     levels = c("0 Hours", 
                                                "15 Minutes", 
                                                "30 Minutes", 
                                                "1 Hour", 
                                                "3 Hours", 
                                                "24 Hours"))
  
  
  # definir maximos en y
  maximo_y <- max(big_data_frame$FC_Average + big_data_frame$SE_FC)
  minimo_y <- min(big_data_frame$FC_Average - big_data_frame$SE_FC)
  
  
  # graficar
  FC_bar_combinado <- ggplot(big_data_frame, 
                             aes(x = Treatment, 
                                 y = FC_Average, 
                                 fill = Sample)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = FC_Average - SE_FC, 
                      ymax = FC_Average + SE_FC), 
                  width = 0.2, 
                  position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = colores_fijos) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      linewidth = 1.2)) +
    ylim(minimo_y * -0.5 , maximo_y * 1.2) +
    labs(title = paste(condicion, ": ", gen_interes,
                       " (", gsub("_", " & ", gen_referencia), ")",
                       sep = ""),
         y = "Fold Change",
         x = "",
         fill = "")
  
  # imprimir en R
  print(FC_bar_combinado)
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    # crear y guardar los heatmaps
    match.fun(i)(paste(ruta_tablas, 
                       paste("/FC_barplot_agrupado_", 
                             condicion, "_", gen_interes, "_", gen_referencia, sep = ""), 
                       ".",i, 
                       sep = ""),
                 res = resolucion,
                 width = fig_ancho,
                 height = fig_alto)
    
    # crear y guardar el heatmpat euclidean
    print(
      FC_bar_combinado
    )
    dev.off()
  }
  
  assign(paste("FC_barplot_agrupado_", condicion, "_", gen_interes, "_", gen_referencia, sep = ""),
         FC_bar_combinado, 
         envir = .GlobalEnv)
  
  ##############################################################################
  
  # definir maximos en y
  maximo_y <- max(big_data_frame$log2FC_Average + big_data_frame$SE_log2FC)
  minimo_y <- min(big_data_frame$log2FC_Average + big_data_frame$SE_log2FC)
  
  # Graficar
  log2FC_bar_combinado <- ggplot(big_data_frame, 
                                 aes(x = Treatment, 
                                     y = log2FC_Average, 
                                     fill = Sample)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = log2FC_Average - SE_log2FC, 
                      ymax = log2FC_Average + SE_log2FC), 
                  width = 0.2, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = colores_fijos) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      linewidth = 1.2)) +
    ylim(ifelse(minimo_y < 0, minimo_y*1.1, minimo_y * -0.5), maximo_y * 1.2) +
    labs(title = paste(condicion, ": ", gen_interes,
                       " (", gsub("_", " & ", gen_referencia), ")",
                       sep = ""),
         y = "log2(Fold Change)",
         x = "",
         fill = "")
  
  # imprimir en R
  print(log2FC_bar_combinado)
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    # crear y guardar los heatmaps
    match.fun(i)(paste(ruta_tablas, 
                       paste("/log2FC_barplot_agrupado_", 
                             condicion, "_", gen_interes, "_", gen_referencia, sep = ""), 
                       ".",i, 
                       sep = ""),
                 res = resolucion,
                 width = fig_ancho,
                 height = fig_alto)
    
    # crear y guardar el heatmpat euclidean
    print(
      log2FC_bar_combinado
    )
    dev.off()
  }
  
  assign(paste("log2FC_barplot_agrupado_", condicion, "_", gen_interes, "_", gen_referencia, sep = ""),
         log2FC_bar_combinado, 
         envir = .GlobalEnv)
  
  ################################# BOXPLOT ######################################
  
  
  # listar documentos en directorio de windows en base a metodo
  archivos_filtrados <- list.files(ruta_tablas, pattern = paste("^Tabla_Fold_Change_Data_", sep = ""), full.names = TRUE)
  
  # determinar longitud de vector de genes de referencia, para distinguir 1 solo gen versus combinaciones
  n_ref_genes <- length(gen_referencia)
  
  if (n_ref_genes == 1) {
    # Solo un gen de referencia: buscar archivos que terminen exactamente con ese nombre
    patron <- paste0(gen_interes, "_", gen_referencia, "\\.csv$")
  } else {
    # Varios genes de referencia: buscar archivos que contengan todos los nombres combinados con "_"
    # definir variable
    gen_referencia <- paste(gen_referencia, collapse = "_")
    patron <- paste(gen_interes, gen_referencia, sep = "_")
  }
  
  
  # filtrar documentos de interes por tratamiento y gen de interes
  archivos_interes <- archivos_filtrados[grepl(patron, 
                                               archivos_filtrados)]
  
  # crear objetos de R: varios data frames unidos en formato de lista
  data_frames_list <- lapply(archivos_interes, read.csv)
  
  # convertir a data frame
  big_data_frame <- as.data.frame(do.call(rbind, data_frames_list))
  
  assign(paste("big_df_box_", condicion, "_", gen_interes, "_", gen_referencia, sep = ""),
         big_data_frame, 
         envir = .GlobalEnv)
  
  # guardar data frame
  write.csv(x = big_data_frame,
            file = paste(paste(ruta_tablas, 
                               paste("/Tabla_Fold_Change_agrupado_", 
                                     condicion, "_", gen_interes, "_", gen_referencia, sep = ""), 
                               sep = ""),
                         ".csv", 
                         sep = ""),
            row.names = FALSE)
  
  #############################################
  
  # volver factor para ordenar eje X
  big_data_frame$Treatment <- factor(big_data_frame$Treatment,
                                     levels = c("0 Hours", 
                                                "15 Minutes", 
                                                "30 Minutes", 
                                                "1 Hour", 
                                                "3 Hours", 
                                                "24 Hours"))
  
  # definir maximos en y
  maximo_y <- max(big_data_frame$log2FC)
  minimo_y <- min(big_data_frame$log2FC)
  
  # Boxplot para log2FC
  log2FC_box_combinado <- ggplot(big_data_frame,
                                 aes(x = Treatment, 
                                     y = log2FC, 
                                     fill = Sample)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black") +
    geom_jitter(width = 0.2, 
                size = 2,
                alpha = 0.5) +
    scale_fill_manual(values = colores_fijos) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          panel.border = element_rect(color = "black", 
                                      fill = NA, 
                                      linewidth = 1.2)) +
    ylim(ifelse(minimo_y < 0, minimo_y*1.2, minimo_y * -0.5), maximo_y * 1.2) +
    labs(title = paste(condicion, ": ", gen_interes,
                       " (", gsub("_", " & ", gen_referencia), ")",
                       sep = ""),
         y = "log2(Fold Change)",
         x = "",
         fill = "")
  
  print(log2FC_box_combinado)
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    # crear y guardar los heatmaps
    match.fun(i)(paste(ruta_tablas, 
                       "/log2FC_boxplot_agrupado_", 
                       condicion, "_", gen_interes, "_", gen_referencia,
                       ".",i, sep = ""),
                 res = resolucion,
                 width = fig_ancho,
                 height = fig_alto)
    
    # crear y guardar grafica
    print(
      log2FC_box_combinado
    )
    dev.off()
  }
  
  
  assign(paste("log2FC_boxplot_agrupado_", condicion, "_", gen_interes, "_", gen_referencia, sep = ""),
         log2FC_box_combinado, 
         envir = .GlobalEnv)
  
  #############################################
  
  # definir maximos en y
  maximo_y <- max(big_data_frame$FoldChange)
  minimo_y <- min(big_data_frame$FoldChange)
  
  # Boxplot para log2FC
  FC_box_combinado <- ggplot(big_data_frame,
                             aes(x = Treatment, 
                                 y = FoldChange, 
                                 fill = Sample)) +
    geom_boxplot(outlier.shape = 21,
                 outlier.fill = "white", 
                 outlier.color = "black") +
    geom_jitter(width = 0.2, 
                size = 2, 
                alpha = 0.5) +
    scale_fill_manual(values = colores_fijos) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          panel.border = element_rect(color = "black", 
                                      fill = NA,
                                      linewidth = 1.2)) +
    ylim(ifelse(minimo_y < 0, minimo_y*1.2, minimo_y * -0.5), maximo_y * 1.2) +
    labs(title = paste(condicion, ": ", gen_interes,
                       " (", gsub("_", " & ", gen_referencia), ")",
                       sep = ""),
         y = "Fold Change",
         x = "",
         fill = "")
  
  print(FC_box_combinado)
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    # crear y guardar los heatmaps
    match.fun(i)(paste(ruta_tablas, 
                       "/FC_boxplot_agrupado_", 
                       condicion, "_", gen_interes, "_", gen_referencia,
                       ".",i, sep = ""),
                 res = resolucion,
                 width = fig_ancho,
                 height = fig_alto)
    
    # crear y guardar grafica
    print(
      FC_box_combinado
    )
    dev.off()
  }
  
  
  assign(paste("FC_boxplot_agrupado_", condicion, "_", gen_interes, "_", gen_referencia, sep = ""),
         FC_box_combinado, 
         envir = .GlobalEnv)
  
}