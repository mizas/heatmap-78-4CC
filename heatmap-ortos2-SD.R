# Cargar las librerias

library(pals)
library(circlize)
library(grid)
library(tibble)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(corrplot)# Clustering
library(cluster) 

# Sección 02. Cargar los resultados de ProteinOrtho

file_list <- list.files(path = "03-proteinortho-graph", pattern = "^file\\d+_cut\\.txt$", full.names = TRUE)
#file_list

# Crear una lista para almacenar los datos de cada archivo
data_list <- list()
#data_list

# Iterar sobre cada archivo en la lista
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t", quote = "")  

  # Tomar los nombres de las dos primeras columnas
  especie <- paste(names(data)[1], names(data)[2], sep = "+")
  
  # Reemplazar "ACARD.faa" y el símbolo "+" por cadenas vacías en la columna "especie"
  data$especie <- gsub("ACARD.faa|\\+", "", especie)
  
  # Almacenar los datos en data_list
  data_list[[file]] <- data
}

# Creamos una función para ordenar los dataframes de la lista
ordenar_dataframe <- function(df) {
  # Seleccionamos todas las columnas excepto la primera y la segunda
  columnas_ordenadas <- c("especie", "ACARD.faa", setdiff(names(df), c("especie", "ACARD.faa")))
  
  # Reordenamos el dataframe según las columnas seleccionadas
  df_ordenado <- df[columnas_ordenadas]
  
  # Devolvemos el dataframe ordenado
  return(df_ordenado)
}

# Aplicamos la función a cada dataframe en data_list
data_list_ordenado <- lapply(data_list, ordenar_dataframe)
#data_list_ordenado[70]


# Cambiar el nombre de la tercera columna por "ortologo"
for (i in seq_along(data_list_ordenado)) {
  names(data_list_ordenado[[i]])[3] <- "ortologo"
}

# Revisar los dataframes
#data_list_ordenado[38]

#Integrar todos los dataframes de "data_list" en un solo dataframe
# Unir todos los dataframes en uno solo
ortologos_graph <- bind_rows(data_list_ordenado)

head(ortologos_graph,50)


#Cargar la matriz de ortologos para crear la matriz de conteos

po.data <- read.csv('./02-proteinortho/myproject.proteinortho.tsv', sep = '\t', header = TRUE)
po.data <- subset(po.data,select = -c(1:3)) #Eliminar las primeras 3 columnas

#####################Parche, corregir en los datos###################################

po.data$ACARD.faa <- gsub("mefF", "mef(F)", po.data$ACARD.faa)
po.data$ACARD.faa <- gsub("mefJ", "mef(J)", po.data$ACARD.faa)

#####################################################################################





po.data.sep <- separate_rows(po.data, ACARD.faa, sep = ",") #convertir a multiples filas las filas que tienen mas de 2 genes
po.data.sep[] <- lapply(po.data.sep, function(x) gsub("\\*", "0", x)) #colocar 0 a todas las celdas que tengan *
po.data.sep<-subset(po.data.sep, po.data.sep$ACARD.faa != 0) #eliminar los ortogolos que no estan en CARD


po.data.sep <- lapply(po.data.sep, as.character)
po.data.sep <- as.data.frame(po.data.sep)

# igualar las columnas de los 2 dataframes, dividir las filas si hay "," en ortologos y anadir la columna de conteo
df_preconteo <-po.data.sep %>%
  pivot_longer(cols = -ACARD.faa,names_to = 'especie', values_to = 'ortologo') %>%
  separate_rows(ortologo, sep = ",") %>%
  mutate(conteo=0) %>%
  dplyr::select(especie,ACARD.faa, ortologo, conteo)


# Elaborar la matriz de conteos

dataframegraph<-ortologos_graph

#####################Parche, corregir en los datos###################################

dataframegraph$ACARD.faa <- gsub("mefF", "mef(F)", dataframegraph$ACARD.faa)
dataframegraph$ACARD.faa <- gsub("mefJ", "mef(J)", dataframegraph$ACARD.faa)

#####################################################################################


dataframemat<-df_preconteo

#####################Parche, corregir en los datos###################################

dataframemat$ACARD.faa <- gsub("mefF", "mef(F)", dataframemat$ACARD.faa)
dataframemat$ACARD.faa <- gsub("mefJ", "mef(J)", dataframemat$ACARD.faa)

#####################################################################################


# La columna ortologo debe ser una cadena de texto, no numérica
dataframegraph$ortologo <- as.character(dataframegraph$ortologo)

# Añadir una columna de conteo al dataframe mat, aunque ya se anadio en el paso anterior
#dataframemat$conteo <- 0

# Con el bucle a continuacion se Itera slobre cada fila de dataframegraph
for (i in 1:nrow(dataframegraph)) {
  # Crear una cadena de texto que representa la fila actual de dataframegraph
  current_row <- paste(dataframegraph[i, "especie"], dataframegraph[i, "ACARD.faa"], dataframegraph[i, "ortologo"], sep = "_")
  # Comprobar si la cadena de texto existe en dataframemat y actualizar el conteo si es así
  if (current_row %in% paste(dataframemat$especie, dataframemat$ACARD.faa, dataframemat$ortologo, sep = "_")) {
    dataframemat[dataframemat$especie == dataframegraph[i, "especie"] & 
                    dataframemat$ACARD.faa == dataframegraph[i, "ACARD.faa"] &
                    dataframemat$ortologo == dataframegraph[i, "ortologo"], "conteo"] <- 1
  }
}




#Arreglar el nombre de las especies
dataframemat$especie <- substring(dataframemat$especie, 1, nchar(dataframemat$especie) - 3)


#Convertir el dataframe a matriz de especies, genes y conteos
dataframemat <- dataframemat[order(dataframemat$especie), ]

dataframemat_wider <- dataframemat %>%
  dplyr::select(-ortologo) %>%
  pivot_wider(names_from = ACARD.faa, values_from = conteo, values_fn = sum, values_fill = 0) 

nombre_filas <- dataframemat_wider$especie

# Convertir el dataframe pivoteado a una matriz
matriz_conteos <- as.matrix(dataframemat_wider[,-1])
rownames(matriz_conteos)<-nombre_filas






### Agregar la informacion de CARD

CARD_db <- read.csv('CARD_database/aro_index.tsv', sep='\t', header = TRUE)
#View(CARD_db)
CARD_db %>%
  dplyr::select(AMR.Gene.Family,CARD.Short.Name,Drug.Class)




#Agregar al dataframe de conteos las familias de genes

dataframemat_fam<-dataframemat
diccionario.1 <- setNames(CARD_db$AMR.Gene.Family,CARD_db$CARD.Short.Name)
dataframemat_fam$family <- diccionario.1[dataframemat$ACARD.faa]

names(dataframemat_fam %>%
  dplyr::select(family, ACARD.faa, especie, ortologo,family,conteo))
dataframemat_fam<-dataframemat_fam %>%
  dplyr::select(family, ACARD.faa, especie, ortologo,family,conteo)

#Convertir el dataframe a matriz de especies, fam_genes y conteos

dataframemat_wider_fam <- dataframemat_fam %>%
  dplyr::select(-ortologo, -ACARD.faa) %>%
  pivot_wider(names_from = family, values_from = conteo, values_fn = sum, values_fill = 0) 

nombre_filas_fam <- dataframemat_wider_fam$especie

# Convertir el dataframe pivoteado a una matriz
matriz_dataframemat_fam <- as.matrix(dataframemat_wider_fam[,-1])
rownames(matriz_dataframemat_fam)<-nombre_filas_fam


#Crear el directorio si no existe y guardar la matriz

# Directory path you want to create
dir_path <- "04-tablas"
# Check if the directory exists
if (!dir.exists(dir_path)) {
  # If it doesn't exist, create it
  dir.create(dir_path)
  cat("Directory", dir_path, "created.\n")
} else {
  cat("Directory", dir_path, "already exists.\n")
}

write.csv(dataframemat, file = "./04-tablas/dataframemat_completa.csv", row.names = TRUE)
write.csv(matriz_dataframemat_fam, file = "./04-tablas/matriz_dataframemat_fam_completa.csv", row.names = TRUE)
write.csv(dataframemat_fam, file = "./04-tablas/dataframemat_fam.csv", row.names = FALSE)


## Eliminar columnas que contengan solo 1

# Función para buscar columnas que tengan solo numeros 1
funcion_unos <- function(columna) {
  length(unique(columna)) == 1}


# Aplicar la función a cada columna de la matriz
columnas_constantes <- apply(matriz_dataframemat_fam, 2, funcion_unos)

# Seleccionar las columnas que cumplen con la condición
columnas_eliminar <- matriz_dataframemat_fam[, columnas_constantes, drop = FALSE]
as.data.frame(columnas_eliminar)
# Generar una nueva matriz para estandarizar
matriz_dataframemat_fam_SD_pre1 <- matriz_dataframemat_fam[, !columnas_constantes]
matriz_dataframemat_fam_SD_pre2 <- scale(matriz_dataframemat_fam_SD_pre1)
matriz_dataframemat_fam_SD_zero <- apply(matriz_dataframemat_fam_SD_pre2, 2, function(column) {column + abs(min(column))
})
matriz_dataframemat_fam_SD_pre2
matriz_dataframemat_fam_SD_zero

write.csv(matriz_dataframemat_fam_SD_zero, file = "./04-tablas/matriz_dataframemat_fam_SD_zero.csv", row.names = TRUE)



# Estimar el número de clusters con Elbow Method, exportar la matriz (no pude instalar factoextra, entonces lo hice en R google colab)
#Elaborar el Heatmap sin anotaciones

heatmap_color<-colorRampPalette(c(  "#46658F", "#183f73","white", "firebrick3","#bf1111"))(50)

#UMPGA por familias de genes

pheatmap_object.upmg <- ComplexHeatmap::pheatmap(
    matriz_dataframemat_fam_SD_zero, #Verificar que matriz se esta cargando
    row_km = 10, 
    column_km = 5,
    clustering_distance_rows = 'euclidean',
    clustering_method = "average",
    color = heatmap_color) 


# Create directory if it doesn't exist
if (!dir.exists("05-figuras")) {
    dir.create("05-figuras")
}

# Save as PNG
png("05-figuras/Heatmap_SA_SD_zero.png", width = 18, height = 13, units = "in", res = 300)
draw(pheatmap_object.upmg)
dev.off()


#Elaborar anotaciones
###Elaborar la anotacion de las filas


dataframemat_fam2 <-dataframemat_fam %>%
  mutate(Organism = especie) %>%
  dplyr::select(-especie)


annotation_row <- dataframemat_fam2 %>% 
  mutate(Phylogeneticgroup = case_when(
    startsWith(Organism, "Bacillus_pumilus_") ~ "Subtilis Clade",
    startsWith(Organism, "Bacillus_atrophaeus_") ~ "Subtilis Clade",
    startsWith(Organism, "Bacillus_altitudinis_") ~ "Subtilis Clade",
    startsWith(Organism, "Bacillus_inaquosorum_") ~ "Subtilis Clade",
    startsWith(Organism, "Bacillus_swezeyi") ~ "Subtilis Clade",
    startsWith(Organism, "Bacillus_thuringiensis") ~ "Core Cereus Clade",
    startsWith(Organism, "Bacillus_cereus") ~ "Core Cereus Clade",
    startsWith(Organism, "Bacillus_mobilis") ~  "Core Cereus Clade",
    startsWith(Organism, "Bacillus_w") ~ "Core Cereus Clade",
    startsWith(Organism, "Bacillus_subtilis") ~ "Subtilis Clade",
    startsWith(Organism, "Bacillus_aquimaris") ~ "Aquimaris Clade",
    startsWith(Organism, "Rossellomorea") ~ "Aquimaris Clade",
    startsWith(Organism, "Jeotgalibacillus") ~ "Jeotgalibacillus sp",
    startsWith(Organism, "Bacillus_infantis") ~ "Bacillus infantis",
    startsWith(Organism, "Cytobacillus") ~ "Cytobacillus gottheilii",
    startsWith(Organism, "Metabacillus") ~ "Metabacillus spp",
    startsWith(Organism, "Sutcli") ~ "Cohnii Clade",
    startsWith(Organism, "Fictibacillus") ~ "Fictibacillus nanhaiensis",
    startsWith(Organism, "Priestia") ~ "Megaterium Clade",
    startsWith(Organism, "Staphylococcus") ~ "Staphylococcus spp",
    startsWith(Organism, "Citricoccus") ~ "Actinomycetales",
    startsWith(Organism, "Micrococcus") ~ "Actinomycetales",
    startsWith(Organism, "Kocuria") ~ "Actinomycetales",
    startsWith(Organism, "Corynebacterium") ~ "Actinomycetales",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(Organism, Phylogeneticgroup) %>%
  unique()

annotation_row_df<-as.data.frame(annotation_row)
rownames(annotation_row_df) <-rownames(matriz_dataframemat_fam)
annotation_row_df <- annotation_row_df %>%
  dplyr::select(Phylogeneticgroup)

### Configurar el color de annotation_row

colors.row <- c("Actinomycetales" = "#e4597a",
               "Aquimaris Clade" = "#4caab6", 
               "Bacillus infantis" = "#ee89ab", 
               "Cohnii Clade" = "#9ccda0",
               "Core Cereus Clade" = "#2095ca",
               "Cytobacillus gottheilii" = "#999999",
               "Fictibacillus nanhaiensis" = "#ba598e",  
               "Jeotgalibacillus sp" = "#381f80", 
               "Megaterium Clade" = "#872154",
               "Metabacillus spp" = "#c3a9d0",
               "Staphylococcus spp" = "#f9d46c", 
               "Subtilis Clade" = "#f08565")

#colorRampPalette(c("gray", "azure" ,"blueviolet","brown","green","blue","red","black"))(12)

#colors()

#UMPGA por familias de genes
pheatmap_object.nuevo <- ComplexHeatmap::pheatmap(
    matriz_dataframemat_fam_SD_zero, 
    row_km = 10, 
    column_km = 5,
    annotation_row = annotation_row_df,
    clustering_distance_rows = 'euclidean',
    clustering_method = "average",
    color = heatmap_color,
    annotation_colors = list(
    Phylogeneticgroup = colors.row)
    ) 

# Save the plot as a PNG file
png("05-figuras/heatmap_AF_SD_zero.png", width = 3000, height = 2000, res = 100)
print(pheatmap_object.nuevo)  # Print the plot
dev.off()




















#Crear diccionario con base de datos de CARD (Genes y Familias)
#Crear tabla de Genes y Familias

annotation_col<- dataframemat_wider %>%
  pivot_longer(cols=-especie, names_to = "Gen", values_to = "conteo",) %>%
  dplyr::select(-conteo, -especie) %>%
  unique()


diccionario.1 <- setNames(CARD_db$AMR.Gene.Family,CARD_db$CARD.Short.Name)
annotation_col$family <- diccionario.1[annotation_col$Gen]

diccionario.2 <- setNames(CARD_db$Resistance.Mechanism,CARD_db$CARD.Short.Name)

annotation_col$RM <- diccionario.2[annotation_col$Gen]
annotation_col <- annotation_col %>%
  separate(RM, into = c("RM", "RM2"), sep = ";", extra = "merge", fill = "right") %>%
  unique()


diccionario.3 <- setNames(CARD_db$Drug.Class,CARD_db$CARD.Short.Name)
annotation_col$DC <- diccionario.3[annotation_col$Gen]

# opcion 1 Crear varias columnas
#annotation_col <-annotation_col %>%
#  separate(DC, into = c("DC","DC2"), sep = ";", extra = "merge", fill = "right")


# Opcion 2 enves de varias columnas utilizar multiresistencia
annotation_col <-annotation_col %>%
  mutate(DC = case_when(
    grepl(";", DC) ~ "Multiresistance",
    TRUE ~ DC
  ))

###################################  Modificar para Estandarizar  #######################################3

#Eliminar la familia que no permite realizar la estandarización, utilizar únicamente para la matriz estandarizada
annotation_col <- subset(annotation_col, family != "sulfonamide resistant sul") #Utilizar solo con estandarizacion

# Elaborar la codificacion de los colores

#Asignar colores para las columnas


colors.RM <- c(
  "antibiotic efflux" = "#DEB877",
  "antibiotic inactivation" = "#ABDD3F",
  "antibiotic target alteration" = "#86F502",
  "antibiotic target protection" = "#EE7621",
  "antibiotic target replacement" = "#A37B6A",
  "antibiotic target alteration_antibiotic efflux" = "#889182",  
  "reduced permeability to antibiotic" = "#00FFFF"
)



colors.RM2 <- c(
  "antibiotic target alteration" = "#86F502",
  "antibiotic target replacement" = "#A37B6A",
  "reduced permeability to antibiotic" = "#00FFFF",
  "NA" = "#FFFFFF"
)

annotation_col %>%
    arrange(DC) %>%
    distinct(DC) %>%
    arrange() %>%
    nrow()


colors.DC <- c(
  "aminocoumarin antibiotic" = "#5050FFFF",
  "aminoglycoside antibiotic" = "#CE3D32FF",
  "antibacterial free fatty acids" = "#749B58FF",
  "bicyclomycin-like antibiotic" = "#F0E685FF",
  "carbapenem" = "#466983FF",
  "cephalosporin" = "#BA6338FF",
  "cephamycin" = "#5DB1DDFF",
  "diaminopyrimidine antibiotic" = "#0A47FFFF",
  "disinfecting agents and antiseptics" = "#6BD76BFF",
  "elfamycin antibiotic" = "#D595A7FF",
  "fluoroquinolone antibiotic" = "#924822FF",
  "fusidane antibiotic" = "#837B8DFF",
  "Multiresistance" = "#C75127FF",
  "glycopeptide antibiotic" = "#D58F5CFF", 
  "lincosamide antibiotic" = "#7A65A5FF", 
  "macrolide antibiotic" = "#991A00FF",
  "mupirocin-like antibiotic" = "#14FFB1FF",
  "nitroimidazole antibiotic" = "#CDDEB7FF",
  "nucleoside antibiotic" = "#AE1F63FF",
  "penam" = "#E7C76FFF",
  "peptide antibiotic" = "#5A655EFF",
  "phenicol antibiotic" = "#1A0099FF",
  "phosphonic acid antibiotic" = "#99CC00FF",
  "pleuromutilin antibiotic"= "#480D60",
  "rifamycin antibiotic" = "#7F9265",
  "sulfonamide antibiotic" = "#A9A9A9FF",
  "tetracycline antibiotic" = "#00991AFF",
  "NA" = "#000000"
)


# Ordenar las filas de annotation_col_agg_df4 con base al orden de col_names_matrix para log
matriz_dataframemat_fam_SD_zero<-as.data.frame(matriz_dataframemat_fam_SD_zero)
col_names_matrix <- colnames(matriz_dataframemat_fam_SD_zero)



#Revisar la matriz y el annotation_col

annotation_col<-annotation_col %>%
  dplyr::select(-Gen)


annotation_col_df<-as.data.frame(annotation_col)

annotation_col_agg <- annotation_col_df %>%
  group_by(family) %>%
  summarize(
    RM = list(unique(RM)),
    RM2 = list(unique(RM2)),
    DC = list(unique(DC))
  )

nombre_filas <-as.data.frame(matriz_dataframemat_fam) %>%
  pivot_longer(cols=1:122,names_to="family",values_to="conteo") %>%
  dplyr::select(family) %>%
  distinct()
nombre_filas <-annotation_col_agg %>%
  dplyr::select(family) %>%
  distinct()

########################################################################
as.data.frame(nombre_filas)

annotation_col_agg_df<-as.data.frame(annotation_col_agg$family)
annotation_col_agg_df2<-as.data.frame(annotation_col_agg)
rownames(annotation_col_agg_df2) <- annotation_col_agg_df$annotation_col_agg

annotation_col_agg_df3<-annotation_col_agg_df2 %>%
  dplyr::select(-family)

#Las columnas RM y RM2 tienen listas con dos elementos, si hay listas con varios elementos entonces las combinaremos en un solo elemento
combinar_elementos <- function(x) {
  if (length(x) > 1) {
    return(paste(x, collapse = "_"))
  } else {
    return(x)
  }
}
annotation_col_agg_df3$RM <- lapply(annotation_col_agg_df3$RM, combinar_elementos)
annotation_col_agg_df3$RM2 <- lapply(annotation_col_agg_df3$RM2, combinar_elementos)


# La columna DC tiene listas con listas en su interior, si hay listas multiples las reduciremos a Multiresistance

multiresistencia_reemplazo <- function(x) {
  if (length(x) > 1) {
    return("Multiresistance")
  } else {
    return(x)
  }
}
annotation_col_agg_df3$DC <- lapply(annotation_col_agg_df3$DC, multiresistencia_reemplazo)

RM<-print(unlist(annotation_col_agg_df3$RM[1:125]))
RM<-as.data.frame(RM)
RM2<-print(unlist(annotation_col_agg_df3$RM2[1:125]))
RM2<-as.data.frame(RM2)

DC<-print(unlist(annotation_col_agg_df3$DC[1:125]))
DC<-as.data.frame(DC)
annotation_col_agg_df4 <- cbind(RM, RM2,DC)

#Asignar nombres a las filas
rownames(annotation_col_agg_df4) <- annotation_col_agg_df$annotation_col_agg

#Ordenar las anotaciones con base a la matriz
annotation_col_agg_df4_ord <- annotation_col_agg_df4[col_names_matrix, , drop = FALSE]
annotation_col_agg_df4_ord
#Corregir valores
annotation_col_agg_df4_ord$RM2 <- gsub("NA_antibiotic target alteration", "antibiotic target alteration", annotation_col_agg_df4_ord$RM2)































# Convert row names of annotation_col_agg_df4 to a character vector
row_names_annotation <- as.character(rownames(annotation_col_agg_df4))

# Extract column names of matriz_dataframemat_fam_log_df
#matriz_dataframemat_fam_log_df<-as.data.frame(matriz_dataframemat_fam_log)
#col_names_matrix <- colnames(matriz_dataframemat_fam_log_df)

# Identify the differences
diff_rownames <- setdiff(row_names_annotation, col_names_matrix)
diff_colnames <- setdiff(col_names_matrix, row_names_annotation)

# Print the differences
print(diff_rownames)
print(diff_colnames)





# Crear la paleta de colores incluyendo negro para cero
#heatmap_color <- colorRampPalette(c("#ffefcd","#fcde9c","#d73027"))(3)
#heatmap_color<-hcl.colors(3,"Zissou 1")
#heatmap_color<-hcl.colors(3,"Plasma")
heatmap_color<-hcl.colors(3,"Geyser")




# Calcular los puntos de corte
min_positive_value <- min(matriz_dataframemat_fam_SD_zero[matriz_dataframemat_fam_SD_zero > 0], na.rm = TRUE)
max_value <- max(matriz_dataframemat_fam_SD_zero, na.rm = TRUE)
intermediate_value <- (min_positive_value + max_value) / 2
min_positive_value
max_value

# Crear la función de color con un valor intermedio
col_fun <- colorRamp2(
  c(0, min_positive_value, intermediate_value, max_value),
  c("#d0d0d0", heatmap_color[1], heatmap_color[2], heatmap_color[3])
)

# Crear anotaciones de fila y columna
row_ha <- rowAnnotation(df = annotation_row_df, col = list(Phylogeneticgroup = colors.row, RM2 = colors.RM2, RM = colors.RM, DC = colors.DC))
col_ha <- HeatmapAnnotation(df = annotation_col_agg_df4_ord, col = list(Phylogeneticgroup = colors.row, RM2 = colors.RM2, RM = colors.RM, DC = colors.DC))

# Crear el heatmap usando ComplexHeatmap con k-means clustering

set.seed(130)
# Convertir a matriz si no lo está
matriz_dataframemat_fam_SD_zero <- as.matrix(matriz_dataframemat_fam_SD_zero)

# Crear el heatmap con su propia leyenda
heatmap <- Heatmap(
  matriz_dataframemat_fam_SD_zero,
  name = "Value",
  row_km = 10,
  column_km = 5,
  col = col_fun,
  top_annotation = col_ha,
  left_annotation = row_ha,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "Value",
    at = c(0, min(matriz_dataframemat_fam_SD_zero[matriz_dataframemat_fam_SD_zero > 0], na.rm = TRUE), max(matriz_dataframemat_fam_SD_zero, na.rm = TRUE)),
    labels = c("Zero", "Min Non-Zero", "Max")
  ),
  row_names_gp = gpar(fontsize = 10),  # Ajusta el tamaño de la letra de los textos de las filas
  column_names_gp = gpar(fontsize = 10)  # Ajusta el tamaño de la letra de los textos de las columnas
)

# Exportar el heatmap a un archivo PNG
png("./05-figuras/heatmap_AFC_SD_zero.png", width = 10500, height = 4500, res = 300)
draw(
  heatmap, 
  heatmap_legend_side = "left",  # Posición de la leyenda del heatmap
  annotation_legend_side = "left",  # Posición de las leyendas de anotaciones
  legend_gap = unit(10, "mm")  # Ajusta la distancia entre la leyenda y el heatmap
)
dev.off()




































## Sección Histograma

dataframemat_fam_hist<-dataframemat_fam
diccionario.2 <- setNames(CARD_db$Resistance.Mechanism,CARD_db$CARD.Short.Name)

dataframemat_fam_hist$RM <- diccionario.2[dataframemat_fam_hist$ACARD.faa]



#dataframemat_fam_hist <- dataframemat_fam_hist %>%
#  separate(RM, into = c("RM", "RM2"), sep = ";", extra = "merge", fill = "right") %>%
#  unique()

diccionario.3 <- setNames(CARD_db$Drug.Class,CARD_db$CARD.Short.Name)
dataframemat_fam_hist$DC <- diccionario.3[dataframemat_fam_hist$ACARD.faa]


dataframemat_fam_hist <-dataframemat_fam_hist %>%
  mutate(DC = case_when(
    grepl(";", DC) ~ "Multiresistance",
    TRUE ~ DC
  ))

dataframemat_fam_hist <- dataframemat_fam_hist %>%
  filter(conteo == 1)

dataframe_fam_hist<-dataframemat_fam_hist %>%
  as.data.frame() %>%
  filter(conteo == 1) %>%
  mutate(Phylogeneticgroup = case_when(
    startsWith(especie, "Bacillus_pumilus_") ~ "Sc", #Subtilis Clade
    startsWith(especie, "Bacillus_atrophaeus_") ~ "Sc",
    startsWith(especie, "Bacillus_altitudinis_") ~ "Sc",
    startsWith(especie, "Bacillus_inaquosorum_") ~ "Sc",
    startsWith(especie, "Bacillus_swezeyi") ~ "Sc",
    startsWith(especie, "Bacillus_thuringiensis") ~ "Ce", #Core Cereus Clade
    startsWith(especie, "Bacillus_cereus") ~ "Ce",
    startsWith(especie, "Bacillus_mobilis") ~  "Ce",
    startsWith(especie, "Bacillus_w") ~ "Ce",
    startsWith(especie, "Bacillus_subtilis") ~ "Sc",
    startsWith(especie, "Bacillus_aquimaris") ~ "Ac",
    startsWith(especie, "Bacillus_endophyticus") ~ "Mc",
    startsWith(especie, "Bacillus_megaterium") ~ "Mc",
    startsWith(especie, "Bacillus_vietnamensis") ~ "Ac",
    startsWith(especie, "Bacillus_marisflavi") ~ "Ac",
    startsWith(especie, "Bacillus_infantis") ~ "Bi", #Bacillus infantis
    startsWith(especie, "Rossellomorea") ~ "Ac", #Aquimaris Clade
    startsWith(especie, "Jeotgalibacillus") ~ "Je", #Jeotgalibacillus sp
    startsWith(especie, "Cytobacillus") ~ "Cg", #Cytobacillus gottheilii
    startsWith(especie, "Metabacillus") ~ "Me", #Metabacillus spp
    startsWith(especie, "Sut") ~ "Cc", #Cohnii Clade
    startsWith(especie, "Fictibacillus") ~ "Fn", #Fictibacillus nanhaiensis
    startsWith(especie, "Priestia") ~ "Mc", #Megaterium Clade
    startsWith(especie, "Staphylococcus") ~ "St",#Staphylococcus spp
    startsWith(especie, "Citricoccus") ~ "At", #Actinomycetales
    startsWith(especie, "Micrococcus") ~ "At",
    startsWith(especie, "Kocuria") ~ "At",
    startsWith(especie, "Corynebacterium") ~ "At",
    TRUE ~ NA_character_
  )) 





dataframe_fam_hist_DC <- dataframe_fam_hist %>%
  dplyr::select(Phylogeneticgroup,DC,conteo) %>%
  group_by(DC,Phylogeneticgroup) %>%
  summarise(total_conteo = sum(conteo)) %>%
  ungroup()

dataframe_fam_hist_RM <- dataframe_fam_hist %>%
  dplyr::select(Phylogeneticgroup,RM,conteo) %>%
  group_by(RM,Phylogeneticgroup) %>%
  summarise(total_conteo = sum(conteo))


### Barplots stacked
# Stacked barplot with multiple groups
barplot_DC <- ggplot(dataframe_fam_hist_DC, aes(x=Phylogeneticgroup, y=total_conteo, fill=DC)) +
#  geom_bar(position="fill", stat="identity") +
  geom_bar(stat="identity")+
  xlab("Phylogenetic group") +
  ylab("Gen frequency") +
  scale_fill_manual(values=c("#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF", "#BA6338FF", "#5DB1DDFF", "#0A47FFFF", "#6BD76BFF", "#D595A7FF", "#924822FF", "#837B8DFF", "#C75127FF", "#D58F5CFF", "#7A65A5FF", "#991A00FF", "#14FFB1FF", "#CDDEB7FF", "#612A79FF", "#AE1F63FF", "#E7C76FFF", "#5A655EFF", "#1A0099FF", "#99CC00FF", "#A9A9A9FF", "#CC9900FF", "#00991AFF"))

barplot_RM <- ggplot(dataframe_fam_hist_RM[,], aes(x=Phylogeneticgroup, y=total_conteo, fill=RM)) +
#  geom_bar(position="fill", stat="identity")+
  geom_bar(stat="identity")+
  xlab("Phylogenetic group") +
  ylab("Gen frequency") +
  scale_fill_manual(values=c("#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF", "#BA6338FF", "#5DB1DDFF", "#0A47FFFF", "#6BD76BFF", "#D595A7FF", "#924822FF", "#837B8DFF", "#C75127FF", "#D58F5CFF", "#7A65A5FF", "#991A00FF", "#14FFB1FF", "#CDDEB7FF", "#612A79FF", "#AE1F63FF", "#E7C76FFF", "#5A655EFF", "#1A0099FF", "#99CC00FF", "#A9A9A9FF", "#CC9900FF", "#00991AFF"))


### Barplots separados en facets
#### pasar los conteos a repeticiones

dataframe_fam_hist_DC_rep <- dataframe_fam_hist_DC %>%
  uncount(total_conteo)
dataframe_fam_hist_RM_rep <- dataframe_fam_hist_RM %>%
  uncount(total_conteo)

barplot_DC<-ggplot(dataframe_fam_hist_DC_rep, aes(x=Phylogeneticgroup, fill=DC)) +
    geom_bar(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), position = "dodge") +
    scale_y_continuous(labels = scales::percent) + 
    facet_wrap(vars(DC), nrow = 1, labeller = label_wrap_gen(width = 8)) + 
    coord_flip() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),  
    axis.text.x = element_text(angle = 270, vjust = 0.5),
    legend.position = "none") +
  labs(
    x = "Phylogenetic group",
    y = "Percentage")


###Agregar etiquetas a Phylogenetic group
####Paso 1


dataframe_fam_conteos<-dataframemat_fam %>%
  as.data.frame() %>%
  filter(conteo == 1) %>%
  mutate(Phylogeneticgroup = case_when(
    startsWith(especie, "Bacillus_pumilus_") ~ "Sc", #Subtilis Clade
    startsWith(especie, "Bacillus_atrophaeus_") ~ "Sc",
    startsWith(especie, "Bacillus_altitudinis_") ~ "Sc",
    startsWith(especie, "Bacillus_inaquosorum_") ~ "Sc",
    startsWith(especie, "Bacillus_swezeyi") ~ "Sc",
    startsWith(especie, "Bacillus_thuringiensis") ~ "Ce", #Core Cereus Clade
    startsWith(especie, "Bacillus_cereus") ~ "Ce",
    startsWith(especie, "Bacillus_mobilis") ~  "Ce",
    startsWith(especie, "Bacillus_w") ~ "Ce",
    startsWith(especie, "Bacillus_subtilis") ~ "Sc",
    startsWith(especie, "Bacillus_aquimaris") ~ "Ac",
    startsWith(especie, "Bacillus_endophyticus") ~ "Mc",
    startsWith(especie, "Bacillus_megaterium") ~ "Mc",
    startsWith(especie, "Bacillus_vietnamensis") ~ "Ac",
    startsWith(especie, "Bacillus_marisflavi") ~ "Ac",
    startsWith(especie, "Bacillus_infantis") ~ "Bi", #Bacillus infantis
    startsWith(especie, "Rossellomorea") ~ "Ac", #Aquimaris Clade
    startsWith(especie, "Jeotgalibacillus") ~ "Je", #Jeotgalibacillus sp
    startsWith(especie, "Cytobacillus") ~ "Cg", #Cytobacillus gottheilii
    startsWith(especie, "Metabacillus") ~ "Me", #Metabacillus spp
    startsWith(especie, "Sut") ~ "Cc", #Cohnii Clade
    startsWith(especie, "Fictibacillus") ~ "Fn", #Fictibacillus nanhaiensis
    startsWith(especie, "Priestia") ~ "Mc", #Megaterium Clade
    startsWith(especie, "Staphylococcus") ~ "St",#Staphylococcus spp
    startsWith(especie, "Citricoccus") ~ "At", #Actinomycetales
    startsWith(especie, "Micrococcus") ~ "At",
    startsWith(especie, "Kocuria") ~ "At",
    startsWith(especie, "Corynebacterium") ~ "At",
    TRUE ~ NA_character_
  )) 

conteo_especies_por_Phylogeneticgroup <- aggregate(especie ~ Phylogeneticgroup, data = dataframe_fam_conteos, FUN = function(x) length(unique(x)))
conteo_ACARD_por_Phylogeneticgroup <- aggregate(ACARD.faa ~ Phylogeneticgroup, data = dataframe_fam_conteos, FUN = function(x) length(unique(x)))
conteo_concat<-cbind(conteo_especies_por_Phylogeneticgroup,conteo_ACARD_por_Phylogeneticgroup)
conteos_juntos <- paste(conteo_especies_por_Phylogeneticgroup$especie, ",", conteo_ACARD_por_Phylogeneticgroup$ACARD.faa, sep = "")

####Paso 2

#### dataframe_fam_hist_RM_rep ####################### 
# Asegurar que haya suficientes etiquetas personalizadas para cada etiqueta de Phylogeneticgroup
# Extraer las etiquetas de Phylogeneticgroup

phylogenetic_labels <- unique(dataframe_fam_hist_RM_rep$Phylogeneticgroup)

#verificar
if (length(conteos_juntos) < length(phylogenetic_labels)) {
  stop("No hay suficientes etiquetas personalizadas para cada etiqueta de Phylogeneticgroup.")
}


# Crear un dataframe para las etiquetas personalizadas
label_df <- data.frame( Phylogeneticgroup = phylogenetic_labels,conteos_juntos)
diccionario.6 <- setNames(conteos_juntos,label_df$Phylogeneticgroup)
dataframe_fam_hist_RM_rep$conteos_juntos <- diccionario.6[dataframe_fam_hist_RM_rep$Phylogeneticgroup]
etiqueta_RM <- paste(dataframe_fam_hist_RM_rep$Phylogeneticgroup, "[", dataframe_fam_hist_RM_rep$conteos_juntos, "]", sep = "")


dataframe_fam_hist_RM_rep <- cbind(dataframe_fam_hist_RM_rep, etiqueta_RM)
colnames(dataframe_fam_hist_RM_rep)[4]<-"etiqueta"

###################################################################################################################################################
#### dataframe_fam_hist_dc_rep ####################### 

dataframe_fam_hist_DC_rep$custom_labels <- diccionario.6[dataframe_fam_hist_DC_rep$Phylogeneticgroup]
etiqueta_DC <- paste(dataframe_fam_hist_DC_rep$Phylogeneticgroup, "[", dataframe_fam_hist_DC_rep$custom_labels, "]", sep = "")
dataframe_fam_hist_DC_rep <- cbind(dataframe_fam_hist_DC_rep, etiqueta_DC)
colnames(dataframe_fam_hist_DC_rep)[4]<-"etiqueta"


#Editar las etiquetas de RM

dataframe_fam_hist_RM_rep <- dataframe_fam_hist_RM_rep %>%
  as.data.frame() %>%
  mutate(RM_label = case_when(
    grepl("^antibiotic efflux$", RM) ~ "AE", #Exactamente 'antibiotic efflux'
    grepl("^antibiotic efflux;reduced permeability to antibiotic$", RM) ~ "AE;RP", #Exactamente 'antibiotic efflux;reduced permeability to antibiotic'
    grepl("^antibiotic target alteration$", RM) ~ "ATA", #Exactamente 'antibiotic target alteration'
    grepl("^antibiotic target protection$", RM) ~ "ATP", #Exactamente 'antibiotic target protection'
    grepl("^reduced permeability to antibiotic$", RM) ~ "RP", #Exactamente 'reduced permeability to antibiotic'
    grepl("^antibiotic efflux;antibiotic target alteration$", RM) ~ "AE;ATA", #Exactamente 'antibiotic efflux;antibiotic target alteration'
    grepl("^antibiotic inactivation$", RM) ~ "AI", #Exactamente 'antibiotic inactivation'
    grepl("^antibiotic target alteration;antibiotic target replacement$", RM) ~ "ATA;ATR", #Exactamente 'antibiotic target alteration;antibiotic target replacement'
    grepl("^antibiotic target replacement$", RM) ~ "ATR", #Exactamente 'antibiotic target replacement'
    TRUE ~ NA_character_
  ))


#Editar las etiquetas de RM

dataframe_fam_hist_DC_rep <- dataframe_fam_hist_DC_rep %>%
  as.data.frame() %>%
  mutate(DC_label = case_when(
    grepl("^Multiresistance$", DC) ~ "MULT", 
    grepl("^antibacterial free fatty acids$", DC) ~ "AFFA", 
    grepl("^cephalosporin$", DC) ~ "CEPH", 
    grepl("^disinfecting agents and antiseptics$", DC) ~ "DAAN", 
    grepl("^fusidane antibiotic$", DC) ~ "FUSD", 
    grepl("^macrolide antibiotic$", DC) ~ "MACR",
    grepl("^nucleoside antibiotic$", DC) ~ "NUCL",
    grepl("^phenicol antibiotic$", DC) ~ "PHEN",
    grepl("^rifamycin antibiotic$", DC) ~ "RIPH",
    grepl("^aminocoumarin antibiotic$", DC) ~ "AMIN",
    grepl("^bicyclomycin-like antibiotic$", DC) ~ "BICY",
    grepl("^cephamycin$", DC) ~ "CEMY", 
    grepl("^elfamycin antibiotic$", DC) ~ "ELFA", 
    grepl("^glycopeptide antibiotic$", DC) ~ "GLYC", 
    grepl("^mupirocin-like antibiotic$", DC) ~ "MUPI", 
    grepl("^penam$", DC) ~ "PENA", 
    grepl("^phosphonic acid antibiotic$", DC) ~ "PHAC",
    grepl("^sulfonamide antibiotic$", DC) ~ "SULF",
    grepl("^aminoglycoside antibiotic$", DC) ~ "AMGL", 
    grepl("^carbapenem$", DC) ~ "CARB",
    grepl("^diaminopyrimidine antibiotic$", DC) ~ "DIAM", 
    grepl("^fluoroquinolone antibiotic$", DC) ~ "FLUO", 
    grepl("^lincosamide antibiotic$", DC) ~ "LINC", 
    grepl("^nitroimidazole antibiotic$", DC) ~ "NIIM",
    grepl("^peptide antibiotic$", DC) ~ "PEPT",
    grepl("^pleuromutilin antibiotic$", DC) ~ "PLEU", 
    grepl("^tetracycline antibiotic$", DC) ~ "TETR",
    TRUE ~ NA_character_
  ))


# Crear el gráfico para RM
barplot_RM.y<- ggplot(dataframe_fam_hist_RM_rep, aes(x = etiqueta, fill = RM)) +
  geom_bar(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), position = "dodge") +
  scale_y_continuous(labels = scales::percent) + 
  facet_wrap(vars(RM_label), nrow = 1, labeller = label_wrap_gen(width = 16)) + 
  coord_flip() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),  
    axis.text.x = element_text(angle = 270, vjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 50, 10, 10)
  ) +
  labs(
    x = "Phylogenetic group",
    y = "Percentage"
  ) 

barplot_RM.x<- ggplot(dataframe_fam_hist_RM_rep, aes(x = RM_label, fill = etiqueta)) +
  geom_bar(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), position = "dodge") +
  scale_y_continuous(labels = scales::percent) + 
  facet_wrap(vars(etiqueta), nrow = 1, labeller = label_wrap_gen(width = 16)) + 
  coord_flip() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),  
    axis.text.x = element_text(angle = 270, vjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 50, 10, 10)
  ) +
  labs(
    x = "Resistance Mechanism",
    y = "Percentage"
  ) 

# Crear el gráfico para DC
barplot_DC.y <- ggplot(dataframe_fam_hist_DC_rep, aes(x = etiqueta, fill = DC)) +
  geom_bar(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), position = "dodge") +
  scale_y_continuous(labels = scales::percent) + 
  facet_wrap(vars(DC_label), nrow = 1, labeller = label_wrap_gen(width = 16)) + 
  coord_flip() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),  
    axis.text.x = element_text(angle = 270, vjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 50, 10, 10)
  ) +
  labs(
    x = "Phylogenetic group",
    y = "Percentage"
  ) 

barplot_DC.x<- ggplot(dataframe_fam_hist_DC_rep, aes(x = DC_label, fill = etiqueta)) +
  geom_bar(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), position = "dodge") +
  scale_y_continuous(labels = scales::percent) + 
  facet_wrap(vars(etiqueta), nrow = 1, labeller = label_wrap_gen(width = 16)) + 
  coord_flip() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),  
    axis.text.x = element_text(angle = 270, vjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 50, 10, 10)
  ) +
  labs(
    x = "Drug Class",
    y = "Percentage"
  )


# Save the plot to a file
ggsave(filename = "05-figuras/barplot_DCy.png", plot = barplot_DC.y, width = 30, height = 18.75, dpi = 300,limitsize = FALSE)
ggsave(filename = "05-figuras/barplot_DCx.png", plot = barplot_DC.x, width = 30, height = 18.75, dpi = 300,limitsize = FALSE)
ggsave(filename = "05-figuras/barplot_RMy.png", plot = barplot_RM.y, width = 24, height = 15, dpi = 300)
ggsave(filename = "05-figuras/barplot_RMx.png", plot = barplot_RM.x, width = 24, height = 15, dpi = 300)


### Estadísticos para la sección de resultados


dataframemat_stats<-dataframemat

diccionario.1 <- setNames(CARD_db$AMR.Gene.Family,CARD_db$CARD.Short.Name)
dataframemat_stats$family <- diccionario.1[dataframemat$ACARD.faa]


diccionario.2 <- setNames(CARD_db$Resistance.Mechanism,CARD_db$CARD.Short.Name)
dataframemat_stats$RM <- diccionario.2[dataframemat$ACARD.faa]

diccionario.3 <- setNames(CARD_db$Drug.Class,CARD_db$CARD.Short.Name)
dataframemat_stats$DC <- diccionario.3[dataframemat$ACARD.faa]

dataframemat_stats$conteo<-as.numeric(dataframemat_stats$conteo)

dataframemat_stats %>%
  dplyr::select(especie, ACARD.faa, family, RM, DC, conteo, ortologo) %>%
  filter(conteo > 0)


dataframemat_stats %>%
  dplyr::select(especie, ACARD.faa, family, RM, DC, conteo, ortologo) %>%
  filter(conteo > 0) %>%
  dplyr::select(RM, conteo) %>%
  group_by(RM) %>%
  summarise(total=sum(conteo)) %>%  
  ungroup() %>%
  mutate(percentage = total / sum(total) * 100)

dataframemat_stats %>%
  dplyr::select(especie, ACARD.faa, family, RM, DC, conteo, ortologo) %>%
  separate_rows(DC, sep = ";") %>%
  filter(conteo > 0) %>%
  dplyr::select(DC, conteo) %>%
  group_by(DC) %>%
  summarise(total=sum(conteo)) %>%  
  ungroup() %>%
  mutate(percentage = total / sum(total) * 100) %>%
  arrange(-percentage)


dataframemat_stats %>%
  dplyr::select(especie, ACARD.faa, family, RM, DC, conteo, ortologo) %>%
  separate_rows(DC, sep = ";") %>%
  filter(conteo > 0) %>%
  group_by(DC, especie) %>%
  summarise(
    RM = paste(unique(RM), collapse = ";"), # Concatenar los valores únicos de RM separados por ;
    total = sum(conteo)
  ) %>%
  ungroup() %>%
  mutate(percentage = total / sum(total) * 100) %>%
  arrange(especie, -percentage)
