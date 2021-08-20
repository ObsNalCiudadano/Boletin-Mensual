####email

###focos rojos

library(tidyverse)
library(lubridate)
library(formattable)
library(png)
library(grid)
library(Hmisc)
library(gridExtra)
library(extrafont)
library(mxmaps)
library(BAMMtools)
library(RColorBrewer)
library(RCurl)

setwd("C:/Users/emcdo/OneDrive/Documents/ONC/ROBOT")

print('starting')

# Usuaria:
usi <- "C:/Users/emcdo"

# Directorios

input <- paste(usi, "Google Drive/big_carpet/Bases Generales/tidy", sep = "/")
input_logos_delitos <- paste(usi, "Google Drive/big_carpet/Disenio/Logos redondos", sep = "/") 
input_im <- paste(usi, "OneDrive/Documents/ONC/ROBOT", sep = "/")
out <- paste(usi, "OneDrive/Documents/ONC/ROBOT", sep = "/")

# Bases

# b_est_t <- read_csv(paste(input, "TasasEstatalMasterTidy.csv", sep = "/"))
# b_est_a <- read_csv(paste(input, "AbsEstatalMasterTidy.csv", sep = "/"))


#################################################################################################
# Creaci?n de objetos para que el mes se seleccione en autom?tico:

if(day(Sys.Date()) >= 20){
  num_mes_boletin <- month(today()) - 1
} else {   num_mes_boletin <- month(today()) - 2 } 

if(month(today()) == 1){
  num_mes_boletin <- 12 + num_mes_boletin
  
}

#temprorary month for testing
 # num_mes_boletin <- 2



if(month(today())== 1){
  
  ao_actual <- year(today())-1
} else{
  ao_actual <- year(today())
}


meses_1 <- c('Enero' = 1, 'Febrero' = 2, 'Marzo' = 3, 'Abril' = 4, 'Mayo' = 5, 'Junio' = 6, 
             'Julio' = 7, 'Agosto' = 8, 'Septiembre' = 9, 'Octubre' = 10, 'Noviembre' = 11, 
             'Diciembre' = 12)


mes_abbr <- c('ene', 'feb', 'mar', 'abr', 'may', 'jun',
              'jul', 'ago', 'sept', 'oct', 'nov', 'dic')



cat_mes <- data.frame(nom_mes = c("Enero", "Febrero", "Marzo",
                                  "Abril", "Mayo", "Junio",
                                  "Julio", "Agosto", "Septiembre",
                                  "Octubre", "Noviembre", "Diciembre"),
                      num_mes = c(1:12))


nom_mes_boletin <- as.character(cat_mes[num_mes_boletin, 1])


# Base codigos delitos
edo_code <- read_csv(paste(input_im, "estados.csv", sep = "/"))


###########subir objetos de obs locales
# C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/obs_locales1.png
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/obs_locales1.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/obs_locales1_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)


Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/obs_locales2.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/obs_locales2_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)


Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/obs_locales3.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/obs_locales3_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)


Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/obs_locales4.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/obs_locales4_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)




img_hd <- readPNG(paste(input_logos_delitos, "HD.png", sep = "/"))
img_hc <- readPNG(paste(input_logos_delitos, "HC.png", sep = "/"))
img_fem <- readPNG(paste(input_logos_delitos, "FEM.png", sep = "/"))
img_ext <- readPNG(paste(input_logos_delitos, "EXT.png", sep = "/"))
img_sec <- readPNG(paste(input_logos_delitos, "SEC.png", sep = "/"))
img_rvi <- readPNG(paste(input_logos_delitos, "RVI.png", sep = "/"))
img_rve <- readPNG(paste(input_logos_delitos, "RVE.png", sep = "/"))
img_rc <- readPNG(paste(input_logos_delitos, "RC.png", sep = "/"))
img_rn <- readPNG(paste(input_logos_delitos, "RN.png", sep = "/"))
img_rt <- readPNG(paste(input_logos_delitos, "RT.png", sep = "/"))
img_rtp <- readPNG(paste(input_logos_delitos, "RTP.png", sep = "/"))
img_vio <- readPNG(paste(input_logos_delitos, "VIO.png", sep = "/"))
img_vf <- readPNG(paste(input_logos_delitos, "VF.png", sep = "/"))
img_tra <- readPNG(paste(input_logos_delitos, "TRA.png", sep = "/"))
img_nar <- readPNG(paste(input_logos_delitos, "NAR.png", sep = "/"))
img_ld <- readPNG(paste(input_logos_delitos, "LD.png", sep = "/"))

img_vacio <- readPNG(paste(input_logos_delitos, "vacio.png", sep = "/"))

hd <- rasterGrob(img_hd)
fem <- rasterGrob(img_fem)
hc <- rasterGrob(img_hc)
sec <- rasterGrob(img_sec)
ext <- rasterGrob(img_ext)
rvi <- rasterGrob(img_rvi)
rve <- rasterGrob(img_rve)
rc <- rasterGrob(img_rc)
rn <- rasterGrob(img_rn)
rt <- rasterGrob(img_rt)
rtp <- rasterGrob(img_rtp)
vio <- rasterGrob(img_vio)
vf <- rasterGrob(img_vf)
tra <- rasterGrob(img_tra)
nar <- rasterGrob(img_nar)
ld <- rasterGrob(img_ld)

vacio <- rasterGrob(img_vacio)




##########subir los demas imagenes
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/ultimo_impacto.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/ultimo_impacto_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)


Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/ultimo_estudio.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/ultimo_estudio_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)


Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/ultimo_francisco.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/ultimo_francisco_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)

Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/ultimo_investigadores.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/ultimo_investigadores_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)

Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/ultimo_podcast.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/podcast_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)

Sys.sleep(1)
ftpUpload("C:/Users/emcdo/Google Drive/Boletin_mensual/Imagenes/ultimo_video.png",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/ultimo_video_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)



#### OJO HAY QUE QUITAR ESTOOOOO PARA CORRERLO:
#nom_mes_boletin <- "Agosto"               #####
###############################################

##### Se crea una carpeta con el nombre del mes para los objetos y se crea un directorio
###dir.create(paste0(out, "/Objetos_", nom_mes_boletin))
#output <- paste0(out, "/Objetos_", nom_mes_boletin)


periodo_actual <- paste0(nom_mes_boletin, " ", ao_actual)
periodo_anterior <- paste0(nom_mes_boletin, " ", ao_actual - 1)


# Delitos
delitos <- c("Homicidio doloso", "Feminicidio", "Homicidio culposo", "Secuestro", "Extorsión",
             "Robo con violencia", "Robo de vehículo", "Robo a casa habitación", "Robo a negocio",
             "Robo a transeúnte total", "Robo en transporte público",
             "Violación","Violencia familiar", "Trata de personas", "Narcomenudeo", "Lesiones dolosas")

# Bases del periodo actual

if(day(Sys.Date()) >20) {
fech_actual <- floor_date(Sys.Date(), '1 month') -1
} else {
  fech_actual <- floor_date(Sys.Date(), '1 month')  %m-% months(1) -1
}


b_est_t <- getONC('est', 'ci_tasa', 1) %>% 
  select(state_code = inegiEntidad, Entidad = inegiEntidadName,
         Crime = typeCrimeName, Fecha = date, Valor = ci_tasa) %>% 
  mutate(Periodo_ = paste(capitalize(as.character(month(Fecha, label = TRUE, abbr = FALSE))),
                          year(Fecha)),
         ao = year(Fecha))

b_est_a <- getONC('est', 'ci_abs', 1) %>% 
  select(state_code = inegiEntidad, Entidad = inegiEntidadName,
         Crime = typeCrimeName, Fecha = date, Valor = ci_abs) %>% 
  mutate(Periodo_ = paste(capitalize(as.character(month(Fecha, label = TRUE, abbr = FALSE))),
                          year(Fecha)),
         ao = year(Fecha))


b_tasa <- b_est_t %>% filter(Fecha == fech_actual) %>% 
  select(state_code, Entidad, Crime, Fecha, Periodo_, Valor)
b_abs <- b_est_a %>% filter(Fecha == fech_actual) %>% 
  select(state_code, Entidad, Crime, Fecha, Periodo_, Valor)


# Base para las gr?ficas sin base:

basesita <- data.frame(x1 = c(1:6),
                       x2 = c(1:6),
                       x3 = c(1:6),
                       x4 = c(1:6),
                       x5 = c(1:6),
                       x6 = c(1:6))



b_cambio_porc <-  b_est_t %>% filter(Periodo_ %in% c(periodo_actual, periodo_anterior), 
                                     state_code != 0) %>% 
  select(Crime, state_code, Periodo_, ao, Valor) %>% 
  ungroup() %>% 
  arrange(Crime, state_code, ao) %>% 
  group_by(Crime, state_code) %>% 
  mutate(cambio_porc = round(((Valor - lag(Valor))/lag(Valor))*100, digits = 2)) %>% 
  filter(!is.na(cambio_porc)) %>% 
  left_join(edo_code, by = c("state_code" = "CVE_INEGI")) %>% 
  ungroup() %>% 
  select(Crime, Estados= state_code.y, Periodo_, ao, Valor, cambio_porc)




base_o <- data.frame(Crime = delitos,
                     Indice = c(1: length(delitos)),
                     importi = c(3, 4, 2, 3, 2, 2, 2, 1, 1, 1, 1, 3, 2, 3, 1, 1)) %>% 
  mutate(Crime_t = as.character(Crime),
         Crime = ifelse(Crime_t == "Robo a transeúnte total", "Robo a transeúnte", Crime_t))


b_fr_edos <- b_cambio_porc %>% 
  filter(!is.infinite(cambio_porc),
         Crime %in% delitos) %>% 
  group_by(Crime) %>% 
  top_n(1) %>% ungroup() %>% 
  left_join(base_o, by = c("Crime" = "Crime_t"))

b_peso <- b_fr_edos %>% 
  count(Estados, wt = importi, sort = T)
b_ocu <- b_fr_edos %>% 
  count(Estados)

b_fr <- b_peso %>% full_join(b_ocu, by = "Estados") %>% 
  mutate(sele = n.x*n.y) %>% 
  arrange(desc(sele))

edo_fr1 <- as.character(b_fr$Estados[1])
edo_fr2 <- as.character(b_fr$Estados[2])
edo_fr3 <- as.character(b_fr$Estados[3])
edo_fr4 <- as.character(b_fr$Estados[4])
edo_fr5 <- as.character(b_fr$Estados[5])

edos_selec <- c(edo_fr1, edo_fr2, edo_fr3, edo_fr4, edo_fr5)


lista_logos <- list(hd, fem, hc, sec, ext, 
                    rvi, rve, rc, rn, rt, rtp,
                    vio, vf, tra, nar, ld)



b_sele_logos <- b_fr_edos %>%
  filter(Estados %in% edos_selec) %>% 
  select(Crime, Estados, importi, Indice) 


l_sele_logo <- list()

for(i in 1:length(edos_selec)) {
  tempo <- b_sele_logos %>% filter(Estados == edos_selec[[i]]) %>% arrange(desc(importi)) 
  
  if(nrow(tempo) >= 3) {
    l_sele_logo[[i]] <- tempo %>% head(n = 3)
  } else {
    tempo_2 <- data.frame(Crime = c(rep(NA, 3 - nrow(tempo))),
                          Estados = c(rep(edos_selec[[i]], 3 - nrow(tempo))),
                          importi = c(rep(NA, 3 - nrow(tempo))),
                          Indice = c(rep(NA, 3 - nrow(tempo))))
    l_sele_logo[[i]] <- tempo %>% bind_rows(tempo_2)
    
  }
}


b_selec <- l_sele_logo %>% reduce(bind_rows) 

b_o <- b_selec %>% mutate(hay = ifelse(is.na(Indice), 0,1)) %>%
  group_by(Estados) %>% 
  summarise(orden = sum(hay))

b_selec_f <- b_selec %>% left_join(b_o, by = "Estados") %>% 
  arrange(desc(orden))


logis <- list()

for(i in 1:15) {
  if(is.na(b_selec_f$Indice[i])) {
    logis[[i]] <- vacio
  } else {
    logis[[i]] <- lista_logos[[as.numeric(b_selec_f$Indice[i])]]
  }
}


edo_o1 <- as.character(b_selec_f$Estados[1])
edo_o2 <- as.character(b_selec_f$Estados[4])
edo_o3 <- as.character(b_selec_f$Estados[7])
edo_o4 <- as.character(b_selec_f$Estados[10])
edo_o5 <- as.character(b_selec_f$Estados[13])


# Elementos gr?ficos
img_tablita <- readPNG(paste(input_im, "tablita.png", sep = "/")) 
tab <- rasterGrob(img_tablita)

size_edo <- 7
x_nom_edo <- 1.2
x_min_izq <- 1.87
size_logo <- 0.35

# Col 1 
xmin1 <- 2.94
xmax1 <- xmin1 + size_logo

ymin1 <- 3.4
ymax1 <- ymin1 + size_logo
ymin2 <- 3.05
ymax2 <- ymin2 + size_logo
ymin3 <- 2.73
ymax3 <- ymin3 + size_logo
ymin4 <- 2.41
ymax4 <- ymin4 + size_logo
ymin5 <- 2.08
ymax5 <- ymin5 + 0.35

# Col 2
xmin2 <- 3.65
xmax2 <- xmin2 + size_logo

# Col 3
xmin3 <- 4.28
xmax3 <- xmin3 + 0.35

ggplot(basesita) +
  annotation_custom(tab, xmin = 1, xmax = 5, ymin = 2, ymax = 4) +
  coord_cartesian(xlim = c(1, 5), ylim = c(2,4)) +
  # Encabezado de la tabla
  geom_text(aes(x = 1.43, y = 3.83, label = "Entidad"), 
            size = 6, family = "Fira Sans Medium", color = "white", fontface = "bold") +
  geom_text(aes(x = 3.1, y = 3.83, label = "Delitos"), 
            size = 6.1, family = "Fira Sans Medium", color = "white", fontface = "bold") +
  # Los nombres de los estados.
  geom_text(aes(x = x_nom_edo, y = 3.55, label = edo_o1), 
            family = "Fira Sans Medium", size = size_edo, hjust = 0) +
  geom_text(aes(x = x_nom_edo, y = 3.23, label = edo_o2), 
            family = "Fira Sans Medium", size = size_edo, hjust = 0) +
  geom_text(aes(x = x_nom_edo, y = 2.9, label = edo_o3), 
            family = "Fira Sans Medium", size = size_edo, hjust = 0) +
  geom_text(aes(x = x_nom_edo, y = 2.55, label = edo_o4), 
            family = "Fira Sans Medium", size = size_edo, hjust = 0) +
  geom_text(aes(x = x_nom_edo, y = 2.25, label = edo_o5), 
            family = "Fira Sans Medium", size = size_edo, hjust = 0) +
  # Los logos de los delitos 
  # Columna 1
  annotation_custom(logis[[1]], xmin = xmin1, xmax = xmax1, ymin = ymin1, ymax = ymax1) +
  annotation_custom(logis[[4]], xmin = xmin1, xmax = xmax1, ymin = ymin2, ymax = ymax2) +
  annotation_custom(logis[[7]], xmin = xmin1, xmax = xmax1, ymin = ymin3, ymax = ymax3) +
  annotation_custom(logis[[10]], xmin = xmin1, xmax = xmax1, ymin = ymin4, ymax = ymax4) +
  annotation_custom(logis[[13]], xmin = xmin1, xmax = xmax1, ymin = ymin5, ymax = ymax5) +
  #Columa 2
  annotation_custom(logis[[2]], xmin = xmin2, xmax = xmax2, ymin = ymin1, ymax = ymax1) +
  annotation_custom(logis[[5]], xmin = xmin2, xmax = xmax2, ymin = ymin2, ymax = ymax2) +
  annotation_custom(logis[[8]], xmin = xmin2, xmax = xmax2, ymin = ymin3, ymax = ymax3) +
  annotation_custom(logis[[11]], xmin = xmin2, xmax = xmax2, ymin = ymin4, ymax = ymax4) +
  annotation_custom(logis[[14]], xmin = xmin2, xmax = xmax2, ymin = ymin5, ymax = ymax5) +
  # Columna 3
  annotation_custom(logis[[3]], xmin = xmin3, xmax = xmax3, ymin = ymin1, ymax = ymax1) +
  annotation_custom(logis[[6]], xmin = xmin3, xmax = xmax3, ymin = ymin2, ymax = ymax2) +
  annotation_custom(logis[[9]], xmin = xmin3, xmax = xmax3, ymin = ymin3, ymax = ymax3) +
  annotation_custom(logis[[12]], xmin = xmin3, xmax = xmax3, ymin = ymin4, ymax = ymax4) +
  annotation_custom(logis[[15]], xmin = xmin3, xmax = xmax3, ymin = ymin5, ymax = ymax5) +
  theme_void() +
  labs(#title = "Focos rojos",
    #subtitle = "Las entidades de México que presentaron aumentos \n importantes* y que podrían convertirse en zonas de riesgo son... \n",
    caption = "Para más información, se puede consultar la Plataforma y la metodología") +
  theme(#plot.title = element_text(family = "Arial Narrow", face = "bold", size = 1, color = 'white'),
    #plot.subtitle = element_text(family = "Arial Narrow", size = 1, color = 'white'),
    #plot.margin = unit(c(rep(1,4)), "cm"),
    plot.caption = element_text(family = "Fira Sans Light", size = 16, color = "gray70", hjust = 0.26)) 

ggsave(paste(out, "focos_rojos.png", sep = "/"), width = 10, height = 6)

ftpUpload("focos_rojos.PNG",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/focos_rojos_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)

print('paso1')


####### Mapas ###########

data("df_mxstate")

polis <- mxstate.map





# Mapa homicidio doloso

b_hd <- b_tasa %>% filter(state_code != 0,
                          Crime == "Homicidio doloso") %>% 
  mutate(id = ifelse(str_length(as.character(state_code)) == 1,
                     paste0("0", state_code), state_code)) 

n_breaks_hd <- getJenksBreaks(b_hd$Valor, 6)

df_breaks <- data.frame(infe_1 = round(c(0, n_breaks_hd[-6]), digits = 2),
                        super_1 = round(c(n_breaks_hd - 0.01), digits = 2)) %>% 
  mutate(infe = ifelse(str_count(infe_1) == 3, paste0(infe_1, "0"), infe_1),
         super = ifelse(str_count(super_1) == 3, paste0(super_1, "0"), super_1))

labels_rango <- c(paste0(df_breaks$infe[1], " - ",df_breaks$super[1]),
                  paste0(df_breaks$infe[2], " - ",df_breaks$super[2]),
                  paste0(df_breaks$infe[3], " - ",df_breaks$super[3]),
                  paste0(df_breaks$infe[4], " - ",df_breaks$super[4]),
                  paste0(df_breaks$infe[5], " - ",df_breaks$super[5]),
                  paste0(df_breaks$infe[6], " - ",df_breaks$super[6]))


#levels_rango <- unique(levels_rango)

b_mapa_hd <-  b_hd  %>%
  mutate(rango = ifelse(Valor < n_breaks_hd[1], labels_rango[1],
                        ifelse(between(Valor,n_breaks_hd[1], n_breaks_hd[2]), labels_rango[2],
                               ifelse(between(Valor, n_breaks_hd[2], n_breaks_hd[3]), labels_rango[3],
                                      ifelse(between(Valor, n_breaks_hd[3], n_breaks_hd[4]), labels_rango[4],
                                             ifelse(between(Valor, n_breaks_hd[4],  n_breaks_hd[5]),
                                                    labels_rango[5], labels_rango[6]))))),
         rango_fac = factor(rango, levels = rev(labels_rango))) %>% 
  right_join(polis, by = "id") %>% 
  arrange(order) 


col_hd0 <- "f1f1f1"
col_hd1 <- "#d7bec6"
col_hd2 <- "#bc8c9c"
col_hd3 <- "#9f5b74"

col_hd <- "#80284f"

col_hd5 <- "#5f213c"
col_hd6 <- "#401a29"


b_mapa_hd %>% 
  ggplot(aes(x = long, y = lat, group = group, fill = rango_fac)) +
  geom_polygon(size = 0.3, color = "gray10") +
  scale_fill_manual(values = rev(c(col_hd0, col_hd1, col_hd2, col_hd3, col_hd, col_hd5))) +
  coord_map() +
  theme_void() +
  labs(title = paste0("Homicidio doloso\n  ", names(meses_1[num_mes_boletin]), ' 2021'),
       subtitle = "\n \n \n",
       fill = "Tasa por cada \n  100,000 \n  habitantes") +
  theme(legend.position = c(-0.055,0.77),
        #legend.background = element_rect(fill = NULL),
        legend.title = element_text(family = "Fira Sans Medium", size = 15, hjust = -0.2),
        legend.text = element_text(family = "Fira Sans Medium", size = 15, color = "gray30"),
        plot.title = element_text(family = "Fira Sans Medium", face = "bold", size = 22, hjust = -0.22, color = 'Black',
                                  vjust = -.3))



ggsave(paste(out, "mapa_hd.png", sep = "/"), width = 10, height = 6)


ftpUpload("mapa_hd.PNG",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/mapa_hd_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)



print('paso2')



# Mapa robo con violencia

b_rvi <- b_tasa %>% filter(state_code != 0,
                           Crime == "Robo con violencia") %>% 
  mutate(id = ifelse(str_length(as.character(state_code)) == 1,
                     paste0("0", state_code), state_code)) 

n_breaks_rvi <- getJenksBreaks(b_rvi$Valor, 6)

df_breaks <- data.frame(infe_1 = round(c(0, n_breaks_rvi[-6]), digits = 2),
                        super_1 = round(c(n_breaks_rvi - 0.01), digits = 2)) %>% 
  mutate(infe = ifelse(str_count(infe_1) == 3, paste0(infe_1, "0"), infe_1),
         super = ifelse(str_count(super_1) == 3, paste0(super_1, "0"), super_1))

labels_rango <- c(paste0(df_breaks$infe[1], " - ",df_breaks$super[1]),
                  paste0(df_breaks$infe[2], " - ",df_breaks$super[2]),
                  paste0(df_breaks$infe[3], " - ",df_breaks$super[3]),
                  paste0(df_breaks$infe[4], " - ",df_breaks$super[4]),
                  paste0(df_breaks$infe[5], " - ",df_breaks$super[5]),
                  paste0(df_breaks$infe[6], " - ",df_breaks$super[6]))




b_mapa_rvi <-  b_rvi  %>%
  mutate(rango = ifelse(Valor < n_breaks_rvi[1], labels_rango[1],
                        ifelse(between(Valor,n_breaks_rvi[1], n_breaks_rvi[2]), labels_rango[2],
                               ifelse(between(Valor, n_breaks_rvi[2], n_breaks_rvi[3]), labels_rango[3],
                                      ifelse(between(Valor, n_breaks_rvi[3], n_breaks_rvi[4]), labels_rango[4],
                                             ifelse(between(Valor, n_breaks_rvi[4],  n_breaks_hd[5]),
                                                    labels_rango[5], labels_rango[6]))))),
         rango_fac = factor(rango, levels = rev(labels_rango))) %>% 
  right_join(polis, by = "id") %>% 
  arrange(order) 


col_rvi1 <- "#b9b9c4"
col_rvi2 <- "#838499"
col_rvi3 <- "#505370"



col_rvi <- "#1E2749"


col_rvi5 <- "#141727"
col_rvi6 <- "#401a29"


b_mapa_rvi %>% 
  ggplot(aes(x = long, y = lat, group = group, fill = rango_fac)) +
  geom_polygon(size = 0.3, color = "gray10") +
  scale_fill_manual(values = rev(c("#f1f1f1",col_rvi1, col_rvi2, col_rvi3, col_rvi, col_rvi5))) +
  coord_map() +
  theme_void() +
  labs(title = paste0("Robo con violencia\n   ", names(meses_1[num_mes_boletin]), ' 2021'),
       subtitle = "\n \n \n",
       fill = "Tasa por cada \n  100,000 \n  habitantes") +
  theme(legend.position = c(-0.055,0.77),
        #legend.background = element_rect(fill = NULL),
        legend.title = element_text(family = "Fira Sans Medium", size = 15, hjust = -0.2),
        legend.text = element_text(family = "Fira Sans Medium", size = 15, color = "gray30"),
        plot.title = element_text(family = "Fira Sans Medium", face = "bold", size = 22, hjust = -0.25, color = 'Black',
                                  vjust = -.3))


ggsave(paste(out, "mapa_rvi.png", sep = "/"), width = 10, height = 6)

ftpUpload("mapa_rvi.PNG",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/mapa_rvi_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)

#############Dia promedio


base_o <- data.frame(Crime = delitos,
                     Indice = c(1: length(delitos)),
                     importi = c(3, 4, 2, 3, 2, 2, 2, 1, 1, 1, 1, 3, 2, 3, 1, 1)) %>% 
  mutate(Crime_t = as.character(Crime),
         Crime = ifelse(Crime_t == "Robo a transeúnte total", "Robo a transeúnte", Crime_t))




base_o_dia <- base_o# %>% filter(Crime != "Robo en transporte público")
delitos_dia <- c("Homicidio doloso", "Feminicidio", "Homicidio culposo", "Secuestro", "Extorsión",
                 "Robo con violencia", "Robo de vehículo", "Robo a casa habitación", "Robo a negocio",
                 "Robo a transeúnte total", 'Robo en transporte público',
                 "Violación","Violencia familiar", "Trata de personas", "Narcomenudeo", "Lesiones dolosas")

b_dia <- b_abs %>% filter(state_code == 0,
                          Crime %in% delitos_dia) %>% 
  mutate(delito_dia = round(Valor/ days_in_month(num_mes_boletin), digits = 0)) %>% 
  mutate(Crime = ifelse(Crime == 'Robo a transeúnte total', 'Robo a transeúnte', Crime)) %>% 
  left_join(base_o_dia, by = "Crime") %>% 
  arrange(Indice) %>% 
  mutate(x_graf = c(rep(.955, 9), rep(1.955, 7)),
         y_graf = c(9.5, 8:1, 8:2),
         Crime = ifelse(Crime == "Robo a transeúnte total", "Robo a transeúnte", Crime),
         etiqueta = paste0(delito_dia, " de ", str_to_lower(Crime)),
         etiqueta = ifelse(Crime == "Homicidio doloso", paste0(delito_dia, " carpetas de investigación de homicidio doloso"), etiqueta)) 

ggplot(b_dia) +
  geom_text(aes(x = x_graf, y = y_graf, label = etiqueta),
            hjust = 0, family = "Fira Sans Medium", size = 6) +
  geom_text(aes(x = x_graf, y = y_graf, label = delito_dia),
            hjust = 0, fontface = "bold", family = "Fira Sans Medium", size = 6) +
  # geom_text(aes(x = 0.8, y = 9.7, hjust = 0, label = "En México en un día promedio se abrieron ... "),
  #           fontface = "bold", size = 11, family = "Arial Narrow") +
  #geom_segment(aes(x = 0.8, xend = 2.8, y = 11, yend = 11), color = "black", size = 1.5) +
  coord_cartesian(xlim = c(0.5, 3), ylim = c(-1, 11.5)) +
  #labs(title = paste0("¿Qué pasó en ", str_to_upper(nom_mes_boletin), " ",ao_actual, "?")) +
  # geom_text(aes(x = 0.8, y = -0.8, label = label_mes), size = 5, family = "Fira Sans Light", hjust = 0) + 
  theme_void() +
  theme(plot.title = element_text(family = "Arial Narrow", size = 32, hjust = 0.5, "bold")) +
  annotation_custom(hd, xmin = 0.68, xmax = 1.08, ymin = 8.9, ymax = 10.1) +
  annotation_custom(rc, xmin = 0.8, xmax = 1, ymin = 1.5, ymax = 2.5) +
  annotation_custom(rve, xmin = 0.8, xmax = 1, ymin = 2.5, ymax = 3.5) +
  annotation_custom(rvi, xmin = 0.8, xmax = 1, ymin = 3.5, ymax = 4.5) +
  annotation_custom(sec, xmin = 0.8, xmax = 1, ymin = 4.5, ymax = 5.5) +
  annotation_custom(ext, xmin = 0.8, xmax = 1, ymin = 5.5, ymax = 6.5) +
  annotation_custom(hc, xmin = 0.8, xmax = 1, ymin = 6.5, ymax = 7.5) +
  annotation_custom(fem, xmin = 0.8, xmax = 1, ymin = 7.5, ymax = 8.5) +
  annotation_custom(ld, xmin = 1.8, xmax = 2, ymin = 1.5, ymax = 2.5) +
  annotation_custom(nar, xmin = 1.8, xmax = 2, ymin = 2.5, ymax = 3.5) +
  annotation_custom(tra, xmin = 1.8, xmax = 2, ymin = 3.5, ymax = 4.5) +
  annotation_custom(vf, xmin = 1.8, xmax = 2, ymin = 4.5, ymax = 5.5) +
  annotation_custom(vio, xmin = 1.8, xmax = 2, ymin = 5.5, ymax = 6.5) +
  annotation_custom(rt, xmin = 1.8, xmax = 2, ymin = 7.5, ymax = 8.5) +
  annotation_custom(rn, xmin = 0.8, xmax = 1, ymin = 0.5, ymax = 1.5) +
  annotation_custom(rtp, xmin = 1.8, xmax = 2, ymin = 6.5, ymax = 7.5) 



ggsave(paste(out, "dia_prom.png", sep = "/"), width = 10, height = 6)


ftpUpload("dia_prom.PNG",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/dia_prom_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)




setwd("C:/Users/emcdo/OneDrive/Documents/ONC/ROBOT")
sink('email_test.txt')

library(blastula)
library(keyring)
library(tidyverse)
library(knitr)
library(rmarkdown)
library(webshot)

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")


render(input = 'tabla_boletin_vic.Rmd')
render(input = 'tabla_ci_boletin.Rmd')

webshot('tabla_boletin_vic.html', 'vic_screenshot.png', vheight = 380)
webshot('tabla_ci_boletin.html', 'ci_screenshot.png')

ftpUpload("vic_screenshot.PNG",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/vic_screenshot_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)

ftpUpload("ci_screenshot.PNG",
          paste0("ftp://uploadservice.admin%40onc.org.mx:sm64T45!@s220884.gridserver.com/imagenes/ci_screenshot_",  
                 mes_abbr[[num_mes_boletin]],
                 ao_actual, ".png")
)






render_email(input = "test1.Rmd")# %>% 
#    smtp_send(
#    to = c(
#"emcdonald@onc.org.mx"
#           ),
#    from = "emcdonald@onc.org.mx",
#    subject = "Boletin Mensual ONC - Octubre 2020",
#    credentials = creds_key("gmail1")
#  )

sink()

Sys.sleep(5)
print("slept 5")

