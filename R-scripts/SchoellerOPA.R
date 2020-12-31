options(digits = 7)
library(tidyverse)
library(shapes)
#### pheno data ####
Names <- c("Allooh_Jonas","Allooh_Moses","Allooh_Noah","Andrews_Doug","Andrews_Ross",
           "Ayer_Carly","Ayer_Lily","Bowen_David","Bowen_Richard","Bowers_Brenda",
           "Braswell_Orlean","Bridges_Carolyn","Burns_Catherine","Burns_McKenzie","Catalano_Joe",
           "Catalano_Sal","Cepeda_Eurides","Cepeda_Ramon","Claucherty_Marion",
           "Croll-Baehre_Emma","Croll-Baehre_Marta","Darren_1", "Darren_2", "Davis_June","Demonet_David","Demonet_Larry",
           "Dwomoh-Piper_Chantelle","Dwomoh-Piper_Danielle","Elder_Marilyn","Estrada_Antonion",
           "Estrada_Jesus","Everingham_Douglas","Everingham_Steve","Furtick_Jason","Furtick_Keith",
           "Furtick_Kevin","Furtick_Victor","Garzilli_Denise","Ginley_Adalynn","Ginley_Breanna",
           "Goldberg_Adam","Goldberg_Peter","Haick_Ava","Haick_Isabel","Hudson_Irene","Humphery_Emma",
           "Humphery_Megan","Humphery_Sarah","Jackson_Sidney","Jacques_Mary","Jensen_Dane","Jensen_Jordan","Karen_1", "Karen_2",
           "Kerner_Sharon","Key_Aidan","Kobylski_Cecelia","Kobylski_Elizabeth","Lipfird_Garry","Lipfird_Larry","Lois_1", "Lois_2",
           "Maassen_Jamie","Martin_Kate","Martin_emily","McClure_Bryan","McClure_Bryce","McLeod_Keisher",
           "McLeod_Teisher","Mielke_Alyssa","Mielke_Emily","Mitchell_Fred","Mitchell_Ned","Moyer_Jean","Mueller_Kirk",
           "Mueller_Nate","Nagel_Jeff","Nagel_Steve","Nelson_Audri","Nelson_Erin","Nick_Skyler",
           "Nick_Spencer","Oliver_David","Oliver_Walter","Parks_Katie","Parks_Sarah","Prijatel","Qualkinbush_Jodie","Rabi_Kate",
           "Ramsdell_Cole","Ramsdell_Seth","Ream_Alyssa","Ream_Rachel","Rowan_Hannah","Rowan_Kaci",
           "Scott-Wolf_Adriana","Serra_Ivory","Serra_Shelter","Smith_Connor","Smith_Kyle","Starr_Ellen",
           "Steeves_Carrie","Steeves_Patty","Stewart_Martha","Stillwell_Diane","Thornhill_Jennie","Timotiwu_Jason",
           "Timotiwu_Jordan","Tynski_Daniel","Tynski_Kristen","Ullman_Helen","Underwood_Helen","Whited_Jackie",
           "Whited_Jessica","Widerka_Amanda","Widerka_Beth")

Co_twins_ID <- c('1','1','1','2','2','3','3','4','4','5','6','7','8','9','10','10','11','11','8','12','12','52','52','13','14','14','15','15','7','16','16','17','17','18','18','18','18','19',
                 '20','20','21','21','22','22','6','23','23','23','9','24','25','25','53','53','26','5','27','27','28','28','54','54','29','30','30','31','31','32','32','33','33','34','34',
                 '13','35','35','36','36','37','37','38','38','39','39','40','40','26','29','41','42','42','43','43','44','44','41','45','45','46','46','47','24','24','8','19','24','48',
                 '48','49','49','8','47','50','50','51','51')

Individual_ID <- c('1A','1B','1C','2A','2B','3A','3B','4A','4B','5A','6A','7A','8A','9A','10A','10B','11A','11B','8B','12A','12B','52A','52B','13A','14A','14B','15A','15B','7B','16A','16B','17A','17B','18A','18B','18C','18D','19A',
                   '20A','20B','21A','21B','22A','22B','6B','23A','23B','23C','9B','24A','25A','25B','53A','53B','26A','5B','27A','27B','28A','28B','54A','54B','29A','30A','30B','31A','31B','32A','32B','33A','33B','34A','34B',
                   '13B','35A','35B','36A','36B','37A','37B','38A','38B','39A','39B','40A','40B','26B','29B','41A','42A','42B','43A','43B','44A','44B','41B','45A','45B','46A','46B','47A','24B','24C','8C','19B','24D','48A',
                   '48B','49A','49B','8D','47B','50A','50B','51A','51B')

Sex <- c('M','M','M','M','M','F','F','M','M','F','F','F','F','F','M','M','M','M','F','F','F','M','M','F','M','M','F','F','F','M','M','M','M','M','M','M','M','F',
         'F','F','M','M','F','F','F','F','F','F','F','F','M','M','F','F','F','M','F','F','M','M','F','F','F','M','M','F','F','F','F','F','F','M','M',
         'F','M','M','M','M','F','F','M','M','M','M','F','F','F','F','F','M','M','F','F','F','F','F','M','M','M','M','F','F','F','F','F','F','M',
         'M','M','F','F','F','F','F','F','F')

Ethnicity <- c('caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','caucasian','caucasian','hispanic','hispanic','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','black','caucasian','caucasian','caucasian','caucasian','caucasian','black','black','black','black','caucasian',
               'black','black','caucasian','caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','black','caucasian','caucasian','black','black',
               'caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','asian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','asian','caucasian','caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','caucasian','caucasian','asian',
               'asian','caucasian','caucasian','caucasian','black','caucasian','caucasian','hispanic','hispanic')

Eye_color <- c('brown','brown','brown','blue','blue','blue','blue','brown','brown','brown','brown','brown','green','brown','brown','brown','brown','brown','green','blue','blue','green','green','blue','blue','blue','brown','brown','brown','brown','brown','brown','brown','brown','brown','brown','brown','blue',
               'brown','brown','blue','blue','blue','blue','brown','blue','blue','blue','brown','blue','blue','blue','blue','blue','blue','brown','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','brown','brown','blue','blue','brown','brown',
               'blue','brown','brown','brown','brown','brown','brown','blue','blue','blue','blue','blue','blue','blue','blue','brown','brown','brown','green','green','blue','blue','brown','brown','brown','green','green','brown','blue','blue','green','blue','blue','brown',
               'brown','brown','brown','brown','brown','blue','blue','brown','brown')


Hair_color <- c('black','black','black','bald','bald','blond','blond','grey','grey','brown','brown','brown','red','black','black','black','black','black','black','red','red','brown','brown','grey','grey','grey','black','black','brown','black','black','grey','grey','black','black','black','black','brown',
                'black','black','bald','bald','brown','brown','brown','blond','blond','blond','blond','blond','blond','blond','blond','blond','blond','grey','brown','brown','grey','grey','blond','blond','brown','blond','blond','blond','blond','black','black','blond','blond','black','black',
                'grey','brown','brown','brown','brown','brown','brown','brown','brown','grey','grey','blond','blond','blond','brown','brown','brown','brown','red','red','blond','blond','brown','black','black','brown','brown','grey','blond','blond','black','brown','blond','black',
                'black','brown','brown','black','grey','blond','blond','brown','brown')

Birth_year <- c('1985','1985','1985','1965','1965','2006','2006','1957','1957','1964','1945','1944','1958','1991','1943','1943','1977','1977','1958','1996','1996','1992','1992','1939','1940','1940','1987','1987','1944','1988','1988','1947','1947','1990','1990','1990','1990','1953',
                '2005','2005','1970','1970','2007','2007','1945','1998','1998','1998','1991','1982','1986','1986','1972','1972','1949','1964','1996','1996','1954','1954','1952','1952','1981','2002','2002','2000','2000','1977','1977','2000','2000','1950','1950',
                '1939','1984','1984','1967','1967','1977','1977','1992','1992','1936','1936','2001','2001','1949','1981','1982','1988','1988','1997','1997','1992','1992','1982','1972','1972','2001','2001','1942','1982','1982','1958','1953','1982','1994',
                '1994','1984','1984','1958','1942','1991','1991','2005','2005')


pheno_data <- data.frame(Names, Individual_ID, Co_twins_ID, Sex, Ethnicity, Birth_year, Eye_color, Hair_color)
pheno_data[,1:8] <- lapply(pheno_data[,1:8] , function(x) as.factor(as.character(x)))
pheno_data$Co_twins_ID <- as.integer(as.character(pheno_data$Co_twins_ID))
pheno_data$Birth_year = as.numeric(as.character(pheno_data$Birth_year))
pheno_data <- pheno_data %>% mutate(., Age = 2012-pheno_data$Birth_year)

rm(list = ls()[!ls() %in% c("pheno_data")])

#### wrangling #####
FACEMESH<- 
  read_csv("~/github/Tensorflowjs_Projects/FACEMESH_schoellerParsed.txt", col_names = "results") %>%
  mutate(results = str_remove(results, " "))
images <- 
  read_csv("~/PycharmProjects/dlib-facial-landmarks/FACEMESH_schoeller_imgnames.csv", col_names = "names")  %>% 
  mutate(names = str_remove(names, " "))

FACEMESH_restructured <- data.frame(matrix(NA, nrow = nrow(images), ncol = 468 * 3))

FACEMESH_names <-
  FACEMESH %>%
  filter(results %in% images$names) %>%
  mutate(results = str_remove(results, ".jpg")) %>%
  rename("names"="results")
FACEMESH_landmarks <-
  FACEMESH %>%
  filter(str_detect(results, "landmark")) %>%
  unique() %>%
  mutate(results_x = paste(results, "_x", sep=""), results_y = paste(results, "_y", sep=""), results_z = paste(results, "_z", sep="")) %>%
  select(-results)
FACEMESH_coords <-
  FACEMESH %>%
  filter(!(results %in% images$names), !(str_detect(results, "landmark")))

landmark_counter = 1
#populate column names
for (i in 1:nrow(FACEMESH_landmarks)){
  
  for (j in 1:ncol(FACEMESH_landmarks)){
    
    colnames(FACEMESH_restructured)[landmark_counter] <- FACEMESH_landmarks[i, j]
    landmark_counter = landmark_counter + 1
  }
  
}
#populate dataframe with coordinates
for (i in 1:ncol(FACEMESH_restructured)){
  
  for (j in 1:nrow(FACEMESH_restructured)){
    
    FACEMESH_restructured[j,i] <- FACEMESH_coords[i+j,]
    
  }
  
}

FACEMESH_restructured <- 
  FACEMESH_restructured %>% 
  mutate_if(is.character, as.numeric)

FACEMESH_schoeller <- 
  cbind(FACEMESH_names, FACEMESH_restructured) %>%
  as.data.frame() %>%
  rename("Names"="names") %>%
  mutate(Names = str_replace(Names, pattern = "Nate_2", replace = "Nate")) %>%
  mutate(Names = str_replace(Names, pattern = "Kirk_2", replace = "Kirk")) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = "Names") %>%
  mutate(Co_twins_ID = as.character(Co_twins_ID)) %>%
  select(Names, Co_twins_ID, everything())

Names_trip <- rep(FACEMESH_schoeller$Names, each = 3)
Co_twins_ID_trip <- rep(FACEMESH_schoeller$Co_twins_ID, each = 3)

FACEMESH_landmarks_unique <-
  colnames(FACEMESH_schoeller)[3:ncol(FACEMESH_schoeller)] %>%
  str_remove(pattern = "_x|_y|_z") %>% 
  unique()

FACEMESH_schoeller_coordinates <- FACEMESH_schoeller[,3:ncol(FACEMESH_schoeller)]
FACEMESH_schoeller_xcoord_indices <- seq(1,ncol(FACEMESH_schoeller_coordinates)-2,3)
FACEMESH_schoeller_ycoord_indices <- seq(2,ncol(FACEMESH_schoeller_coordinates)-1,3)
FACEMESH_schoeller_zcoord_indices <- seq(3,ncol(FACEMESH_schoeller_coordinates)-0,3)
FACEMESH_schoeller_xcoord <- FACEMESH_schoeller_coordinates[,FACEMESH_schoeller_xcoord_indices]
FACEMESH_schoeller_ycoord <- FACEMESH_schoeller_coordinates[,FACEMESH_schoeller_ycoord_indices]
FACEMESH_schoeller_zcoord <- FACEMESH_schoeller_coordinates[,FACEMESH_schoeller_zcoord_indices]

FACEMESH_schoeller_coords_reshaped <- data.frame(matrix(NA, nrow = nrow(FACEMESH_schoeller_coordinates) * 3, ncol = length(FACEMESH_landmarks_unique)))
counter_id <- 0
for (j in seq(0, (nrow(FACEMESH_schoeller_coordinates) * 3) - 3, by = 3)){
  counter_id<-counter_id+1
  for (i in 1:length(FACEMESH_landmarks_unique)){
    
    FACEMESH_schoeller_coords_reshaped[j+1,i]<-(FACEMESH_schoeller_xcoord[counter_id,i])
    FACEMESH_schoeller_coords_reshaped[j+2,i]<-(FACEMESH_schoeller_ycoord[counter_id,i])
    FACEMESH_schoeller_coords_reshaped[j+3,i]<-(FACEMESH_schoeller_zcoord[counter_id,i])
    
  } 
}
colnames(FACEMESH_schoeller_coords_reshaped) <- FACEMESH_landmarks_unique
FACEMESH_schoeller_coords_reshaped <-
  cbind(Names_trip, Co_twins_ID_trip, FACEMESH_schoeller_coords_reshaped) %>%
  rename("Names"="Names_trip","Co_twins_ID"="Co_twins_ID_trip")

#### Ordinary Procrustes Analysis #####
FACEMESH_schoeller_twins <- FACEMESH_schoeller_coords_reshaped[FACEMESH_schoeller_coords_reshaped$Co_twins_ID %in% names(which(table(FACEMESH_schoeller_coords_reshaped$Co_twins_ID) == 6)), ]
FACEMESH_schoeller_triplets <- FACEMESH_schoeller_coords_reshaped[FACEMESH_schoeller_coords_reshaped$Co_twins_ID %in% names(which(table(FACEMESH_schoeller_coords_reshaped$Co_twins_ID) == 9)), ]
FACEMESH_schoeller_quadruplets <- FACEMESH_schoeller_coords_reshaped[FACEMESH_schoeller_coords_reshaped$Co_twins_ID %in% names(which(table(FACEMESH_schoeller_coords_reshaped$Co_twins_ID) == 12)), ]

#twins
L.name = list()
counter = 1
counter2 = 4
twins_ID <- c()
OPA_rmsd_twins <- c()
OPA_twins_OSS <- c()
OPA_df_twins <- data.frame()
FACEMESH_schoeller_twins <- FACEMESH_schoeller_twins %>% arrange(Co_twins_ID)
FACEMESH_schoeller_twins_matrix <- FACEMESH_schoeller_twins[,3:ncol(FACEMESH_schoeller_twins)] %>% t() %>% as.matrix() 
for (j in seq(6,ncol(FACEMESH_schoeller_twins_matrix),by = 6)){
  
  nam = paste0(FACEMESH_schoeller_twins[counter,1],"__",FACEMESH_schoeller_twins[counter2,1])
  L.name[counter.names]<-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_twins_matrix[,counter:(counter+2)], 
                      B = FACEMESH_schoeller_twins_matrix[,counter2:(counter2+2)], 
                      scale = TRUE)

  cat(paste(FACEMESH_schoeller_twins[counter,which(names(FACEMESH_schoeller_twins) %in% c("Names"))], 
            FACEMESH_schoeller_twins[counter2,which(names(FACEMESH_schoeller_twins) %in% c("Names"))],
            sep = "__"), "\n")
  
  OPA_rmsd_twins <- c(OPA_rmsd_twins,OPA_temp[["rmsd"]])
  OPA_twins_OSS <- c(OPA_twins_OSS,OPA_temp[["OSS"]])
  OPA_twins_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_twins_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_twins <- rbind(OPA_df_twins,OPA_twins_Ahat,OPA_twins_Bhat)
  
  names_twins <- do.call(rbind,L.name)
  names_twins <- as.data.frame(names_twins)
  twins_ID <- rbind(twins_ID, names_twins)
  
  counter = 6 + counter
  counter2 = 6 + counter2
  
}
colnames(OPA_df_twins) <- names(FACEMESH_schoeller_twins[,3:ncol(FACEMESH_schoeller_twins)])
twins_ID_temp <- strsplit(sub('__', '\\1 \\2', twins_ID$V1), ' ')
twins_ID_temp <- flatten(twins_ID_temp)
twins_ID_temp <- do.call(rbind,twins_ID_temp)
Names_xyz <- rep(twins_ID_temp[,1], each = 3)
FACEMESH_OPA_twins <- 
  cbind(Names_xyz, OPA_df_twins) %>%
  rename("Names"="Names_xyz") %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = "Names") %>%
  select(Names, Co_twins_ID, everything())


OPA_twins_OSS <- 
  cbind(twins_ID,OPA_twins_OSS)

OPA_twins_OSS <- 
  OPA_twins_OSS %>% 
  as.data.frame() %>% 
  rename("OSS" = "OPA_twins_OSS", "twins_ID"="V1")

OPA_twins_OSS_cotwins <-
  OPA_twins_OSS$twins_ID %>% 
  str_split(pattern = "__") %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rename("cotwin_1" = "V1") %>% 
  rename("cotwin_2" = 'V2') %>% 
  unnest(cols = c(cotwin_1, cotwin_2)) 

FACEMESH_OPA_twins_OSS <- 
  cbind(OPA_twins_OSS_cotwins,OPA_twins_OSS$OSS) %>%
  rename("OSS"="OPA_twins_OSS$OSS")

FACEMESH_OPA_twins_RMSD <- 
  cbind(OPA_twins_OSS_cotwins,OPA_rmsd_twins) %>% 
  rename("RMSD" = "OPA_rmsd_twins")

#triplets
L.name = list()
names_trip_all <- data.frame()
counter.names = 1
counter = 1
counter2 = 4
counter3 = 7
OPA_triplets_OSS <- c()
OPA_rmsd_triplets <- c()
OPA_df_triplets <- data.frame()
FACEMESH_schoeller_triplets <- FACEMESH_schoeller_triplets %>% arrange(Co_twins_ID)
FACEMESH_schoeller_triplets_matrix <- FACEMESH_schoeller_triplets[,3:ncol(FACEMESH_schoeller_triplets)] %>% t() %>% as.matrix() 
for (j in seq(6, ncol(FACEMESH_schoeller_triplets_matrix)-6, by = 6)){
  
  nam = paste0(FACEMESH_schoeller_triplets[counter,1],"__",FACEMESH_schoeller_triplets[counter2,1])
  L.name[counter.names]<-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_triplets_matrix[,counter:(counter+2)], 
                      B = FACEMESH_schoeller_triplets_matrix[,counter2:(counter2+2)], 
                      scale = TRUE)
  
  OPA_rmsd_triplets <- c(OPA_rmsd_triplets, OPA_temp[["rmsd"]])
  OPA_triplets_OSS <- c(OPA_triplets_OSS, OPA_temp[["OSS"]])
  OPA_triplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_triplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_triplets <- rbind(OPA_df_triplets,OPA_triplets_Ahat,OPA_triplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_triplets[counter,1],"__",FACEMESH_schoeller_triplets[counter3,1])
  L.name[counter.names+1] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_triplets_matrix[,counter:(counter+2)], 
                      B = FACEMESH_schoeller_triplets_matrix[,counter3:(counter3+2)], 
                      scale = TRUE)
  
  OPA_rmsd_triplets <- c(OPA_rmsd_triplets,OPA_temp[["rmsd"]])
  OPA_triplets_OSS <- c(OPA_triplets_OSS, OPA_temp[["OSS"]])
  OPA_triplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_triplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_triplets <- rbind(OPA_df_triplets,OPA_triplets_Ahat,OPA_triplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_triplets[counter2,1],"__",FACEMESH_schoeller_triplets[counter3,1])
  L.name[counter.names+2] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_triplets_matrix[,counter2:(counter2+2)], 
                      B = FACEMESH_schoeller_triplets_matrix[,counter3:(counter3+2)], 
                      scale = TRUE)
  
  OPA_rmsd_triplets <- c(OPA_rmsd_triplets,OPA_temp[["rmsd"]])
  OPA_triplets_OSS <- c(OPA_triplets_OSS, OPA_temp[["OSS"]])
  OPA_triplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_triplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_triplets <- rbind(OPA_df_triplets,OPA_triplets_Ahat,OPA_triplets_Bhat)
  
  names_trip <- do.call(rbind,L.name)
  names_trip <- as.data.frame(names_trip)
  names_trip_all <- rbind(names_trip_all, names_trip)
  
  counter = 9 + counter
  counter2 = 9 + counter2
  counter3 = 9 + counter3
  
}
names_trip_all_temp <- strsplit(sub('__', '\\1 \\2', names_trip_all$V1), ' ')
names_trip_all_temp <- flatten(names_trip_all_temp)
names_trip_all_temp <- do.call(rbind,names_trip_all_temp)
Names_xyz <- rep(names_trip_all_temp[,1], each = 3)
FACEMESH_OPA_triplets <- 
  cbind(Names_xyz, OPA_df_triplets) %>%
  rename("Names"="Names_xyz") %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = "Names") %>%
  select(Names, Co_twins_ID, everything())
colnames(FACEMESH_OPA_triplets)[3:ncol(FACEMESH_OPA_triplets)] <- names(FACEMESH_schoeller_triplets[,3:ncol(FACEMESH_schoeller_triplets)])

cotwins_trip <-
  names_trip_all$V1 %>%
  str_split("__") %>%
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rename("cotwin_1" = "V1") %>% 
  rename("cotwin_2" = 'V2') %>% 
  unnest(cols = c(cotwin_1, cotwin_2)) 

FACEMESH_OPA_triplets_OSS <-
  cbind(cotwins_trip, OPA_triplets_OSS) %>% 
  rename("OSS"="OPA_triplets_OSS")

FACEMESH_OPA_triplets_RMSD <- 
  cbind(cotwins_trip,OPA_rmsd_triplets) %>%
  rename("RMSD"="OPA_rmsd_triplets")

#quadruplets
L.name = list()
names_quad_all <- data.frame()
counter.names = 1
counter = 1
counter2 = 4
counter3 = 7
counter4 = 10
OPA_rmsd_quadruplets <- c()
OPA_quad_OSS <- c()
OPA_df_quadruplets <- data.frame()
FACEMESH_schoeller_quadruplets <- FACEMESH_schoeller_quadruplets %>% arrange(Co_twins_ID)
FACEMESH_schoeller_quadruplets_matrix <- FACEMESH_schoeller_quadruplets[,3:ncol(FACEMESH_schoeller_quadruplets)] %>% t() %>% as.matrix() 
for (j in seq(12,ncol(FACEMESH_schoeller_quadruplets_matrix),by = 12)){
  
  nam = paste0(FACEMESH_schoeller_quadruplets[counter,1],"__",FACEMESH_schoeller_quadruplets[counter2,1])
  L.name[counter.names]<-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_quadruplets_matrix[,counter:(counter+2)], 
                      B = FACEMESH_schoeller_quadruplets_matrix[,counter2:(counter2+2)], 
                      scale = TRUE)
  
  OPA_rmsd_quadruplets <- c(OPA_rmsd_quadruplets,OPA_temp[["rmsd"]])
  OPA_quad_OSS <- c(OPA_quad_OSS,OPA_temp[["OSS"]])
  OPA_quadruplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_quadruplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_quadruplets <- rbind(OPA_df_quadruplets,OPA_quadruplets_Ahat,OPA_quadruplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_quadruplets[counter,1],"__",FACEMESH_schoeller_quadruplets[counter3,1])
  L.name[counter.names+1] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_quadruplets_matrix[,counter:(counter+2)], 
                      B = FACEMESH_schoeller_quadruplets_matrix[,counter3:(counter3+2)], 
                      scale = TRUE)
  
  OPA_rmsd_quadruplets <- c(OPA_rmsd_quadruplets,OPA_temp[["rmsd"]])
  OPA_quad_OSS <- c(OPA_quad_OSS,OPA_temp[["OSS"]])
  OPA_quadruplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_quadruplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_quadruplets <- rbind(OPA_df_quadruplets,OPA_quadruplets_Ahat,OPA_quadruplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_quadruplets[counter,1],"__",FACEMESH_schoeller_quadruplets[counter4,1])
  L.name[counter.names+2] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_quadruplets_matrix[,counter:(counter+2)], 
                      B = FACEMESH_schoeller_quadruplets_matrix[,counter4:(counter4+2)], 
                      scale = TRUE)
  
  OPA_rmsd_quadruplets <- c(OPA_rmsd_quadruplets,OPA_temp[["rmsd"]])
  OPA_quad_OSS <- c(OPA_quad_OSS,OPA_temp[["OSS"]])
  OPA_quadruplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_quadruplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_quadruplets <- rbind(OPA_df_quadruplets,OPA_quadruplets_Ahat,OPA_quadruplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_quadruplets[counter2,1],"__",FACEMESH_schoeller_quadruplets[counter3,1])
  L.name[counter.names+3] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_quadruplets_matrix[,counter2:(counter2+2)], 
                      B = FACEMESH_schoeller_quadruplets_matrix[,counter3:(counter3+2)], 
                      scale = TRUE)
  
  OPA_rmsd_quadruplets <- c(OPA_rmsd_quadruplets,OPA_temp[["rmsd"]])
  OPA_quad_OSS <- c(OPA_quad_OSS,OPA_temp[["OSS"]])
  OPA_quadruplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_quadruplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_quadruplets <- rbind(OPA_df_quadruplets,OPA_quadruplets_Ahat,OPA_quadruplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_quadruplets[counter2,1],"__",FACEMESH_schoeller_quadruplets[counter4,1])
  L.name[counter.names+4] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_quadruplets_matrix[,counter2:(counter2+2)], 
                      B = FACEMESH_schoeller_quadruplets_matrix[,counter4:(counter4+2)], 
                      scale = TRUE)
  
  OPA_rmsd_quadruplets <- c(OPA_rmsd_quadruplets,OPA_temp[["rmsd"]])
  OPA_quad_OSS <- c(OPA_quad_OSS,OPA_temp[["OSS"]])
  OPA_quadruplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_quadruplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_quadruplets <- rbind(OPA_df_quadruplets,OPA_quadruplets_Ahat,OPA_quadruplets_Bhat)
  
  nam = paste0(FACEMESH_schoeller_quadruplets[counter3,1],"__",FACEMESH_schoeller_quadruplets[counter4,1])
  L.name[counter.names+5] <-nam
  OPA_temp <- procOPA(A = FACEMESH_schoeller_quadruplets_matrix[,counter3:(counter3+2)], 
                      B = FACEMESH_schoeller_quadruplets_matrix[,counter4:(counter4+2)], 
                      scale = TRUE)
  
  OPA_rmsd_quadruplets <- c(OPA_rmsd_quadruplets,OPA_temp[["rmsd"]])
  OPA_quad_OSS <- c(OPA_quad_OSS,OPA_temp[["OSS"]])
  OPA_quadruplets_Ahat <- OPA_temp[["Ahat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  OPA_quadruplets_Bhat <- OPA_temp[["Bhat"]] %>% as.matrix() %>% t() %>% as.data.frame()
  
  OPA_df_quadruplets <- rbind(OPA_df_quadruplets,OPA_quadruplets_Ahat,OPA_quadruplets_Bhat)
  
  names_quad <- do.call(rbind,L.name)
  names_quad <- as.data.frame(names_quad)
  names_quad_all <- rbind(names_quad_all, names_quad)
  
  counter = 12 + counter
  counter2 = 12 + counter2
  counter3 = 12 + counter3
  counter4 = 12 + counter4
  
}
names_quad_all_temp <- strsplit(sub('__', '\\1 \\2', names_quad_all$V1), ' ')
names_quad_all_temp <- flatten(names_quad_all_temp)
names_quad_all_temp <- do.call(rbind,names_quad_all_temp)
Names_xyz <- rep(names_quad_all_temp[,1], each = 3)
FACEMESH_OPA_quadruplets <- 
  cbind(Names_xyz, OPA_df_quadruplets) %>%
  rename("Names"="Names_xyz") %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names", "Co_twins_ID"))], by = "Names") %>%
  select(Names, Co_twins_ID, everything())
colnames(FACEMESH_OPA_quadruplets)[3:ncol(FACEMESH_OPA_quadruplets)] <- names(FACEMESH_schoeller_quadruplets[,3:ncol(FACEMESH_schoeller_quadruplets)])

cotwins_trip <-
  names_quad_all$V1 %>%
  str_split("__") %>%
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rename("cotwin_1" = "V1") %>% 
  rename("cotwin_2" = 'V2') %>% 
  unnest(cols = c(cotwin_1, cotwin_2)) 

FACEMESH_OPA_quad_OSS <-
  cbind(cotwins_trip, OPA_triplets_OSS) %>% 
  rename("OSS"="OPA_triplets_OSS")

FACEMESH_OPA_quad_RMSD <- 
  cbind(cotwins_trip,OPA_rmsd_triplets) %>%
  rename("RMSD"="OPA_rmsd_triplets")

FACEMESH_OPA <- rbind(FACEMESH_OPA_twins, FACEMESH_OPA_triplets, FACEMESH_OPA_quadruplets)

FACEMESH_OPA_OSS <- 
  rbind(FACEMESH_OPA_twins_OSS, FACEMESH_OPA_triplets_OSS, FACEMESH_OPA_quad_OSS) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names", "Co_twins_ID"))], by = c("cotwin_1" = "Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, OSS)
FACEMESH_OPA_RMSD <- 
  rbind(FACEMESH_OPA_twins_RMSD, FACEMESH_OPA_triplets_RMSD, FACEMESH_OPA_quad_RMSD) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names", "Co_twins_ID"))], by = c("cotwin_1" = "Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, RMSD)

rm(list = ls()[!ls() %in% c("FACEMESH_OPA","FACEMESH_OPA_OSS","FACEMESH_OPA_RMSD")])
