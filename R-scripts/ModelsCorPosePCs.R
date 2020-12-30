options(digits = 7)
library(tidyverse)
library(factoextra)

#### Schoeller Pheno data ####
#names of Schoeller twins in alphabetical order (as that is how they are sorted in the S3 bucket originally)
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
           "McLeod_Teisher","Mielke_Alyssa","Mielke_Emily","Mitchell_Fred","Mitchell_Ned","Moyer_Jean","Mueller_Kirk_2",
           "Mueller_Nate_2","Nagel_Jeff","Nagel_Steve","Nelson_Audri","Nelson_Erin","Nick_Skyler",
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

#### Pose ####
aws_pose <- read_csv(file = "~/github/facial-polyphenism/CSVs/schoellerFifth_raw_nodiv_pose.csv", col_names = TRUE) %>% 
  t() %>% 
  as.data.frame()

#remove rownames
rownames(aws_pose) <- c()

#make first row the column names
colnames(aws_pose) <- aws_pose[1,]

#remove first row
aws_pose <- slice_tail(aws_pose,n=(nrow(aws_pose)-1))

#clean first column
aws_pose$ID <- aws_pose$ID %>% str_sub(start = 1L, end = -5L)

#remove brackets from string vectors of numeric pose variables
aws_pose$Pitch <- aws_pose$Pitch %>% str_sub(start = 2L, end = -2L)
aws_pose$Yaw <- aws_pose$Yaw %>% str_sub(start = 2L, end = -2L)
aws_pose$Roll <- aws_pose$Roll %>% str_sub(start = 2L, end = -2L)

#coerce the string vectors into numeric vectors
aws_pose[,which(names(aws_pose) %in% c("Pitch","Yaw","Roll"))] <- sapply(aws_pose[,which(names(aws_pose) %in% c("Pitch","Yaw","Roll"))], as.numeric)

rm(list = ls()[!ls() %in% c("pheno_data","aws_pose")])

#### functions ####

#by default, returns pw distance dataframe with co-twin pairs arranged by the mean distance of 
#the variables, which are the top 100 loadings of PC1; if return_pw == FALSE and return_arr == TRUE, 
#then the mean of the top 100 variables of PC1 is returned; if both are false, the variable loadings
#and their contribution to PC1 are returned
meanTopLoadings <- function(df, id, columns=df, return_pw=TRUE, return_arr=FALSE){
  
  require(tidyverse)
  
  scaled_matrix <- scale(df[,columns], scale = FALSE)
  scaled_df <- as.data.frame(scaled_matrix)
  
  res.pca <- prcomp(scaled_df)
  
  PC1_loadings <- 
    res.pca[["rotation"]] %>% 
    as.data.frame()
  
  PC1_top100_loadingNames <- 
    PC1_loadings %>% 
    rownames_to_column(var = "variables") %>%
    select(variables, PC1) %>% 
    mutate(PC1_contr = ((PC1)**2)*100) %>%
    select(variables, PC1_contr) %>%
    arrange(desc(PC1_contr)) %>%
    head(100) %>%
    select(variables)
  
  df_pca_top100Loadings <- df[,which(names(df) %in% PC1_top100_loadingNames$variables)]
  meanTopLoadingsTemp <- rowMeans(df_pca_top100Loadings)
  
  schoeller_dlib194_norm_adp <- 
    df %>%
    cbind(., meanTopLoadingsTemp) %>%
    arrange(Co_twins_ID, meanTopLoadingsTemp) %>%
    select(-meanTopLoadingsTemp)
  
  if(return_arr==TRUE){
    
    data.frame(id, meanTopLoadingsTemp)
    
  } else if (return_pw==FALSE){
    
    PC1_loadings %>% 
      rownames_to_column(var = "variables") %>%
      select(variables, PC1) %>% 
      mutate(PC1_contr = ((PC1)**2)*100) %>%
      select(variables, PC1_contr) %>%
      arrange(desc(PC1_contr)) %>%
      head(100) %>%
      select(variables, PC1_contr)
    
    
  } else{
    
    df %>%
      cbind(., meanTopLoadingsTemp) %>%
      arrange(Co_twins_ID, meanTopLoadingsTemp) %>%
      select(-meanTopLoadingsTemp)
    
  }
}

#calculates the difference in the pitch between co-twin pairs
pitch_diff <- function(df, columns, pose_table=aws_pose){
  
  df_Names <- df$Names
  df_loadings <- 
    meanTopLoadings(df, id=df_Names, columns = columns, return_pw = FALSE, return_arr=TRUE)
  
  pitch_arr <- 
    aws_pose %>%
    left_join(., df_loadings, by = c("ID"="id")) %>%
    drop_na(meanTopLoadingsTemp) %>%
    left_join(., df[,which(names(df) %in% c("Names","Co_twins_ID"))], by = c("ID"="Names")) %>%
    select(ID, Co_twins_ID, meanTopLoadingsTemp, Pitch) %>%
    arrange(Co_twins_ID, meanTopLoadingsTemp) %>%
    select(-meanTopLoadingsTemp)
  
  #pitch delta
  CoIds= pitch_arr$Co_twins_ID
  
  UniqIds=unique(CoIds)
  
  count=1
  L.Pitch<-list()
  
  for (UniqId in UniqIds){
    
    idx<-as.numeric(which(CoIds==UniqId))
    
    for (v in 1:(length(idx)-1)){
      for (w in (v+1):length(idx)){
        i=idx[v]
        j=idx[w]
        nam= paste0(pitch_arr[i,1],"__",pitch_arr[j,1])
        pitch<-as.numeric(pitch_arr[j,which(names(pitch_arr) %in% c("Pitch"))] - pitch_arr[i,which(names(pitch_arr) %in% c("Pitch"))])
        L.Pitch[[nam]]<-pitch
        cat(nam,"\n")
      }
    }
  }
  
  pitch_arr_delta <- 
    do.call(rbind,L.Pitch) %>% 
    as.data.frame()
  colnames(pitch_arr_delta)<-"Pitch"
  
  pitch_arr_delta_temp <- pitch_arr_delta %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cotwins") 
  
  pitch_arr_delta_cotwins <- pitch_arr_delta_temp$cotwins %>% 
    str_split(pattern = "__") %>% 
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    rename("cotwin_1" = "V1") %>% 
    rename("cotwin_2" = 'V2') %>% 
    unnest(cols = c(cotwin_1, cotwin_2)) 
  
  pitch_arr_delta <- cbind(pitch_arr_delta_cotwins,pitch_arr_delta)
  rownames(pitch_arr_delta) <- c()
  pitch_arr_delta
  
}

#calculates the difference in the yaw between co-twin pairs
yaw_diff <- function(df, columns, pose_table=aws_pose){
  
  df_loadings <- 
    meanTopLoadings(df, id=df$Names, columns = columns, return_pw = FALSE, return_arr=TRUE)
  
  yaw_arr <- 
    aws_pose %>%
    left_join(., df_loadings, by = c("ID"="id")) %>%
    drop_na(meanTopLoadingsTemp) %>%
    left_join(., df[,which(names(df) %in% c("Names","Co_twins_ID"))], by = c("ID"="Names")) %>%
    select(ID, Co_twins_ID, meanTopLoadingsTemp, Yaw) %>%
    arrange(Co_twins_ID, meanTopLoadingsTemp) %>%
    select(-meanTopLoadingsTemp)

  #yaw delta
  CoIds= yaw_arr$Co_twins_ID
  
  UniqIds=unique(CoIds)
  
  count=1
  L.Yaw<-list()
  
  for (UniqId in UniqIds){
    
    idx<-as.numeric(which(CoIds==UniqId))
    
    for (v in 1:(length(idx)-1)){
      for (w in (v+1):length(idx)){
        i=idx[v]
        j=idx[w]
        nam= paste0(yaw_arr[i,1],"__",yaw_arr[j,1])
        Yaw<-as.numeric(yaw_arr[j,3] - yaw_arr[i,3])
        L.Yaw[[nam]]<-Yaw
        cat(nam,"\n")
      }
    }
  }
  
  yaw_arr_delta <- do.call(rbind,L.Yaw) %>% as.data.frame()
  colnames(yaw_arr_delta)<-colnames(yaw_arr)[3]
  
  yaw_arr_delta_temp <- yaw_arr_delta %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cotwins") 
  
  yaw_arr_delta_cotwins <- yaw_arr_delta_temp$cotwins %>% 
    str_split(pattern = "__") %>% 
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    rename("cotwin_1" = "V1") %>% 
    rename("cotwin_2" = 'V2') %>% 
    unnest(cols = c(cotwin_1, cotwin_2)) 
  
  yaw_arr_delta <- cbind(yaw_arr_delta_cotwins,yaw_arr_delta)
  rownames(yaw_arr_delta) <- c()
  yaw_arr_delta
}

#calculates the correlation of each PC with the Pitch and Yaw of the co-twin pairs
corPosePC <- function(df, columns = 1:ncol(df), returnPitch=TRUE){
  
  adp_pcs_pitch_cor <- c()
  adp_pcs_yaw_cor <- c()
  counter = 0
  for (i in columns){
    
    if (str_detect(names(df)[i], pattern = "PC") == TRUE){
      
      counter = counter + 1
      temp_pc_pitch_value <- cor.test(df[,i], schoeller_adp_pcs_pose[,which(names(df) %in% "Pitch")])
      temp_pc_yaw_value <-cor.test(df[,i], schoeller_adp_pcs_pose[,which(names(df) %in% "Yaw")])
      
      temp_pc_pitch_value <- cbind(temp_pc_pitch_value$estimate, temp_pc_pitch_value$p.value)
      temp_pc_yaw_value <- cbind(temp_pc_yaw_value$estimate, temp_pc_yaw_value$p.value)
      
      adp_pcs_pitch_cor <- rbind(adp_pcs_pitch_cor, temp_pc_pitch_value)
      adp_pcs_yaw_cor <- rbind(adp_pcs_yaw_cor, temp_pc_yaw_value)
      
      rownames(adp_pcs_pitch_cor)[counter] <- names(df)[i]
      rownames(adp_pcs_yaw_cor)[counter] <- names(df)[i]
      
    }
    
  }
  colnames(adp_pcs_pitch_cor) <- c("estimate","pvalue")
  colnames(adp_pcs_yaw_cor) <-c("estimate","pvalue")
  
  if (returnPitch){
    
    adp_pcs_pitch_cor
    
  } else {
    
    adp_pcs_yaw_cor
    
  }
  
}

#removes any infinite values in the dataset
rmInf <- function(df){
  
  columns2remove <- colnames(df)[colSums(is.infinite(as.matrix(df))) > 0]
  
  if (is.null(columns2remove) == FALSE){
    
    df_removedInf <- df[,-which(names(df)%in%columns2remove)]
    
  }
  return(df_removedInf)
  
}

#comparison of co-twins; fc = fold change, delta = raw difference
fc_diff <- function(df, individual_index, cotwin_index, starting_index){
  
  count=1
  L<-list()
  
  CoIds= df[,cotwin_index]
  
  UniqIds=unique(CoIds)
  
  for (UniqId in UniqIds){
    
    idx<-as.numeric(which(CoIds==UniqId))
    
    for (v in 1:(length(idx)-1)){
      for (w in (v+1):length(idx)){
        i=idx[v]
        j=idx[w]
        
        nam=paste0(df[i,individual_index],"__",df[j,individual_index])
        ratio<-as.numeric(df[j,starting_index:dim(df)[2]]) / as.numeric(df[i,starting_index:dim(df)[2]])
        L[[nam]]<-ratio
        cat(nam,"\n")
      }
    }
  }
  
  fc_df<-do.call(rbind,L)
  colnames(fc_df)<-colnames(df)[starting_index:ncol(df)]
  
  fc_df_temp <- fc_df %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cotwins") %>%
    select(cotwins)
  
  fc_df_cotwins <- fc_df_temp$cotwins %>% 
    str_split(pattern = "__") %>% 
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    rename("cotwin_1" = "V1") %>% 
    rename("cotwin_2" = 'V2') %>% 
    unnest(cols = c()) 
  
  fc_df <- cbind(fc_df_cotwins,fc_df)
  rownames(fc_df) <- c()
  fc_df
  
}
delta_diff <- function(df, individual_index, cotwin_index, starting_index){
  
  count=1
  L<-list()
  
  CoIds= df[,cotwin_index]
  
  UniqIds=unique(CoIds)
  
  for (UniqId in UniqIds){
    
    idx<-as.numeric(which(CoIds==UniqId))
    
    for (v in 1:(length(idx)-1)){
      for (w in (v+1):length(idx)){
        i=idx[v]
        j=idx[w]
        
        nam=paste0(df[i,individual_index],"__",df[j,individual_index])
        delta<-as.numeric(df[j,starting_index:dim(df)[2]])  - as.numeric(df[i,starting_index:dim(df)[2]])
        L[[nam]]<-delta
        cat(nam,"\n")
      }
    }
  }
  
  delta_df<-do.call(rbind,L)
  colnames(delta_df)<-colnames(df)[starting_index:ncol(df)]
  
  delta_df_temp <- delta_df %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cotwins") %>%
    select(cotwins)
  
  delta_df_cotwins <- delta_df_temp$cotwins %>% 
    str_split(pattern = "__") %>% 
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    rename("cotwin_1" = "V1") %>% 
    rename("cotwin_2" = 'V2') %>% 
    unnest(cols = c()) 
  
  delta_df <- cbind(delta_df_cotwins,delta_df)
  rownames(delta_df) <- c()
  delta_df
  
}

#normalization of PW distance data by argued IOD pw distance
IODnorm <- function(df, IODcolumn){
  
  dist_colnames <- colnames(df)
  dist_rownames <- rownames(df)
  df_norm <- data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
  counter <- 0
  for (i in 1:ncol(df)){
    
    counter <- counter + 1
    df_norm[,i] <- df[,i] / df[,which(names(df) %in% IODcolumn)]
    cat(paste(round(counter/ncol(df), 6)*100, "% done"),"\n")
    
  }
  colnames(df_norm) <- dist_colnames
  rownames(df_norm) <- dist_rownames
  df_norm
}

#function defining Euclidean distance for two landmarks of xy coordinates
eud<-function(x1,y1,x2,y2){
  sqrt((x2-x1)**2+(y2-y1)**2)
}

#nested-for loop to calculate all non-redundant PW distances; choose(194,2) = 18721
pw_dist_calc <- function(df, landmark_names){
  colnames_df <- colnames(df)
  L = list()
  for (i in seq(1,ncol(df)-2,2)){
    for (j in seq(i+2,ncol(df),2)){
      n1=landmark_names[i,1]
      n2=landmark_names[j,1]
      r=as.numeric()
      for (v in 1:nrow(df)){
        d1= eud(df[v,i],df[v,i+1],df[v,j],df[v,j+1])
        r=c(r,d1)
      }
      name=paste(n1,n2,sep="_")
      L[[name]]=r
      cat(name, "\n")
    }
  }
  M=do.call(rbind,L)
  M %>% t() %>% as.data.frame() 
}

##### dlib194 cascade depth 15 tree depth 4 ####
# all: pw distances
load("~/Desktop/OneDrive-Van-Andel-Institute/schoeller_fifth_pwDistance_dataframes.RData")

#PCA calculation
all_dist_pheno <- schoeller_dlib194_norm_adp

#center pairwise data (this is automatically done when calling prcomp() but I find it reassuring)
scaled_all_dist_pheno <- scale(all_dist_pheno[,10:ncol(all_dist_pheno)], scale = FALSE)
scaled_all_dist_pheno <- as.data.frame(scaled_all_dist_pheno)
scaled_all_dist_pheno <- cbind(schoeller_dlib194_norm_adp[1:9],scaled_all_dist_pheno)

#PCA of centered pw distances
res.pca <- prcomp(scaled_all_dist_pheno[,10:ncol(scaled_all_dist_pheno)])

adp_loadings <- res.pca[["rotation"]] 

schoeller_adp_pcs <- 
  res.pca[["x"]] %>% 
  as.data.frame()

schoeller_adp_pcs <- cbind(schoeller_dlib194_norm_adp[,1:9], schoeller_adp_pcs)

eigenvalues <- get_eig(res.pca)
eigenvalues_adp <- 
  eigenvalues %>% 
  select(variance.percent, cumulative.variance.percent)

schoeller_adp_pcs_pose <- 
  schoeller_adp_pcs %>% 
  left_join(., aws_pose[,which(names(aws_pose) %in% c("ID","Pitch","Yaw"))], by=c("Names"="ID")) 


adp_pcs_pitch_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = TRUE)
adp_pcs_yaw_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/Plots/Annotation-plots/dlib194_pop_norm_pcs_pitch.pdf")
adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .8, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/Plots/Annotation-plots/dlib194_pop_norm_pcs_yaw.pdf")
adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .8, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20)) 
dev.off()

adp_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
dlib194_pop_cumvar <-
  adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)


schoeller_dlib194_norm_adp_arr <- meanTopLoadings(df=schoeller_dlib194_norm_adp, columns = 10:ncol(schoeller_dlib194_norm_adp))
dlib194_norm_adp_arr_fc <- fc_diff(df = schoeller_dlib194_norm_adp_arr, individual_index = 1, cotwin_index = 2, starting_index = 10)
dlib194_norm_adp_arr_delta <- delta_diff(df = schoeller_dlib194_norm_adp_arr, individual_index = 1, cotwin_index = 2, starting_index = 10)

dlib194_pitch_diff <- pitch_diff(df = schoeller_dlib194_norm_adp, columns = 10:ncol(schoeller_dlib194_norm_adp), pose_table = aws_pose)
dlib194_yaw_diff <- yaw_diff(df = schoeller_dlib194_norm_adp, columns = 10:ncol(schoeller_dlib194_norm_adp), pose_table = aws_pose)


dlib194_pitch_diff <- 
  pitch_diff(df = schoeller_dlib194_norm_adp_arr, columns = 10:ncol(schoeller_dlib194_norm_adp_arr), pose_table = aws_pose) %>%
  left_join(., schoeller_dlib194_norm_adp_arr[,which(names(schoeller_dlib194_norm_adp_arr) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
dlib194_yaw_diff <- 
  yaw_diff(df = schoeller_dlib194_norm_adp_arr, columns = 10:ncol(schoeller_dlib194_norm_adp_arr), pose_table = aws_pose) %>%
  left_join(., schoeller_dlib194_norm_adp_arr[,which(names(schoeller_dlib194_norm_adp_arr) %in% c("Names","Co_twins_ID"))], 
            by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
dlib194_pose_diff <- 
  dlib194_pitch_diff %>%
  cbind(., dlib194_yaw_diff$Yaw) %>%
  rename("Yaw"="dlib194_yaw_diff$Yaw")

#fc
dlib194_pop_norm_adp_arr_fc_scaled <- 
  scale(dlib194_norm_adp_arr_fc[,3:ncol(dlib194_norm_adp_arr_fc)], scale = FALSE) %>%
  as.data.frame()
columns2remove <- colnames(dlib194_pop_norm_adp_arr_fc_scaled)[colSums(is.infinite(as.matrix(dlib194_pop_norm_adp_arr_fc_scaled))) > 0]
if (is.null(columns2remove) == FALSE){
  
  dlib194_pop_norm_adp_arr_fc_scaled <- dlib194_pop_norm_adp_arr_fc_scaled[,-which(names(dlib194_pop_norm_adp_arr_fc_scaled)%in%columns2remove)]
  
}
dlib194_pop_norm_adp_arr_fc_scaled <- rmInf(dlib194_pop_norm_adp_arr_fc_scaled)
res.pca <- prcomp(dlib194_pop_norm_adp_arr_fc_scaled, scale = FALSE, center = FALSE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between  PCs and pose
fc_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_fc_pcs_pose <- 
  fc_norm_pcs %>% 
  cbind(dlib194_norm_adp_arr_fc[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., dlib194_pose_diff[,which(names(dlib194_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

fc_pcs_pitch_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = TRUE)
fc_pcs_yaw_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = FALSE)


pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib194_schoeller_fc_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib194_schoeller_fc_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
fc_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

fc_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
dlib194_fc_pcs_cor_sigFiltered <- 
  fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

dlib194_fc_pcs_cor_sigFiltered
max(dlib194_fc_pcs_cor_sigFiltered$cumulative.variance.percent)

#delta
res.pca <- prcomp(dlib194_norm_adp_arr_delta[,3:ncol(dlib194_norm_adp_arr_delta)], scale = FALSE, center = TRUE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between PCs and pose
delta_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_delta_pcs_pose <- 
  delta_norm_pcs %>% 
  cbind(dlib194_norm_adp_arr_delta[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., dlib194_pose_diff[,which(names(dlib194_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

delta_pcs_pitch_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = TRUE)
delta_pcs_yaw_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib194_schoeller_delta_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib194_schoeller_delta_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
delta_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

delta_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
dlib194_delta_pcs_cor_sigFiltered <- 
  delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

dlib194_delta_pcs_cor_sigFiltered
max(dlib194_delta_pcs_cor_sigFiltered$cumulative.variance.percent)
#### AWS ####

AWS_schoeller <- read_csv("/Users/Warren.Sink/github/facial-polyphenism/CSVs/schoeller_raw_for_comparison.csv",col_names = FALSE, skip = 1)
AWS_schoeller <- as.data.frame(t(AWS_schoeller))
rownames(AWS_schoeller) <- c()
AWS_schoeller$V1 <- AWS_schoeller$V1 %>% str_remove(pattern = ".jpg")
colnames(AWS_schoeller)<-AWS_schoeller[1,]
AWS_schoeller<-AWS_schoeller[-1,]
rownames(AWS_schoeller) <- c()
AWS_schoeller <- AWS_schoeller %>% as.data.frame() %>% column_to_rownames(var = "ID")
rownames_AWS_schoeller <- rownames(AWS_schoeller)
#remove string character parts for turning into numeric vectors 
for (i in (1:ncol(AWS_schoeller))){
  for (j in (1:nrow(AWS_schoeller))){
    AWS_schoeller[j,i]<-str_remove(string = AWS_schoeller[j,i], pattern = "\\[|\\]")
    AWS_schoeller[j,i]<-str_remove(string = AWS_schoeller[j,i], pattern = "\\]|\\[")
  }
}

#first, we'll look at the difference in quality in capturing the 
AWS_schoeller <- splitstackshape::cSplit(AWS_schoeller, names(AWS_schoeller)) 
AWS_schoeller <- cbind(rownames_AWS_schoeller, AWS_schoeller)
AWS_schoeller <- AWS_schoeller %>% column_to_rownames(var = "rownames_AWS_schoeller")

landmark_names <- 
  names(AWS_schoeller) %>%
  as.data.frame() %>%
  rename("names"=".") %>%
  mutate(names = str_remove(names, "_x|_y|_z"))

AWS_schoeller_pw_dist_calc <- pw_dist_calc(df = AWS_schoeller, landmark_names = landmark_names)
rownames(AWS_schoeller_pw_dist_calc) <- rownames(AWS_schoeller)
saveRDS(AWS_schoeller_pw_dist_calc, "~/github/facial-polyphenism/AWS_schoeller_pw_dist_calc.rds")

AWS_schoeller_pw_dist_calc <- readRDS("~/github/facial-polyphenism/AWS_schoeller_pw_dist_calc.rds")
AWS_schoeller_pw_dist_calc <-
  AWS_schoeller_pw_dist_calc %>%
  rownames_to_column(var = "Names") %>%
  #demonets aren't spelled consistently; have added "_2" because they were copied twice
  mutate(Names = str_replace(Names, pattern = "David_2", replacement = "David")) %>%
  mutate(Names = str_replace(Names, pattern = "Larry_2", replacement = "Larry")) %>%
  column_to_rownames(var = "Names")

AWS_pw_norm <- IODnorm(df = AWS_schoeller_pw_dist_calc, IODcolumn = c("LeftEyeLeft_1_RightEyeRight_1"))

res.pca <- prcomp(AWS_pw_norm, scale = FALSE, center = TRUE)

adp_loadings <- res.pca[["rotation"]] 

schoeller_adp_pcs <- 
  res.pca[["x"]] %>% 
  as.data.frame()

rownames(schoeller_adp_pcs) <- rownames(AWS_pw_norm)
schoeller_adp_pcs <- 
  schoeller_adp_pcs %>%
  rownames_to_column(var = "Names")
eigenvalues <- get_eig(res.pca)
eigenvalues_adp <- 
  eigenvalues %>% 
  select(variance.percent, cumulative.variance.percent)

schoeller_adp_pcs_pose <- 
  schoeller_adp_pcs %>% 
  left_join(., aws_pose[,which(names(aws_pose) %in% c("ID","Pitch","Yaw"))], by=c("Names"="ID")) 

adp_pcs_pitch_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = TRUE)
adp_pcs_yaw_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/Plots/Annotation-plots/AWS_pop_norm_pcs_pitch.pdf")
adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .8, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/Plots/Annotation-plots/AWS_pop_norm_pcs_yaw.pdf")
adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .8, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20)) 
dev.off()

adp_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
AWS_pop_cumvar <-
  adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

AWS_pw_norm <- 
  AWS_pw_norm %>%
  rownames_to_column(var = "Names") %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = "Names") %>%
  select(Names, Co_twins_ID, everything())
  
AWS_pw_norm_adp_arr <- meanTopLoadings(df=AWS_pw_norm, columns = 3:ncol(AWS_pw_norm)) 
AWS_pw_norm_adp_arr_fc <- fc_diff(df = AWS_pw_norm_adp_arr, individual_index = 1, cotwin_index = 2, starting_index = 3)
AWS_pw_norm_adp_arr_delta <- delta_diff(df = AWS_pw_norm_adp_arr, individual_index = 1, cotwin_index = 2, starting_index = 3)


AWS_pitch_diff <- pitch_diff(df = AWS_pw_norm, columns = 3:ncol(AWS_pw_norm), pose_table = aws_pose)
AWS_yaw_diff <- yaw_diff(df = AWS_pw_norm, columns = 3:ncol(AWS_pw_norm), pose_table = aws_pose)

AWS_pitch_diff <- 
  pitch_diff(df = dlib68_pop_norm_adp_arr, columns = 3:ncol(dlib68_pop_norm_adp_arr), pose_table = aws_pose) %>%
  left_join(., dlib68_pop_norm_adp_arr[,which(names(dlib68_pop_norm_adp_arr) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
AWS_yaw_diff <- 
  yaw_diff(df = dlib68_pop_norm_adp_arr, columns = 3:ncol(dlib68_pop_norm_adp_arr), pose_table = aws_pose) %>%
  left_join(., dlib68_pop_norm_adp_arr[,which(names(dlib68_pop_norm_adp_arr) %in% c("Names","Co_twins_ID"))], 
            by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
AWS_pose_diff <- 
  AWS_pitch_diff %>%
  cbind(., AWS_yaw_diff$Yaw) %>%
  rename("Yaw"="AWS_yaw_diff$Yaw")

#fc
AWS_pop_norm_adp_arr_fc_scaled <- 
  scale(AWS_pw_norm_adp_arr_fc[,3:ncol(AWS_pw_norm_adp_arr_fc)], scale = FALSE) %>%
  as.data.frame()
AWS_pop_norm_adp_arr_fc_scaled <- rmInf(AWS_pop_norm_adp_arr_fc_scaled)
res.pca <- prcomp(AWS_pop_norm_adp_arr_fc_scaled, scale = FALSE, center = FALSE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between  PCs and pose
fc_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_fc_pcs_pose <- 
  fc_norm_pcs %>% 
  cbind(AWS_pw_norm_adp_arr_fc[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., AWS_pose_diff[,which(names(AWS_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

fc_pcs_pitch_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = TRUE)
fc_pcs_yaw_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/AWS_schoeller_fc_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/AWS_schoeller_fc_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
fc_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

fc_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
AWS_fc_pcs_cor_sigFiltered <- 
  fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

AWS_fc_pcs_cor_sigFiltered
max(AWS_fc_pcs_cor_sigFiltered$cumulative.variance.percent)

#delta
res.pca <- prcomp(AWS_pop_norm_adp_arr_delta[,3:ncol(AWS_pop_norm_adp_arr_delta)], scale = FALSE, center = TRUE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between PCs and pose
delta_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_delta_pcs_pose <- 
  delta_norm_pcs %>% 
  cbind(AWS_pop_norm_adp_arr_delta[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., AWS_pose_diff[,which(names(AWS_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

delta_pcs_pitch_cor <- corPosePC(delta_pcs_pitch_cor, returnPitch = TRUE)
delta_pcs_yaw_cor <- corPosePC(delta_pcs_yaw_cor, returnPitch = FALSE)


pdf("~/github/facial-polyphenism/CSVs/FACEMESH/AWS_schoeller_delta_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/AWS_schoeller_delta_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
delta_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

delta_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
AWS_delta_pcs_cor_sigFiltered <- 
  delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

AWS_delta_pcs_cor_sigFiltered
max(AWS_delta_pcs_cor_sigFiltered$cumulative.variance.percent)


#### dlib68 ####
shape_predictor_68_face_landmarks_coord_median <- read_csv(file = "/Users/Warren.Sink/PycharmProjects/dlib-facial-landmarks/shape_predictor_68_face_landmarks_schoeller_coord_33.csv", col_names = FALSE) %>% 
  as_tibble()

shape_predictor_68_face_landmarks_names_median <- read_csv(file = "/Users/Warren.Sink/PycharmProjects/dlib-facial-landmarks/shape_predictor_68_face_landmarks_schoeller_names_33.csv", col_names = FALSE) %>%
  as_tibble()

#get list of image files
setwd("~/PycharmProjects/dlib-facial-landmarks/194_landmark_trained_model_cascade15_tdepth4/resized_faces_median")
img_file_list <- list.files(pattern = "\\.jpg$")

#get image dimensions
library(magick)
jpeg_dim_median <- data.frame()
for (i in 1:length(img_file_list)){
  temp_data <- image_read(img_file_list[i])
  img <- image_info(temp_data)
  img_attributes <- img %>% select(format, width, height)
  jpeg_dim_median <- rbind(jpeg_dim_median, img_attributes)
}
cotwins_names <- rep(shape_predictor_68_face_landmarks_names_median$X1, each=2)
cotwins_names <- str_remove(cotwins_names, " copy")

#landmark names: numbers indicate going from lateral to median or toward the center
landmark_names <- c("LeftJaw_8","LeftJaw_7","LeftJaw_6","LeftJaw_5","LeftJaw_4","LeftJaw_3",
                    "LeftJaw_2","LeftJaw_1","ChinBottom","RightJaw_1","RightJaw_2","RightJaw_3",
                    "RightJaw_4","RightJaw_5","RightJaw_6","RightJaw_7","RightJaw_8","LeftEyebrow_5",
                    "LeftEyebrow_4","LeftEyebrow_3","LeftEyebrow_2","LeftEyebrow_1","RightEyebrow_1",
                    "RightEyebrow_2","RightEyebrow_3","RightEyebrow_4","RightEyebrow_5","MiddleNose_5",
                    "MiddleNose_4","MiddleNose_3","MiddleNose_2","LeftNose_2","LeftNose_1","MiddleNose_1",
                    "RightNose_1","RightNose_2","LeftEyeLeft","LeftEyeTop_2","LeftEyeTop_1","LeftEyeRight",
                    "LeftEyeBottom_2","LeftEyeBottom_1","RightEyeLeft","RightEyeTop_1","RightEyeTop_2",
                    "RightEyeRight","RightEyeBottom_2","RightEyeBottom_1","LeftMouthCorner","UpperLipTopLeft_2",
                    "UpperLipTopLeft_1","UpperLipCenterTop","UpperLipTopRight_1","UpperLipTopRight_2",
                    "RightMouthCorner","LowerLipBottomRight_2","LowerLipBottomRight_1","LowerLipCenterBottom",
                    "LowerLipBottomLeft_1","LowerLipBottomLeft_2","UpperLipBottomLeft_2","UpperLipBottomLeft_1",
                    "UpperLipCenterBottom","UpperLipBottomRight_1",'UpperLipBottomRight_2',"LowerLipTopRight_1",
                    "LowerLipTopCenter","LowerLipTopLeft_1")

#replicating each landmark name twice
landmark_names_rep <- rep(landmark_names, each = 2)

#assigning x and y coordinates for each landmark
for (i in 1:(length(landmark_names_rep)-1)){
  if (landmark_names_rep[i] == landmark_names_rep[i+1]){
    landmark_names_rep[i] = paste(landmark_names_rep[i], "x",sep = "_")
    landmark_names_rep[i+1] = paste(landmark_names_rep[i+1], "y",sep = "_")
  }
}

dlib_facial_landmarks_median <- as.data.frame(t(shape_predictor_68_face_landmarks_coord_median))

dlib_facial_landmarks_median <- cbind(cotwins_names,dlib_facial_landmarks_median)
colnames(dlib_facial_landmarks_median)[2:ncol(dlib_facial_landmarks_median)] <- landmark_names
dlib_facial_landmarks_xcoord_indices <- seq(1,nrow(dlib_facial_landmarks_median)-1,2)
dlib_facial_landmarks_ycoord_indices <- seq(2,nrow(dlib_facial_landmarks_median),2)
dlib_facial_landmarks_median_xcoord <- dlib_facial_landmarks_median[dlib_facial_landmarks_xcoord_indices,]
dlib_facial_landmarks_median_ycoord <- dlib_facial_landmarks_median[dlib_facial_landmarks_ycoord_indices,]
rownames(dlib_facial_landmarks_median_xcoord) <- c()
rownames(dlib_facial_landmarks_median_ycoord) <- c()
dlib_facial_landmarks_median_xcoord <- dlib_facial_landmarks_median_xcoord %>% arrange(cotwins_names) %>% column_to_rownames(.data = .,var = "cotwins_names")
dlib_facial_landmarks_median_ycoord <- dlib_facial_landmarks_median_ycoord %>% arrange(cotwins_names) %>% column_to_rownames(.data = .,var = "cotwins_names")

dlib_facial_landmarks_xcoord_twinnames <- 
  dlib_facial_landmarks_median_xcoord %>% 
  rownames_to_column(.,var = "rowname") %>%
  select(rowname)

testdf <- data.frame(matrix(NA, nrow = ncol(dlib_facial_landmarks_median_xcoord)*2, ncol = nrow(dlib_facial_landmarks_median_xcoord)))
counter = 1
for (i in seq(from = 1, to = ncol(dlib_facial_landmarks_median_xcoord)*2, by = 1)){
  if (i %% 2 > 0){
    testdf[i,]<-(dlib_facial_landmarks_median_xcoord[,((i+1)/2)])
  } else {
    testdf[i,]<-(dlib_facial_landmarks_median_ycoord[,(i/2)])
  }
}
temp_facial_landmarks_median <- as.data.frame(testdf)
rownames(temp_facial_landmarks_median) <- c()
colnames(temp_facial_landmarks_median) <- c()

for (i in 1:ncol(temp_facial_landmarks_median)){
  colnames(temp_facial_landmarks_median)[i] <- dlib_facial_landmarks_xcoord_twinnames[i,1]
}
for (i in 1:nrow(temp_facial_landmarks_median)){
  rownames(temp_facial_landmarks_median)[i] <- landmark_names_rep[i]
}
temp_facial_landmarks_median <- temp_facial_landmarks_median %>%
  rownames_to_column(var = "rowname")

dlib68_facial_landmarks_median <- temp_facial_landmarks_median
dlib68_facial_landmarks_median <- dlib68_facial_landmarks_median %>% column_to_rownames(var = "rowname") 
dlib68_facial_landmarks_median <- as.data.frame(t(dlib68_facial_landmarks_median))
rownames_dlib68_facial_landmarks_median <- rownames(dlib68_facial_landmarks_median)
colnames_dlib68_facial_landmarks_median <- colnames(dlib68_facial_landmarks_median)
#compute the actual coordinates for each image
testdf <- data.frame()
for (i in (1:(nrow(dlib68_facial_landmarks_median)))){
  for (j in (1:(ncol(dlib68_facial_landmarks_median)))){
    if (j %% 2 > 0){
      testdf[i,j]<-(dlib68_facial_landmarks_median[i,j]/jpeg_dim_median[i,2])
    } else {
      testdf[i,j]<-((1-(dlib68_facial_landmarks_median[i,j])/jpeg_dim_median[i,3]))
    }
  }
}
dlib68_facial_landmarks_median <- testdf
colnames(dlib68_facial_landmarks_median) <- c()
colnames(dlib68_facial_landmarks_median) <- colnames_dlib68_facial_landmarks_median
rownames(dlib68_facial_landmarks_median) <- rownames_dlib68_facial_landmarks_median

#assigning image names to columns
colnames(temp_facial_landmarks) <- shape_predictor_68_face_landmarks_names_median$X1

landmark_names <- 
  names(dlib68_facial_landmarks_median) %>%
  as.data.frame() %>%
  rename("names"=".") %>%
  mutate(names = str_remove(names, "_x|_y|_z"))

dlib68_facial_landmarks_median_pw_dist_calc <- pw_dist_calc(df = dlib68_facial_landmarks_median, landmark_names = landmark_names)
rownames(dlib68_facial_landmarks_median_pw_dist_calc) <- rownames(dlib68_facial_landmarks_median)
saveRDS(dlib68_facial_landmarks_median_pw_dist_calc, "~/github/facial-polyphenism/dlib68_facial_landmarks_median_pw_dist_calc.rds")

dlib68_pop_norm <- IODnorm(df = dlib68_facial_landmarks_median_pw_dist_calc, IODcolumn = c("LeftEyeLeft_RightEyeRight"))

res.pca <- prcomp(dlib68_pop_norm, scale = FALSE, center = TRUE)

adp_loadings <- res.pca[["rotation"]] 

schoeller_adp_pcs <- 
  res.pca[["x"]] %>% 
  as.data.frame()

eigenvalues <- get_eig(res.pca)
eigenvalues_adp <- 
  eigenvalues %>% 
  select(variance.percent, cumulative.variance.percent)

rownames(schoeller_adp_pcs) <- rownames(dlib68_pop_norm)
schoeller_adp_pcs <- 
  schoeller_adp_pcs %>%
  rownames_to_column(var = "Names")
schoeller_adp_pcs_pose <- 
  schoeller_adp_pcs %>% 
  left_join(., aws_pose[,which(names(aws_pose) %in% c("ID","Pitch","Yaw"))], by=c("Names"="ID")) 

adp_pcs_pitch_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = TRUE)
adp_pcs_yaw_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = FALSE)


pdf("~/github/facial-polyphenism/Plots/Annotation-plots/dlib68_pop_norm_pcs_pitch.pdf")
adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .8, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/Plots/Annotation-plots/dlib68_pop_norm_pcs_yaw.pdf")
adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .8, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20)) 
dev.off()

adp_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
dlib68_pop_cumvar <-
  adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues_adp) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

dlib68_pop_norm <-
  dlib68_pop_norm %>%
  rownames_to_column(var = "Names") %>%
  #demonets aren't spelled consistently; have added "_2" because they were copied twice
  mutate(Names = str_replace(Names, pattern = "David_2", replacement = "David")) %>%
  mutate(Names = str_replace(Names, pattern = "Larry_2", replacement = "Larry")) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = "Names") %>%
  select(Names, Co_twins_ID, everything())
dlib68_pop_norm_adp_arr <- meanTopLoadings(df=dlib68_pop_norm, columns = 3:ncol(dlib68_pop_norm)) 
dlib68_pop_norm_adp_arr_fc <- fc_diff(df = dlib68_pop_norm_adp_arr, individual_index = 1, cotwin_index = 2, starting_index = 3)
dlib68_pop_norm_adp_arr_delta <- delta_diff(df = dlib68_pop_norm_adp_arr, individual_index = 1, cotwin_index = 2, starting_index = 3)


dlib68_pitch_diff <- pitch_diff(df = dlib68_pop_norm, columns = 3:ncol(dlib68_pop_norm), pose_table = aws_pose)
dlib68_yaw_diff <- yaw_diff(df = dlib68_pop_norm, columns = 3:ncol(dlib68_pop_norm), pose_table = aws_pose)

dlib68_pitch_diff <- 
  pitch_diff(df = dlib68_pop_norm_adp_arr, columns = 3:ncol(dlib68_pop_norm_adp_arr), pose_table = aws_pose) %>%
  left_join(., dlib68_pop_norm_adp_arr[,which(names(dlib68_pop_norm_adp_arr) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
dlib68_yaw_diff <- 
  yaw_diff(df = dlib68_pop_norm_adp_arr, columns = 3:ncol(dlib68_pop_norm_adp_arr), pose_table = aws_pose) %>%
  left_join(., dlib68_pop_norm_adp_arr[,which(names(dlib68_pop_norm_adp_arr) %in% c("Names","Co_twins_ID"))], 
            by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
dlib68_pose_diff <- 
  dlib68_pitch_diff %>%
  cbind(., dlib68_yaw_diff$Yaw) %>%
  rename("Yaw"="dlib68_yaw_diff$Yaw")

#fc
dlib68_pop_norm_adp_arr_fc_scaled <- 
  scale(dlib68_pop_norm_adp_arr_fc[,3:ncol(dlib68_pop_norm_adp_arr_fc)], scale = FALSE) %>%
  as.data.frame()
dlib68_pop_norm_adp_arr_fc_scaled <- rmInf(dlib68_pop_norm_adp_arr_fc_scaled)
res.pca <- prcomp(dlib68_pop_norm_adp_arr_fc_scaled, scale = FALSE, center = FALSE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between  PCs and pose
fc_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_fc_pcs_pose <- 
  fc_norm_pcs %>% 
  cbind(dlib68_pop_norm_adp_arr_fc[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., dlib68_pose_diff[,which(names(dlib68_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

fc_pcs_pitch_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = TRUE)
fc_pcs_yaw_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = FALSE)


pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib68_schoeller_fc_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib68_schoeller_fc_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
fc_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

fc_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
dlib68_fc_pcs_cor_sigFiltered <- 
  fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

dlib68_fc_pcs_cor_sigFiltered
max(dlib68_fc_pcs_cor_sigFiltered$cumulative.variance.percent)

#delta
res.pca <- prcomp(dlib68_pop_norm_adp_arr_delta[,3:ncol(dlib68_pop_norm_adp_arr_delta)], scale = FALSE, center = TRUE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between PCs and pose
delta_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_delta_pcs_pose <- 
  delta_norm_pcs %>% 
  cbind(dlib68_pop_norm_adp_arr_delta[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., dlib68_pose_diff[,which(names(dlib68_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

delta_pcs_pitch_cor <- corPosePC(schoeller_delta_pcs_pose, returnPitch = TRUE)
delta_pcs_yaw_cor <- corPosePC(schoeller_delta_pcs_pose, returnPitch = FALSE)


pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib68_schoeller_delta_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/dlib68_schoeller_delta_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
delta_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

delta_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
dlib68_delta_pcs_cor_sigFiltered <- 
  delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

dlib68_delta_pcs_cor_sigFiltered
max(dlib68_delta_pcs_cor_sigFiltered$cumulative.variance.percent)


#### FACEMESH ####
#import
load("~/github/facial-polyphenism/FACEMESH_schoeller_PWdist_dfs.RData")

res.pca <- prcomp(FACEMESH_schoeller_pw_norm_arranged[,10:ncol(FACEMESH_schoeller_pw_norm_arranged)], scale = FALSE, center = TRUE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between PCs and pose
pop_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_adp_pcs_pose <- 
  pop_norm_pcs %>% 
  cbind(FACEMESH_schoeller_pw_norm_arranged[,c(1,3)], .) %>%
  left_join(., aws_pose[,which(names(aws_pose) %in% c("ID","Pitch","Yaw"))], by=c("Names"="ID")) 

adp_pcs_pitch_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = TRUE)
adp_pcs_yaw_cor <- corPosePC(schoeller_adp_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/FACEMESH_schoeller_pop_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/FACEMESH_schoeller_pop_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

adp_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
FACEMESH_pcs_cor_sigFiltered <- 
  adp_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

FACEMESH_pcs_cor_sigFiltered
max(FACEMESH_pcs_cor_sigFiltered$cumulative.variance.percent)

FACEMESH_pop_norm_adp_arr_fc <- fc_diff(df = FACEMESH_schoeller_pw_norm_arranged, individual_index = 1, cotwin_index = 3, starting_index = 10)
FACEMESH_pop_norm_adp_arr_delta <- delta_diff(df = FACEMESH_schoeller_pw_norm_arranged, individual_index = 1, cotwin_index = 3, starting_index = 10)


FACEMESH_pitch_diff <- 
  pitch_diff(df = FACEMESH_schoeller_pw_norm_arranged, columns = 10:ncol(FACEMESH_schoeller_pw_norm_arranged), pose_table = aws_pose) %>%
  left_join(., FACEMESH_schoeller_pw_norm_arranged[,which(names(FACEMESH_schoeller_pw_norm_arranged) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
FACEMESH_yaw_diff <- 
  yaw_diff(df = FACEMESH_schoeller_pw_norm_arranged, columns = 10:ncol(FACEMESH_schoeller_pw_norm_arranged), pose_table = aws_pose) %>%
  left_join(., FACEMESH_schoeller_pw_norm_arranged[,which(names(FACEMESH_schoeller_pw_norm_arranged) %in% c("Names","Co_twins_ID"))], 
            by = c("cotwin_1"="Names")) %>%
  select(cotwin_1, cotwin_2, Co_twins_ID, everything())
FACEMESH_pose_diff <- 
  FACEMESH_pitch_diff %>%
  cbind(., FACEMESH_yaw_diff$Yaw) %>%
  rename("Yaw"="FACEMESH_yaw_diff$Yaw")

#fc
res.pca <- prcomp(FACEMESH_pop_norm_adp_arr_fc[,3:ncol(FACEMESH_pop_norm_adp_arr_fc)], scale = FALSE, center = TRUE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between  PCs and pose
fc_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_fc_pcs_pose <- 
  fc_norm_pcs %>% 
  cbind(FACEMESH_pop_norm_adp_arr_fc[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., FACEMESH_pose_diff[,which(names(FACEMESH_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

fc_pcs_pitch_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = TRUE)
fc_pcs_yaw_cor <- corPosePC(schoeller_fc_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/FACEMESH_schoeller_fc_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/FACEMESH_schoeller_fc_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
fc_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

fc_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
FACEMESH_fc_pcs_cor_sigFiltered <- 
  fc_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

FACEMESH_fc_pcs_cor_sigFiltered
max(FACEMESH_fc_pcs_cor_sigFiltered$cumulative.variance.percent)

#delta
res.pca <- prcomp(FACEMESH_pop_norm_adp_arr_delta[,3:ncol(FACEMESH_pop_norm_adp_arr_delta)], scale = FALSE, center = TRUE) 
eigenvalues <- get_eig(res.pca)
#testing correlation between PCs and pose
delta_norm_pcs <- 
  res.pca[["x"]] %>%
  as.data.frame()
schoeller_delta_pcs_pose <- 
  delta_norm_pcs %>% 
  cbind(FACEMESH_pop_norm_adp_arr_delta[,1:2], .) %>%
  left_join(., pheno_data[,which(names(pheno_data) %in% c("Names","Co_twins_ID"))], by = c("cotwin_1"="Names")) %>%
  left_join(., FACEMESH_pose_diff[,which(names(FACEMESH_pose_diff) %in% c("cotwin_1","cotwin_2","Pitch","Yaw"))], 
            by=c("cotwin_1","cotwin_2"))

delta_pcs_pitch_cor <- corPosePC(schoeller_delta_pcs_pose, returnPitch = TRUE)
delta_pcs_yaw_cor <- corPosePC(schoeller_delta_pcs_pose, returnPitch = FALSE)

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/FACEMESH_schoeller_delta_IODnorm_mf_corPitchPCsSig.pdf",width = 15,height = 20)
delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

pdf("~/github/facial-polyphenism/CSVs/FACEMESH/FACEMESH_schoeller_delta_IODnorm_mf_corYawPCsSig",width = 15,height = 20)
delta_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  mutate(PCnum = str_sub(PC, 3L, str_length(PC)), PCnum = as.integer(PCnum),
         PCnum = case_when(PCnum < 6 ~ "PCs 01 to 05", 
                           PCnum < 11 ~ "PCs 06 to 10", 
                           PCnum < 16 ~ "PCs 11 to 15", 
                           PCnum < 21 ~ "PCs 16 to 20", 
                           PCnum < 26 ~ "PCs 21 to 25", 
                           PCnum < 31 ~ "PCs 26 to 30",
                           TRUE ~ "PCs greater than 30")) %>%
  select(PC, PCnum, pvalue, variance.percent) %>%
  arrange(pvalue, desc()) %>%
  head(20) %>%
  mutate(PC = as.factor(PC), pvalue = -log10(pvalue), PC = fct_reorder(PC, pvalue)) %>%
  ggplot(aes(x=PC, y=pvalue, fill=PCnum)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = ifelse(pvalue > 3, "*", "")), vjust = .7, size = 20) +
  coord_flip() +
  labs(x="", y = "-log10(pvalue)") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title = element_text(size=20),
        legend.text=element_text(size=20))  
dev.off()

delta_pcs_yaw_cor_temp <-  
  adp_pcs_yaw_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue) %>%
  rename("pvalue_yaw"="pvalue")
FACEMESH_delta_pcs_cor_sigFiltered <- 
  delta_pcs_pitch_cor %>% 
  as.data.frame() %>%
  cbind(., eigenvalues) %>%
  rownames_to_column(var = "PC") %>%
  select(PC, pvalue, variance.percent) %>%
  left_join(., adp_pcs_yaw_cor_temp, by = "PC") %>%
  filter(pvalue > 0.001, pvalue_yaw > 0.001) %>%
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>%
  select(PC, cumulative.variance.percent)

FACEMESH_delta_pcs_cor_sigFiltered
max(FACEMESH_delta_pcs_cor_sigFiltered$cumulative.variance.percent)

