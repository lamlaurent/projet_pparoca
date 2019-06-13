###################################
######PACKAGES#####################
###################################


pacman::p_load(tidyverse, 
               DescTools,
               readxl, 
               janitor, 
               lubridate,
               Hmisc,
               questionr,
               hydroTSM,
               prettyR,
               tableone,
               skimr,
               mice, 
               magrittr,
               DataExplorer,
               dataMaid,
               tidyr,
               dplyr,
               psych,
               questionr,
               pastecs,
               summarytools,
               magrittr,
               naniar,
               epiDisplay,
               GGally,
               reshape,
               lme4,
               nlme,
               gplots,
               lmerTest,
               ggplot2,
               splitstackshape)

select <- dplyr::select

######################## importation du dataframe et création des fichiers long et wide ##################################

dl <- read.csv2("base.csv",na.strings=c("•","N/C",""))


### reformatage du nom des variables
dl <- dl %>% 
  clean_names() 

names(dl)

# suppression des colonnes vides
dl <- dl[1:41]
dl <- dl[1:84,]

#repérage et recodage des dates de naissance variant en fonction des consultations 

dl <- dl %>% group_by(id_patient) %>% mutate(date_naissance=date_naissance[1]) %>% ungroup() %>% as.data.frame() 

# recodage du type de variables
dl<- dl %>% mutate_at(c("id_patient"),funs(factor(.)))

# création d'un fichier wide
 
dw <- reshape(dl,direction="wide",timevar ="visite" ,idvar="id_patient",v.names = c("date_visite","temps","intro_traitement",
                            "audc","audc_dose_mg_j",
                           "fibrate","type_de_fibrate","fibrate_dose_mg_j","ocaliva",            
                           "ocaliva_dose_mg_j","prurit","prurit_echelle_0_10","asat_ui_l","asat_x_n","alat_ui_l","alat_x_n",
                            "ggt_ui_l","ggt_x_n","pal_ui_l","pal_x_n","bilirubine_µmol_l","tp","albumine_g_l","plaquettes_µ_l",     
                            "creatinine_µmol_l","cholesterol_g_l","ig_m_g_l","myalgies","date_fibroscan","fibroscan"))  


## représentations graphiques ###########

# taux de PAL (xN) en fonction des visites (ensemble des patients)

means <- aggregate(pal_ui_l ~  visite, dl, mean)
ggplot(dl) +
  aes(x = visite, y = pal_ui_l) +
  geom_boxplot() +
  xlab("numéro de visite") +
  ylab("PAL (UI/L)") +
  ggtitle("PAL en UI/L en fonction des visites")+
  ylim(0,500)+
  geom_text(data = means, aes(label = round(pal_ui_l,digits = 0), y = pal_ui_l + 0))
  

 
############################ TEST T APPARIE ###############################################################################

### taux de PAL entre V2 et V4 

t.test(dw$pal_ui_l.V2,dw$pal_ui_l.V4,paired = TRUE) # p-value = 1.962e-05 / difference de 262.1481 entre V2 et V4

### taux de PAL entre V2 et V3
t.test(dw$pal_ui_l.V2,dw$pal_ui_l.V3,paired = TRUE) # p-value = 0.0003433 / difference de 172.3704 entre V2 et V3

### taux de PAL entre V3 et V4
t.test(dw$pal_ui_l.V3,dw$pal_ui_l.V4,paired = TRUE) # p-value = 3.868e-06 / difference de 86.57143 entre V3 et V4


# évolution du taux moyen entre V2, V3 et V4 
mean(dw$pal_ui_l.V2,na.rm = TRUE) # taux moyen de 400 à V2
mean(dw$pal_ui_l.V3,na.rm = TRUE) # taux moyen de 221 à V3
mean(dw$pal_ui_l.V4,na.rm = TRUE) # taux moyen de 135 à V4






############## regression linéaire  ##############################################

## entre v2, V3 et V4
summary(lmer(pal_ui_l~temps + (1|id_patient), data=dl))


## entre V2 et V4
dl24 <- dl %>% select (id_patient, visite,pal_ui_l) %>% filter (visite=="V2"|visite=="V4") 

lm24 <- lm(pal_ui_l~visite + id_patient, data=dl24)
summary(lm24)
summary(lmer(pal_ui_l~visite + (1|id_patient), data=dl24)) # même résultat quand on met les sujets en facteur aléatoire
drop1(lm24,.~.,test = "F") # p=1.962e-05

## entre V2 et V3
dl23 <- dl %>% select (id_patient, visite,pal_ui_l) %>% filter (visite=="V2"|visite=="V3") 
lm23 <- lm(pal_ui_l~visite + id_patient, data=dl23)
summary(lm23)
drop1(lm23,.~.,test = "F") #p= 0.0003433

## entre V3 et V4
dl34 <- dl %>% select (id_patient, visite,pal_ui_l) %>% filter (visite=="V3"|visite=="V4") 
lm34 <- lm(pal_ui_l~visite + id_patient, data=dl34)
summary(lm34)
drop1(lm34,.~.,test = "F") #p=3.868e-06

###########################################################
###################### Modèle mixte ######################
############################################################

# le taux de PAL(xN) en fonction du temps (ensemble des patients) ####################################################

div1<-matrix(data=c(1), ncol=1, nrow=1)
layout(div1)
plotmeans(dl$pal_ui_l~dl$visite,pch=16,barcol="black", xlab="numéro de la
visite", ylab="Moyenne du taux de PAL en UI/L")


# spaghetti plot

p <- ggplot(data = dl, aes(x = temps, y = pal_ui_l, group = id_patient))
p + geom_line() # donne la ligne évolutive de Hamilton pour chaque patient
p + geom_line() + facet_grid(. ~ groupe) #va nous donner deux graphes séparés pour chaque groupe
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun.y = mean,shape =17, size = 3, col="#CC0066") + facet_grid(. ~ groupe) # même graphe que précédemment mais avec la moyenne

# modele mixte sans varexplicative trithérapie

mod_pal_ui <- lmer(pal_ui_l~groupe+temps+(groupe*temps)+(1|id_patient),data=dl)
summary(mod_pal_ui)


# création d'une variable trithérapie

dl <- dl %>% mutate(tritherapie=ifelse(visite=="V3"|visite=="V4",1,0))


 
# modele mixte avec trithérapie      ----------------- QUE FAIRE ?



mod_pal_ui_tt <- lmer(pal_ui_l~temps+tritherapie+(1|id_patient),data=dl) # avec PAL UI
summary(mod_pal_ui_tt ) 


mod_pal_n <- lmer(pal_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) # avec PAL xN
summary(mod_pal_n)

#conditions de validité
par(mfrow=c(1,2))
hist(resid(mod_pal_ui_tt),col="cornflowerblue",main="Distribution du bruit",xlab="Résidus")
qqnorm(resid(mod_pal_ui_tt))
qqline(resid(mod_pal_ui_tt))



# les taux de GGT (xN)
mod_ggt_n <- lmer(ggt_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_ggt_n)

# ASAT (xN)
mod_asat_n <- lmer(asat_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_asat_n)

# ALAT (xN)
mod_alat_n <- lmer(alat_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_alat_n)

# Bilirubine totale
mod_bil <- lmer(bilirubine_µmol_l~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_bil)

# Echelle du prurit
mod_prurit <- lmer(prurit_echelle_0_10~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_prurit)

# Albumine
mod_albumine <- lmer(albumine_g_l~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_albumine)

# IgM 
mod_igm <- lmer(ig_m_g_l~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_igm)

# Creatinine
mod_creat <- lmer(creatinine_µmol_l~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_creat)


#TP
mod_tp <- lmer(ggt_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_tp)

#Plaquettes
mod_plaq <- lmer(ggt_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_plaq)

#Cholestérol 
mod_cho <- lmer(ggt_x_n~groupe+temps+(groupe*temps)+(1|id_patient)+tritherapie,data=dl) 
summary(mod_cho) 

names(dlong)
######################################################################################################################
##### création d'une base avec 2  groupes bithérapie et trithérapie débutant tous les deux à t0 pour la comparaison ##
######################################################################################################################
dl_bi <- dl %>% filter(visite=="V2"|visite=="V3")
dl_tri <- dl %>% filter(visite=="V3"|visite=="V4")

dl_tri <- dl_tri %>% group_by(id_patient) %>% mutate(temps_m = temps - lag(temps, default = first(temps))) %>% ungroup() %>% as.data.frame()
dl_bi$temps_m <- dl_bi$temps

# changement des ids des sujets du dataframe dl_tri (multiplié par 10)

dl_bi$id_patient_m <- dl_bi$id_patient
dl_bi$id_patient_m <- as.numeric(dl_bi$id_patient) 
dl_tri <- dl_tri %>% mutate(id_patient_m=as.numeric(id_patient)*100) 
dl_tri$id_patient_m <- as.numeric(dl_tri$id_patient_m)
dl_tri$tritherapie_m <- 1
dl_bi$tritherapie_m <- 0

# on regroupe les deux bases en une base dl_bi_tri pour faire nos analyses
dl_bi_tri <-  bind_rows(dl_bi,dl_tri) 
dl_bi_tri$id_patient_m <- as.factor(dl_bi_tri$id_patient_m)
dl_bi_tri$tritherapie_m<- as.factor(dl_bi_tri$tritherapie_m)

######## y a t'il une différence de pente entre le groupe trithérapie et le groupe bithérapie - modele mixte

##### PAL ####

mod_pal_ui_tt <- lmer(pal_ui_l~tritherapie_m+temps_m+tritherapie_m*temps+(1|id_patient),data=dl_bi_tri) # avec PAL UI
summary(mod_pal_ui_tt )

        