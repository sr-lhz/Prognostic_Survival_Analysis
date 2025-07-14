
############
### (1) IMPORT Input dataframe:

### Input dataframe:
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Clinical-Survival_Data")

dat.mv <- read.csv2("Survival-Multivariate-R_mCRC_200pts-QT_26May2020.csv",
                    header=TRUE,sep=";",stringsAsFactors=FALSE)
dim(dat.mv)
#[1] 200  35



#
#setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Clinical-Survival_Data/New_Databases-Adela_26May")

#dat.mv2 <- read.csv2("Survival-Multivariate-R_mCRC_200pts-QT_26May2020_AR.csv",
#                    header=TRUE,sep=";",stringsAsFactors=FALSE)
#dim(dat.mv2)
#[1] 200  35

#dat.mv2$Dcho_1.izq_2.recto_3 <- dat.mv$Location_Right1_Left2_Transverse3

#
#dat.mv$PFS_days <- dat.mv2$PFS_days
#dat.mv$OS_days <- dat.mv2$OS_days



############
### (2) FIX Input dataframe & Clinical variables:

### In Mutation-Status columns, convert 1/0 per "Mut"/"No_Mut":
dat.mv[,colnames(dat.mv) %in% c("KRAS","BRAF","FBXW7","NRAS","PIK3CA","PTEN","SMAD4","TP53")] <- ifelse(dat.mv[,colnames(dat.mv) %in% c("KRAS","BRAF","FBXW7","NRAS","PIK3CA","PTEN","SMAD4","TP53")]==1,"Mut","No_Mut")


### Create column "RAS" (which includes KRAS/NRAS):
dat.mv$RAS <- "No_Mut"

for(i in 1:nrow(dat.mv))
{
  if(dat.mv$KRAS[i]=="Mut" | dat.mv$NRAS[i]=="Mut")
  { dat.mv$RAS[i] <- "Mut" }
}


### CLASS of column-variables:
sapply(dat.mv,class)

dat.mv$Sex <- as.character(dat.mv$Sex)
dat.mv$Sex <- ifelse(dat.mv$Sex=="H","Man","Woman")

# We want this variable "Num_Lesions_Hep_1is1to3_2is4to9_3isMore9" to be CATEGORICAL:
dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9 <- as.character(dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9)
dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9 <- gsub("9","3",dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9)

dat.mv$Number_Organs <- as.character(dat.mv$Number_Organs)

dat.mv$Number_Organs_2 <- dat.mv$Number_Organs
dat.mv$Number_Organs_2 <- gsub("0","1",dat.mv$Number_Organs_2)
dat.mv$Number_Organs_2 <- ifelse(dat.mv$Number_Organs_2=="1","1","More")

dat.mv$Number_Organs <- as.numeric(dat.mv$Number_Organs)

#dat.mv$ECOG_PS_DxM1 <- as.character(dat.mv$ECOG_PS_DxM1)

dat.mv$Alk_Phosph_.116 <- ifelse(dat.mv$Alk_Phosph>=116,1,0)

dat.mv$LEUCOS_.11100 <- ifelse(dat.mv$LEUCOS>=11100,1,0)

dat.mv$PCR <- as.numeric(as.character(dat.mv$PCR))
dat.mv$PCR_.1 <- as.numeric(dat.mv$PCR_.1)
dat.mv$PCR_.1 <- ifelse(dat.mv$PCR>=1,1,0)

dat.mv$CEA <- as.numeric(as.character(dat.mv$CEA))
dat.mv$CEA_.5 <- ifelse(dat.mv$CEA>=5,1,0)

dat.mv$MSS0_MSI1 <- as.character(dat.mv$MSS0_MSI1)
dat.mv$MSS0_MSI1 <- ifelse(dat.mv$MSS0_MSI1=="0","MSS","MSI")

dat.mv$GEMCAD <- as.character(dat.mv$GEMCAD)

dat.mv$Location_Right1_Left2_Transverse3 <- as.character(dat.mv$Location_Right1_Left2_Transverse3)
dat.mv$Location_Right1_Left2_Transverse3 <- gsub("3","1",dat.mv$Location_Right1_Left2_Transverse3)
colnames(dat.mv)[34] <- "Location_Right1_Left2"


## FIX "Treatment" variable:
levels(factor(dat.mv$Treatment))
#[1] "CAPOX"        "FOLFIRI"      "FOLFIRI+BEV"  "FOLFIRI+CET" 
#[5] "FOLFIRI+PANI" "FOLFOX"       "FOLFOX+BEV"   "FOLFOX+CET"  
#[9] "FOLFOX+PANI"  "XELOX"        "XELOX+BEV"

for(i in 1:nrow(dat.mv))
{
  if(dat.mv$Treatment[i] %in% c("CAPOX","FOLFIRI","FOLFOX","XELOX"))
  { dat.mv$Treatment_Group[i] <- "Doublet" }
  
  else 
  { dat.mv$Treatment_Group[i] <- "Doublet_plus_Targeted-Agent" }
}

#
nrow(dat.mv[dat.mv$Treatment_Group=="Doublet",])
#[1] 139
nrow(dat.mv[dat.mv$Treatment_Group=="Doublet_plus_Targeted-Agent",])
#[1] 61



## COLUMNS (for PFS / OS):
#Variable 
#N events (%)
#Median survival (mo)

# Univariate:
#HR
#95% CI
#P-val (Univariate) 

# Multivariate:
#HR
#95% CI
#P-val (multivariate) 


# FUNCTION for obtaining MEDIAN-TIMES of Survival: 
time.unv.sex <- survfit(Surv(PFS_days,Progression) ~ Sex, 
                        data = dat.mv)


############
### ASSOCIATIONS of Mutations:

nrow(dat.mv[dat.mv$PIK3CA=="Mut" & dat.mv$FBXW7=="Mut",])
#[1] 4
nrow(dat.mv[dat.mv$SMAD4=="Mut" & dat.mv$FBXW7=="Mut",])
#[1] 4
nrow(dat.mv[dat.mv$PIK3CA=="Mut" & dat.mv$SMAD4=="Mut",])
#[1] 5


    ### 3 COMBINATIONS:
dat.mv$Combo1 <- ifelse(dat.mv$PIK3CA=="Mut" & dat.mv$FBXW7=="Mut","Double_PIK3CA-FBXW7","Else")
dat.mv$Combo2 <- ifelse(dat.mv$PIK3CA=="Mut" & dat.mv$SMAD4=="Mut","Double_PIK3CA-SMAD4","Else")
dat.mv$Combo3 <- ifelse(dat.mv$SMAD4=="Mut" & dat.mv$FBXW7=="Mut","Double_SMAD4-FBXW7","Else")

# COX Model for PFS:
pfs.combo1 <- coxph(Surv(PFS_days,Progression) ~ Combo1, 
                    data = dat.mv)
pfs.combo2 <- coxph(Surv(PFS_days,Progression) ~ Combo2, 
                    data = dat.mv)
pfs.combo3 <- coxph(Surv(PFS_days,Progression) ~ Combo3, 
                    data = dat.mv)

# COX Model for OS:
os.combo1 <- coxph(Surv(OS_days,Exitus) ~ Combo1, 
                    data = dat.mv)
os.combo2 <- coxph(Surv(OS_days,Exitus) ~ Combo2, 
                    data = dat.mv)
os.combo3 <- coxph(Surv(OS_days,Exitus) ~ Combo3, 
                    data = dat.mv)


# GLM (Logit) Model for SEX:
glm.sex <- glm(OS_days ~ Sex, data = dat.mv)

summary(glm.sex)



####################################################################################
####################################################################################
### (4) COMBINATIONS of GENES:
####################################################################################

  # FUNCTION for obtaining MEDIAN-TIMES of Survival:
time.unv.sex <- survfit(Surv(PFS_days,Progression) ~ Sex, 
                        data = dat.mv)

############
### (4.1) COMBO 1 [TP53+/SMAD4+ vs. TP53+/SMAD4-wt]:

dat.tp53 <- dat.mv[dat.mv$TP53=="Mut",]

#
pfs.smad4.tp53 <- survfit(Surv(PFS_days,Progression) ~ SMAD4, 
                          data = dat.tp53)

os.smad4.tp53 <- survfit(Surv(OS_days,Exitus) ~ SMAD4, 
                         data = dat.tp53)


### PFS-Univariate:
cox.unv.smad4.2 <- coxph(Surv(PFS_days,Progression) ~ SMAD4, 
                            data = dat.tp53)
unv.smad4.2 <- cox_as_data_frame(
  cox.unv.smad4.2,
  factor_id_sep = ": "
)
unv.smad4.2$Analysis <- "PFS_UNV_TP53-SMAD4"


### PFS-Multivariate:
cox.mv.smad4.2 <- coxph(Surv(PFS_days,Progression) ~ SMAD4+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                        data = dat.tp53)
mv.smad4.2 <- cox_as_data_frame(
  cox.mv.smad4.2,
  factor_id_sep = ": "
)
mv.smad4.2$Analysis <- "PFS_MV_TP53-SMAD4"


### OS-Univariate:
cox.unv.os.smad4.2 <- coxph(Surv(OS_days,Exitus) ~ SMAD4, 
                           data = dat.tp53)
unv.os.smad4.2 <- cox_as_data_frame(
  cox.unv.os.smad4.2,
  factor_id_sep = ": "
)
unv.os.smad4.2$Analysis <- "OS_UNV_TP53-SMAD4"


### OS-Multivariate:
cox.mv.os.smad4.2 <- coxph(Surv(OS_days,Exitus) ~ SMAD4+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.tp53)
mv.os.smad4.2 <- cox_as_data_frame(
  cox.mv.os.smad4.2,
  factor_id_sep = ": "
)
mv.os.smad4.2$Analysis <- "OS_MV_TP53-SMAD4"


### JOIN:
combo1 <- rbind(unv.smad4.2,mv.smad4.2[1,])
combo1 <- rbind(combo1,unv.os.smad4.2)
combo1 <- rbind(combo1,mv.os.smad4.2[1,])



############
### (4.2) COMBO 2 [TP53+/FBXW7+ vs. TP53+/FBXW7-wt]:

#
pfs.fb7.tp53 <- survfit(Surv(PFS_days,Progression) ~ FBXW7, 
                          data = dat.tp53)

os.fb7.tp53 <- survfit(Surv(OS_days,Exitus) ~ FBXW7, 
                         data = dat.tp53)

### PFS-Univariate:
cox.unv.fb7.2 <- coxph(Surv(PFS_days,Progression) ~ FBXW7, 
                         data = dat.tp53)
unv.fb7.2 <- cox_as_data_frame(
  cox.unv.fb7.2,
  factor_id_sep = ": "
)
unv.fb7.2$Analysis <- "PFS_UNV_TP53-FBXW7"


### PFS-Multivariate:
cox.mv.fb7.2 <- coxph(Surv(PFS_days,Progression) ~ FBXW7+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                        data = dat.tp53)
mv.fb7.2 <- cox_as_data_frame(
  cox.mv.fb7.2,
  factor_id_sep = ": "
)
mv.fb7.2$Analysis <- "PFS_MV_TP53-FBXW7"


### OS-Univariate:
cox.unv.os.fb7.2 <- coxph(Surv(OS_days,Exitus) ~ FBXW7, 
                            data = dat.tp53)
unv.os.fb7.2 <- cox_as_data_frame(
  cox.unv.os.fb7.2,
  factor_id_sep = ": "
)
unv.os.fb7.2$Analysis <- "OS_UNV_TP53-FBXW7"


### OS-Multivariate:
cox.mv.os.fb7.2 <- coxph(Surv(OS_days,Exitus) ~ FBXW7+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                           data = dat.tp53)
mv.os.fb7.2 <- cox_as_data_frame(
  cox.mv.os.fb7.2,
  factor_id_sep = ": "
)
mv.os.fb7.2$Analysis <- "OS_MV_TP53-FBXW7"


### JOIN:
combo2 <- rbind(unv.fb7.2,mv.fb7.2[1,])
combo2 <- rbind(combo2,unv.os.fb7.2)
combo2 <- rbind(combo2,mv.os.fb7.2[1,])


############
### (4.3) COMBO 3 [TP53+/SMAD4+ or TP53+/FBXW7+ vs. TP53+/Else]:

#
for(i in 1:nrow(dat.tp53)){
  if(dat.tp53$SMAD4[i]=="Mut" | dat.tp53$FBXW7[i]=="Mut")
  { dat.tp53$Groups_Combo3[i] <- "TP53_SMAD4-or-FBXW7" }
  
  else
  { dat.tp53$Groups_Combo3[i] <- "TP53_Else" }
}

  ## 30 August 2020 (for Poster ESMO 2020):
dat.tp53 <- dat.tp53[dat.tp53$OS_days<1400,]


#
pfs.combo3 <- survfit(Surv(PFS_days,Progression) ~ Groups_Combo3, 
                        data = dat.tp53)

os.combo3 <- survfit(Surv(OS_days,Exitus) ~ Groups_Combo3, 
                       data = dat.tp53)

### PFS-Univariate:
cox.unv.combo3 <- coxph(Surv(PFS_days,Progression) ~ Groups_Combo3, 
                       data = dat.tp53)
unv.combo3 <- cox_as_data_frame(
  cox.unv.combo3,
  factor_id_sep = ": "
)
unv.combo3$Analysis <- "PFS_UNV_Combo3"


### PFS-Multivariate:
cox.mv.combo3 <- coxph(Surv(PFS_days,Progression) ~ Groups_Combo3+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                      data = dat.tp53)
mv.combo3 <- cox_as_data_frame(
  cox.mv.combo3,
  factor_id_sep = ": "
)
mv.combo3$Analysis <- "PFS_MV_Combo3"


### OS-Univariate:
cox.unv.os.combo3 <- coxph(Surv(OS_days,Exitus) ~ Groups_Combo3, 
                          data = dat.tp53)
unv.os.combo3 <- cox_as_data_frame(
  cox.unv.os.combo3,
  factor_id_sep = ": "
)
unv.os.combo3$Analysis <- "OS_UNV_Combo3"

# Survfit object for ggsurvplot:
survfit.combo3 <- survfit(Surv(OS_days,Exitus) ~ Groups_Combo3, 
                           data = dat.tp53)

### OS-Multivariate:
cox.mv.os.combo3 <- coxph(Surv(OS_days,Exitus) ~ Groups_Combo3+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                         data = dat.tp53)
mv.os.combo3 <- cox_as_data_frame(
  cox.mv.os.combo3,
  factor_id_sep = ": "
)
mv.os.combo3$Analysis <- "OS_MV_Combo3"


### JOIN:
combo3 <- rbind(unv.combo3,mv.combo3[1,])
combo3 <- rbind(combo3,unv.os.combo3)
combo3 <- rbind(combo3,mv.os.combo3[1,])

#
ggsurvplot(survfit.combo3, #pval = TRUE,
           #risk.table = TRUE,   #risk.table="percentage",
           conf.int=FALSE,   #risk.table.row=c(0,200,400,600,800,1000,2000,3000,4000,5000),
           palette="lancet", #title= "Progression-Free Survival (PFS) Probability",
           break.time.by=350,
           xlab="Time elapsed (days)",ylab="Accumulated overall survival")#,
           #surv.median.line=c("h","v"))#,
#legend="right",legend.labs=c("QT-Doublet alone","QT-Doublet plus Targeted Agent")#,
#)




############
### (4.4) COMBO 4 [SMAD4+ or FBXW7+ vs. Else]:

#
for(i in 1:nrow(dat.mv)){
  if(dat.mv$SMAD4[i]=="Mut" | dat.mv$FBXW7[i]=="Mut")
  { dat.mv$Groups_Combo4[i] <- "SMAD4-or-FBXW7" }
  
  else
  { dat.mv$Groups_Combo4[i] <- "Else" }
}

#
pfs.combo4 <- survfit(Surv(PFS_days,Progression) ~ Groups_Combo4, 
                      data = dat.mv)

os.combo4 <- survfit(Surv(OS_days,Exitus) ~ Groups_Combo4, 
                     data = dat.mv)

### PFS-Univariate:
cox.unv.combo4 <- coxph(Surv(PFS_days,Progression) ~ Groups_Combo4, 
                        data = dat.mv)
unv.combo4 <- cox_as_data_frame(
  cox.unv.combo4,
  factor_id_sep = ": "
)
unv.combo4$Analysis <- "PFS_UNV_Combo4"


### PFS-Multivariate:
cox.mv.combo4 <- coxph(Surv(PFS_days,Progression) ~ Groups_Combo4+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                       data = dat.mv)
mv.combo4 <- cox_as_data_frame(
  cox.mv.combo4,
  factor_id_sep = ": "
)
mv.combo4$Analysis <- "PFS_MV_Combo4"


### OS-Univariate:
cox.unv.os.combo4 <- coxph(Surv(OS_days,Exitus) ~ Groups_Combo4, 
                           data = dat.mv)
unv.os.combo4 <- cox_as_data_frame(
  cox.unv.os.combo4,
  factor_id_sep = ": "
)
unv.os.combo4$Analysis <- "OS_UNV_Combo4"


### OS-Multivariate:
cox.mv.os.combo4 <- coxph(Surv(OS_days,Exitus) ~ Groups_Combo4+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                          data = dat.mv)
mv.os.combo4 <- cox_as_data_frame(
  cox.mv.os.combo4,
  factor_id_sep = ": "
)
mv.os.combo4$Analysis <- "OS_MV_Combo4"


### JOIN:
combo4 <- rbind(unv.combo4,mv.combo4[1,])
combo4 <- rbind(combo4,unv.os.combo4)
combo4 <- rbind(combo4,mv.os.combo4[1,])



############
### (4.5) COMBO 5 [TP53+/SMAD4+ or TP53+/FBXW7+ vs. Else]:

#
for(i in 1:nrow(dat.mv)){
  if(dat.mv$TP53[i]=="Mut" & dat.mv$SMAD4[i]=="Mut")
  { dat.mv$Groups_Combo5[i] <- "Double_Mut" }
  
  else if(dat.mv$TP53[i]=="Mut" & dat.mv$FBXW7[i]=="Mut")
  { dat.mv$Groups_Combo5[i] <- "Double_Mut" }
  
  else
  { dat.mv$Groups_Combo5[i] <- "Else" }
}

#
pfs.combo5 <- survfit(Surv(PFS_days,Progression) ~ Groups_Combo5, 
                      data = dat.mv)

os.combo5 <- survfit(Surv(OS_days,Exitus) ~ Groups_Combo5, 
                     data = dat.mv)

### PFS-Univariate:
cox.unv.combo5 <- coxph(Surv(PFS_days,Progression) ~ Groups_Combo5, 
                        data = dat.mv)
unv.combo5 <- cox_as_data_frame(
  cox.unv.combo5,
  factor_id_sep = ": "
)
unv.combo5$Analysis <- "PFS_UNV_Combo5"


### PFS-Multivariate:
cox.mv.combo5 <- coxph(Surv(PFS_days,Progression) ~ Groups_Combo5+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                       data = dat.mv)
mv.combo5 <- cox_as_data_frame(
  cox.mv.combo5,
  factor_id_sep = ": "
)
mv.combo5$Analysis <- "PFS_MV_Combo5"


### OS-Univariate:
cox.unv.os.combo5 <- coxph(Surv(OS_days,Exitus) ~ Groups_Combo5, 
                           data = dat.mv)
unv.os.combo5 <- cox_as_data_frame(
  cox.unv.os.combo5,
  factor_id_sep = ": "
)
unv.os.combo5$Analysis <- "OS_UNV_Combo5"


### OS-Multivariate:
cox.mv.os.combo5 <- coxph(Surv(OS_days,Exitus) ~ Groups_Combo5+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                          data = dat.mv)
mv.os.combo5 <- cox_as_data_frame(
  cox.mv.os.combo5,
  factor_id_sep = ": "
)
mv.os.combo5$Analysis <- "OS_MV_Combo5"


### JOIN:
combo5 <- rbind(unv.combo5,mv.combo5[1,])
combo5 <- rbind(combo5,unv.os.combo5)
combo5 <- rbind(combo5,mv.os.combo5[1,])



##################
### FINAL JOINING:

combos.all <- rbind(combo1,combo2)
combos.all <- rbind(combos.all,combo3)
combos.all <- rbind(combos.all,combo4)
combos.all <- rbind(combos.all,combo5)

#
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Survival_Analyses/Outputs/14June2020_Poster_Gene-COMBINATIONS")

WriteXLS(combos.all,"Table-Survival_mCRC_Uni-Multi-Variate_14June2020_Lahoz.xls",
         col.names=TRUE,row.names=FALSE)



####################################################################################
####################################################################################
### (3) INDIVIDUAL VARIABLES (Individual Genes + Clinical):
####################################################################################

############
### (3.1) UNIVARIATE-PFS:

#Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+MSS0_MSI1+Location_Right1_Left2

### Sex:
cox.unv.sex <- coxph(Surv(PFS_days,Progression) ~ Sex, 
                     data = dat.mv)
unv.sex <- cox_as_data_frame(
  cox.unv.sex,
  factor_id_sep = ": "
)
unv.sex$Analysis <- "PFS_UNV_Sex"


### Age:
cox.unv.age <- coxph(Surv(PFS_days,Progression) ~ Age_Dx.M1, 
                     data = dat.mv)
unv.age <- cox_as_data_frame(
  cox.unv.age,
  factor_id_sep = ": "
)
unv.age$Analysis <- "PFS_UNV_Age"


### Number of Organs:
cox.unv.norg <- coxph(Surv(PFS_days,Progression) ~ Number_Organs_2, 
                      data = dat.mv)
unv.norg <- cox_as_data_frame(
  cox.unv.norg,
  factor_id_sep = ": "
)
unv.norg$Analysis <- "PFS_UNV_Number-Organs"


### ECOG PS:
cox.unv.ecog <- coxph(Surv(PFS_days,Progression) ~ ECOG_PS_DxM1, 
                      data = dat.mv)
unv.ecog <- cox_as_data_frame(
  cox.unv.ecog,
  factor_id_sep = ": "
)
unv.ecog$Analysis <- "PFS_UNV_ECOG"


### LDH: 
cox.unv.ldh <- coxph(Surv(PFS_days,Progression) ~ LDH, 
                     data = dat.mv)
unv.ldh <- cox_as_data_frame(
  cox.unv.ldh,
  factor_id_sep = ": "
)
unv.ldh$Analysis <- "PFS_UNV_LDH"


### Alk_Phosph:
cox.unv.alp <- coxph(Surv(PFS_days,Progression) ~ Alk_Phosph, 
                     data = dat.mv)
unv.alp <- cox_as_data_frame(
  cox.unv.alp,
  factor_id_sep = ": "
)
unv.alp$Analysis <- "PFS_UNV_ALP"


### LEUCOS:
cox.unv.leucos <- coxph(Surv(PFS_days,Progression) ~ LEUCOS, 
                        data = dat.mv)
unv.leucos <- cox_as_data_frame(
  cox.unv.leucos,
  factor_id_sep = ": "
)
unv.leucos$Analysis <- "PFS_UNV_Leucos"


### PCR:
cox.unv.pcr <- coxph(Surv(PFS_days,Progression) ~ PCR, 
                     data = dat.mv)
unv.pcr <- cox_as_data_frame(
  cox.unv.pcr,
  factor_id_sep = ": "
)
unv.pcr$Analysis <- "PFS_UNV_PCR"


### MSI Status:
cox.unv.msi <- coxph(Surv(PFS_days,Progression) ~ MSS0_MSI1, 
                     data = dat.mv)
unv.msi <- cox_as_data_frame(
  cox.unv.msi,
  factor_id_sep = ": "
)
unv.msi$Analysis <- "PFS_UNV_MSI"


### Location:
cox.unv.loc <- coxph(Surv(PFS_days,Progression) ~ Location_Right1_Left2, 
                     data = dat.mv)
unv.loc <- cox_as_data_frame(
  cox.unv.loc,
  factor_id_sep = ": "
)
unv.loc$Analysis <- "PFS_UNV_Location"


### TP53:
cox.unv.tp53 <- coxph(Surv(PFS_days,Progression) ~ TP53, 
                      data = dat.mv)
unv.tp53 <- cox_as_data_frame(
  cox.unv.tp53,
  factor_id_sep = ": "
)
unv.tp53$Analysis <- "PFS_UNV_TP53"


### RAS:
cox.unv.ras <- coxph(Surv(PFS_days,Progression) ~ RAS, 
                     data = dat.mv)
unv.ras <- cox_as_data_frame(
  cox.unv.ras,
  factor_id_sep = ": "
)
unv.ras$Analysis <- "PFS_UNV_RAS"


### PIK3CA:
cox.unv.pik3 <- coxph(Surv(PFS_days,Progression) ~ PIK3CA, 
                      data = dat.mv)
unv.pik3 <- cox_as_data_frame(
  cox.unv.pik3,
  factor_id_sep = ": "
)
unv.pik3$Analysis <- "PFS_UNV_PIK3CA"


### SMAD4:
cox.unv.smad4 <- coxph(Surv(PFS_days,Progression) ~ SMAD4, 
                       data = dat.mv)
unv.smad4 <- cox_as_data_frame(
  cox.unv.smad4,
  factor_id_sep = ": "
)
unv.smad4$Analysis <- "PFS_UNV_SMAD4"


### FBXW7:
cox.unv.fb7 <- coxph(Surv(PFS_days,Progression) ~ FBXW7, 
                     data = dat.mv)
unv.fb7 <- cox_as_data_frame(
  cox.unv.fb7,
  factor_id_sep = ": "
)
unv.fb7$Analysis <- "PFS_UNV_FBXW7"


### BRAF:
cox.unv.braf <- coxph(Surv(PFS_days,Progression) ~ BRAF, 
                      data = dat.mv)
unv.braf <- cox_as_data_frame(
  cox.unv.braf,
  factor_id_sep = ": "
)
unv.braf$Analysis <- "PFS_UNV_BRAF"


###
unv.pfs.all <- rbind(unv.sex,unv.age)
unv.pfs.all <- rbind(unv.pfs.all,unv.norg)
unv.pfs.all <- rbind(unv.pfs.all,unv.ecog)
unv.pfs.all <- rbind(unv.pfs.all,unv.ldh)
unv.pfs.all <- rbind(unv.pfs.all,unv.alp) 
unv.pfs.all <- rbind(unv.pfs.all,unv.leucos)
unv.pfs.all <- rbind(unv.pfs.all,unv.pcr)
unv.pfs.all <- rbind(unv.pfs.all,unv.msi)
unv.pfs.all <- rbind(unv.pfs.all,unv.loc)
unv.pfs.all <- rbind(unv.pfs.all,unv.tp53)
unv.pfs.all <- rbind(unv.pfs.all,unv.ras)
unv.pfs.all <- rbind(unv.pfs.all,unv.pik3)
unv.pfs.all <- rbind(unv.pfs.all,unv.smad4)
unv.pfs.all <- rbind(unv.pfs.all,unv.fb7)
unv.pfs.all <- rbind(unv.pfs.all,unv.braf)



############
### (3.2) UNIVARIATE-OS:

### Sex:
cox.unv.sex <- coxph(Surv(OS_days,Exitus) ~ Sex, 
                     data = dat.mv)
os.unv.sex <- cox_as_data_frame(
  cox.unv.sex,
  factor_id_sep = ": "
)
os.unv.sex$Analysis <- "OS_UNV_Sex"


### Age:
cox.unv.age <- coxph(Surv(OS_days,Exitus) ~ Age_Dx.M1, 
                     data = dat.mv)
os.unv.age <- cox_as_data_frame(
  cox.unv.age,
  factor_id_sep = ": "
)
os.unv.age$Analysis <- "OS_UNV_Age"


### Number of Organs:
cox.unv.norg <- coxph(Surv(OS_days,Exitus) ~ Number_Organs_2, 
                      data = dat.mv)
os.unv.norg <- cox_as_data_frame(
  cox.unv.norg,
  factor_id_sep = ": "
)
os.unv.norg$Analysis <- "OS_UNV_Number-Organs"


### ECOG PS:
cox.unv.ecog <- coxph(Surv(OS_days,Exitus) ~ ECOG_PS_DxM1, 
                      data = dat.mv)
os.unv.ecog <- cox_as_data_frame(
  cox.unv.ecog,
  factor_id_sep = ": "
)
os.unv.ecog$Analysis <- "OS_UNV_ECOG"


### LDH:
cox.unv.ldh <- coxph(Surv(OS_days,Exitus) ~ LDH, 
                     data = dat.mv)
os.unv.ldh <- cox_as_data_frame(
  cox.unv.ldh,
  factor_id_sep = ": "
)
os.unv.ldh$Analysis <- "OS_UNV_LDH"


### Alk_Phosph:
cox.unv.alp <- coxph(Surv(OS_days,Exitus) ~ Alk_Phosph, 
                     data = dat.mv)
os.unv.alp <- cox_as_data_frame(
  cox.unv.alp,
  factor_id_sep = ": "
)
os.unv.alp$Analysis <- "OS_UNV_ALP"


### LEUCOS:
cox.unv.leucos <- coxph(Surv(OS_days,Exitus) ~ LEUCOS, 
                        data = dat.mv)
os.unv.leucos <- cox_as_data_frame(
  cox.unv.leucos,
  factor_id_sep = ": "
)
os.unv.leucos$Analysis <- "OS_UNV_Leucos"


### PCR:
cox.unv.pcr <- coxph(Surv(OS_days,Exitus) ~ PCR, 
                     data = dat.mv)
os.unv.pcr <- cox_as_data_frame(
  cox.unv.pcr,
  factor_id_sep = ": "
)
os.unv.pcr$Analysis <- "OS_UNV_PCR"


### MSI Status:
cox.unv.msi <- coxph(Surv(OS_days,Exitus) ~ MSS0_MSI1, 
                     data = dat.mv)
os.unv.msi <- cox_as_data_frame(
  cox.unv.msi,
  factor_id_sep = ": "
)
os.unv.msi$Analysis <- "OS_UNV_MSI"


### Location:
cox.unv.loc <- coxph(Surv(OS_days,Exitus) ~ Location_Right1_Left2, 
                     data = dat.mv)
os.unv.loc <- cox_as_data_frame(
  cox.unv.loc,
  factor_id_sep = ": "
)
os.unv.loc$Analysis <- "OS_UNV_Location"


### TP53:
cox.unv.tp53 <- coxph(Surv(OS_days,Exitus) ~ TP53, 
                      data = dat.mv)
os.unv.tp53 <- cox_as_data_frame(
  cox.unv.tp53,
  factor_id_sep = ": "
)
os.unv.tp53$Analysis <- "OS_UNV_TP53"


### RAS:
cox.unv.ras <- coxph(Surv(OS_days,Exitus) ~ RAS, 
                     data = dat.mv)
os.unv.ras <- cox_as_data_frame(
  cox.unv.ras,
  factor_id_sep = ": "
)
os.unv.ras$Analysis <- "OS_UNV_RAS"


### PIK3CA:
cox.unv.pik3 <- coxph(Surv(OS_days,Exitus) ~ PIK3CA, 
                      data = dat.mv)
os.unv.pik3 <- cox_as_data_frame(
  cox.unv.pik3,
  factor_id_sep = ": "
)
os.unv.pik3$Analysis <- "OS_UNV_PIK3CA"


### SMAD4:
cox.unv.smad4 <- coxph(Surv(OS_days,Exitus) ~ SMAD4, 
                       data = dat.mv)
os.unv.smad4 <- cox_as_data_frame(
  cox.unv.smad4,
  factor_id_sep = ": "
)
os.unv.smad4$Analysis <- "OS_UNV_SMAD4"


### FBXW7:
cox.unv.fb7 <- coxph(Surv(OS_days,Exitus) ~ FBXW7, 
                     data = dat.mv)
os.unv.fb7 <- cox_as_data_frame(
  cox.unv.fb7,
  factor_id_sep = ": "
)
os.unv.fb7$Analysis <- "OS_UNV_FBXW7"


### BRAF:
cox.unv.braf <- coxph(Surv(OS_days,Exitus) ~ BRAF, 
                      data = dat.mv)
os.unv.braf <- cox_as_data_frame(
  cox.unv.braf,
  factor_id_sep = ": "
)
os.unv.braf$Analysis <- "OS_UNV_BRAF"


###
unv.os.all <- rbind(os.unv.sex,os.unv.age)
unv.os.all <- rbind(unv.os.all,os.unv.norg)
unv.os.all <- rbind(unv.os.all,os.unv.ecog)
unv.os.all <- rbind(unv.os.all,os.unv.ldh)
unv.os.all <- rbind(unv.os.all,os.unv.alp) 
unv.os.all <- rbind(unv.os.all,os.unv.leucos)
unv.os.all <- rbind(unv.os.all,os.unv.pcr)
unv.os.all <- rbind(unv.os.all,os.unv.msi)
unv.os.all <- rbind(unv.os.all,os.unv.loc)
unv.os.all <- rbind(unv.os.all,os.unv.tp53)
unv.os.all <- rbind(unv.os.all,os.unv.ras)
unv.os.all <- rbind(unv.os.all,os.unv.pik3)
unv.os.all <- rbind(unv.os.all,os.unv.smad4)
unv.os.all <- rbind(unv.os.all,os.unv.fb7)
unv.os.all <- rbind(unv.os.all,os.unv.braf)



############
### (4.1) MULTIVARIATE-PFS:

#Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+MSS0_MSI1+Location_Right1_Left2

### Sex:
cox.mv.clin <- coxph(Surv(PFS_days,Progression) ~ Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
mv.clin <- cox_as_data_frame(
  cox.mv.clin,
  factor_id_sep = ": "
)
mv.clin$Analysis <- "PFS_MV_Clinical"


### TP53:
cox.mv.tp53 <- coxph(Surv(PFS_days,Progression) ~ TP53+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
mv.tp53 <- cox_as_data_frame(
  cox.mv.tp53,
  factor_id_sep = ": "
)
mv.tp53$Analysis <- "PFS_MV_TP53"


### RAS:
cox.mv.ras <- coxph(Surv(PFS_days,Progression) ~ RAS+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                    data = dat.mv)
mv.ras <- cox_as_data_frame(
  cox.mv.ras,
  factor_id_sep = ": "
)
mv.ras$Analysis <- "PFS_MV_RAS"


### PIK3CA:
cox.mv.pik3 <- coxph(Surv(PFS_days,Progression) ~ PIK3CA+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
mv.pik3 <- cox_as_data_frame(
  cox.mv.pik3,
  factor_id_sep = ": "
)
mv.pik3$Analysis <- "PFS_MV_PIK3CA"


### SMAD4:
cox.mv.smad4 <- coxph(Surv(PFS_days,Progression) ~ SMAD4+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                      data = dat.mv)
mv.smad4 <- cox_as_data_frame(
  cox.mv.smad4,
  factor_id_sep = ": "
)
mv.smad4$Analysis <- "PFS_MV_SMAD4"


### FBXW7:
cox.mv.fb7 <- coxph(Surv(PFS_days,Progression) ~ FBXW7+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                    data = dat.mv)
mv.fb7 <- cox_as_data_frame(
  cox.mv.fb7,
  factor_id_sep = ": "
)
mv.fb7$Analysis <- "PFS_MV_FBXW7"


### BRAF:
cox.mv.braf <- coxph(Surv(PFS_days,Progression) ~ BRAF+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
mv.braf <- cox_as_data_frame(
  cox.mv.braf,
  factor_id_sep = ": "
)
mv.braf$Analysis <- "PFS_MV_BRAF"


###
mv.pfs.all <- rbind(mv.clin,mv.tp53[1,])
mv.pfs.all <- rbind(mv.pfs.all,mv.ras[1,])
mv.pfs.all <- rbind(mv.pfs.all,mv.pik3[1,])
mv.pfs.all <- rbind(mv.pfs.all,mv.smad4[1,])
mv.pfs.all <- rbind(mv.pfs.all,mv.fb7[1,])
mv.pfs.all <- rbind(mv.pfs.all,mv.braf[1,])



############
### (4.2) MULTIVARIATE-OS:

### CLINICAL Variables:
cox.mv.clin <- coxph(Surv(OS_days,Exitus) ~ Sex+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
os.mv.clin <- cox_as_data_frame(
  cox.mv.clin,
  factor_id_sep = ": "
)
os.mv.clin$Analysis <- "OS_MV_Clinical"


### TP53:
cox.mv.tp53 <- coxph(Surv(OS_days,Exitus) ~ TP53+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
os.mv.tp53 <- cox_as_data_frame(
  cox.mv.tp53,
  factor_id_sep = ": "
)
os.mv.tp53$Analysis <- "OS_MV_TP53"


### RAS:
cox.mv.ras <- coxph(Surv(OS_days,Exitus) ~ RAS+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                    data = dat.mv)
os.mv.ras <- cox_as_data_frame(
  cox.mv.ras,
  factor_id_sep = ": "
)
os.mv.ras$Analysis <- "OS_MV_RAS"


### PIK3CA:
cox.mv.pik3 <- coxph(Surv(OS_days,Exitus) ~ PIK3CA+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
os.mv.pik3 <- cox_as_data_frame(
  cox.mv.pik3,
  factor_id_sep = ": "
)
os.mv.pik3$Analysis <- "OS_MV_PIK3CA"


### SMAD4:
cox.mv.smad4 <- coxph(Surv(OS_days,Exitus) ~ SMAD4+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                      data = dat.mv)
os.mv.smad4 <- cox_as_data_frame(
  cox.mv.smad4,
  factor_id_sep = ": "
)
os.mv.smad4$Analysis <- "OS_MV_SMAD4"


### FBXW7:
cox.mv.fb7 <- coxph(Surv(OS_days,Exitus) ~ FBXW7+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                    data = dat.mv)
os.mv.fb7 <- cox_as_data_frame(
  cox.mv.fb7,
  factor_id_sep = ": "
)
os.mv.fb7$Analysis <- "OS_MV_FBXW7"


### BRAF:
cox.mv.braf <- coxph(Surv(OS_days,Exitus) ~ BRAF+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                     data = dat.mv)
os.mv.braf <- cox_as_data_frame(
  cox.mv.braf,
  factor_id_sep = ": "
)
os.mv.braf$Analysis <- "OS_MV_BRAF"


###
mv.os.all <- rbind(os.mv.clin,os.mv.tp53[1,])
mv.os.all <- rbind(mv.os.all,os.mv.ras[1,])
mv.os.all <- rbind(mv.os.all,os.mv.pik3[1,])
mv.os.all <- rbind(mv.os.all,os.mv.smad4[1,])
mv.os.all <- rbind(mv.os.all,os.mv.fb7[1,])
mv.os.all <- rbind(mv.os.all,os.mv.braf[1,])



###### ###### ###### ###### ###### ###### 
### FINAL-JOINING of 4 main dataframes:
surv.all <- rbind(unv.pfs.all,unv.os.all)
surv.all <- rbind(surv.all,mv.pfs.all)
surv.all <- rbind(surv.all,mv.os.all)

dim(surv.all)
#[1] 64 11


#
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Survival_Analyses/Outputs/2June2020_Table-Poster-Adela")

WriteXLS(surv.all,"Table-Survival_mCRC_Uni-Multi-Variate_3June2020_Lahoz.xls",     
         col.names=TRUE,row.names=FALSE)







