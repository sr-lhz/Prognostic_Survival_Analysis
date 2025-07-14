############
### (5.1) UNIVARIATE-PFS --- TP53 (baseline) w/ other 5 genes:

### TP53 & RAS (Univariate):
pfs.unv.tp53.ras <- dat.mv[dat.mv$TP53=="Mut",]
pfs.unv.tp53.ras$Groups <- ifelse(pfs.unv.tp53.ras$RAS=="Mut","TP53+/RAS+","TP53+")

cox.pfs.unv.tp53.ras <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                              data = pfs.unv.tp53.ras)
df.pfs.unv.tp53.ras <- cox_as_data_frame(
  cox.pfs.unv.tp53.ras,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.tp53.ras$Analysis <- "PFS_UNV_TP53_RAS"

#
surv.sub.unv <- survfit(Surv(PFS_days,Progression) ~ Groups, data = pfs.unv.tp53.ras)
ggsurvplot(surv.sub.unv, pval = TRUE,
           #risk.table = TRUE,   # risk.table = TRUE, #risk.table="percentage",
           conf.int=FALSE,   #risk.table.row=c(0,200,400,600,800,1000,2000,3000,4000,5000),
           palette="lancet", title= "Progression-Free Survival (PFS) Probability",
           break.time.by=500,xlab="Time elapsed (days)",ylab="Accumulated survival",
           legend="right",legend.labs=c("TP53+/RAS-wt","TP53+/RAS+"))#,
#pval.method=TRUE)



### TP53 & SMAD4 (Univariate):
pfs.unv.tp53.smad4 <- dat.mv[dat.mv$TP53=="Mut",]
pfs.unv.tp53.smad4$Groups <- ifelse(pfs.unv.tp53.smad4$SMAD4=="Mut","TP53+/SMAD4+","TP53+")

cox.pfs.unv.tp53.smad4 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                                data = pfs.unv.tp53.smad4)
df.pfs.unv.tp53.smad4 <- cox_as_data_frame(
  cox.pfs.unv.tp53.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.tp53.smad4$Analysis <- "PFS_UNV_TP53_SMAD4"


### TP53 & FBXW7 (Univariate):
pfs.unv.tp53.fb7 <- dat.mv[dat.mv$TP53=="Mut",]
pfs.unv.tp53.fb7$Groups <- ifelse(pfs.unv.tp53.fb7$FBXW7=="Mut","TP53+/FBXW7+","TP53+")

cox.pfs.unv.tp53.fb7 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                              data = pfs.unv.tp53.fb7)
df.pfs.unv.tp53.fb7 <- cox_as_data_frame(
  cox.pfs.unv.tp53.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.tp53.fb7$Analysis <- "PFS_UNV_TP53_FBXW7"


### TP53 & PIK3CA (Univariate):
pfs.unv.tp53.pik3 <- dat.mv[dat.mv$TP53=="Mut",]
pfs.unv.tp53.pik3$Groups <- ifelse(pfs.unv.tp53.pik3$PIK3CA=="Mut","TP53+/PIK3CA+","TP53+")

cox.pfs.unv.tp53.pik3 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                               data = pfs.unv.tp53.pik3)
df.pfs.unv.tp53.pik3 <- cox_as_data_frame(
  cox.pfs.unv.tp53.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.tp53.pik3$Analysis <- "PFS_UNV_TP53_PIK3CA"


### TP53 & BRAF (Univariate):
pfs.unv.tp53.braf <- dat.mv[dat.mv$TP53=="Mut",]
pfs.unv.tp53.braf$Groups <- ifelse(pfs.unv.tp53.braf$BRAF=="Mut","TP53+/BRAF+","TP53+")

cox.pfs.unv.tp53.braf <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                               data = pfs.unv.tp53.braf)
df.pfs.unv.tp53.braf <- cox_as_data_frame(
  cox.pfs.unv.tp53.braf,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.tp53.braf$Analysis <- "PFS_UNV_TP53_BRAF"


### JOINT dataframe "unv.pfs.tp53":
unv.pfs.tp53 <- rbind(df.pfs.unv.tp53.ras,df.pfs.unv.tp53.smad4)
unv.pfs.tp53 <- rbind(unv.pfs.tp53,df.pfs.unv.tp53.fb7)
unv.pfs.tp53 <- rbind(unv.pfs.tp53,df.pfs.unv.tp53.pik3)
unv.pfs.tp53 <- rbind(unv.pfs.tp53,df.pfs.unv.tp53.braf)



############
### (5.2) UNIVARIATE-OS --- TP53 (baseline) w/ other 5 genes:

### TP53 & RAS (Univariate):
os.unv.tp53.ras <- dat.mv[dat.mv$TP53=="Mut",]
os.unv.tp53.ras$Groups <- ifelse(os.unv.tp53.ras$RAS=="Mut","TP53+/RAS+","TP53+")

cox.os.unv.tp53.ras <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                             data = os.unv.tp53.ras)
df.os.unv.tp53.ras <- cox_as_data_frame(
  cox.os.unv.tp53.ras,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.tp53.ras$Analysis <- "OS_UNV_TP53_RAS"


### TP53 & SMAD4 (Univariate):
os.unv.tp53.smad4 <- dat.mv[dat.mv$TP53=="Mut",]
os.unv.tp53.smad4$Groups <- ifelse(os.unv.tp53.smad4$SMAD4=="Mut","TP53+/SMAD4+","TP53+")

cox.os.unv.tp53.smad4 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                               data = os.unv.tp53.smad4)
df.os.unv.tp53.smad4 <- cox_as_data_frame(
  cox.os.unv.tp53.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.tp53.smad4$Analysis <- "OS_UNV_TP53_SMAD4"


### TP53 & FBXW7 (Univariate):
os.unv.tp53.fb7 <- dat.mv[dat.mv$TP53=="Mut",]
os.unv.tp53.fb7$Groups <- ifelse(os.unv.tp53.fb7$FBXW7=="Mut","TP53+/FBXW7+","TP53+")

cox.os.unv.tp53.fb7 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                             data = os.unv.tp53.fb7)
df.os.unv.tp53.fb7 <- cox_as_data_frame(
  cox.os.unv.tp53.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.tp53.fb7$Analysis <- "OS_UNV_TP53_FBXW7"


### TP53 & PIK3CA (Univariate):
os.unv.tp53.pik3 <- dat.mv[dat.mv$TP53=="Mut",]
os.unv.tp53.pik3$Groups <- ifelse(os.unv.tp53.pik3$PIK3CA=="Mut","TP53+/PIK3CA+","TP53+")

cox.os.unv.tp53.pik3 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                              data = os.unv.tp53.pik3)
df.os.unv.tp53.pik3 <- cox_as_data_frame(
  cox.os.unv.tp53.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.tp53.pik3$Analysis <- "OS_UNV_TP53_PIK3CA"


### TP53 & BRAF (Univariate):
os.unv.tp53.braf <- dat.mv[dat.mv$TP53=="Mut",]
os.unv.tp53.braf$Groups <- ifelse(os.unv.tp53.braf$BRAF=="Mut","TP53+/BRAF+","TP53+")

cox.os.unv.tp53.braf <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                              data = os.unv.tp53.braf)
df.os.unv.tp53.braf <- cox_as_data_frame(
  cox.os.unv.tp53.braf,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.tp53.braf$Analysis <- "OS_UNV_TP53_BRAF"


### JOINT dataframe "unv.os.tp53":
unv.os.tp53 <- rbind(df.os.unv.tp53.ras,df.os.unv.tp53.smad4)
unv.os.tp53 <- rbind(unv.os.tp53,df.os.unv.tp53.fb7)
unv.os.tp53 <- rbind(unv.os.tp53,df.os.unv.tp53.pik3)
unv.os.tp53 <- rbind(unv.os.tp53,df.os.unv.tp53.braf)



############
### (5.3) UNIVARIATE-PFS --- RAS (baseline) w/ other 5 genes:

### RAS & TP53 (Univariate):
pfs.unv.ras.tp53 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.unv.ras.tp53$Groups <- ifelse(pfs.unv.ras.tp53$TP53=="Mut","RAS+/TP53+","RAS+")

cox.pfs.unv.ras.tp53 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                              data = pfs.unv.ras.tp53)
df.pfs.unv.ras.tp53 <- cox_as_data_frame(
  cox.pfs.unv.ras.tp53,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.ras.tp53$Analysis <- "PFS_UNV_RAS_TP53"


### RAS & SMAD4 (Univariate):
pfs.unv.ras.smad4 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.unv.ras.smad4$Groups <- ifelse(pfs.unv.ras.smad4$SMAD4=="Mut","RAS+/SMAD4+","RAS+")

cox.pfs.unv.ras.smad4 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                               data = pfs.unv.ras.smad4)
df.pfs.unv.ras.smad4 <- cox_as_data_frame(
  cox.pfs.unv.ras.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.ras.smad4$Analysis <- "PFS_UNV_RAS_SMAD4"


### RAS & FBXW7 (Univariate):
pfs.unv.ras.fb7 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.unv.ras.fb7$Groups <- ifelse(pfs.unv.ras.fb7$FBXW7=="Mut","RAS+/FBXW7+","RAS+")

cox.pfs.unv.ras.fb7 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                             data = pfs.unv.ras.fb7)
df.pfs.unv.ras.fb7 <- cox_as_data_frame(
  cox.pfs.unv.ras.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.ras.fb7$Analysis <- "PFS_UNV_RAS_FBXW7"


### RAS & PIK3CA (Univariate):
pfs.unv.ras.pik3 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.unv.ras.pik3$Groups <- ifelse(pfs.unv.ras.pik3$PIK3CA=="Mut","RAS+/PIK3CA+","RAS+")

cox.pfs.unv.ras.pik3 <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                              data = pfs.unv.ras.pik3)
df.pfs.unv.ras.pik3 <- cox_as_data_frame(
  cox.pfs.unv.ras.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.ras.pik3$Analysis <- "PFS_UNV_RAS_PIK3CA"


### RAS & BRAF (Univariate):
pfs.unv.ras.braf <- dat.mv[dat.mv$RAS=="Mut",]
pfs.unv.ras.braf$Groups <- ifelse(pfs.unv.ras.braf$BRAF=="Mut","RAS+/BRAF+","RAS+")

cox.pfs.unv.ras.braf <- coxph(Surv(PFS_days,Progression) ~ Groups, 
                              data = pfs.unv.ras.braf)
df.pfs.unv.ras.braf <- cox_as_data_frame(
  cox.pfs.unv.ras.braf,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.unv.ras.braf$Analysis <- "PFS_UNV_RAS_BRAF"


### JOINT dataframe "unv.pfs.ras":
unv.pfs.ras <- rbind(df.pfs.unv.ras.tp53,df.pfs.unv.ras.smad4)
unv.pfs.ras <- rbind(unv.pfs.ras,df.pfs.unv.ras.fb7)
unv.pfs.ras <- rbind(unv.pfs.ras,df.pfs.unv.ras.pik3)
unv.pfs.ras <- rbind(unv.pfs.ras,df.pfs.unv.ras.braf)



############
### (5.4) UNIVARIATE-OS --- RAS (baseline) w/ other 5 genes:

### RAS & TP53 (Univariate):
os.unv.ras.tp53 <- dat.mv[dat.mv$RAS=="Mut",]
os.unv.ras.tp53$Groups <- ifelse(os.unv.ras.tp53$TP53=="Mut","RAS+/TP53+","RAS+")

cox.os.unv.ras.tp53 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                             data = os.unv.ras.tp53)
df.os.unv.ras.tp53 <- cox_as_data_frame(
  cox.os.unv.ras.tp53,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.ras.tp53$Analysis <- "OS_UNV_RAS_TP53"


### RAS & SMAD4 (Univariate):
os.unv.ras.smad4 <- dat.mv[dat.mv$RAS=="Mut",]
os.unv.ras.smad4$Groups <- ifelse(os.unv.ras.smad4$SMAD4=="Mut","RAS+/SMAD4+","RAS+")

cox.os.unv.ras.smad4 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                              data = os.unv.ras.smad4)
df.os.unv.ras.smad4 <- cox_as_data_frame(
  cox.os.unv.ras.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.ras.smad4$Analysis <- "OS_UNV_RAS_SMAD4"


### RAS & FBXW7 (Univariate):
os.unv.ras.fb7 <- dat.mv[dat.mv$RAS=="Mut",]
os.unv.ras.fb7$Groups <- ifelse(os.unv.ras.fb7$FBXW7=="Mut","RAS+/FBXW7+","RAS+")

cox.os.unv.ras.fb7 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                            data = os.unv.ras.fb7)
df.os.unv.ras.fb7 <- cox_as_data_frame(
  cox.os.unv.ras.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.ras.fb7$Analysis <- "OS_UNV_RAS_FBXW7"


### RAS & PIK3CA (Univariate):
os.unv.ras.pik3 <- dat.mv[dat.mv$RAS=="Mut",]
os.unv.ras.pik3$Groups <- ifelse(os.unv.ras.pik3$PIK3CA=="Mut","RAS+/PIK3CA+","RAS+")

cox.os.unv.ras.pik3 <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                             data = os.unv.ras.pik3)
df.os.unv.ras.pik3 <- cox_as_data_frame(
  cox.os.unv.ras.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.ras.pik3$Analysis <- "OS_UNV_RAS_PIK3CA"


### RAS & BRAF (Univariate) ---> Only 1 patient having RAS+/BRAF+!
os.unv.ras.braf <- dat.mv[dat.mv$RAS=="Mut",]
os.unv.ras.braf$Groups <- ifelse(os.unv.ras.braf$BRAF=="Mut","RAS+/BRAF+","RAS+")

cox.os.unv.ras.braf <- coxph(Surv(OS_days,Exitus) ~ Groups, 
                             data = os.unv.ras.braf)
df.os.unv.ras.braf <- cox_as_data_frame(
  cox.os.unv.ras.braf,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.unv.ras.braf$Analysis <- "OS_UNV_RAS_BRAF"


### JOINT dataframe "unv.os.ras":
unv.os.ras <- rbind(df.os.unv.ras.tp53,df.os.unv.ras.smad4)
unv.os.ras <- rbind(unv.os.ras,df.os.unv.ras.fb7)
unv.os.ras <- rbind(unv.os.ras,df.os.unv.ras.pik3)
unv.os.ras <- rbind(unv.os.ras,df.os.unv.ras.braf)



############
### (6.1) MULTIVARIATE-PFS --- TP53 (baseline) w/ other 5 genes:

### TP53 & RAS (Multivariate):
pfs.mv.tp53.ras <- dat.mv[dat.mv$TP53=="Mut",]
pfs.mv.tp53.ras$Groups <- ifelse(pfs.mv.tp53.ras$RAS=="Mut","TP53+/RAS+","TP53+")

cox.pfs.mv.tp53.ras <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                              data = pfs.mv.tp53.ras)
df.pfs.mv.tp53.ras <- cox_as_data_frame(
  cox.pfs.mv.tp53.ras,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.tp53.ras$Analysis <- "PFS_MV-1_TP53_RAS"

#
surv.sub.unv <- survfit(Surv(PFS_days,Progression) ~ Groups, data = pfs.unv.tp53.ras)
ggsurvplot(surv.sub.unv, pval = TRUE,
           #risk.table = TRUE,   # risk.table = TRUE, #risk.table="percentage",
           conf.int=FALSE,   #risk.table.row=c(0,200,400,600,800,1000,2000,3000,4000,5000),
           palette="lancet", title= "Progression-Free Survival (PFS) Probability",
           break.time.by=500,xlab="Time elapsed (days)",ylab="Accumulated survival",
           legend="right",legend.labs=c("TP53+/RAS-wt","TP53+/RAS+"))#,
#pval.method=TRUE)



### TP53 & SMAD4 (Multivariate):
pfs.mv.tp53.smad4 <- dat.mv[dat.mv$TP53=="Mut",]
pfs.mv.tp53.smad4$Groups <- ifelse(pfs.mv.tp53.smad4$SMAD4=="Mut","TP53+/SMAD4+","TP53+")

cox.pfs.mv.tp53.smad4 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, , 
                                data = pfs.mv.tp53.smad4)
df.pfs.mv.tp53.smad4 <- cox_as_data_frame(
  cox.pfs.mv.tp53.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.tp53.smad4$Analysis <- "PFS_MV-1_TP53_SMAD4"


### TP53 & FBXW7 (Multivariate):
pfs.mv.tp53.fb7 <- dat.mv[dat.mv$TP53=="Mut",]
pfs.mv.tp53.fb7$Groups <- ifelse(pfs.mv.tp53.fb7$FBXW7=="Mut","TP53+/FBXW7+","TP53+")

cox.pfs.mv.tp53.fb7 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                              data = pfs.mv.tp53.fb7)
df.pfs.mv.tp53.fb7 <- cox_as_data_frame(
  cox.pfs.mv.tp53.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.tp53.fb7$Analysis <- "PFS_MV-1_TP53_FBXW7"


### TP53 & PIK3CA (Multivariate):
pfs.mv.tp53.pik3 <- dat.mv[dat.mv$TP53=="Mut",]
pfs.mv.tp53.pik3$Groups <- ifelse(pfs.mv.tp53.pik3$PIK3CA=="Mut","TP53+/PIK3CA+","TP53+")

cox.pfs.mv.tp53.pik3 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                               data = pfs.mv.tp53.pik3)
df.pfs.mv.tp53.pik3 <- cox_as_data_frame(
  cox.pfs.mv.tp53.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.tp53.pik3$Analysis <- "PFS_MV-1_TP53_PIK3CA"


### TP53 & BRAF (Multivariate):
pfs.mv.tp53.braf <- dat.mv[dat.mv$TP53=="Mut",]
pfs.mv.tp53.braf$Groups <- ifelse(pfs.mv.tp53.braf$BRAF=="Mut","TP53+/BRAF+","TP53+")

cox.pfs.mv.tp53.braf <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                               data = pfs.mv.tp53.braf)
df.pfs.mv.tp53.braf <- cox_as_data_frame(
  cox.pfs.mv.tp53.braf,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.tp53.braf$Analysis <- "PFS_MV-1_TP53_BRAF"


### JOINT dataframe "mv.pfs.tp53":
mv.pfs.tp53 <- rbind(df.pfs.mv.tp53.ras[df.pfs.mv.tp53.ras$factor.name=="Groups",],
                     df.pfs.mv.tp53.smad4[df.pfs.mv.tp53.smad4$factor.name=="Groups",])
mv.pfs.tp53 <- rbind(mv.pfs.tp53,df.pfs.mv.tp53.fb7[df.pfs.mv.tp53.fb7$factor.name=="Groups",])
mv.pfs.tp53 <- rbind(mv.pfs.tp53,df.pfs.mv.tp53.pik3[df.pfs.mv.tp53.pik3$factor.name=="Groups",])
mv.pfs.tp53 <- rbind(mv.pfs.tp53,df.pfs.mv.tp53.braf[df.pfs.mv.tp53.braf$factor.name=="Groups",])



############
### (6.2) MULTIVARIATE-OS --- TP53 (baseline) w/ other 5 genes:

### TP53 & RAS (Multivariate):
os.mv.tp53.ras <- dat.mv[dat.mv$TP53=="Mut",]
os.mv.tp53.ras$Groups <- ifelse(os.mv.tp53.ras$RAS=="Mut","TP53+/RAS+","TP53+")

cox.os.mv.tp53.ras <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                             data = os.mv.tp53.ras)
df.os.mv.tp53.ras <- cox_as_data_frame(
  cox.os.mv.tp53.ras,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.tp53.ras$Analysis <- "OS_MV-1_TP53_RAS"


### TP53 & SMAD4 (Multivariate):
os.mv.tp53.smad4 <- dat.mv[dat.mv$TP53=="Mut",]
os.mv.tp53.smad4$Groups <- ifelse(os.mv.tp53.smad4$SMAD4=="Mut","TP53+/SMAD4+","TP53+")

cox.os.mv.tp53.smad4 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                               data = os.mv.tp53.smad4)
df.os.mv.tp53.smad4 <- cox_as_data_frame(
  cox.os.mv.tp53.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.tp53.smad4$Analysis <- "OS_MV-1_TP53_SMAD4"


### TP53 & FBXW7 (Multivariate):
os.mv.tp53.fb7 <- dat.mv[dat.mv$TP53=="Mut",]
os.mv.tp53.fb7$Groups <- ifelse(os.mv.tp53.fb7$FBXW7=="Mut","TP53+/FBXW7+","TP53+")

cox.os.mv.tp53.fb7 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                             data = os.mv.tp53.fb7)
df.os.mv.tp53.fb7 <- cox_as_data_frame(
  cox.os.mv.tp53.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.tp53.fb7$Analysis <- "OS_MV-1_TP53_FBXW7"


### TP53 & PIK3CA (Multivariate):
os.mv.tp53.pik3 <- dat.mv[dat.mv$TP53=="Mut",]
os.mv.tp53.pik3$Groups <- ifelse(os.mv.tp53.pik3$PIK3CA=="Mut","TP53+/PIK3CA+","TP53+")

cox.os.mv.tp53.pik3 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                              data = os.mv.tp53.pik3)
df.os.mv.tp53.pik3 <- cox_as_data_frame(
  cox.os.mv.tp53.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.tp53.pik3$Analysis <- "OS_MV-1_TP53_PIK3CA"


### TP53 & BRAF (Multivariate):
os.mv.tp53.braf <- dat.mv[dat.mv$TP53=="Mut",]
os.mv.tp53.braf$Groups <- ifelse(os.mv.tp53.braf$BRAF=="Mut","TP53+/BRAF+","TP53+")

cox.os.mv.tp53.braf <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
                              data = os.mv.tp53.braf)
df.os.mv.tp53.braf <- cox_as_data_frame(
  cox.os.mv.tp53.braf,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.tp53.braf$Analysis <- "OS_MV-1_TP53_BRAF"


### JOINT dataframe "mv.os.tp53":
mv.os.tp53 <- rbind(df.os.mv.tp53.ras[df.os.mv.tp53.ras$factor.name=="Groups",],
                    df.os.mv.tp53.smad4[df.os.mv.tp53.smad4$factor.name=="Groups",])
mv.os.tp53 <- rbind(mv.os.tp53,df.os.mv.tp53.fb7[df.os.mv.tp53.fb7$factor.name=="Groups",])
mv.os.tp53 <- rbind(mv.os.tp53,df.os.mv.tp53.pik3[df.os.mv.tp53.pik3$factor.name=="Groups",])
mv.os.tp53 <- rbind(mv.os.tp53,df.os.mv.tp53.braf[df.os.mv.tp53.braf$factor.name=="Groups",])



############
### (6.3) MULTIVARIATE-PFS --- RAS (baseline) w/ other 5 genes:

### RAS & TP53 (Multivariate):
pfs.mv.ras.tp53 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.mv.ras.tp53$Groups <- ifelse(pfs.mv.ras.tp53$TP53=="Mut","RAS+/TP53+","RAS+")

cox.pfs.mv.ras.tp53 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                              data = pfs.mv.ras.tp53)
df.pfs.mv.ras.tp53 <- cox_as_data_frame(
  cox.pfs.mv.ras.tp53,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.ras.tp53$Analysis <- "PFS_MV-1_RAS_TP53"


### RAS & SMAD4 (Multivariate):
pfs.mv.ras.smad4 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.mv.ras.smad4$Groups <- ifelse(pfs.mv.ras.smad4$SMAD4=="Mut","RAS+/SMAD4+","RAS+")

cox.pfs.mv.ras.smad4 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                               data = pfs.mv.ras.smad4)
df.pfs.mv.ras.smad4 <- cox_as_data_frame(
  cox.pfs.mv.ras.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.ras.smad4$Analysis <- "PFS_MV-1_RAS_SMAD4"


### RAS & FBXW7 (Multivariate):
pfs.mv.ras.fb7 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.mv.ras.fb7$Groups <- ifelse(pfs.mv.ras.fb7$FBXW7=="Mut","RAS+/FBXW7+","RAS+")

cox.pfs.mv.ras.fb7 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                             data = pfs.mv.ras.fb7)
df.pfs.mv.ras.fb7 <- cox_as_data_frame(
  cox.pfs.mv.ras.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.ras.fb7$Analysis <- "PFS_MV-1_RAS_FBXW7"


### RAS & PIK3CA (Multivariate):
pfs.mv.ras.pik3 <- dat.mv[dat.mv$RAS=="Mut",]
pfs.mv.ras.pik3$Groups <- ifelse(pfs.mv.ras.pik3$PIK3CA=="Mut","RAS+/PIK3CA+","RAS+")

cox.pfs.mv.ras.pik3 <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                              data = pfs.mv.ras.pik3)
df.pfs.mv.ras.pik3 <- cox_as_data_frame(
  cox.pfs.mv.ras.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.pfs.mv.ras.pik3$Analysis <- "PFS_MV-1_RAS_PIK3CA"


### RAS & BRAF (Multivariate) ---> Only 1 case that is RAS+/BRAF+!!!
#pfs.mv.ras.braf <- dat.mv[dat.mv$RAS=="Mut",]
#pfs.mv.ras.braf$Groups <- ifelse(pfs.mv.ras.braf$BRAF=="Mut","RAS+/BRAF+","RAS+")

#cox.pfs.mv.ras.braf <- coxph(Surv(PFS_days,Progression) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
#                              data = pfs.mv.ras.braf)
#df.pfs.mv.ras.braf <- cox_as_data_frame(
#  cox.pfs.mv.ras.braf,
#  unmangle_dict = NULL,
#  factor_id_sep = ": ",
#  sort_by = NULL
#)
#df.pfs.mv.ras.braf$Analysis <- "PFS_MV-1_RAS_BRAF"


### JOINT dataframe "mv.pfs.ras":
mv.pfs.ras <- rbind(df.pfs.mv.ras.tp53[df.pfs.mv.ras.tp53$factor.name=="Groups",],
                    df.pfs.mv.ras.smad4[df.pfs.mv.ras.smad4$factor.name=="Groups",])
mv.pfs.ras <- rbind(mv.pfs.ras,df.pfs.mv.ras.fb7[df.pfs.mv.ras.fb7$factor.name=="Groups",])
mv.pfs.ras <- rbind(mv.pfs.ras,df.pfs.mv.ras.pik3[df.pfs.mv.ras.pik3$factor.name=="Groups",])
#mv.pfs.ras <- rbind(mv.pfs.ras,df.pfs.mv.ras.braf)



############
### (6.4) MULTIVARIATE-OS --- RAS (baseline) w/ other 5 genes:

### RAS & TP53 (Multivariate):
os.mv.ras.tp53 <- dat.mv[dat.mv$RAS=="Mut",]
os.mv.ras.tp53$Groups <- ifelse(os.mv.ras.tp53$TP53=="Mut","RAS+/TP53+","RAS+")

cox.os.mv.ras.tp53 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                             data = os.mv.ras.tp53)
df.os.mv.ras.tp53 <- cox_as_data_frame(
  cox.os.mv.ras.tp53,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.ras.tp53$Analysis <- "OS_MV-1_RAS_TP53"


### RAS & SMAD4 (Multivariate):
os.mv.ras.smad4 <- dat.mv[dat.mv$RAS=="Mut",]
os.mv.ras.smad4$Groups <- ifelse(os.mv.ras.smad4$SMAD4=="Mut","RAS+/SMAD4+","RAS+")

cox.os.mv.ras.smad4 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                              data = os.mv.ras.smad4)
df.os.mv.ras.smad4 <- cox_as_data_frame(
  cox.os.mv.ras.smad4,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.ras.smad4$Analysis <- "OS_MV-1_RAS_SMAD4"


### RAS & FBXW7 (Multivariate):
os.mv.ras.fb7 <- dat.mv[dat.mv$RAS=="Mut",]
os.mv.ras.fb7$Groups <- ifelse(os.mv.ras.fb7$FBXW7=="Mut","RAS+/FBXW7+","RAS+")

cox.os.mv.ras.fb7 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                            data = os.mv.ras.fb7)
df.os.mv.ras.fb7 <- cox_as_data_frame(
  cox.os.mv.ras.fb7,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.ras.fb7$Analysis <- "OS_MV-1_RAS_FBXW7"


### RAS & PIK3CA (Multivariate):
os.mv.ras.pik3 <- dat.mv[dat.mv$RAS=="Mut",]
os.mv.ras.pik3$Groups <- ifelse(os.mv.ras.pik3$PIK3CA=="Mut","RAS+/PIK3CA+","RAS+")

cox.os.mv.ras.pik3 <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+Location_Right1_Left2, 
                             data = os.mv.ras.pik3)
df.os.mv.ras.pik3 <- cox_as_data_frame(
  cox.os.mv.ras.pik3,
  unmangle_dict = NULL,
  factor_id_sep = ": ",
  sort_by = NULL
)
df.os.mv.ras.pik3$Analysis <- "OS_MV-1_RAS_PIK3CA"




### RAS & BRAF (Multivariate) ---> Onlye 1 patient that is RAS+/BRAF+:
#os.mv.ras.braf <- dat.mv[dat.mv$RAS=="Mut",]
#os.mv.ras.braf$Groups <- ifelse(os.mv.ras.braf$BRAF=="Mut","RAS+/BRAF+","RAS+")
#
#cox.os.mv.ras.braf <- coxph(Surv(OS_days,Exitus) ~ Groups+Sex+Age_Dx.M1+Number_Organs_2+ECOG_PS_DxM1+LDH+Alk_Phosph+LEUCOS+PCR+MSS0_MSI1+Location_Right1_Left2, 
#                             data = os.mv.ras.braf)
#df.os.mv.ras.braf <- cox_as_data_frame(
#  cox.os.mv.ras.braf,
#  unmangle_dict = NULL,
#  factor_id_sep = ": ",
#  sort_by = NULL
#)
#df.os.mv.ras.braf$Analysis <- "OS_MV-1_RAS_BRAF"


### JOINT dataframe "mv.os.ras":
mv.os.ras <- rbind(df.os.mv.ras.tp53[df.os.mv.ras.tp53$factor.name=="Groups",],
                   df.os.mv.ras.smad4[df.os.mv.ras.smad4$factor.name=="Groups",])
mv.os.ras <- rbind(mv.os.ras,df.os.mv.ras.fb7[df.os.mv.ras.fb7$factor.name=="Groups",])
mv.os.ras <- rbind(mv.os.ras,df.os.mv.ras.pik3[df.os.mv.ras.pik3$factor.name=="Groups",])
#mv.os.ras <- rbind(mv.os.ras,df.os.mv.ras.braf)



###############################################
###############################################
### FUSION dataframe for all Analyses:

# (a) UNIVARIATE:
unv.all <- rbind(unv.pfs.tp53,unv.pfs.ras)
unv.all <- rbind(unv.all,unv.os.tp53)
unv.all <- rbind(unv.all,unv.os.ras)

dim(unv.all)
#[1] 20 11

unv.all2 <- unv.all

# (b) MULTIVARIATE:
mv.all <- rbind(mv.pfs.tp53,mv.pfs.ras)
mv.all <- rbind(mv.all,mv.os.tp53)
mv.all <- rbind(mv.all,mv.os.ras)

dim(mv.all)
#[1] 20 11

mv.all2 <- mv.all


### (c) ALL - Univariate + Multivariate:

#unv.all2$Factor_value <- "Mut vs No-Mut"
unv.all2 <- unv.all2[,c("factor.name","HR","Lower_CI","Upper_CI","p","Analysis")]

colnames(unv.all2) <- c("Factor","HR","Lower_HR_CI","Upper_HR_CI","p-val","Analysis")
unv.all2 <- unv.all2[,c("Factor","Analysis","p-val","HR","Lower_HR_CI","Upper_HR_CI")]

#
#mv.all2$Factor_value <- "Mut vs No-Mut"
mv.all2 <- mv.all2[,c("factor.name","HR","Lower_CI","Upper_CI","p","Analysis")]

colnames(mv.all2) <- c("Factor","HR","Lower_HR_CI","Upper_HR_CI","p-val","Analysis")
mv.all2 <- mv.all2[,c("Factor","Analysis","p-val","HR","Lower_HR_CI","Upper_HR_CI")]



### FINAL-JOIN:

surv.all <- rbind(unv.all2,mv.all2)
dim(surv.all)
#[1] 38  6

surv.all

#
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Survival_Analyses/Outputs")

#write.table(surv.all,"Survival-mCRC_COMBINATIONS-Genes_Uni-Multi-Variate_15May2020_Lahoz.txt",
#            col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
#WriteXLS(surv.all,"Survival-mCRC_COMBINATIONS-Genes_Uni-Multi-Variate_15May2020_Lahoz.xls",     
#         col.names=TRUE,row.names=FALSE)









