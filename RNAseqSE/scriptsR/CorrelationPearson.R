library("psych")
###############
#CD4 (row 11 ignored)
##########
ignore<-11
#Biometric1: air temperature(3)
ct_CD4_air_temp<-corr.test(x=CD4_air_temp_df[-ignore,], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric2: calories(4)
ct_CD4_calories<-corr.test(x=CD4_calories_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric3: gsr(4)
ct_CD4_gsr<-corr.test(x=CD4_gsr_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric4: heart_rate(3)
ct_CD4_heart_rate<-corr.test(x=CD4_heart_rate_df[-ignore,], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric5: skin_temp(3)
ct_CD4_skin_temp<-corr.test(x=CD4_skin_temp_df[-ignore,], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric6: steps(4)
ct_CD4_steps<-corr.test(x=CD4_steps_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric7: sleep_calories(4)
ct_CD4_sleep_calories<-corr.test(x=CD4_sleep_calories_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric8: sleep_deep_mins(4)
ct_CD4_sleep_deep_mins<-corr.test(x=CD4_sleep_deep_mins_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric9: sleep_duration(4)
ct_CD4_sleep_duration<-corr.test(x=CD4_sleep_duration_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_heart_rate(3)
ct_CD4_sleep_heart_rate<-corr.test(x=CD4_sleep_heart_rate_df[-ignore,], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometricignore: sleep_light_mins(4)
ct_CD4_sleep_light_mins<-corr.test(x=CD4_sleep_light_mins_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_rem_mins(4)
ct_CD4_sleep_rem_mins<-corr.test(x=CD4_sleep_rem_mins_df[-ignore,-5], y=CD4_gene_expresion[-ignore,], method="pearson", adjust="none")


###############
#CD8 (row 10 ignored)
##########
ignore<-10
#Biometric1: air temperature(3)
ct_CD8_air_temp<-corr.test(x=CD8_air_temp_df[-ignore,], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric2: calories(4)
ct_CD8_calories<-corr.test(x=CD8_calories_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric3: gsr(4)
ct_CD8_gsr<-corr.test(x=CD8_gsr_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric8: heart_rate(3)
ct_CD8_heart_rate<-corr.test(x=CD8_heart_rate_df[-ignore,], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric5: skin_temp(3)
ct_CD8_skin_temp<-corr.test(x=CD8_skin_temp_df[-ignore,], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric6: steps(4)
ct_CD8_steps<-corr.test(x=CD8_steps_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric7: sleep_calories(4)
ct_CD8_sleep_calories<-corr.test(x=CD8_sleep_calories_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric8: sleep_deep_mins(4)
ct_CD8_sleep_deep_mins<-corr.test(x=CD8_sleep_deep_mins_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric9: sleep_duration(4)
ct_CD8_sleep_duration<-corr.test(x=CD8_sleep_duration_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_heart_rate(3)
ct_CD8_sleep_heart_rate<-corr.test(x=CD8_sleep_heart_rate_df[-ignore,], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometricignore: sleep_light_mins(4)
ct_CD8_sleep_light_mins<-corr.test(x=CD8_sleep_light_mins_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_rem_mins(4)
ct_CD8_sleep_rem_mins<-corr.test(x=CD8_sleep_rem_mins_df[-ignore,-5], y=CD8_gene_expresion[-ignore,], method="pearson", adjust="none")


###############
#CD14 (row 8 ignored)
##########
ignore<-8
#Biometric1: air temperature(3)
ct_CD14_air_temp<-corr.test(x=CD14_air_temp_df[-ignore,], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric2: calories(4)
ct_CD14_calories<-corr.test(x=CD14_calories_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric3: gsr(4)
ct_CD14_gsr<-corr.test(x=CD14_gsr_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric14: heart_rate(3)
ct_CD14_heart_rate<-corr.test(x=CD14_heart_rate_df[-ignore,], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric5: skin_temp(3)
ct_CD14_skin_temp<-corr.test(x=CD14_skin_temp_df[-ignore,], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric6: steps(4)
ct_CD14_steps<-corr.test(x=CD14_steps_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric7: sleep_calories(4)
ct_CD14_sleep_calories<-corr.test(x=CD14_sleep_calories_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric14: sleep_deep_mins(4)
ct_CD14_sleep_deep_mins<-corr.test(x=CD14_sleep_deep_mins_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric9: sleep_duration(4)
ct_CD14_sleep_duration<-corr.test(x=CD14_sleep_duration_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_heart_rate(3)
ct_CD14_sleep_heart_rate<-corr.test(x=CD14_sleep_heart_rate_df[-ignore,], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometricignore: sleep_light_mins(4)
ct_CD14_sleep_light_mins<-corr.test(x=CD14_sleep_light_mins_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_rem_mins(4)
ct_CD14_sleep_rem_mins<-corr.test(x=CD14_sleep_rem_mins_df[-ignore,-5], y=CD14_gene_expresion[-ignore,], method="pearson", adjust="none")



###############
#CD19 (row 7 ignored)
##########
ignore<-7
#Biometric1: air temperature(3)
ct_CD19_air_temp<-corr.test(x=CD19_air_temp_df[-ignore,], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric2: calories(4)
ct_CD19_calories<-corr.test(x=CD19_calories_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric3: gsr(4)
ct_CD19_gsr<-corr.test(x=CD19_gsr_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric19: heart_rate(3)
ct_CD19_heart_rate<-corr.test(x=CD19_heart_rate_df[-ignore,], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric5: skin_temp(3)
ct_CD19_skin_temp<-corr.test(x=CD19_skin_temp_df[-ignore,], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric6: steps(4)
ct_CD19_steps<-corr.test(x=CD19_steps_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric7: sleep_calories(4)
ct_CD19_sleep_calories<-corr.test(x=CD19_sleep_calories_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric19: sleep_deep_mins(4)
ct_CD19_sleep_deep_mins<-corr.test(x=CD19_sleep_deep_mins_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric9: sleep_duration(4)
ct_CD19_sleep_duration<-corr.test(x=CD19_sleep_duration_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_heart_rate(3)
ct_CD19_sleep_heart_rate<-corr.test(x=CD19_sleep_heart_rate_df[-ignore,], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometricignore: sleep_light_mins(4)
ct_CD19_sleep_light_mins<-corr.test(x=CD19_sleep_light_mins_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_rem_mins(4)
ct_CD19_sleep_rem_mins<-corr.test(x=CD19_sleep_rem_mins_df[-ignore,-5], y=CD19_gene_expresion[-ignore,], method="pearson", adjust="none")


###############
#PBMC (row 13 ignored)
##########
ignore<-13
#Biometric1: air temperature(3)
ct_PBMC_air_temp<-corr.test(x=PBMC_air_temp_df[-ignore,], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric2: calories(4)
ct_PBMC_calories<-corr.test(x=PBMC_calories_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric3: gsr(4)
ct_PBMC_gsr<-corr.test(x=PBMC_gsr_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric4: heart_rate(3)
ct_PBMC_heart_rate<-corr.test(x=PBMC_heart_rate_df[-ignore,], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric5: skin_temp(3)
ct_PBMC_skin_temp<-corr.test(x=PBMC_skin_temp_df[-ignore,], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric6: steps(4)
ct_PBMC_steps<-corr.test(x=PBMC_steps_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric7: sleep_calories(4)
ct_PBMC_sleep_calories<-corr.test(x=PBMC_sleep_calories_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric4: sleep_deep_mins(4)
ct_PBMC_sleep_deep_mins<-corr.test(x=PBMC_sleep_deep_mins_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric9: sleep_duration(4)
ct_PBMC_sleep_duration<-corr.test(x=PBMC_sleep_duration_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_heart_rate(3)
ct_PBMC_sleep_heart_rate<-corr.test(x=PBMC_sleep_heart_rate_df[-ignore,], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometricignore: sleep_light_mins(4)
ct_PBMC_sleep_light_mins<-corr.test(x=PBMC_sleep_light_mins_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")
#Biometric10: sleep_rem_mins(4)
ct_PBMC_sleep_rem_mins<-corr.test(x=PBMC_sleep_rem_mins_df[-ignore,-5], y=PBMC_gene_expresion[-ignore,], method="pearson", adjust="none")

####PLOTS(HEATMAPS)####
#############
library(gplots)
####CD4
CD4<-cbind(t(t(ct_CD4_air_temp$r[1,])), t(t(ct_CD4_calories$r[1,])), t(t(ct_CD4_gsr$r[1,])), t(t(ct_CD4_heart_rate$r[1,])),
           t(t(ct_CD4_skin_temp$r[1,])), t(t(ct_CD4_steps$r[1,])), t(t(ct_CD4_sleep_calories$r[1,])), t(t(ct_CD4_sleep_deep_mins$r[1,])),
           t(t(ct_CD4_sleep_duration$r[1,])), t(t(ct_CD4_sleep_heart_rate$r[1,])), t(t(ct_CD4_sleep_light_mins$r[1,])), t(t(ct_CD4_sleep_rem_mins$r[1,])))
colnames(CD4)<-c("air_temp", "calories", "gsr", "heart_rate", "skin", "steps", "sleep_calories", "sleep_deep_mins", "sleep_duration", "sleep_heart_rate", "sleep_light_mins", "sleep_rem_mins") 

heatmap.2(CD4, col=bluered(245), key=TRUE, key.xlab ="", scale="none", trace="none",
          xlab= "Factors", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="CD4 pearson")

####CD8
CD8<-cbind(t(t(ct_CD8_air_temp$r[1,])), t(t(ct_CD8_calories$r[1,])), t(t(ct_CD8_gsr$r[1,])), t(t(ct_CD8_heart_rate$r[1,])),
           t(t(ct_CD8_skin_temp$r[1,])), t(t(ct_CD8_steps$r[1,])), t(t(ct_CD8_sleep_calories$r[1,])), t(t(ct_CD8_sleep_deep_mins$r[1,])),
           t(t(ct_CD8_sleep_duration$r[1,])), t(t(ct_CD8_sleep_heart_rate$r[1,])), t(t(ct_CD8_sleep_light_mins$r[1,])), t(t(ct_CD8_sleep_rem_mins$r[1,])))
colnames(CD8)<-c("air_temp", "calories", "gsr", "heart_rate", "skin", "steps", "sleep_calories", "sleep_deep_mins", "sleep_duration", "sleep_heart_rate", "sleep_light_mins", "sleep_rem_mins") 

heatmap.2(CD8, col=bluered(245), key=TRUE, key.xlab ="", scale="none", trace="none",
          xlab= "Factors", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="CD8 pearson")


####CD14
CD14<-cbind(t(t(ct_CD14_air_temp$r[1,])), t(t(ct_CD14_calories$r[1,])), t(t(ct_CD14_gsr$r[1,])), t(t(ct_CD14_heart_rate$r[1,])),
            t(t(ct_CD14_skin_temp$r[1,])), t(t(ct_CD14_steps$r[1,])), t(t(ct_CD14_sleep_calories$r[1,])), t(t(ct_CD14_sleep_deep_mins$r[1,])),
            t(t(ct_CD14_sleep_duration$r[1,])), t(t(ct_CD14_sleep_heart_rate$r[1,])), t(t(ct_CD14_sleep_light_mins$r[1,])), t(t(ct_CD14_sleep_rem_mins$r[1,])))
colnames(CD14)<-c("air_temp", "calories", "gsr", "heart_rate", "skin", "steps", "sleep_calories", "sleep_deep_mins", "sleep_duration", "sleep_heart_rate", "sleep_light_mins", "sleep_rem_mins") 

heatmap.2(CD14, col=bluered(245), key=TRUE, key.xlab ="", scale="none", trace="none",
          xlab= "Factors", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="CD14 pearson")




#CD19
CD19<-cbind(t(t(ct_CD19_air_temp$r[1,])), t(t(ct_CD19_calories$r[1,])), t(t(ct_CD19_gsr$r[1,])), t(t(ct_CD19_heart_rate$r[1,])),
            t(t(ct_CD19_skin_temp$r[1,])), t(t(ct_CD19_steps$r[1,])), t(t(ct_CD19_sleep_calories$r[1,])), t(t(ct_CD19_sleep_deep_mins$r[1,])),
            t(t(ct_CD19_sleep_duration$r[1,])), t(t(ct_CD19_sleep_heart_rate$r[1,])), t(t(ct_CD19_sleep_light_mins$r[1,])), t(t(ct_CD19_sleep_rem_mins$r[1,])))
colnames(CD19)<-c("air_temp", "calories", "gsr", "heart_rate", "skin", "steps", "sleep_calories", "sleep_deep_mins", "sleep_duration", "sleep_heart_rate", "sleep_light_mins", "sleep_rem_mins") 

heatmap.2(CD19, col=bluered(245), key=TRUE, key.xlab ="", scale="none", trace="none",
          xlab= "Factors", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="CD19 pearson")


#PBMC
PBMC<-cbind(t(t(ct_PBMC_air_temp$r[1,])), t(t(ct_PBMC_calories$r[1,])), t(t(ct_PBMC_gsr$r[1,])), t(t(ct_PBMC_heart_rate$r[1,])),
            t(t(ct_PBMC_skin_temp$r[1,])), t(t(ct_PBMC_steps$r[1,])), t(t(ct_PBMC_sleep_calories$r[1,])), t(t(ct_PBMC_sleep_deep_mins$r[1,])),
            t(t(ct_PBMC_sleep_duration$r[1,])), t(t(ct_PBMC_sleep_heart_rate$r[1,])), t(t(ct_PBMC_sleep_light_mins$r[1,])), t(t(ct_PBMC_sleep_rem_mins$r[1,])))
colnames(PBMC)<-c("air_temp", "calories", "gsr", "heart_rate", "skin", "steps", "sleep_calories", "sleep_deep_mins", "sleep_duration", "sleep_heart_rate", "sleep_light_mins", "sleep_rem_mins") 

heatmap.2(PBMC, col=bluered(245), key=TRUE, key.xlab ="", scale="none", trace="none",
          xlab= "Factors", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="PBMC pearson")

