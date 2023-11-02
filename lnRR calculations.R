warmdata <- read.csv("Survival Project Data/Survival project effect size reftemp FF warm.csv")

library(metafor)


# Create a list of column names for each set of columns
column_sets <- list(warmdata[c("Temp1", "Mean1", "SD1", "N1", "ani_1")],
                    warmdata[c("Temp2", "Mean2", "SD2", "N2", "ani_2")],
                    warmdata[c("Temp3", "Mean3", "SD3", "N3", "ani_3")],
                    warmdata[c("Temp4", "Mean4", "SD4", "N4", "ani_4")],
                    warmdata[c("Temp5", "Mean5", "SD5", "N5", "ani_5")],
                    warmdata[c("Temp6", "Mean6", "SD6", "N6", "ani_6")],
                    warmdata[c("Temp7", "Mean7", "SD7", "N7", "ani_7")],
                    warmdata[c("Temp8", "Mean8", "SD8", "N8", "ani_8")],
                    warmdata[c("Temp9", "Mean9", "SD9", "N9", "ani_9")])


### make sure escalc uses number of animals not number of 
#replicates when i've added these to spreadhseet
ref_temp <- warmdata$ref_temp
ref_mean <- warmdata$ref_mean
ref_sd <- warmdata$ref_sd
ref_ani <- warmdata$ref_ani
ref_N <- warmdata$ref_N



#calculate_smd_old <- function(rtemp1, rmean1, rsd1, rani1, rn1, temp2, mean2, sd2, ani2, n2){
#  if (is.na(temp2)) {
#    smd <- NA  # Initialize as NA (numeric)
#    smd_v <- NA
#  } else if (rtemp1 > temp2) {
#    smd <- escalc(measure = "SMD", m1i=rmean1, sd1i=rsd1, n1i=rn1, m2i=mean2, sd2i=sd2, n2i=n2)[1]
#    smd_v <- escalc(measure = "SMD", m1i=rmean1, sd1i=rsd1, n1i=rani1, m2i=mean2, sd2i=sd2, n2i=ani2)[2]
#  } else if (rtemp1 < temp2) {
#    smd <- escalc(measure = "SMD", m2i=rmean1, sd2i=rsd1, n2i=rn1, m1i=mean2, sd1i=sd2, n1i=n2)[1]
#    smd_v <- escalc(measure = "SMD", m2i=rmean1, sd2i=rsd1, n2i=rani1, m1i=mean2, sd1i=sd2, n1i=ani2)[2]
#  } else {
#    smd <- 0  
#    smd_v <- 0# rtemp1 equals temp2, so smd is set to 0
#  }
#  Smd <- cbind(smd, smd_v, rtemp1-temp2, rtemp1, temp2)
#  return(Smd)
#}
### calculate pairwise effect sizes for each treatment compared with ref temp


calculate_smd <- function(rtemp1, rmean1, rsd1, rani1, rn1, temp2, mean2, sd2, ani2, n2){
  if (is.na(temp2)) {
    smd <- NA  # Initialize as NA (numeric)
    smd_v <- NA
  } else if (rtemp1 != temp2) {
    smd <- escalc(measure = "ROM", m2i=rmean1, sd2i=rsd1, n2i=rn1, m1i=mean2, sd1i=sd2, n1i=n2)[1]
    smd_v <- escalc(measure = "ROM", m2i=rmean1, sd2i=rsd1, n2i=rani1, m1i=mean2, sd1i=sd2, n1i=ani2)[2]
  } else {
    smd <- 0  
    smd_v <- 0 # rtemp1 equals temp2, so smd is set to 0
  }
  Smd <- cbind(smd, smd_v, temp2-rtemp1, rtemp1, temp2, n2)
  return(Smd)
}


esdata <- c()
for(i in 1:nrow(warmdata)){
  row <- calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[1]]$Temp1[i], column_sets[[1]]$Mean1[i], column_sets[[1]]$SD1[i], column_sets[[1]]$ani_1[i], column_sets[[1]]$N1[i])
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[2]]$Temp2[i], column_sets[[2]]$Mean2[i], column_sets[[2]]$SD2[i], column_sets[[2]]$ani_2[i], column_sets[[2]]$N2[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[3]]$Temp3[i], column_sets[[3]]$Mean3[i], column_sets[[3]]$SD3[i], column_sets[[3]]$ani_3[i], column_sets[[3]]$N3[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[4]]$Temp4[i], column_sets[[4]]$Mean4[i], column_sets[[4]]$SD4[i], column_sets[[4]]$ani_4[i], column_sets[[4]]$N4[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[5]]$Temp5[i], column_sets[[5]]$Mean5[i], column_sets[[5]]$SD5[i], column_sets[[5]]$ani_5[i], column_sets[[5]]$N5[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[6]]$Temp6[i], column_sets[[6]]$Mean6[i], column_sets[[6]]$SD6[i], column_sets[[6]]$ani_6[i], column_sets[[6]]$N6[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[7]]$Temp7[i], column_sets[[7]]$Mean7[i], column_sets[[7]]$SD7[i], column_sets[[7]]$ani_7[i], column_sets[[7]]$N7[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[8]]$Temp8[i], column_sets[[8]]$Mean8[i], column_sets[[8]]$SD8[i], column_sets[[8]]$ani_8[i], column_sets[[8]]$N8[i]))
  row <- cbind(row, calculate_smd(ref_temp[i], ref_mean[i], ref_sd[i], ref_ani[i], ref_N[i], column_sets[[9]]$Temp9[i], column_sets[[9]]$Mean9[i], column_sets[[9]]$SD9[i], column_sets[[9]]$ani_9[i], column_sets[[9]]$N9[i]))
  
  colnames(row) <- c("es.1", "v.1", "diff", "reftemp", "treattemp", "n",
                     "es.2", "v.2", "diff", "reftemp", "treattemp", "n",
                     "es.3", "v.3", "diff", "reftemp", "treattemp", "n",
                     "es.4", "v.4", "diff", "reftemp", "treattemp", "n",
                     "es.5", "v.5", "diff", "reftemp", "treattemp", "n",
                     "es.6", "v.6", "diff", "reftemp", "treattemp", "n",
                     "es.7", "v.7", "diff", "reftemp", "treattemp", "n",
                     "es.8", "v.8", "diff", "reftemp", "treattemp", "n",
                     "es.9", "v.9", "diff", "reftemp", "treattemp", "n")
  esdata <- rbind(esdata, row)
  
}
### creating error when N=1 because of SMD formula. dont use for now. 

esdata <- as.data.frame(esdata)


pairwise.test <- esdata
studydat <- warmdata[,1:34]

colnames(pairwise.test) <- c("es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n",
                             "es", "v", "diff", "reftemp", "treattemp", "n")






test.1 <- cbind(studydat, pairwise.test[,1:6])
test.2 <- cbind(studydat, pairwise.test[,7:12])
test.3 <- cbind(studydat, pairwise.test[,13:18])
test.4 <- cbind(studydat, pairwise.test[,19:24])
test.5 <- cbind(studydat, pairwise.test[,25:30])
test.6 <- cbind(studydat, pairwise.test[,31:36])
test.7 <- cbind(studydat, pairwise.test[,37:42])
test.8 <- cbind(studydat, pairwise.test[,43:48])
test.9 <- cbind(studydat, pairwise.test[,49:54])


total.test <- rbind(test.1, test.2, test.3, test.4, test.5, test.6, test.7, test.8, test.9)

cleaned_df <- total.test[!is.na(total.test$es),]

#write.csv(cleaned_df, "pairwise.props.es.csv")

cleaned_df$gtest <- (cleaned_df$es / sqrt(cleaned_df$v)) * sqrt(cleaned_df$n)

#(pairwise.es, "pairwise.es.csv")

