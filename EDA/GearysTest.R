## Geary's test (mean* sqrt(n) )/ SD

warmdata <- read.csv("EDA/Survival project effect size reftemp FF warm.csv")


# Create a list of column names for each treatment
column_sets <- list(warmdata[c("Temp1", "Mean1", "SD1", "N1", "ani_1")],
                    warmdata[c("Temp2", "Mean2", "SD2", "N2", "ani_2")],
                    warmdata[c("Temp3", "Mean3", "SD3", "N3", "ani_3")],
                    warmdata[c("Temp4", "Mean4", "SD4", "N4", "ani_4")],
                    warmdata[c("Temp5", "Mean5", "SD5", "N5", "ani_5")],
                    warmdata[c("Temp6", "Mean6", "SD6", "N6", "ani_6")],
                    warmdata[c("Temp7", "Mean7", "SD7", "N7", "ani_7")],
                    warmdata[c("Temp8", "Mean8", "SD8", "N8", "ani_8")],
                    warmdata[c("Temp9", "Mean9", "SD9", "N9", "ani_9")])


#Gearys <- NA

calculate_g <- function(temp, mean, sd, n){
  if (is.na(temp)){
    g <- NA  # Initialize as NA (numeric)
  } else {
    g <- (mean*sqrt(n))/sd
  }
  
  return(g)
}


 ################ I  can combine this into the paiwise calc so I can easily remove specific es.
 ref_temp <- warmdata$ref_temp
 ref_mean <- warmdata$ref_mean
 ref_sd <- warmdata$ref_sd
 ref_ani <- warmdata$ref_ani
 ref_N <- warmdata$ref_N
 
 
 
 calculate_smd <- function(rtemp1, rmean1, rsd1, rani1, rn1, temp2, mean2, sd2, ani2, n2){
   if (is.na(temp2)) {
     smd <- NA  # Initialize as NA (numeric)
     smd_v <- NA
     gtest <- NA
   } else if (rtemp1 != temp2) {
     smd <- escalc(measure = "SMD", m2i=rmean1, sd2i=rsd1, n2i=rn1, m1i=mean2, sd1i=sd2, n1i=n2)[1]
     smd_v <- escalc(measure = "SMD", m2i=rmean1, sd2i=rsd1, n2i=rani1, m1i=mean2, sd1i=sd2, n1i=ani2)[2]
   } else {
     smd <- 0  
     smd_v <- 0 # rtemp1 equals temp2, so smd is set to 0
   }
   Smd <- cbind(smd, smd_v, temp2-rtemp1, rtemp1, temp2, calculate_g(temp2, mean2, sd2, n2))
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
   
   colnames(row) <- c("es.1", "v.1", "??", "reftemp", "treattemp", "gtest",
                      "es.2", "v.2", "??", "reftemp", "treattemp", "gtest", 
                      "es.3", "v.3", "??", "reftemp", "treattemp", "gtest",
                      "es.4", "v.4", "??", "reftemp", "treattemp", "gtest",
                      "es.5", "v.5", "??", "reftemp", "treattemp", "gtest",
                      "es.6", "v.6", "??", "reftemp", "treattemp", "gtest",
                      "es.7", "v.7", "??", "reftemp", "treattemp", "gtest",
                      "es.8", "v.8", "??", "reftemp", "treattemp", "gtest",
                      "es.9", "v.9", "??", "reftemp", "treattemp", "gtest")
   esdata <- rbind(esdata, row)
   
 }
 
 
 
 pairwise.test <- esdata
 studydat <- warmdata[,1:34]
 
 colnames(pairwise.test) <- c("es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest",
                              "es", "v", "diff", "reftemp", "treattemp", "gtest")
 
 
 
 
 
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
 
 
 
 
 ##  Count how many values in total are < 3 i.e. how many effect sizes we'd end up removing
table(cleaned_df$gtest < 3) #= 256
 
 