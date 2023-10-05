## Geary's test (mean* sqrt(n) )/ SD

warmdata <- read.csv("Survival project effect size reftemp FF warm.csv")


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


gdata <- c()
for(i in 1:nrow(warmdata)){
  row <- calculate_g(column_sets[[1]]$Temp1[i], column_sets[[1]]$Mean1[i], column_sets[[1]]$SD1[i], column_sets[[1]]$N1[i])
  row <- cbind(row, calculate_g(column_sets[[2]]$Temp2[i], column_sets[[2]]$Mean2[i], column_sets[[2]]$SD2[i], column_sets[[2]]$N2[i]))
  row <- cbind(row, calculate_g(column_sets[[3]]$Temp3[i], column_sets[[3]]$Mean3[i], column_sets[[3]]$SD3[i], column_sets[[3]]$N3[i]))
  row <- cbind(row, calculate_g(column_sets[[4]]$Temp4[i], column_sets[[4]]$Mean4[i], column_sets[[4]]$SD4[i], column_sets[[4]]$N4[i]))
  row <- cbind(row, calculate_g(column_sets[[5]]$Temp5[i], column_sets[[5]]$Mean5[i], column_sets[[5]]$SD5[i], column_sets[[5]]$N5[i]))
  row <- cbind(row, calculate_g(column_sets[[6]]$Temp6[i], column_sets[[6]]$Mean6[i], column_sets[[6]]$SD6[i], column_sets[[6]]$N6[i]))
  row <- cbind(row, calculate_g(column_sets[[7]]$Temp7[i], column_sets[[7]]$Mean7[i], column_sets[[7]]$SD7[i], column_sets[[7]]$N7[i]))
  row <- cbind(row, calculate_g(column_sets[[8]]$Temp8[i], column_sets[[8]]$Mean8[i], column_sets[[8]]$SD8[i], column_sets[[8]]$N8[i]))
  row <- cbind(row, calculate_g(column_sets[[9]]$Temp9[i], column_sets[[9]]$Mean9[i], column_sets[[9]]$SD9[i], column_sets[[9]]$N9[i]))

  gdata <- rbind(gdata, row)
  
}
### creating error when N=1 because of SMD formula. dont use for now. 

gdata <- as.data.frame(gdata)

new_column <- apply(gdata, 1, function(row) any(!is.na(row) & row < 3))
gearys_test <- cbind(warmdata$Paper.code, gdata, new_column)



# Use apply to count how many values in total are < 3 i.e. how many effect sizes we'd end up removing
result <- apply(gdata, 2, function(col) sum(!is.na(col) & col < 3))

# Sum the results to get the total count = 256 .... worth removing. 
total_count <- sum(result)




#### Indentifying which effect sizes would be removed by importing treattemp

 newdata <- cbind(warmdata$Paper.code, 
                  column_sets[[1]][1], gearys_test$row,
                  column_sets[[2]][1], gearys_test$V2,
                  column_sets[[3]][1], gearys_test$V3,
                  column_sets[[4]][1], gearys_test$V4,
                  column_sets[[5]][1], gearys_test$V5,
                  column_sets[[6]][1], gearys_test$V6,
                  column_sets[[7]][1], gearys_test$V7,
                  column_sets[[8]][1], gearys_test$V8,
                  column_sets[[9]][1], gearys_test$V9,
                  gearys_test$new_column)

 colnames(newdata) <- c("Paper.code", 
                        "Temp1", "Gtest1",
                        "Temp2", "Gtest2",
                        "Temp3", "Gtest3",
                        "Temp4", "Gtest4",
                        "Temp5", "Gtest5",
                        "Temp6", "Gtest6",
                        "Temp7", "Gtest7",
                        "Temp8", "Gtest8",
                        "Temp9", "Gtest9",
                        "Outlier")
 
 outlier_temps <- NA
 
 #### need to amend code as sometimes there are more than one treatments per row which have a gtest < 3
 for(i in 1:nrow(gdata)){
   if(gearys_test$new_column[i]){
 ot <- which(gdata[i,] < 3)
 outlier_temps[i] <- ot  
   }
  }
 
 outlier_temps <- as.data.frame(outlier_temps)
 