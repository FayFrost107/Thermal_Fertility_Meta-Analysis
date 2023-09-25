library(metafor)

calculate_meta_for_different_CORFACTOR <- function(rdata, CORFACTOR_values) {
  results <- data.frame(CORFACTOR = numeric(0), Estimate = numeric(0), CI_lower = numeric(0), CI_upper = numeric(0),
                        I2_1 = numeric(0), I2_2 = numeric(0), I2_3 = numeric(0), I2_4 = numeric(0), I2_5 = numeric(0), AIC = numeric(0))
  
  for (CORFACTOR in CORFACTOR_values) {
    cat("Calculating for CORFACTOR =", CORFACTOR, "\n")
    
    VCV_shared <- matrix(0, nrow = dim(rdata)[1], ncol = dim(rdata)[1])
    rownames(VCV_shared) <- rdata[, "es"]
    colnames(VCV_shared) <- rdata[, "es"]
    
    shared_coord <- which(rdata[, "shared_control"] %in% rdata[duplicated(rdata[, "shared_control"]), "shared_control"] == TRUE) 
    
    combinations <- do.call("rbind", tapply(shared_coord, rdata[shared_coord,  "shared_control"], function(x) t(utils::combn(x, 2)))) 
    
    for (i in 1:dim(combinations)[1]) {
      p1 <- combinations[i, 1]
      p2 <- combinations[i, 2]
      p1_p2_cov <- CORFACTOR * sqrt(rdata[p1, "v"]) * sqrt(rdata[p2, "v"])
      VCV_shared[p1, p2] <- p1_p2_cov
      VCV_shared[p2, p1] <- p1_p2_cov
    }
    
    diag(VCV_shared) <- rdata[, "v"]
    
    meta_model <- rma.mv(es, VCV_shared, random= list(~ 1|species, ~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
    
    AIC_value <- AIC(meta_model)
    estimate <- coef(summary(meta_model))$estimate[1]
    ci_lower <- meta_model$ci.lb
    ci_upper <- meta_model$ci.ub
    
    i2 <- i2_ml(meta_model)  # Calculate I²
    
    results <- rbind(results, data.frame(CORFACTOR = CORFACTOR, Estimate = estimate, CI_lower = ci_lower, CI_upper = ci_upper,
                                         I2_1 = i2[1], I2_2 = i2[2], I2_3 = i2[3], I2_4 = i2[4], I2_5 = i2[5], AIC = AIC_value))
  }
  return(results)
}

# Example usage
CORFACTOR<- seq(0.1, 0.9, by = 0.1)
results <- calculate_meta_for_different_CORFACTOR(rdata, CORFACTOR)

colnames(results) <- c("CORFACTOR", "Estimate",  "CI lower",  "CI upper", 
                       "I^2 Total", "I^2 species","I^2 studycode", "I^2 shared control", "I^2 obs", "AIC")  

results