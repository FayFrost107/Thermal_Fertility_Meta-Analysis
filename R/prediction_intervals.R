# Read in data and models

## Longevity 
mv_long<- readRDS(here("output", "models", "meta_longevity_3.rds"))
long_data <- readRDS(here("output", "Output data", "long_data.rds"))

## Reproduction
mv_rep <- readRDS(here("output", "models", "meta_rep_2.rds")) 
rep_data <- readRDS(here("output", "Output data", "rep_data.rds"))

# Make predicitons
preds.long <- predict(mv_long,addx=TRUE)
preds.rep <- predict(mv_rep, addx=TRUE)

# Add predicitons and prediciton intervals into data
long_data$pred <-preds.long$pred  
long_data$pred.lb <- preds.long$pi.lb
long_data$pred.ub <- preds.long$pi.ub

rep_data$pred <- preds.rep$pred  
rep_data$pred.lb <- preds.rep$pi.lb
rep_data$pred.ub <- preds.rep$pi.ub


### Plot curves and prediciton intervals
# Plot for the first dataset (rdata)
plot_ldata <- ggplot(data = long_data, aes(x = c_treattemp, y = pred)) +
  geom_ribbon(aes(ymin = pred.lb, ymax = pred.ub), fill = "lightblue", alpha = 0.5) +
  geom_line() +
  geom_point() +
  theme_bw()

# Add the plot for the second dataset (ldata) to the first plot
final_plot <- plot_ldata +
  geom_line(data = rep_data, aes(x = c_treattemp, y = pred), color = "red") +
  geom_point(data = rep_data, aes(x = c_treattemp, y=pred), col = "red")+
  geom_ribbon(data = rep_data, aes(x = c_treattemp, ymin = pred.lb, ymax = pred.ub), fill = "lightcoral", alpha = 0.5) +
  labs(title = "Confidence Interval Plot",
       x = "X-axis Label",
       y = "Y-axis Label")
