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
rep_data$c.lb <- preds.rep$ci.lb
rep_data$c.ub <- preds.rep$ci.ub


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


### Multivariate models
mv_mlma_4 <- readRDS(here("output", "models", "mv_mlma_4.rds"))
preds.mv <- predict(mv_mlma_4, addx=TRUE)


mvs_data <- readRDS(here("output", "Output data", "mv_data.rds"))
mvs_data$pred <- preds.mv$pred  
mvs_data$pred.lb <- preds.mv$pi.lb
mvs_data$pred.ub <- preds.mv$pi.ub
mvs_data$c.lb <- preds.mv$ci.lb
mvs_data$c.ub <- preds.mv$ci.ub

ggplot(data = mvs_data, aes(x = c_treattemp, y = pred, col= outcome)) +
  geom_ribbon(aes(ymin = pred.lb, ymax = pred.ub, fill=outcome), alpha=0.25) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  labs(title = "Prediciton Interval Plot",
       x = "X-axis Label",
       y = "Y-axis Label")

ggplot(data = mvs_data, aes(x = c_treattemp, y = pred, col= outcome)) +
  geom_ribbon(aes(ymin = c.lb, ymax = c.ub, fill=outcome), alpha=0.25) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  labs(title = "Confidence Interval Plot",
       x = "X-axis Label",
       y = "Y-axis Label")


### prediciton intervals for sex exposed (reproduction)
sex_rep_meta <- readRDS(here("output", "models", "meta_treat_rep_sex.rds")) 
sex_rep_data <- readRDS(here("output", "Output data", "data_rep_sex.rds"))

preds.rep.sex <- predict(sex_rep_meta, addx=TRUE)


sex_rep_data$pred <- preds.rep.sex$pred  
sex_rep_data$pred.lb <- preds.rep.sex$pi.lb
sex_rep_data$pred.ub <- preds.rep.sex$pi.ub
sex_rep_data$c.lb <- preds.rep.sex$ci.lb
sex_rep_data$c.ub <- preds.rep.sex$ci.ub

library(ggplot2)

ggplot(data = sex_rep_data, aes(x = c_treattemp, y = pred, col = Sex.exposed, linetype = Sex.exposed)) +
  geom_ribbon(aes(ymin = c.lb, ymax = c.ub, fill = Sex.exposed), alpha = 0.15) +
  geom_line() +
  scale_color_manual(values = c("purple", "orange")) +  # Set colors for the points and lines
  scale_fill_manual(values = c("purple", "orange")) +    # Set colors for the ribbons
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme_bw() + 
  labs(title = "Temperature effect on reproduciton",
       x = "Temperature devaition from 25C",
       y = "Effect size",
       color = "Sex exposed",  # Set the title of the color legend
       fill = "Sex exposed",
       linetype = "Sex exposed")   # Set the title of the fill legend (for ribbons)



### prediciton intervals for life stage (reproduction)
sex_rep_meta <- readRDS(here("output", "models", "meta_treat_rep_ls.rds")) 
sex_rep_data <- readRDS(here("output", "Output data", "data_rep_ls.rds"))

preds.rep.sex <- predict(sex_rep_meta, addx=TRUE)


sex_rep_data$pred <- preds.rep.sex$pred  
sex_rep_data$pred.lb <- preds.rep.sex$pi.lb
sex_rep_data$pred.ub <- preds.rep.sex$pi.ub
sex_rep_data$c.lb <- preds.rep.sex$ci.lb
sex_rep_data$c.ub <- preds.rep.sex$ci.ub

library(ggplot2)

ggplot(data = sex_rep_data, aes(x = c_treattemp, y = pred, col = Life.stage.of.animal, linetype = Life.stage.of.animal,)) +
  geom_ribbon(aes(ymin = c.lb, ymax = c.ub, fill = Life.stage.of.animal,), alpha = 0.15) +
  geom_line() +
  scale_color_manual(values = c("purple", "orange")) +  # Set colors for the points and lines
  scale_fill_manual(values = c("purple", "orange")) +    # Set colors for the ribbons
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme_bw() + 
  labs(title = "Temperature effect on reproduction",
       x = "Temperature devaition from 25C",
       y = "Effect size",
       color = "Life stage",  # Set the title of the color legend
       fill = "Life stage",
       linetype = "Life stage")   # Set the title of the fill legend (for ribbons)



################################ longevity plots ####################

## choose which to plot from below
### prediciton intervals for life stage (longevity)
sex_rep_meta <- readRDS(here("output", "models", "meta_treat_long_ls.rds")) 
sex_rep_data <- readRDS(here("output", "Output data", "data_long_ls.rds"))

### prediciton intervals for sex exposed (longevity)
sex_rep_meta <- readRDS(here("output", "models", "meta_treat_long_sex.rds")) 
sex_rep_data <- readRDS(here("output", "Output data", "data_long_sex.rds"))


preds.rep.sex <- predict(sex_rep_meta, addx=TRUE)


sex_rep_data$pred <- preds.rep.sex$pred  
sex_rep_data$pred.lb <- preds.rep.sex$pi.lb
sex_rep_data$pred.ub <- preds.rep.sex$pi.ub
sex_rep_data$c.lb <- preds.rep.sex$ci.lb
sex_rep_data$c.ub <- preds.rep.sex$ci.ub

library(ggplot2)

ggplot(data = sex_rep_data, aes(x = c_treattemp, y = pred, col = Life.stage.of.animal, linetype = Life.stage.of.animal,)) +
  geom_ribbon(aes(ymin = c.lb, ymax = c.ub, fill = Life.stage.of.animal,), alpha = 0.15) +
  geom_line() +
  scale_color_manual(values = c("purple", "orange")) +  # Set colors for the points and lines
  scale_fill_manual(values = c("purple", "orange")) +    # Set colors for the ribbons
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme_bw() + 
  labs(title = "Temperature effect on Longevity",
       x = "Temperature devaition from 25C",
       y = "Effect size",
       color = "Life stage",  # Set the title of the color legend
       fill = "Life stage",
       linetype = "Life stage")   # Set the title of the fill legend (for ribbons)



ggplot(data = sex_rep_data, aes(x = c_treattemp, y = pred, col = Sex.exposed, linetype = Sex.exposed)) +
  geom_ribbon(aes(ymin = c.lb, ymax = c.ub, fill = Sex.exposed), alpha = 0.15) +
  geom_line() +
  scale_color_manual(values = c("purple", "orange")) +  # Set colors for the points and lines
  scale_fill_manual(values = c("purple", "orange")) +    # Set colors for the ribbons
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme_bw() + 
  labs(title = "Temperature effect on Longevity",
       x = "Temperature devaition from 25C",
       y = "Effect size",
       color = "Sex.exposed",  # Set the title of the color legend
       fill = "Sex.exposed",
       linetype = "Sex.exposed")   # Set the title of the fill legend (for ribbons)


### can also predict more data points with this code
xs <- seq(min(rdata$c_treattemp), max(rdata$c_treattemp), length.out=100)
preds <- data.frame(predict(mv_long, newmods=unname(poly(xs, degree=3, raw=TRUE)), addx=TRUE))
