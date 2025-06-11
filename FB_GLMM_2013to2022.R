library(glmmTMB)
library(readODS)
library(tidyverse)
library(deming)
library(brms)
library(lme4)
library(car)
library(loo)
library(cv)
library(multcomp)
library(conflicted)
library(bayesplot)
library(bayestestR)
options(mc.cores = parallel::detectCores())

conflicted::conflict_prefer("%--%", "lubridate")
conflicted::conflict_prefer("accumulate", "purrr")
conflicted::conflict_prefer("ar", "brms")
conflicted::conflict_prefer("as_data_frame", "dplyr")
conflicted::conflict_prefer("combine", "dplyr")
conflicted::conflict_prefer("compare", "loo")
conflicted::conflict_prefer("compose", "purrr")
conflicted::conflict_prefer("crossing", "tidyr")
conflicted::conflict_prefer("decompose", "stats")
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflict_prefer("extract", "tidyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("geyser", "MASS")
conflicted::conflict_prefer("kidney", "survival")
conflicted::conflict_prefer("lag", "dplyr")
conflicted::conflict_prefer("lognormal", "glmmTMB")
conflicted::conflict_prefer("loo", "loo")
conflicted::conflict_prefer("ngrps", "lme4")
conflicted::conflict_prefer("pack", "tidyr")
conflicted::conflict_prefer("recode", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("simplify", "purrr")
conflicted::conflict_prefer("some", "purrr")
conflicted::conflict_prefer("spectrum", "stats")
conflicted::conflict_prefer("union", "base")
conflicted::conflict_prefer("unpack", "tidyr")
conflicted::conflict_prefer("when", "purrr")
conflicted::conflict_prefer("smiths", "tidyr")
conflicted::conflict_prefer("rhat", "bayesplot")
conflicted::conflict_prefer("layout", "plotly")
conflicted::conflict_prefer("parnames", "brms")
conflicted::conflict_prefer("slice", "plotly")


# Check for additional conflicts
conflicted::conflict_scout()
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


# Set working directory
setwd("ENTER WORKING DIRECTORY HERE")
#-------------------------------------
#-------------------------------------
# 2024

# Function to read specific columns from LibreOffice Calc file by sheet and
# column name.
read_libreoffice_file_columns <-
  function(file_path, sheet_name, columns) {
    data <- read_ods(file_path, sheet = sheet_name)
    names(data) <- gsub(" ", "", names(data))
    data <- dplyr::select(data, one_of(columns))
    return(data)
  }

# Define the folder path
file_path <- "ENTER PATH TO GLMData_Entiat2013to2022.ods"

# Specify the names of the data sheet included in the .ods files that you want to read
sheet <- 'Sheet'


# Specify the column names you want to read in each sheet and file
selected_columns <- c(
  "Segment",
  "Treatment",
  "ReachName",
  "Year",
  "Area",
  "Chinook",
  "Steelhead",
  "RestorationAge"
)

# Read ods file
data <- read_libreoffice_file_columns(file_path, sheet_name = sheet, columns = selected_columns)

# Change treatment names to be more descriptive
data$Treatment[data$Treatment == "RestoredUntreated"] <- "Untreated/Restored"
data$Treatment[data$Treatment == "UnrestoredUntreated"] <- "Untreated/Unrestored"
data$Treatment[data$Treatment == "RestoredTreated"] <- "Treated/Restored"

# Create data subset from Upper Valley Segment
data_upper <- data %>% filter(Segment == "Upper")
data_upper$Treatment <- relevel(factor(data_upper$Treatment), ref = "Untreated/Restored")
# Fit the model ('restored only') for Upper Valley Segment Chinook (multi-year)
model_upper_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)) + (1|Year),
    data = data_upper,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
# Get the Bayes R2 value for the model
bayes_R2(model_upper_chinook)
# Add the LOO criterion to the model
model_upper_chinook <- add_criterion(model_upper_chinook, "loo")
prior_summary(model_upper_chinook)
y_upper <- data_upper$Chinook
y_upper_pred <- posterior_predict(model_upper_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_pred_matrix <- matrix(y_upper_pred[1, 1:341], nrow = 1)


# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper, y_upper_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_chinook)
#save the plot
# ggsave("model_upper_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_plot <- conditional_effects(model_upper_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_chinook <- upper_plot$Treatment
# Relevel the Treatment factor for plotting
data_upper$Treatment <- factor(
  data_upper$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_upper, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  geom_violin(data = data_upper, aes(x = Treatment, y = Chinook/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = upper_data_chinook, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_chinook, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Chinook / m"^2)) +
  scale_y_continuous(limits = c(0, 2.5, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_chinook.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

# Get posterior samples from the model and calculate treatment effects
model_upper_chinook_draws <- posterior_samples(model_upper_chinook)
effect_upper_chinook <- model_upper_chinook_draws %>%
  transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_chinook <- mean(effect_upper_chinook$treatment_effect)

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Calculate mode and HDI of treatment effect
mode_upper_chinook <- getmode(effect_upper_chinook$treatment_effect)
hdi_upper_chinook <- hdi(effect_upper_chinook$treatment_effect, prob = 0.95)

# Repeat the process for Steelhead in the Upper Valley Segment
data_upper$Treatment <- relevel(factor(data_upper$Treatment), ref = "Untreated/Restored")
model_upper_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)) + (1|Year),
    data = data_upper,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_steelhead)
model_upper_steelhead <- add_criterion(model_upper_steelhead, "loo")
y_upper_steelhead <- data_upper$Steelhead
y_upper_steelhead_pred <- posterior_predict(model_upper_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_steelhead_pred_matrix <- matrix(y_upper_steelhead_pred[1, 1:341], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_steelhead, y_upper_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_steelhead)
#save the plot
# ggsave("model_upper_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_steelhead_plot <- conditional_effects(model_upper_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_steelhead <- upper_steelhead_plot$Treatment
data_upper$Treatment <- factor(
  data_upper$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_upper, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  geom_violin(data = data_upper, aes(x = Treatment, y = Steelhead/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = upper_data_steelhead, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_steelhead, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Steelhead / m"^2)) +
  scale_y_continuous(limits = c(0, 2, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_steelhead.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_steelhead_draws <- posterior_samples(model_upper_steelhead)
effect_upper_steelhead <- model_upper_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_steelhead <- mean(effect_upper_steelhead$treatment_effect)
mode_upper_steelhead <- getmode(effect_upper_steelhead$treatment_effect)
hdi_upper_steelhead <- hdi(effect_upper_steelhead$treatment_effect, prob = 0.95)

# Create data subset from Lower Valley Segment
data_lower <- data %>% filter(Segment == "Lower")
data_lower$Treatment <- relevel(factor(data_lower$Treatment), ref = "Untreated/Restored")
# Fit the model ('restored/control') for Lower Valley Segment Chinook (multi-year)
model_lower_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)) + (1|Year),
    data = data_lower,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
# Get the Bayes R2 value for the model
bayes_R2(model_lower_chinook)
# Add the LOO criterion to the model
model_lower_chinook <- add_criterion(model_lower_chinook, "loo")
y_lower <- data_lower$Chinook
y_lower_pred <- posterior_predict(model_lower_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_pred_matrix <- matrix(y_lower_pred[1, 1:495], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower, y_lower_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_chinook)
#save the plot
# ggsave("model_lower_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_plot <- conditional_effects(model_lower_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_chinook <- lower_plot$Treatment
# Relevel the Treatment factor for plotting
data_lower$Treatment <- factor(
data_lower$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored", "Untreated/Unrestored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_lower, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  geom_violin(data = data_lower, aes(x = Treatment, y = Chinook/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = lower_data_chinook, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_chinook, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Chinook / m"^2)) +
  scale_y_continuous(limits = c(0, 2.5, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_chinook.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

# Get posterior samples from the model and calculate treatment effects
model_lower_chinook_draws <- posterior_samples(model_lower_chinook)
effect_lower_chinook_restun <- model_lower_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_chinook_unrrestun <- model_lower_chinook_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
# Calculate mean, mode, and HDI of treatment effects
mean_lower_chinook_restun <- mean(effect_lower_chinook_restun$treatment_effect)
mean_lower_chinook_unrrestun <- mean(effect_lower_chinook_unrrestun$treatment_effect)
mode_lower_chinook_restun <- getmode(effect_lower_chinook_restun$treatment_effect)
mode_lower_chinook_unrrestun <- getmode(effect_lower_chinook_unrrestun$treatment_effect)
hdi_lower_chinook_restun <- hdi(effect_lower_chinook_restun$treatment_effect, prob = 0.95)
hdi_lower_chinook_unrrestun <- hdi(effect_lower_chinook_unrrestun$treatment_effect, prob = 0.95)

# Repeat the process for Steelhead in the Lower Valley Segment
data_lower$Treatment <- relevel(factor(data_lower$Treatment), ref = "Untreated/Restored")
model_lower_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)) + (1|Year),
    data = data_lower,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_steelhead)
model_lower_steelhead <- add_criterion(model_lower_steelhead, "loo")
y_lower_steelhead <- data_lower$Steelhead
y_lower_steelhead_pred <- posterior_predict(model_lower_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_steelhead_pred_matrix <- matrix(y_lower_steelhead_pred[1, 1:495], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_steelhead, y_lower_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_steelhead)
#save the plot
# ggsave("model_lower_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_steelhead_plot <- conditional_effects(model_lower_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_steelhead <- lower_steelhead_plot$Treatment
data_lower$Treatment <- factor(
data_lower$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored", "Untreated/Unrestored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_lower, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  geom_violin(data = data_lower, aes(x = Treatment, y = Steelhead/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = lower_data_steelhead, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_steelhead, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Steelhead / m"^2)) +
  scale_y_continuous(limits = c(0, 2, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_steelhead.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_steelhead_draws <- posterior_samples(model_lower_steelhead)
effect_lower_steelhead_restun <- model_lower_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_steelhead_unrrestun <- model_lower_steelhead_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_steelhead_restun <- mean(effect_lower_steelhead_restun$treatment_effect)
mean_lower_steelhead_unrrestun <- mean(effect_lower_steelhead_unrrestun$treatment_effect)
mode_lower_steelhead_restun <- getmode(effect_lower_steelhead_restun$treatment_effect)
mode_lower_steelhead_unrrestun <- getmode(effect_lower_steelhead_unrrestun$treatment_effect)
hdi_lower_steelhead_restun <- hdi(effect_lower_steelhead_restun$treatment_effect, prob = 0.95)
hdi_lower_steelhead_unrrestun <- hdi(effect_lower_steelhead_unrrestun$treatment_effect, prob = 0.95)

# Fit the models for the year 2022 only for both Upper and Lower Valley Segments and both species (single-year)
data_upper_2022 <- data_upper %>% filter(Year == 2022)
data_upper_2022$Treatment <- relevel(factor(data_upper_2022$Treatment), ref = "Untreated/Restored")
model_upper_2022_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_upper_2022,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_2022_chinook)
model_upper_2022_chinook <- add_criterion(model_upper_2022_chinook, "loo")
y_upper_2022_chinook <- data_upper_2022$Chinook
y_upper_2022_chinook_pred <- posterior_predict(model_upper_2022_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_2022_chinook_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_2022_chinook_pred_matrix <- matrix(y_upper_2022_chinook_pred[1, 1:184], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_2022_chinook, y_upper_2022_chinook_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_2022_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_2022_chinook)
#save the plot
# ggsave("model_upper_2022_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_2022_chinook_plot <- conditional_effects(model_upper_2022_chinook, conditions = data.frame(Area = 1))  
# ggsave(
#   filename = "upper_chinook_2022_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_chinook_2022 <- upper_2022_chinook_plot$Treatment
data_upper_2022$Treatment <- factor(
  data_upper_2022$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_upper_2022, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  geom_violin(data = data_upper_2022, aes(x = Treatment, y = Chinook/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = upper_data_chinook_2022, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_chinook_2022, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Chinook / m"^2)) +
   scale_y_continuous(limits = c(0, 2.5, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_2022_chinook.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_2022_chinook_draws <- posterior_samples(model_upper_2022_chinook)
effect_upper_2022_chinook <- model_upper_2022_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_2022_chinook <- mean(effect_upper_2022_chinook$treatment_effect)
mode_upper_2022_chinook <- getmode(effect_upper_2022_chinook$treatment_effect)
hdi_upper_2022_chinook <- hdi(effect_upper_2022_chinook$treatment_effect, prob = 0.95)

data_upper_2022$Treatment <- relevel(factor(data_upper_2022$Treatment), ref = "Untreated/Restored")
model_upper_2022_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_upper_2022,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_2022_steelhead)
model_upper_2022_steelhead <- add_criterion(model_upper_2022_steelhead, "loo")
y_upper_2022_steelhead <- data_upper_2022$Steelhead
y_upper_2022_steelhead_pred <- posterior_predict(model_upper_2022_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_2022_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_2022_steelhead_pred_matrix <- matrix(y_upper_2022_steelhead_pred[1, 1:184], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_2022_steelhead, y_upper_2022_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_2022_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_2022_steelhead)
#save the plot
# ggsave("model_upper_2022_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_2022_steelhead_plot <- conditional_effects(model_upper_2022_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_2022_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_steelhead_2022 <- upper_2022_steelhead_plot$Treatment
data_upper_2022$Treatment <- factor(
  data_upper_2022$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_upper_2022, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  geom_violin(data = data_upper_2022, aes(x = Treatment, y = Steelhead/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = upper_data_steelhead_2022, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_steelhead_2022, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Steelhead / m"^2)) +
  scale_y_continuous(limits = c(0, 2, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_2022_steelhead.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_2022_steelhead_draws <- posterior_samples(model_upper_2022_steelhead)
effect_upper_2022_steelhead <- model_upper_2022_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_2022_steelhead <- mean(effect_upper_2022_steelhead$treatment_effect)
mode_upper_2022_steelhead <- getmode(effect_upper_2022_steelhead$treatment_effect)
hdi_upper_2022_steelhead <- hdi(effect_upper_2022_steelhead$treatment_effect, prob = 0.95)

data_lower_2022 <- data_lower %>% filter(Year == 2022)
data_lower_2022$Treatment <- relevel(factor(data_lower_2022$Treatment), ref = "Untreated/Restored")
model_lower_2022_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_2022,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_2022_chinook)
# model_lower_2022_chinook <- add_criterion(model_lower_2022_chinook, "loo")
y_lower_2022_chinook <- data_lower_2022$Chinook
y_lower_2022_chinook_pred <- posterior_predict(model_lower_2022_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_2022_chinook_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_2022_chinook_pred_matrix <- matrix(y_lower_2022_chinook_pred[1, 1:167], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_2022_chinook, y_lower_2022_chinook_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_2022_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_2022_chinook)
#save the plot
# ggsave("model_lower_2022_chinook_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_2022_chinook_plot <- conditional_effects(model_lower_2022_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_chinook_2022_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_chinook_2022 <- lower_2022_chinook_plot$Treatment
data_lower_2022$Treatment <- factor(
  data_lower_2022$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored", "Untreated/Unrestored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_lower_2022, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  geom_violin(data = data_lower_2022, aes(x = Treatment, y = Chinook/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = lower_data_chinook_2022, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_chinook_2022, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Chinook / m"^2)) +
  scale_y_continuous(limits = c(0, 2.5, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_2022_chinook.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_2022_chinook_draws <- posterior_samples(model_lower_2022_chinook)
effect_lower_2022_chinook_restun <- model_lower_2022_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_2022_chinook_unrrestun <- model_lower_2022_chinook_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_2022_chinook_restun <- mean(effect_lower_2022_chinook_restun$treatment_effect)
mean_lower_2022_chinook_unrrestun <- mean(effect_lower_2022_chinook_unrrestun$treatment_effect)
mode_lower_2022_chinook_restun <- getmode(effect_lower_2022_chinook_restun$treatment_effect)
mode_lower_2022_chinook_unrrestun <- getmode(effect_lower_2022_chinook_unrrestun$treatment_effect)
hdi_lower_2022_chinook_restun <- hdi(effect_lower_2022_chinook_restun$treatment_effect, prob = 0.95)
hdi_lower_2022_chinook_unrrestun <- hdi(effect_lower_2022_chinook_unrrestun$treatment_effect, prob = 0.95)

data_lower_2022$Treatment <- relevel(factor(data_lower_2022$Treatment), ref = "Untreated/Restored")
model_lower_2022_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_2022,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_2022_steelhead)
model_lower_2022_steelhead <- add_criterion(model_lower_2022_steelhead, "loo")
y_lower_2022_steelhead <- data_lower_2022$Steelhead
y_lower_2022_steelhead_pred <- posterior_predict(model_lower_2022_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_2022_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_2022_steelhead_pred_matrix <- matrix(y_lower_2022_steelhead_pred[1, 1:167], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_2022_steelhead, y_lower_2022_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_2022_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_2022_steelhead)
#save the plot
# ggsave("model_lower_2022_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_2022_steelhead_plot <- conditional_effects(model_lower_2022_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_2022_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_steelhead_2022 <- lower_2022_steelhead_plot$Treatment
data_lower_2022$Treatment <- factor(
  data_lower_2022$Treatment,
  levels = c("Treated/Restored", "Untreated/Restored", "Untreated/Unrestored")
)
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_lower_2022, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  geom_violin(data = data_lower_2022, aes(x = Treatment, y = Steelhead/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = lower_data_steelhead_2022, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_steelhead_2022, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Steelhead / m"^2)) +
  scale_y_continuous(limits = c(0, 2, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_2022_steelhead.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_2022_steelhead_draws <- posterior_samples(model_lower_2022_steelhead)
effect_lower_2022_steelhead_restun <- model_lower_2022_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_2022_steelhead_unrrestun <- model_lower_2022_steelhead_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_2022_steelhead_restun <- mean(effect_lower_2022_steelhead_restun$treatment_effect)
mean_lower_2022_steelhead_unrrestun <- mean(effect_lower_2022_steelhead_unrrestun$treatment_effect)
mode_lower_2022_steelhead_restun <- getmode(effect_lower_2022_steelhead_restun$treatment_effect)
mode_lower_2022_steelhead_unrrestun <- getmode(effect_lower_2022_steelhead_unrrestun$treatment_effect)
hdi_lower_2022_steelhead_restun <- hdi(effect_lower_2022_steelhead_restun$treatment_effect, prob = 0.95)
hdi_lower_2022_steelhead_unrrestun <- hdi(effect_lower_2022_steelhead_unrrestun$treatment_effect, prob = 0.95)

# Fit the models for the year 2016 only for both Upper and Lower Valley Segments and both species (single-year)
data_upper_2016 <- data_upper %>% filter(Year == 2016)
data_upper_2016$Treatment <- relevel(factor(data_upper_2016$Treatment), ref = "Untreated/Restored")
model_upper_2016_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_upper_2016,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_2016_chinook)
model_upper_2016_chinook <- add_criterion(model_upper_2016_chinook, "loo")
y_upper_2016 <- data_upper_2016$Chinook
y_upper_2016_pred <- posterior_predict(model_upper_2016_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_2016_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_2016_pred_matrix <- matrix(y_upper_2016_pred[1, 1:40], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_2016, y_upper_2016_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_2016_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_2016_chinook)
#save the plot
# ggsave("model_upper_2016_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_2016_plot <- conditional_effects(model_upper_2016_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_2016_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_chinook_2016 <- upper_2016_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_upper_2016, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = upper_data_chinook_2016, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_chinook_2016, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Chinook / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_2016_chinook.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_2016_chinook_draws <- posterior_samples(model_upper_2016_chinook)
effect_upper_2016_chinook <- model_upper_2016_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_2016_chinook <- mean(effect_upper_2016_chinook$treatment_effect)
mode_upper_2016_chinook <- getmode(effect_upper_2016_chinook$treatment_effect)
hdi_upper_2016_chinook <- hdi(effect_upper_2016_chinook$treatment_effect, prob = 0.95)

model_upper_2016_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_upper_2016,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_2016_steelhead)
model_upper_2016_steelhead <- add_criterion(model_upper_2016_steelhead, "loo")
y_upper_2016_steelhead <- data_upper_2016$Steelhead
y_upper_2016_steelhead_pred <- posterior_predict(model_upper_2016_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_2016_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_2016_steelhead_pred_matrix <- matrix(y_upper_2016_steelhead_pred[1, 1:40], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_2016_steelhead, y_upper_2016_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_2016_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_2016_steelhead)
#save the plot
# ggsave("model_upper_2016_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_2016_steelhead_plot <- conditional_effects(model_upper_2016_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_2016_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_steelhead_2016 <- upper_2016_steelhead_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_upper_2016, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = upper_data_steelhead_2016, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_steelhead_2016, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Steelhead / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_2016_steelhead.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_2016_steelhead_draws <- posterior_samples(model_upper_2016_steelhead)
effect_upper_2016_steelhead <- model_upper_2016_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_2016_steelhead <- mean(effect_upper_2016_steelhead$treatment_effect)
mode_upper_2016_steelhead <- getmode(effect_upper_2016_steelhead$treatment_effect)
hdi_upper_2016_steelhead <- hdi(effect_upper_2016_steelhead$treatment_effect, prob = 0.95)

data_lower_2016 <- data_lower %>% filter(Year == 2016)
data_lower_2016$Treatment <- relevel(factor(data_lower_2016$Treatment), ref = "Untreated/Restored")
model_lower_2016_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_2016,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_2016_chinook)
model_lower_2016_chinook <- add_criterion(model_lower_2016_chinook, "loo")
y_lower_2016 <- data_lower_2016$Chinook
y_lower_2016_pred <- posterior_predict(model_lower_2016_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_2016_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_2016_pred_matrix <- matrix(y_lower_2016_pred[1, 1:165], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_2016, y_lower_2016_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_2016_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_2016_chinook)
#save the plot
# ggsave("model_lower_2016_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_2016_plot <- conditional_effects(model_lower_2016_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_2016_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_chinook_2016 <- lower_2016_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_lower_2016, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = lower_data_chinook_2016, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_chinook_2016, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Chinook / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_2016_chinook.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_2016_chinook_draws <- posterior_samples(model_lower_2016_chinook)
effect_lower_2016_chinook_restun <- model_lower_2016_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_2016_chinook_unrrestun <- model_lower_2016_chinook_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_2016_chinook_restun <- mean(effect_lower_2016_chinook_restun$treatment_effect)
mean_lower_2016_chinook_unrrestun <- mean(effect_lower_2016_chinook_unrrestun$treatment_effect)
mode_lower_2016_chinook_restun <- getmode(effect_lower_2016_chinook_restun$treatment_effect)
mode_lower_2016_chinook_unrrestun <- getmode(effect_lower_2016_chinook_unrrestun$treatment_effect)
hdi_lower_2016_chinook_restun <- hdi(effect_lower_2016_chinook_restun$treatment_effect, prob = 0.95)
hdi_lower_2016_chinook_unrrestun <- hdi(effect_lower_2016_chinook_unrrestun$treatment_effect, prob = 0.95)

model_lower_2016_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_2016,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_2016_steelhead)
model_lower_2016_steelhead <- add_criterion(model_lower_2016_steelhead, "loo")
y_lower_2016_steelhead <- data_lower_2016$Steelhead
y_lower_2016_steelhead_pred <- posterior_predict(model_lower_2016_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_2016_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_2016_steelhead_pred_matrix <- matrix(y_lower_2016_steelhead_pred[1, 1:165], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_2016_steelhead, y_lower_2016_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_2016_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_2016_steelhead)
#save the plot
# ggsave("model_lower_2016_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_2016_steelhead_plot <- conditional_effects(model_lower_2016_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_2016_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_steelhead_2016 <- lower_2016_steelhead_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_lower_2016, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = lower_data_steelhead_2016, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_steelhead_2016, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Steelhead / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_2016_steelhead.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_2016_steelhead_draws <- posterior_samples(model_lower_2016_steelhead)
effect_lower_2016_steelhead_restun <- model_lower_2016_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_2016_steelhead_unrrestun <- model_lower_2016_steelhead_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_2016_steelhead_restun <- mean(effect_lower_2016_steelhead_restun$treatment_effect)
mean_lower_2016_steelhead_unrrestun <- mean(effect_lower_2016_steelhead_unrrestun$treatment_effect)
mode_lower_2016_steelhead_restun <- getmode(effect_lower_2016_steelhead_restun$treatment_effect)
mode_lower_2016_steelhead_unrrestun <- getmode(effect_lower_2016_steelhead_unrrestun$treatment_effect)
hdi_lower_2016_steelhead_restun <- hdi(effect_lower_2016_steelhead_restun$treatment_effect, prob = 0.95)
hdi_lower_2016_steelhead_unrrestun <- hdi(effect_lower_2016_steelhead_unrrestun$treatment_effect, prob = 0.95)

# Fit the models for the year 2016 only for both Upper and Lower Valley Segments and both species (single-year)
data_upper_2015 <- data_upper %>% filter(Year == 2015)
data_upper_2015$Treatment <- relevel(factor(data_upper_2015$Treatment), ref = "Untreated/Restored")
model_upper_2015_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_upper_2015,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_2015_chinook)
model_upper_2015_chinook <- add_criterion(model_upper_2015_chinook, "loo")
y_upper_2015 <- data_upper_2015$Chinook
y_upper_2015_pred <- posterior_predict(model_upper_2015_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_2015_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_2015_pred_matrix <- matrix(y_upper_2015_pred[1, 1:39], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_2015, y_upper_2015_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_2015_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_2015_chinook)
#save the plot
# ggsave("model_upper_2015_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_2015_plot <- conditional_effects(model_upper_2015_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_2015_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_chinook_2015 <- upper_2015_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_upper_2015, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = upper_data_chinook_2015, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_chinook_2015, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Chinook / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_2015_chinook.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_2015_chinook_draws <- posterior_samples(model_upper_2015_chinook)
effect_upper_2015_chinook <- model_upper_2015_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_2015_chinook <- mean(effect_upper_2015_chinook$treatment_effect)
mode_upper_2015_chinook <- getmode(effect_upper_2015_chinook$treatment_effect)
hdi_upper_2015_chinook <- hdi(effect_upper_2015_chinook$treatment_effect, prob = 0.95)

model_upper_2015_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_upper_2015,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_2015_steelhead)
model_upper_2015_steelhead <- add_criterion(model_upper_2015_steelhead, "loo")
y_upper_2015_steelhead <- data_upper_2015$Steelhead
y_upper_2015_steelhead_pred <- posterior_predict(model_upper_2015_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_2015_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_upper_2015_steelhead_pred_matrix <- matrix(y_upper_2015_steelhead_pred[1, 1:39], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_2015_steelhead, y_upper_2015_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_upper_2015_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_upper_2015_steelhead)
#save the plot
# ggsave("model_upper_2015_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

upper_2015_steelhead_plot <- conditional_effects(model_upper_2015_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "upper_2015_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
upper_data_steelhead_2015 <- upper_2015_steelhead_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_upper_2015, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = upper_data_steelhead_2015, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_steelhead_2015, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Steelhead / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_2015_steelhead.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_upper_2015_steelhead_draws <- posterior_samples(model_upper_2015_steelhead)
effect_upper_2015_steelhead <- model_upper_2015_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_2015_steelhead <- mean(effect_upper_2015_steelhead$treatment_effect)
mode_upper_2015_steelhead <- getmode(effect_upper_2015_steelhead$treatment_effect)
hdi_upper_2015_steelhead <- hdi(effect_upper_2015_steelhead$treatment_effect, prob = 0.95)

data_lower_2015 <- data_lower %>% filter(Year == 2015)
data_lower_2015$Treatment <- relevel(factor(data_lower_2015$Treatment), ref = "Untreated/Restored")
model_lower_2015_chinook <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_2015,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_2015_chinook)
model_lower_2015_chinook <- add_criterion(model_lower_2015_chinook, "loo")
y_lower_2015 <- data_lower_2015$Chinook
y_lower_2015_pred <- posterior_predict(model_lower_2015_chinook, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_2015_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_2015_pred_matrix <- matrix(y_lower_2015_pred[1, 1:163], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_2015, y_lower_2015_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_2015_chinook.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_2015_chinook)
#save the plot
# ggsave("model_lower_2015_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_2015_plot <- conditional_effects(model_lower_2015_chinook, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_2015_conditional_effects_plot_chinook.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_chinook_2015 <- lower_2015_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_lower_2015, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = lower_data_chinook_2015, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_chinook_2015, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Chinook / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_2015_chinook.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_2015_chinook_draws <- posterior_samples(model_lower_2015_chinook)
effect_lower_2015_chinook_restun <- model_lower_2015_chinook_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_2015_chinook_unrrestun <- model_lower_2015_chinook_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_2015_chinook_restun <- mean(effect_lower_2015_chinook_restun$treatment_effect)
mean_lower_2015_chinook_unrrestun <- mean(effect_lower_2015_chinook_unrrestun$treatment_effect)
mode_lower_2015_chinook_restun <- getmode(effect_lower_2015_chinook_restun$treatment_effect)
mode_lower_2015_chinook_unrrestun <- getmode(effect_lower_2015_chinook_unrrestun$treatment_effect)
hdi_lower_2015_chinook_restun <- hdi(effect_lower_2015_chinook_restun$treatment_effect, prob = 0.95)
hdi_lower_2015_chinook_unrrestun <- hdi(effect_lower_2015_chinook_unrrestun$treatment_effect, prob = 0.95)

model_lower_2015_steelhead <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_2015,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_lower_2015_steelhead)
model_lower_2015_steelhead <- add_criterion(model_lower_2015_steelhead, "loo")
y_lower_2015_steelhead <- data_lower_2015$Steelhead
y_lower_2015_steelhead_pred <- posterior_predict(model_lower_2015_steelhead, draw = 500)

# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_2015_steelhead_pred)
color_scheme_set("brightblue")

# Convert the vector to a matrix
y_lower_2015_steelhead_pred_matrix <- matrix(y_lower_2015_steelhead_pred[1, 1:163], nrow = 1)

# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_2015_steelhead, y_lower_2015_steelhead_pred_matrix) + xlim(0, 50)
# save the plot
# ggsave("ppc_dens_overlay_lower_2015_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

# plot(model_lower_2015_steelhead)
#save the plot
# ggsave("model_lower_2015_steelhead.png",
# width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )

lower_2015_steelhead_plot <- conditional_effects(model_lower_2015_steelhead, conditions = data.frame(Area = 1))
# ggsave(
#   filename = "lower_2015_steelhead_conditional_effects_plot_steelhead.png",
#   width = 10,           # Width in inches
#   height = 8,           # Height in inches
#   units = "in",         # Units for width and height
#   dpi = 300             # Resolution in dots per inch
# )
lower_data_steelhead_2015 <- lower_2015_steelhead_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_lower_2015, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = lower_data_steelhead_2015, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_steelhead_2015, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Steelhead / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_2015_steelhead.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)

model_lower_2015_steelhead_draws <- posterior_samples(model_lower_2015_steelhead)
effect_lower_2015_steelhead_restun <- model_lower_2015_steelhead_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
effect_lower_2015_steelhead_unrrestun <- model_lower_2015_steelhead_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_2015_steelhead_restun <- mean(effect_lower_2015_steelhead_restun$treatment_effect)
mean_lower_2015_steelhead_unrrestun <- mean(effect_lower_2015_steelhead_unrrestun$treatment_effect)
mode_lower_2015_steelhead_restun <- getmode(effect_lower_2015_steelhead_restun$treatment_effect)
mode_lower_2015_steelhead_unrrestun <- getmode(effect_lower_2015_steelhead_unrrestun$treatment_effect)
hdi_lower_2015_steelhead_restun <- hdi(effect_lower_2015_steelhead_restun$treatment_effect, prob = 0.95)
hdi_lower_2015_steelhead_unrrestun <- hdi(effect_lower_2015_steelhead_unrrestun$treatment_effect, prob = 0.95)


##########################################################################################################

# Introduce 'restoration age' models for Upper Valley Segment (multi-year) and both species
data_upper$Treatment <- relevel(factor(data_upper$Treatment), ref = "Untreated/Restored")
model_upper_chinook_age <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)) + RestorationAge + Treatment:RestorationAge,
    data = data_upper,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4)
bayes_R2(model_upper_chinook_age)
model_upper_chinook_age <- add_criterion(model_upper_chinook_age, "loo")
y_upper_chinook_age <- data_upper$Chinook
y_upper_chinook_age_pred <- posterior_predict(model_upper_chinook_age, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_chinook_age_pred)
color_scheme_set("brightblue")
# Convert the vector to a matrix
y_upper_chinook_age_pred_matrix <- matrix(y_upper_chinook_age_pred[1, 1:341], nrow = 1)
# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_chinook_age, y_upper_chinook_age_pred_matrix) + xlim(0, 50)
# plot(model_upper_chinook_age)
upper_chinook_age_plot <- conditional_effects(model_upper_chinook_age, conditions = data.frame(Area = 1))
upper_data_chinook_age <- upper_chinook_age_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_upper, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  geom_violin(data = data_upper, aes(x = Treatment, y = Chinook/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = upper_data_chinook_age, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_chinook_age, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Chinook / m"^2)) +
   scale_y_continuous(limits = c(0, 2.5, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_chinook_age.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)
model_upper_chinook_age_draws <- posterior_samples(model_upper_chinook_age)
effect_upper_chinook_age <- model_upper_chinook_age_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_chinook_age <- mean(effect_upper_chinook_age$treatment_effect)
mode_upper_chinook_age <- getmode(effect_upper_chinook_age$treatment_effect)
hdi_upper_chinook_age <- hdi(effect_upper_chinook_age$treatment_effect, prob = 0.95)



model_upper_steelhead_age <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)) + RestorationAge + Treatment:RestorationAge,
    data = data_upper,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
bayes_R2(model_upper_steelhead_age)
model_upper_steelhead_age <- add_criterion(model_upper_steelhead_age, "loo")
y_upper_steelhead_age <- data_upper$Steelhead
y_upper_steelhead_age_pred <- posterior_predict(model_upper_steelhead_age, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_upper_steelhead_age_pred)print(summary(model_upper_chinook), digits = 3)print(summary(model_upper_chinook), digits = 3)
color_scheme_set("brightblue")
# Convert the vector to a matrix
y_upper_steelhead_age_pred_matrix <- matrix(y_upper_steelhead_age_pred[1, 1:341], nrow = 1)
# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_upper_steelhead_age, y_upper_steelhead_age_pred_matrix) + xlim(0, 50)
# plot(model_upper_steelhead_age)
upper_steelhead_age_plot <- conditional_effects(model_upper_steelhead_age, conditions = data.frame(Area = 1))
upper_data_steelhead_age <- upper_steelhead_age_plot$Treatment
data_upper$Treatment <- factor(data_upper$Treatment, levels = c("Treated/Restored", "Untreated/Restored"))
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_upper, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  geom_violin(data = data_upper, aes(x = Treatment, y = Steelhead/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = upper_data_steelhead_age, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = upper_data_steelhead_age, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Steelhead / m"^2)) +
    scale_y_continuous(limits = c(0, 2, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_upper_steelhead_age.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)
model_upper_steelhead_age_draws <- posterior_samples(model_upper_steelhead_age)
effect_upper_steelhead_age <- model_upper_steelhead_age_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_upper_steelhead_age <- mean(effect_upper_steelhead_age$treatment_effect)
mode_upper_steelhead_age <- getmode(effect_upper_steelhead_age$treatment_effect)
hdi_upper_steelhead_age <- hdi(effect_upper_steelhead_age$treatment_effect, prob = 0.95)

# # from data_lower remove all rows where Treatment is UnrestoredUntreated
data_lower_restored <- data_lower %>% filter(Treatment != "Untreated/Unrestored")
data_lower_restored$Treatment <- relevel(factor(data_lower_restored$Treatment), ref = "Untreated/Restored")

# Introduce 'restoration age' models for Lower Valley Segment (multi-year) and both species
model_lower_chinook_age <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)) + RestorationAge + Treatment:RestorationAge,
    data = data_lower_restored,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
bayes_R2(model_lower_chinook_age)
model_lower_chinook_age <- add_criterion(model_lower_chinook_age, "loo")
model_lower_chinook_treated <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_restored,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
model_lower_chinook_treated <- add_criterion(model_lower_chinook_treated, "loo")
y_lower_chinook_age <- data_lower_restored$Chinook
y_lower_chinook_age_pred <- posterior_predict(model_lower_chinook_age, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_chinook_age_pred)
color_scheme_set("brightblue")
# Convert the vector to a matrix
y_lower_chinook_age_pred_matrix <- matrix(y_lower_chinook_age_pred[1, 1:315], nrow = 1)
# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_chinook_age, y_lower_chinook_age_pred_matrix) + xlim(0, 50)
# plot(model_lower_chinook_age)
lower_chinook_age_plot <- conditional_effects(model_lower_chinook_age, conditions = data.frame(Area = 1))
lower_data_chinook_age <- lower_chinook_age_plot$Treatment
data_lower_restored$Treatment <- factor(data_lower_restored$Treatment, levels = c("Treated/Restored", "Untreated/Restored"))
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_lower_restored, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  geom_violin(data = data_lower_restored, aes(x = Treatment, y = Chinook/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = lower_data_chinook_age, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_chinook_age, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Chinook / m"^2)) +
    scale_y_continuous(limits = c(0, 2.5, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_chinook_age.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)
model_lower_chinook_age_draws <- posterior_samples(model_lower_chinook_age)
effect_lower_chinook_age <- model_lower_chinook_age_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_lower_chinook_age <- mean(effect_lower_chinook_age$treatment_effect)
mode_lower_chinook_age <- getmode(effect_lower_chinook_age$treatment_effect)
hdi_lower_chinook_age <- hdi(effect_lower_chinook_age$treatment_effect, prob = 0.95)

data_lower_restored$Treatment <- relevel(factor(data_lower_restored$Treatment), ref = "Untreated/Restored")
model_lower_steelhead_age <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)) + RestorationAge + Treatment:RestorationAge,
    data = data_lower_restored,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
bayes_R2(model_lower_steelhead_age)
model_lower_steelhead_age <- add_criterion(model_lower_steelhead_age, "loo")  
model_lower_steelhead_treated <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_restored,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
model_lower_steelhead_treated <- add_criterion(model_lower_steelhead_treated, "loo")
y_lower_steelhead_age <- data_lower_restored$Steelhead
y_lower_steelhead_age_pred <- posterior_predict(model_lower_steelhead_age, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_steelhead_age_pred)
color_scheme_set("brightblue")
# Convert the vector to a matrix
y_lower_steelhead_age_pred_matrix <- matrix(y_lower_steelhead_age_pred[1, 1:315], nrow = 1)
# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_steelhead_age, y_lower_steelhead_age_pred_matrix) + xlim(0, 50)
# plot(model_lower_steelhead_age)
lower_steelhead_age_plot <- conditional_effects(model_lower_steelhead_age, conditions = data.frame(Area = 1))
lower_data_steelhead_age <- lower_steelhead_age_plot$Treatment
data_lower_restored$Treatment <- factor(data_lower_restored$Treatment, levels = c("Treated/Restored", "Untreated/Restored"))
ggplot() +
  # Add the actual treated numbers first to send them to the back
  # geom_jitter(data = data_lower_restored, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  geom_violin(data = data_lower_restored, aes(x = Treatment, y = Steelhead/Area), fill = "gray", alpha = 0.5) +
  # Add the conditional effects layers
  geom_point(data = lower_data_steelhead_age, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_steelhead_age, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Steelhead / m"^2)) +
    scale_y_continuous(limits = c(0, 2, by = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_steelhead_age.png",
  width = 10,           # Width in inches
  height = 10,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)
model_lower_steelhead_age_draws <- posterior_samples(model_lower_steelhead_age)
effect_lower_steelhead_age <- model_lower_steelhead_age_draws %>%
transmute(treatment_effect = exp(b_Intercept + b_TreatmentTreatedDRestored)-exp(b_Intercept))
mean_lower_steelhead_age <- mean(effect_lower_steelhead_age$treatment_effect)
mode_lower_steelhead_age <- getmode(effect_lower_steelhead_age$treatment_effect)
hdi_lower_steelhead_age <- hdi(effect_lower_steelhead_age$treatment_effect, prob = 0.95)

# Compare 'restored only' and 'restoration age' models for Upper and Lower Valley Segments and both species
loo_compare(model_upper_chinook, model_upper_chinook_age)
loo_compare(model_upper_steelhead, model_upper_steelhead_age)
loo_compare(model_lower_chinook_treated, model_lower_chinook_age)
loo_compare(model_lower_steelhead_treated, model_lower_steelhead_age)

# Get untreated data for 'untreated only' models 
data_lower_unrestored <- data_lower %>% filter(Treatment != "Treated/Restored")
data_lower_unrestored$Treatment <- relevel(factor(data_lower_unrestored$Treatment), ref = "Untreated/Restored")

# Introduce 'unrestored' models for Lower Valley Segment (single year) and both species 
model_lower_chinook_unrestored <- brm(Chinook ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_unrestored,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
bayes_R2(model_lower_chinook_unrestored)
# model_lower_chinook_unrestored <- add_criterion(model_lower_chinook_unrestored, "loo")
y_lower_chinook_unrestored <- data_lower_unrestored$Chinook
y_lower_chinook_unrestored_pred <- posterior_predict(model_lower_chinook_unrestored, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_chinook_unrestored_pred)
color_scheme_set("brightblue")
# Convert the vector to a matrix
y_lower_chinook_unrestored_pred_matrix <- matrix(y_lower_chinook_unrestored_pred[1, 1:405], nrow = 1)
# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_chinook_unrestored, y_lower_chinook_unrestored_pred_matrix) + xlim(0, 50)
# plot(model_lower_chinook_unrestored)
lower_chinook_unrestored_plot <- conditional_effects(model_lower_chinook_unrestored, conditions = data.frame(Area = 1))
lower_data_chinook_unrestored <- lower_chinook_unrestored_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_lower_unrestored, aes(x = Treatment, y = Chinook/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = lower_data_chinook_unrestored, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_chinook_unrestored, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Chinook / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_chinook_unrestored.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)
model_lower_chinook_unrestored_draws <- posterior_samples(model_lower_chinook_unrestored)
effect_lower_chinook_unrestored <- model_lower_chinook_unrestored_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_chinook_unrestored <- mean(effect_lower_chinook_unrestored$treatment_effect)
mode_lower_chinook_unrestored <- getmode(effect_lower_chinook_unrestored$treatment_effect)
hdi_lower_chinook_unrestored <- hdi(effect_lower_chinook_unrestored$treatment_effect, prob = 0.95)

model_lower_steelhead_unrestored <- brm(Steelhead ~ Treatment + (1|ReachName) + offset(log(Area)),
    data = data_lower_unrestored,
    family = zero_inflated_negbinomial,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    save_pars = save_pars(all = TRUE))
bayes_R2(model_lower_steelhead_unrestored)
# model_lower_steelhead_unrestored <- add_criterion(model_lower_steelhead_unrestored, "loo")
y_lower_steelhead_unrestored <- data_lower_unrestored$Steelhead
y_lower_steelhead_unrestored_pred <- posterior_predict(model_lower_steelhead_unrestored, draw = 500)
# use ggplot2 to plot the mean of predicted vlaued and the 95% credible interval from the model
dim(y_lower_steelhead_unrestored_pred)
color_scheme_set("brightblue")
# Convert the vector to a matrix
y_lower_steelhead_unrestored_pred_matrix <- matrix(y_lower_steelhead_unrestored_pred[1, 1:405], nrow = 1)
# Use ppc_dens_overlay with the matrix
ppc_dens_overlay(y_lower_steelhead_unrestored, y_lower_steelhead_unrestored_pred_matrix) + xlim(0, 50)
# plot(model_lower_steelhead_unrestored)
lower_steelhead_unrestored_plot <- conditional_effects(model_lower_steelhead_unrestored, conditions = data.frame(Area = 1))
lower_data_steelhead_unrestored <- lower_steelhead_unrestored_plot$Treatment
ggplot() +
  # Add the actual treated numbers first to send them to the back
  geom_point(data = data_lower_unrestored, aes(x = Treatment, y = Steelhead/Area), size = 3, color = "gray") +
  # Add the conditional effects layers
  geom_point(data = lower_data_steelhead_unrestored, aes(x = Treatment, y = estimate__), size = 3.5, color = "black") +
  geom_errorbar(data = lower_data_steelhead_unrestored, aes(x = Treatment, ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, color = "black") +
   labs(x = "Treatment", y = expression("Number Steelhead / m"^2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ggsave(
  filename = "glmm_lower_steelhead_unrestored.png",
  width = 10,           # Width in inches
  height = 8,           # Height in inches
  units = "in",         # Units for width and height
  dpi = 300             # Resolution in dots per inch
)
model_lower_steelhead_unrestored_draws <- posterior_samples(model_lower_steelhead_unrestored)
effect_lower_steelhead_unrestored <- model_lower_steelhead_unrestored_draws %>%
  transmute(treatment_effect = exp(b_Intercept)-exp(b_Intercept + b_TreatmentUntreatedDUnrestored))
mean_lower_steelhead_unrestored <- mean(effect_lower_steelhead_unrestored$treatment_effect)
mode_lower_steelhead_unrestored <- getmode(effect_lower_steelhead_unrestored$treatment_effect)
hdi_lower_steelhead_unrestored <- hdi(effect_lower_steelhead_unrestored$treatment_effect, prob = 0.95)