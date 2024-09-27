library(tidyverse)
library(psych)
library(readr)
library(ggplot2)
library(tidyr)
library(DESeq2)
library(apeglm)
library(biomaRt)
library(gplots)
library(gtools)
library(ggrepel)
library(Hmisc)
library(plotly)
library(ggpubr)
library(lme4)
library(MASS) # load it only when needed

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                      dataset = "hsapiens_gene_ensembl", 
                      host = "www.ensembl.org")


# Get data
thesis_data_strict <- read_csv("thesis_data_strict2.csv")

tdf <- as.data.frame(thesis_data_strict)

#Descriptive statistics
summary(tdf$AtoI_TOTALmatches)
summary(tdf$AtoI)
summary(tdf$CtoU)

#Histograms to visualise distributions

ggplot(tdf, aes(x=CtoU)) +
  geom_histogram(binwidth=40, color="black", fill="cyan") +
  labs(x="Count of C to U edits", y="Frequency") +
  theme_minimal() +
  ggtitle("Histogram of counts of C to U edits")


ggplot(tdf, aes(x=AtoI_TOTALmatches)) +
  geom_histogram(binwidth=80, color="black", fill="lightgreen") +
  labs(x="Count of A to I edits per sample", y="Frequency") +
  theme_minimal() +
  ggtitle("Histogram of counts of A to I edits")

ggplot(tdf, aes(x=AtoI)) +
  geom_histogram(binwidth=200, color="black", fill="lightgreen") +
  labs(x="Count of A-to-I edits per sample", y="Frequency") +
  theme_dark() +
  ggtitle("Histogram of counts of A-to-I edits")

#trying facetwrap

tdf_long <- tdf %>%
  pivot_longer(cols = c(AtoI, CtoU), 
               names_to = "Edit_Type", 
               values_to = "Count")

# Combined plot using facet_wrap and the dark theme
ggplot(tdf_long, aes(x=Count)) +
  geom_histogram(data = subset(tdf_long, Edit_Type == "AtoI"), 
                 binwidth = 300, color = "black", fill = "lightgreen") +
  geom_histogram(data = subset(tdf_long, Edit_Type == "CtoU"), 
                 binwidth = 80, color = "black", fill = "lightblue") +
  facet_wrap(~Edit_Type, scales = "free_x", 
             labeller = as_labeller(c(AtoI = "A-to-I", CtoU = "C-to-U"))) +
  labs(x = "Count of RNA edits per sample", y = "Frequency") +
  ggtitle("Histograms of counts of RNA edits") +
  theme_dark() +
  theme(plot.title = element_text(hjust = 0.5))

#Checking for normality and equality of variance of cytokine data
shapiro.test(tdf$LPS_IL_10)


#Scatter plot to see correlation between atoi and ctou
cor_result <- cor.test(tdf$AtoI, tdf$CtoU, method = "spearman")

cor_result
# Extract the correlation coefficient and p-value
spearman_corr <- round(cor_result$estimate, 2)
p_value <- round(cor_result$p.value, 10)
cor_result$p.value
# Create the plot with the Spearman correlation coefficient and p-value
# Create the scatter plot with one global regression line and Spearman correlation
ggscatter(tdf, x = "AtoI", y = "CtoU", 
          add = "reg.line", conf.int = TRUE,  # Global regression line with confidence interval
          cor.coef = TRUE, cor.method = "spearman",  # Display Spearman correlation and p-value
          xlab = "A-to-I count", ylab = "C-to-U count",
          color = "Time_point",  # Color points by Time_point but use a global trendline
          legend.title = "Group", add.params = list(color = "black")) +
  geom_point(size = 3, alpha = 0.7, shape = 21, stroke = 0.5, aes(fill = as.factor(Time_point)), color = "black") +  # Black outline for points
  ggtitle("Correlation of Counts of RNA Edits") +
  theme_dark() +  # Apply dark theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust x-axis text
        legend.position = "bottom")




#Find the mean variant counts
mean_var_counts <- tdf %>% 
  select(Time_point, CtoU, AtoI) %>% 
  group_by(Time_point) %>% 
  summarise("C to U" = mean(CtoU), "A to I" = mean(AtoI))

mean_var_counts

# Descriptive statistics of each group
descriptive_ctou_acute <- tdf %>%
  select(Time_point,CtoU) %>% 
  filter(Time_point == "acute") %>% 
  summary()

descriptive_ctou_recovery <- tdf %>%
  select(Time_point,CtoU) %>% 
  filter(Time_point == "recovery") %>% 
  summary()

descriptive_ctou_control <- tdf %>%
  select(Time_point,CtoU) %>% 
  filter(Time_point == "control") %>% 
  summary()


descriptive_ctou_acute
descriptive_ctou_recovery
descriptive_ctou_control

descriptive_ctou_acute <- tdf %>%
  select(Time_point,CtoU) %>% 
  filter(Time_point == "acute") %>%
  pull(CtoU)

descriptive_ctou_recovery <- tdf %>%
  select(Time_point,CtoU) %>% 
  filter(Time_point == "recovery") %>%
  pull(CtoU)

descriptive_ctou_control <- tdf %>%
  select(Time_point,CtoU) %>% 
  filter(Time_point == "control") %>%
  pull(CtoU)

hist(descriptive_ctou_acute$CtoU)
hist(descriptive_ctou_recovery$CtoU)
hist(descriptive_ctou_control$CtoU)

descriptive_atoi_acute <- tdf %>%
  select(Time_point,AtoI) %>% 
  filter(Time_point == "acute") %>%
  pull(AtoI)

descriptive_atoi_recovery <- tdf %>%
  select(Time_point,AtoI) %>% 
  filter(Time_point == "recovery") %>%
  pull(AtoI)

descriptive_atoi_control <- tdf %>%
  select(Time_point,AtoI) %>% 
  filter(Time_point == "control") %>%
  pull(AtoI)

wilcox.test(descriptive_atoi_acute, descriptive_atoi_control) # not significant
wilcox.test(descriptive_ctou_acute, descriptive_ctou_control) # not significant

# negative binomial atoi control vs acute
combined_data <- data.frame(
  count = c(descriptive_atoi_acute, descriptive_atoi_control),
  group = factor(c(rep("acute", length(descriptive_atoi_acute)), 
                   rep("control", length(descriptive_atoi_control))))
)

combined_data$group <- relevel(combined_data$group, ref = "control")


nb_model <- glm.nb(count ~ group, data = combined_data)
summary(nb_model)

# negative binomial atoi control vs recovery
combined_data <- data.frame(
  count = c(descriptive_atoi_recovery, descriptive_atoi_control),
  group = factor(c(rep("recovery", length(descriptive_atoi_recovery)), 
                   rep("control", length(descriptive_atoi_control))))
)
combined_data
nb_model <- glm.nb(count ~ group, data = combined_data)
nb_model
# negative binomial ctou control vs acute
combined_data <- data.frame(
  count = c(descriptive_ctou_acute, descriptive_ctou_control),
  group = factor(c(rep("acute", length(descriptive_ctou_acute)), 
                   rep("control", length(descriptive_ctou_control))))
)

nb_model <- glm.nb(count ~ group, data = combined_data)
summary (nb_model)

# negative binomial ctou control vs recovery
combined_data <- data.frame(
  count = c(descriptive_ctou_recovery, descriptive_ctou_control),
  group = factor(c(rep("recovery", length(descriptive_ctou_recovery)), 
                   rep("control", length(descriptive_ctou_control))))
)

nb_model <- glm.nb(count ~ group, data = combined_data)

# Print the model summary
summary(nb_model)

stat<-describe(descriptive_atoi_acute$AtoI)
summary(descriptive_atoi_recovery$AtoI)
summary(descriptive_atoi_control$AtoI)

hist(descriptive_atoi_acute$AtoI)
hist(descriptive_atoi_recovery$AtoI)
hist(descriptive_atoi_control$AtoI)

# trying to plot a 3d histogram

fig <- plot_ly() %>%
  # Acute group
  add_trace(x = h1$mids, y = rep(1, length(h1$mids)), z = h1$counts,
            type = 'scatter3d', mode = 'lines+markers',
            marker = list(size = 5, color = 'red'),
            line = list(color = 'red', width = 5),
            name = 'Acute') %>%
  # Recovery group
  add_trace(x = h2$mids, y = rep(2, length(h2$mids)), z = h2$counts,
            type = 'scatter3d', mode = 'lines+markers',
            marker = list(size = 5, color = 'blue'),
            line = list(color = 'blue', width = 5),
            name = 'Recovery') %>%
  # Control group
  add_trace(x = h3$mids, y = rep(3, length(h3$mids)), z = h3$counts,
            type = 'scatter3d', mode = 'lines+markers',
            marker = list(size = 5, color = 'green'),
            line = list(color = 'green', width = 5),
            name = 'Control') %>%
  layout(scene = list(
    xaxis = list(title = 'RNA-editing Count (C-to-U)'),
    yaxis = list(title = 'Group', tickvals = c(1, 2, 3), ticktext = c('Acute', 'Recovery', 'Control')),
    zaxis = list(title = 'Count'),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
  ))

fig

# Create the 2D overlay of histograms
data_combined <- rbind(
  data.frame(CtoU = descriptive_ctou_control$CtoU, Group = 'Control'),
  data.frame(CtoU = descriptive_ctou_recovery$CtoU, Group = 'Recovery'),
  data.frame(CtoU = descriptive_ctou_acute$CtoU, Group = 'Acute')
)

data_combined$Group <- factor(data_combined$Group, levels = c('Control', 'Recovery', 'Acute'))

# Create a single ggplot with faceting
ggplot(data_combined, aes(x = CtoU)) +
  geom_histogram(aes(fill = Group), color = "black", alpha = 0.6, bins = 30) +
  scale_fill_manual(values = c('Control' = 'lightgreen', 'Recovery' = 'lightblue', 'Acute' = 'pink')) +
  facet_wrap(~ Group, scales = 'free_y') +
  labs(title = 'Histograms of C-to-U RNA-editing Counts by Group', 
       x = 'Count of A-to-I editing sites', 
       y = 'Count') +
  theme_dark() +
  theme(legend.position = 'none')

# Same as above but A to I

# Create the 2D overlay of histograms
data_combined_atoi <- rbind(
  data.frame(CtoU = descriptive_atoi_control$AtoI, Group = 'Control'),
  data.frame(CtoU = descriptive_atoi_recovery$AtoI, Group = 'Recovery'),
  data.frame(CtoU = descriptive_atoi_acute$AtoI, Group = 'Acute')
)

data_combined_atoi$Group <- factor(data_combined$Group, levels = c('Control', 'Recovery', 'Acute'))

# Create a single ggplot with faceting
ggplot(data_combined_atoi, aes(x = CtoU)) +
  geom_histogram(aes(fill = Group), color = "black", alpha = 0.6, bins = 30) +
  scale_fill_manual(values = c('Control' = 'lightgreen', 'Recovery' = 'lightblue', 'Acute' = 'pink')) +
  facet_wrap(~ Group, scales = 'free_y') +
  labs(title = 'Histograms of A-to-I RNA-editing Counts by Group', 
       x = 'Count of A-to-I editing sites', 
       y = 'Count') +
  theme_dark() +
  theme(legend.position = 'none')


# Display the plot
print(plot)

# Reshape data to long format
mean_var_counts_long <- pivot_longer(mean_var_counts, cols = c("C to U", "A to I"),
                                     names_to = "Variable", values_to = "Mean")

# Plot mean C to U and A to I in the same plot
ggplot(mean_var_counts_long, aes(x = Time_point, y = Mean, fill = Time_point)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.2)) +
  labs(title = "Mean RNA editing counts",
       x = "Time Point",
       y = "Mean count of variants") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  guides(fill = guide_legend(title = "Time point"))

#box plots showing C to U and A to I variants by Time_point
ggplot(tdf, aes(x = Time_point, y = CtoU, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Distribution of C to U edit counts",
       x = "Time point",
       y = "count of C-to-U variants",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(tdf, aes(x = Time_point, y = AtoI, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Distribution of all A to I edit counts",
       x = "Time point",
       y = "count of A to I variants", color = "Time point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# facetwrap to combine them

tdf_long <- tdf %>%
  pivot_longer(cols = c(CtoU, AtoI),
               names_to = "Edit_Type",
               values_to = "Count")

# Create the combined plot
ggplot(tdf_long, aes(x = Time_point, y = Count, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA, aes(fill = Time_point),
               color = "black") +  # Set boxplot outline color to black
  geom_point(position = position_jitter(width = 0.1), size = 1, color = "black") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.75) +  # Add black median lines
  labs(title = "Distribution of Edit Counts",
       x = "Group",
       y = "Count of Variants") +
  facet_wrap(~ Edit_Type, 
             scales = "free_y",
             labeller = labeller(Edit_Type = c(CtoU = "C-to-U", AtoI = "A-to-I"))) +
  theme_dark() +  # Apply dark theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggplot(tdf, aes(x = Time_point, y = AtoIandCtoU, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Distribution of A to I edit counts",
       x = "Time point",
       y = "count of A to I variants", color = "Time point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Q-Q plot for CtoU
qqnorm(tdf$CtoU)
qqline(tdf$CtoU)

# Q-Q plot for AtoI
qqnorm(tdf$AtoI_TOTALmatches)
qqline(tdf$AtoI_TOTALmatches)


# Convert log counts to counts
acute_count <- exp(intercept)
control_count <- exp(intercept + control_estimate)
recovery_count <- exp(intercept + recovery_estimate)

# Calculate percentage changes
percentage_change_control <- ((acute_count - control_count) / acute_count) * 100
percentage_change_recovery <- ((acute_count - recovery_count) / acute_count) * 100

percentage_change_control
percentage_change_recovery

#Negative binomial test results
tdf$Time_point <- factor(tdf$Time_point, levels = c("acute", "control", "recovery"))
nb_test_result_CtoU <- glm.nb(CtoU ~ Time_point, data = tdf)
nb_test_result_AtoI <- glm.nb(AtoI ~ Time_point, data = tdf)
summary(nb_test_result_CtoU)
summary(nb_test_result_AtoI)

rnaedit_no_acute <- tdf %>% 
  select(AtoI, CtoU, Time_point) %>% 
  filter(Time_point != "acute")

nb_test_result_no_acute <- glm.nb(AtoI ~ Time_point, data = rnaedit_no_acute)
summary(nb_test_result_no_acute)

#Kruskal-wallis test
kruskal.test(AtoI ~ Time_point, data = tdf)
kruskal.test(CtoU ~ Time_point, data = tdf)



# C to U and A to I count differences in controls, paired acute and paired rec
paired_editcounts <- tdf %>% 
  select(Participant, Paired_single, Time_point, CtoU, AtoI_TOTALmatches, AtoIandCtoU, AtoI) %>% 
  filter(Paired_single =="paired")# %>% 
  #filter(!Participant_ID =="1090") %>% #outlier
  #filter(!Participant_ID == "3008") #outlier
paired_editcounts

paired_editcounts$Participant <- as.factor(paired_editcounts$Participant)

# negative binomial with mixed effects for paired data:

library(glmmTMB)

# Fit the negative binomial mixed-effects model
nb_model <- glmmTMB(CtoU ~ Time_point + (1 | Participant), 
                    data = paired_editcounts, 
                    family = nbinom2)  # Use nbinom1 if appropriate for your dispersion

# Print the model summary
summary(nb_model)


# Paired bar plots and box plots of paired Acute vs Recovery 
ggplot(paired_editcounts, aes(x = Participant, y = CtoU, fill = Time_point)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Paired Bar Plot of C to U Counts for Acute and Recovery",
       x = "Participant ID",
       y = "C to U Counts",
       fill = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(paired_editcounts, aes(x = Time_point, y = CtoU, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Distribution of C to U edit counts - paired samples",
       x = "Time point",
       y = "count of C to U variants", color = "Time point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(paired_editcounts, aes(x = Time_point, y = AtoI, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Distribution of all A to I edit counts - paired samples",
       x = "Time point",
       y = "count of A to I variants", color = "Time point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(paired_editcounts, aes(x = Participant, y = AtoI_TOTALmatches, fill = Time_point)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Paired Bar Plot of A to I Counts for Acute and Recovery",
       x = "Participant ID",
       y = "A to I Counts",
       fill = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# negative binomial test for paired data acute vs recovery
nb_test_result_CtoU_paired <- glm.nb(CtoU ~ Time_point, data = paired_editcounts)
summary(nb_test_result_CtoU_paired)

nb_test_result_AtoI_paired <- glm.nb(AtoI_TOTALmatches ~ Time_point, data = paired_data_controls)
summary(nb_test_result_AtoI_paired)

nb_test_result_AtoI_paired <- glm.nb(AtoI ~ Time_point, data = paired_editcounts)
summary(nb_test_result_AtoI_paired)

# Wilcoxon Signed Ranks tests
# A to I
paired_editcounts_atoi_acute <- paired_editcounts %>% 
  select(Participant, Time_point, AtoI) %>% 
  filter(Time_point == "acute")

paired_editcounts_atoi_recovery <- paired_editcounts %>% 
  select(Participant, Time_point, AtoI) %>% 
  filter(Time_point == "recovery")

paired_atoi <- merge(paired_editcounts_atoi_acute, paired_editcounts_atoi_recovery,
                     by = "Participant")

wilcox.test(paired_atoi$AtoI.x, paired_atoi$AtoI.y, paired=TRUE)

# C to U
paired_editcounts_ctou_acute <- paired_editcounts %>% 
  select(Participant, Time_point, CtoU) %>% 
  filter(Time_point == "acute")

paired_editcounts_ctou_recovery <- paired_editcounts %>% 
  select(Participant, Time_point, CtoU) %>% 
  filter(Time_point == "recovery")

paired_ctou <- merge(paired_editcounts_ctou_acute, paired_editcounts_ctou_recovery,
                     by = "Participant")

wilcox.test(paired_ctou$CtoU.y, paired_ctou$CtoU.x, paired=TRUE)

# Age

CtoU_controls <- tdf %>%
  select(age, Time_point, CtoU) %>% 
  filter(Time_point == "control")

AtoI_control <- tdf %>%
  select(age, Time_point, AtoI) %>% 
  filter(Time_point == "control")

CtoU_acute <- tdf %>%
  select(age, Time_point, CtoU) %>%
  filter(Time_point == "acute")

AtoI_acute <- tdf %>%
  select(age, Time_point, AtoI) %>% 
  filter(Time_point == "acute")

CtoU_rec <- tdf %>% 
  select(age, Time_point, CtoU) %>%
  filter(Time_point == "recovery")

AtoI_recovery <- tdf %>%
  select(age, Time_point, AtoI) %>% 
  filter(Time_point == "recovery")

# Define the age breaks
age_breaks <- seq(20, 100, by = 10)

# Create age groups using cut
age_groups <- cut(tdf$age, breaks = c(-Inf, age_breaks, Inf), 
                  labels = FALSE, include.lowest = TRUE)

age_groups
hist(tdf$age)

ggplot(tdf, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(title = "Histogram of ages of study subjects",
       x = "Age",
       y = "Frequency") +
  scale_y_continuous(breaks = seq(0, max(tdf$age), by = 5)) +  # Set y-axis breaks to increments of 5
  scale_x_continuous(breaks = seq(min(tdf$age), max(tdf$age), by = 5)) +  # Set x-axis labels at 5-year intervals
  theme_dark()

# Count the number of individuals in each age group
table(age_groups)
table(tdf$age)

# Scatter plot C to U vs Age
ggplot(tdf, aes(x = age, y = CtoU, color = as.factor(Time_point))) +
  geom_point(size = 3, alpha = 0.7, shape = 21, stroke = 0.5, color = "black", aes(fill = as.factor(Time_point))) +  # Black outline with filled color
  labs(title = "Scatter Plot of C-to-U by Age",
       x = "Age",
       y = "Count of C-to-U Sites",
       fill = "Group") +  # Update legend title
  theme_dark() +  # Apply dark theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Adjust x-axis text for readability
  scale_x_continuous(breaks = seq(min(tdf$age), max(tdf$age), by = 5)) +  # Label x-axis every 5 years
  scale_y_continuous(breaks = seq(0, max(tdf$CtoU), by = 100))

# Scatter plot A to I vs age
ggplot(tdf, aes(x = age, y = AtoI, color = as.factor(Time_point))) +
  geom_point(size = 3, alpha = 0.7, shape = 21, stroke = 0.5, color = "black", aes(fill = as.factor(Time_point))) +  # Black outline with filled color
  labs(title = "Scatter Plot of A-to-I by Age",
       x = "Age",
       y = "Count of A-to-I Sites",
       fill = "Group") +  # Update legend title
  theme_dark() +  # Apply dark theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust x-axis text for readability
        legend.position = "bottom") +  # Move legend to the bottom
  scale_x_continuous(breaks = seq(min(tdf$age), max(tdf$age), by = 5))

#boxplot everyone C to U vs age groups

tdf$age_group <- cut(tdf$age, breaks = seq(20, 100, by = 10),
                     labels = FALSE)
age_groups_names <- c("31-40","41-50","51-60","61-70","71-80","81-90", "91-100")

ggplot(tdf, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in all samples",
       x = "Age Group",
       y = "Count of C to U variants") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")


ggplot(tdf, aes(x = as.factor(age_group), y = AtoI_TOTALmatches, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in all samples",
       x = "Age Group",
       y = "Count of A to I variants") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")


ggplot(AtoI_acute, aes(x = as.factor(age_group), y = AtoI, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in Acute sample",
       x = "Age Group",
       y = "Count of A to I variants") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")


# grouping individuals by age group

CtoU_acute$age_group <- cut(CtoU_acute$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

CtoU_controls$age_group <- cut(CtoU_controls$age, breaks = seq(20, 100, by = 10),
                               labels = FALSE)

CtoU_rec$age_group <- cut(CtoU_rec$age, breaks = seq(20, 100, by = 10),
                          labels = FALSE)


AtoI_acute$age_group <- cut(AtoI_acute$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

AtoI_recovery$age_group <- cut(AtoI_recovery$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

AtoI_control$age_group <- cut(AtoI_control$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

CtoU_controls$age_group <- cut(CtoU_controls$age, breaks = seq(20, 100, by = 10),
                               labels = FALSE)

CtoU_rec$age_group <- cut(CtoU_rec$age, breaks = seq(20, 100, by = 10),
                          labels = FALSE)

# Plotting box plots

ggplot(CtoU_acute, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in Acute samples",
       x = "Age Group",
       y = "Count of C to U variants")  + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

ggplot(AtoI_acute, aes(x = as.factor(age_group), y = AtoI, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in Acute samples",
       x = "Age Group",
       y = "Count of A to I variants")  + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

ggplot(AtoI_recovery, aes(x = as.factor(age_group), y = AtoI, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in recovery samples",
       x = "Age Group",
       y = "Count of A to I variants")  + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

ggplot(AtoI_control, aes(x = as.factor(age_group), y = AtoI, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in control samples",
       x = "Age Group",
       y = "Count of A to I variants")  + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

ggplot(CtoU_controls, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in Control samples",
       x = "Age Group",
       y = "Count of C to U variants")  + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")


ggplot(CtoU_rec, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in Recovery samples",
       x = "Age Group",
       y = "Count of C to U variants") + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

# trying to combine them into a single plot:

# Combine the three data frames into a single data frame with a new 'Group' column
CtoU_acute$Group <- "Acute"
CtoU_controls$Group <- "Control"
CtoU_rec$Group <- "Recovery"

combined_data_ctou <- rbind(CtoU_acute, CtoU_controls, CtoU_rec)

# Plot the combined data using facet_wrap
ggplot(combined_data_ctou, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Boxplot of C-to-U counts by Age Group",
       x = "Age Group",
       y = "Count of C-to-U edits") + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges") +
  theme_dark()

ggplot(CtoU_acute, aes(x = age, y = CtoU, color = age)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear trend line
  labs(title = "Scatter Plot of C to U by Age in Acute Group",
       x = "Age",
       y = "Count of C to U variants") +
  scale_color_gradient(low = "lightblue", high = "darkblue")

lm_ctoU_acute <- lm(CtoU ~ age, data = tdf)
summary(lm_ctoU_acute)

CtoU_acute$age_group <- as.factor(CtoU_acute$age_group)

# combining a to i vs age box plots

AtoI_acute$Group <- "Acute"
AtoI_recovery$Group <- "Recovery"
AtoI_control$Group <- "Control"

# Combine the datasets
combined_data_atoi <- rbind(AtoI_acute, AtoI_recovery, AtoI_control)

# Plot the combined data using facet_wrap
ggplot(combined_data_atoi, aes(x = as.factor(age_group), y = AtoI, fill = as.factor(age_group))) +
  geom_boxplot() +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Boxplot of A-to-I counts by Age Group",
       x = "Age Group",
       y = "Count of A-to-I edits") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names, name = "Age ranges") +
  theme_dark() +
  theme(legend.position = "bottom")

# Kruskal-Wallis test
kruskal.test(CtoU ~ age, data = tdf)


# Sex
summary(as.factor(tdf$sex))

ggplot(tdf, aes(x = as.factor(sex), y = CtoU)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = as.factor(Time_point)), position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of C to U by sex",
       x = "Sex",
       y = "Count of C to U variants")

ggplot(tdf, aes(x = as.factor(sex), y = AtoI)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = as.factor(Time_point)), position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of A to I by sex",
       x = "Sex",
       y = "Count of A to I variants")

# combining the two plots

tdf_long <- tdf %>%
  pivot_longer(cols = c(AtoI, CtoU), 
               names_to = "Edit_Type", 
               values_to = "Count")

# Reorder the levels of Edit_Type so A-to-I is on the left and C-to-U is on the right
tdf_long$Edit_Type <- factor(tdf_long$Edit_Type, levels = c("AtoI", "CtoU"))


# Generate the combined plot
ggplot(tdf_long, aes(x = as.factor(sex), y = Count, fill = as.factor(sex))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = as.factor(Time_point)),
             position = position_jitter(width = 0.1),
             size = 2) +
  facet_wrap(~ Edit_Type, scales = "free_y", labeller = as_labeller(c(AtoI = "A-to-I", CtoU = "C-to-U"))) +
  scale_fill_manual(values = c("Male" = "steelblue", "Female" = "lightpink")) +  # Set custom colors for sexes
  labs(title = "Boxplot of RNA edits by Sex",
       x = "Sex",
       y = "Count of Variants",
       fill = "Sex", color = "Group") +
  theme_dark() +
  theme(legend.position = "bottom")

#Negative binomial test results
nb_test_result_CtoU_sex <- glm.nb(CtoU ~ sex, data = tdf)
nb_test_result_AtoI_sex <- glm.nb(AtoI_TOTALmatches ~ sex, data = tdf)
summary(nb_test_result_CtoU_sex)
summary(nb_test_result_AtoI_sex)

# Mann-whitney U test
wilcox.test(CtoU ~ sex, data = tdf)
wilcox.test(AtoI ~ sex, data = tdf)

# COPD
summary(as.factor(tdf$copd))
table(tdf$Time_point, tdf$copd)

edits_copd <- tdf %>% 
  select(CtoU, AtoI, Time_point, copd)

edits_acute_copd <- tdf %>%
  select(CtoU, AtoI_TOTALmatches, AtoI, Time_point, copd) %>% 
  filter(Time_point == "acute")

ggplot(edits_acute_copd, aes(x = as.factor(copd), y = CtoU)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of C to U by COPD in acute samples",
       x = "COPD",
       y = "Count of C-to-U variants")

ggplot(edits_acute_copd, aes(x = as.factor(copd), y = AtoI)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of A to I by COPD in acute samples",
       x = "COPD",
       y = "Count of A-to-I variants")

# Combining the two plots
edits_copd_long <- edits_copd %>%
  pivot_longer(cols = c(CtoU, AtoI), 
               names_to = "Edit_Type", 
               values_to = "Count")

# Plot the combined data
ggplot(edits_copd_long, aes(x = as.factor(copd), y = Count)) +
  geom_boxplot(aes(fill = as.factor(copd)), outlier.shape = NA) +  # Boxplot filled by COPD
  geom_point(aes(color = as.factor(Time_point)), position = position_jitter(width = 0.1), size = 2, stroke = 0.5) +  # Points with black outline
  facet_wrap(~ Edit_Type, scales = "free_y", labeller = as_labeller(c(CtoU = "C-to-U", AtoI = "A-to-I"))) +
  #scale_fill_manual(values = c("0" = "steelblue", "1" = "lightpink")) +  # Custom colors for COPD
  #scale_color_manual(values = c("0" = "steelblue", "1" = "lightpink")) +  # Custom colors for points
  labs(title = "Boxplot of RNA Edits by COPD",
       x = "COPD",
       y = "Count of Edits",
       fill = "COPD", color = "COPD") +
  theme_dark() +
  theme(legend.position = "bottom")

wilcox.test(CtoU ~ copd, data = edits_copd)
wilcox.test(AtoI ~ copd, data = edits_copd)

# Kruskal tests for chromosomal analysis
# C to U

library(broom)

# Columns to consider
# Columns to consider
ctou_chr_col <- grep("^CtoU_", names(tdf), value = TRUE)

# Initialize an empty data frame to store the results
ctou_kw_results <- data.frame()

# Loop over each chromosome column
for (col in ctou_chr_col) {
  # Create the formula for the test
  formula <- as.formula(paste(col, "~ Time_point"))
  
  # Perform the test
  test_result <- kruskal.test(formula, data = tdf)
  
  # Extract the p-value
  p_value <- tidy(test_result)$p.value
  test_statistic <- tidy(test_result)$statistic
  
  # Add the result to the data frame
  ctou_kw_results <- rbind(ctou_kw_results, data.frame(Column = col, 
                                                       Test_Statistic = test_statistic,
                                                       P_Value = p_value))
}

# Adjust the p-values using the BH method
ctou_kw_results$Adjusted_P_Value <- p.adjust(ctou_kw_results$P_Value, method = "BH")

# Print the results
ctou_kw_results

#post-hoc analysis
library(FSA)
dunnTest(CtoU_chr15 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr14 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr16 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr9 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr4 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr19 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr18 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chrX ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr12 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr1 ~ Time_point, data=tdf, method="bh")
dunnTest(CtoU_chr3 ~ Time_point, data=tdf, method="bh")


chr_atoi_long <- tdf %>%
  select(Time_point, starts_with("AtoI_Ch")) %>%
  select(-ends_with("matches")) %>% 
  pivot_longer(cols = starts_with("AtoI_Ch"), names_to = "Chromosome", values_to = "Count")

chr_ctou_long <- tdf %>%
  select(Time_point, starts_with("CtoU_")) %>% 
  pivot_longer(cols = starts_with("CtoU_"), names_to = "Chromosome", values_to = "Count")

# A to I
# Columns to consider
atoi_chr_col <- grep("^AtoI_Ch", names(tdf), value = TRUE)

# remove those names that have "matches" ... a to i variants that 
# only match the REDI portal database.
not_ending_with_matches <- !grepl("matches$", atoi_chr_col)
atoi_chr_col <- atoi_chr_col[not_ending_with_matches]


# Initialize an empty data frame to store the results
atoi_kw_results <- data.frame()

# Loop over each chromosome column
for (col in atoi_chr_col) {
  # Create the formula for the test
  formula <- as.formula(paste(col, "~ Time_point"))
  
  # Perform the test
  test_result <- kruskal.test(formula, data = tdf)
  
  # Extract the p-value
  p_value <- tidy(test_result)$p.value
  test_statistic <- tidy(test_result)$statistic
  
  # Add the result to the data frame
  atoi_kw_results <- rbind(atoi_kw_results, data.frame(Column = col, 
                                                       Test_Statistic = test_statistic,
                                                       P_Value = p_value))
}

# Adjust the p-values using the Bonferroni method
atoi_kw_results$Adjusted_P_Value <- p.adjust(atoi_kw_results$P_Value, method = "BH")

# Print the results
atoi_kw_results
tdf$AtoI_Ch16matches

#post-hoc analysis
library(FSA)
dunnTest(AtoI_Chr12 ~ Time_point, data=tdf, method="bh")

# Finding the mean count of A to I or C to U variants on each chromosome
# C to U
mean_ctou_recovery <- tdf %>% 
  filter(Time_point == "recovery") %>% 
  select(starts_with("CtoU_"))

mean_ctou_controls <- tdf %>% 
  filter(Time_point == "control") %>% 
  select(starts_with("CtoU_"))

mean_ctou_acute <- tdf %>% 
  filter(Time_point == "acute") %>% 
  select(starts_with("CtoU_"))

mean_c_rec <-colMeans(mean_ctou_recovery)
mean_c_con <-colMeans(mean_ctou_controls)
mean_c_acu <-colMeans(mean_ctou_acute)

mean_ctou_acute <- as.data.frame(mean_c_acu)
mean_ctou_controls <-as.data.frame(mean_c_con)
mean_ctou_recovery <-as.data.frame(mean_c_rec)

mean_ctou <- cbind(mean_ctou_acute, mean_ctou_controls, mean_ctou_recovery)
mean_ctou_long <- melt(mean_ctou)
mean_ctou_long

mean_ctou <- data.frame(
  chromosome = rep(c(1:22, "X", "Y", "MT"), times = 75),
  group = rep(c("acute", "control","recovery"), each = 25),
  mean = mean_ctou_long$value)

mean_ctou <- mean_ctou[1:75,]

# Plot the bar graph
ggplot(mean_ctou,
       aes(x = factor(chromosome, levels = mixedsort(unique(chromosome))),
           y = mean, fill = factor(group))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "Mean count of C-to-U edits by Chromosome and Group",
       x = "Chromosome",
       y = "Mean count", 
       fill = "Group") +
  theme(axis.text.x = element_text(hjust = 1)) +
  theme_dark()

# A to I
mean_atoi_recovery <- tdf %>% 
  filter(Time_point == "recovery") %>% 
  select(starts_with("AtoI_Ch")) %>% 
  select(-ends_with("matches"))

mean_atoi_controls <- tdf %>% 
  filter(Time_point == "control") %>% 
  select(starts_with("AtoI_Ch")) %>% 
  select(-ends_with("matches"))

mean_atoi_acute <- tdf %>% 
  filter(Time_point == "acute") %>% 
  select(starts_with("AtoI_Ch")) %>% 
  select(-ends_with("matches"))

mean_a_rec <-colMeans(mean_atoi_recovery)
mean_a_con <-colMeans(mean_atoi_controls)
mean_a_acu <-colMeans(mean_atoi_acute)

mean_atoi_acute <- as.data.frame(mean_a_acu)
mean_atoi_controls <-as.data.frame(mean_a_con)
mean_atoi_recovery <-as.data.frame(mean_a_rec)

mean_atoi <- cbind(mean_atoi_acute, mean_atoi_controls, mean_atoi_recovery)
mean_atoi_long <- melt(mean_atoi)
mean_atoi_long

mean_atoi <- data.frame(
  chromosome = rep(c(1:22, "X", "Y", "MT"), times = 75),
  group = rep(c("acute", "control","recovery"), each = 25),
  mean = mean_atoi_long$value)

mean_atoi <- mean_atoi[1:75,]

mean_atoi

# Plot the bar graph
ggplot(mean_atoi,
       aes(x = factor(chromosome, levels = mixedsort(unique(chromosome))),
           y = mean, fill = factor(group))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +  # Add black outline
  labs(title = "Mean count of A-to-I edits by Chromosome and Group",
       x = "Chromosome",
       y = "Mean count", 
       fill = "Group") +
  theme(axis.text.x = element_text(hjust = 1)) +
  theme_dark()


# Create a box plot with facets
ggplot(chr_ctou_long, aes(x = as.factor(Time_point), y = Count, fill = as.factor(Time_point))) +
  geom_boxplot() +
  facet_wrap(~Chromosome, scales = "free_y", ncol = 5) +
  labs(title = "Boxplot of C to U counts by Chromosome and Time_point",
       x = "Time_point",
       y = "Count of C to U variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reorder the Chromosome levels
chr_ctou_long$Chromosome <- factor(chr_ctou_long$Chromosome, 
                                   levels = c("CtoU_chr1", "CtoU_chr2", "CtoU_chr3", "CtoU_chr4", "CtoU_chr5", 
                                              "CtoU_chr6", "CtoU_chr7", "CtoU_chr8", "CtoU_chr9", "CtoU_chr10", 
                                              "CtoU_chr11", "CtoU_chr12", "CtoU_chr13", "CtoU_chr14", "CtoU_chr15", 
                                              "CtoU_chr16", "CtoU_chr17", "CtoU_chr18", "CtoU_chr19", "CtoU_chr20", 
                                              "CtoU_chr21", "CtoU_chr22", "CtoU_chrX", "CtoU_chrY", "CtoU_chrMT"))

# Custom labels for each Chromosome
facet_labels <- c("CtoU_chr1" = "Chr 1", "CtoU_chr2" = "Chr 2", "CtoU_chr3" = "Chr 3", 
                  "CtoU_chr4" = "Chr 4", "CtoU_chr5" = "Chr 5", "CtoU_chr6" = "Chr 6",
                  "CtoU_chr7" = "Chr 7", "CtoU_chr8" = "Chr 8", "CtoU_chr9" = "Chr 9",
                  "CtoU_chr10" = "Chr 10", "CtoU_chr11" = "Chr 11", "CtoU_chr12" = "Chr 12",
                  "CtoU_chr13" = "Chr 13", "CtoU_chr14" = "Chr 14", "CtoU_chr15" = "Chr 15",
                  "CtoU_chr16" = "Chr 16", "CtoU_chr17" = "Chr 17", "CtoU_chr18" = "Chr 18",
                  "CtoU_chr19" = "Chr 19", "CtoU_chr20" = "Chr 20", "CtoU_chr21" = "Chr 21",
                  "CtoU_chr22" = "Chr 22", "CtoU_chrX" = "Chr X", "CtoU_chrY" = "Chr Y", "CtoU_chrMT" = "Chr MT")

# Create the plot
ggplot(chr_ctou_long, aes(x = as.factor(Time_point), y = Count, fill = as.factor(Time_point))) +
  geom_boxplot() +
  facet_wrap(~Chromosome, scales = "free_y", ncol = 5, 
             labeller = as_labeller(facet_labels)) +  # Apply custom titles
  labs(title = "Boxplot of C-to-U Counts by Chromosome and Group",
       x = "Time Point",
       y = "Count of C-to-U Variants",
       fill = "Group") +  # Change legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_dark() +  # Apply dark theme
  theme(legend.position = "bottom")


ggplot(chr_atoi_long, aes(x = as.factor(Time_point), y = Count, fill = as.factor(Time_point))) +
  geom_boxplot() +
  facet_wrap(~Chromosome, scales = "free_y", ncol = 5) +
  labs(title = "Boxplot of A to I counts by Chromosome and Time_point",
       x = "Time_point",
       y = "Count of A to I variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reorder the Chromosome levels
chr_atoi_long$Chromosome <- factor(chr_atoi_long$Chromosome, 
                                   levels = c("AtoI_Chr1", "AtoI_Chr2", "AtoI_Chr3", "AtoI_Chr4", "AtoI_Chr5", 
                                              "AtoI_Chr6", "AtoI_Chr7", "AtoI_Chr8", "AtoI_Chr9", "AtoI_Chr10", 
                                              "AtoI_Chr11", "AtoI_Chr12", "AtoI_Chr13", "AtoI_Chr14", "AtoI_Chr15", 
                                              "AtoI_Chr16", "AtoI_Chr17", "AtoI_Chr18", "AtoI_Chr19", "AtoI_Chr20", 
                                              "AtoI_Chr21", "AtoI_Chr22", "AtoI_ChrX", "AtoI_ChrY", "AtoI_ChrMT"))

# Custom labels for each Chromosome
facet_labels <- c("AtoI_Chr1" = "Chr 1", "AtoI_Chr2" = "Chr 2", "AtoI_Chr3" = "Chr 3", 
                  "AtoI_Chr4" = "Chr 4", "AtoI_Chr5" = "Chr 5", "AtoI_Chr6" = "Chr 6",
                  "AtoI_Chr7" = "Chr 7", "AtoI_Chr8" = "Chr 8", "AtoI_Chr9" = "Chr 9",
                  "AtoI_Chr10" = "Chr 10", "AtoI_Chr11" = "Chr 11", "AtoI_Chr12" = "Chr 12",
                  "AtoI_Chr13" = "Chr 13", "AtoI_Chr14" = "Chr 14", "AtoI_Chr15" = "Chr 15",
                  "AtoI_Chr16" = "Chr 16", "AtoI_Chr17" = "Chr 17", "AtoI_Chr18" = "Chr 18",
                  "AtoI_Chr19" = "Chr 19", "AtoI_Chr20" = "Chr 20", "AtoI_Chr21" = "Chr 21",
                  "AtoI_Chr22" = "Chr 22", "AtoI_ChrX" = "Chr X", "AtoI_ChrY" = "Chr Y", "AtoI_ChrMT" = "Chr MT")

# Create the plot
ggplot(chr_atoi_long, aes(x = as.factor(Time_point), y = Count, fill = as.factor(Time_point))) +
  geom_boxplot() +
  facet_wrap(~Chromosome, scales = "free_y", ncol = 5, 
             labeller = as_labeller(facet_labels)) +  # Apply custom titles
  labs(title = "Boxplot of A-to-I Counts by Chromosome and Time Point",
       x = "Time Point",
       y = "Count of A-to-I Variants",
       fill = "Time Point") +  # Change legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_dark()+  # Apply dark theme
  theme(legend.position = "bottom")



# Common C to U variant count on each chromosome in the groups:
# The following block plots the diveristy of C to U variants on
# each chromosome per group
ctou_var_counts <- ctou_var_counts %>%
  rownames_to_column(var = "Chromosome")

# Reshape the data for ggplot
ctou_var_counts_long <- ctou_var_counts %>%
  pivot_longer(cols = -Chromosome, names_to = "Group", values_to = "Count") %>%
  arrange(factor(Chromosome, levels = mixedsort(unique(ctou_var_counts$Chromosome)))) %>% 
  arrange(factor(Group, levels = c("acute", "recovery", "controls")))

# Plot the bar graph
ggplot(ctou_var_counts_long,
       aes(x = factor(Chromosome, levels = mixedsort(unique(ctou_var_counts$Chromosome))),
           y = Count, fill = factor(Group, levels = c("acute", "recovery", "controls")))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Diversity of C to U Variants on Chromosomes",
       x = "Chromosome",
       y = "Count", 
       fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ctou_var_counts <- ctou_var_counts %>% 
  mutate(difacutecontrol = ((acute - controls) / acute) * 100) %>% 
  mutate(difacuterec = ((acute - recovery) / acute) * 100) %>% 
  mutate(difacuterecabs = acute - recovery) %>% 
  mutate(difacutecontrolabs = acute - controls)
ctou_var_counts

max(ctou_var_counts$difacutecontrol)
max(ctou_var_counts$difacuterec)

table(acute_var_unique$chr)
table(acute_varNotInControls$chr)
table(acute_varNotInRec$chr)

#####DID NOT TOUCH THE BELOW 2024 ########
#------------------------#
#Common Variant analysis
hosp_com_var <- read_csv("hosp_com_var.csv")
controls_com_var <- read_csv("controls_com_var.csv")
recovery_com_var <- read_csv("recovery_com_var.csv")

# Renaming a variable and cleaning data from tailing "," from it's records
hosp_com_var <- hosp_com_var %>% 
  rename(found_in = "Found in")

# Remove the last comma from the "found_in" column
hosp_com_var$found_in <- substring(hosp_com_var$found_in, 1, 
                                   nchar(hosp_com_var$found_in) - 1)

controls_com_var <- controls_com_var %>% 
  rename(found_in = "Found in")

controls_com_var$found_in <- substring(controls_com_var$found_in, 1, 
                                       nchar(controls_com_var$found_in) - 1)

recovery_com_var <- recovery_com_var %>% 
  rename(found_in = "Found in")

recovery_com_var$found_in <- substring(recovery_com_var$found_in, 1, 
                                       nchar(recovery_com_var$found_in) - 1)

# Finding variants found only in hospital but not in other groups and show 
# their % occurance.

acute_varNotInControls <- anti_join(hosp_com_var, controls_com_var, by = "unique_var_id")
controls_var_NotInAcute <- anti_join(controls_com_var, hosp_com_var, by = "unique_var_id")
acute_varNotInRec <- anti_join(hosp_com_var, recovery_com_var, by = "unique_var_id")
acute_var_unique <- anti_join(acute_varNotInControls, recovery_com_var, by =  "unique_var_id")

output_file1 <- "acute_varNotInControls.csv"
output_file2 <- "controls_var_NotInAcute.csv"
output_file3 <- "acute_varNotInRec.csv"
output_file4 <- "acute_var_unique.csv"

# Export tables as tab-delimited text files.
write.table(acute_varNotInControls, file = output_file1, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(controls_var_NotInAcute, file = output_file2, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(acute_varNotInRec, file = output_file3, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(acute_var_unique, file = output_file4, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

# Split the comma-delimited "found_in" values into separate rows
hosp_com_var <- hosp_com_var %>%
  separate_rows(found_in, sep = ",")

controls_com_var <- controls_com_var %>%
  separate_rows(found_in, sep = ",")

recovery_com_var <- recovery_com_var %>%
  separate_rows(found_in, sep = ",")

# Merging all variant tables:
all_variants <- bind_rows(hosp_com_var, controls_com_var, recovery_com_var)

all_variants <- all_variants %>% 
  select(-Time_point.y)


# Filter paired, hospital admission and not recovery
paired_variants <- all_variants %>% 
  select(unique_var_id, ref, alt, matches, perc_hosp, Time_point.x, 
         Paired_single, Participant_ID) %>% 
  filter(Time_point.x == "Hospital admission" | Time_point.x == "After 1 month recovery", 
         Paired_single == "paired")

output_file <- "paired_variants.txt"

# Export the dataframe as tab-delimited text
write.table(paired_variants, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


# Perform a left join to get the corresponding data from the main df
hosp_com_var <- hosp_com_var %>%
  left_join(thesis_datav2, by = c("found_in" = "Run"))

controls_com_var <- controls_com_var %>%
  left_join(thesis_datav2, by = c("found_in" = "Run"))

recovery_com_var <- recovery_com_var %>%
  left_join(thesis_datav2, by = c("found_in" = "Run"))


# Creating a df with unique variant ids and their respective % occurance in groups

controls_per <- controls_com_var %>% 
  select(unique_var_id, perc_controls)

hosp_per <- hosp_com_var %>% 
  select(unique_var_id, perc_hosp)

recovery_per <- recovery_com_var %>% 
  select(unique_var_id, perc_rec)

combined_df <- full_join(hosp_per, controls_per, by = "unique_var_id") %>%
  full_join(., recovery_per, by = "unique_var_id")

#Change all NA values in df to "0"
combined_df <- combined_df %>%
  mutate_all(~ifelse(is.na(.), 0, .))

combined_df <- combined_df %>%
  mutate(diff_hosp_controls = abs(perc_hosp - perc_controls),
         diff_hosp_rec = abs(perc_hosp - perc_rec),
         diff_controls_rec = abs(perc_controls - perc_rec))

unique_var_id_with_biggest_diff <- combined_df %>%
  arrange(desc(pmax(diff_hosp_controls, diff_hosp_rec, diff_controls_rec))) %>%
  select(unique_var_id, diff_hosp_controls, diff_hosp_rec, diff_controls_rec)


#####DID NOT TOUCH THE above 2024 ########
#---------------Gene expression analysis--------------------#
#Filter-out 0 counts

gene_counts <- gene_counts[which(rowSums(gene_counts[, -1]) > 0), ]


# getting gene names
BiocManager::install("org.Hs.eg.db")
BiocManager::install("apeglm")
library("org.Hs.eg.db")

# convert gene_counts to regular dataframe
gene_counts <- as.data.frame(gene_counts)
str(gene_counts)

# check if the genes are all different
any(duplicated(gene_counts$Gene))

# Exclude the first column (gene IDs)
gene_counts_without_gene_id <- gene_counts[, -1]

str(gene_counts_without_gene_id)

# setting the row names
row.names(gene_counts_without_gene_id) <- gene_counts$Gene

#maybe read: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

# Creating a list of factors indicating the group of each sample
condition <- factor(c(thesis_datav4$Time_point))
coldata <- data.frame(row.names = colnames(gene_counts_without_gene_id), condition)
coldata

# Construct DESEQDataSet Object
dds2 <- DESeqDataSetFromMatrix(countData = gene_counts_without_gene_id,
                               colData = coldata,
                               design = ~condition)

model.matrix(~0+~condition, data=colData(dds2))

#no. of genes in dds2 object
nrow(dds2)

#check number of genes with expression of less than 10
genes_less_10_reads <- dds2[rowSums(counts(dds2)) < 10, ]
nrow(genes_less_10_reads)

# estimate size factors
dds2_sizefactors <- estimateSizeFactors(dds2)
dds2_sizefactors$sizeFactor

plot(dds2_sizefactors$sizeFactor)

# Running the DEseq function
DEdds2 <- DESeq(dds2)

# Variance stabilizing transformation 
# to stabilize the variance across the mean. 
# In RNA-seq data analysis, the variance of counts tends to increase 
# with the mean.
vsdata2 <- vst(dds2, blind = FALSE)

#PCA plot
plotPCA(vsdata2, intgroup = "condition")

# Dispersion Estimates plot of DEdds2
plotDispEsts(DEdds2)

# results
res2 <- results(DEdds2, contrast = c("condition", "acute", "control"))

write.csv(res2, file="res_acutevscontrol.csv", row.names = TRUE)

res2_rec <- results(DEdds2, contrast = c("condition", "acute", "recovery"))
write.csv(res2_rec, file="res_acutevsrecovery.csv", row.names = TRUE)

res2_recvcon <- results(DEdds2, contrast = c("condition", "recovery", "control"))
write.csv(res2_recvcon, file="res_recvcon.csv", row.names = TRUE)



# significant results
sig_acutevscontrol <-na.omit(res2)
sig_acutevscontrol <- sig_acutevscontrol[sig_acutevscontrol$padj < 0.05,]
sig_acutevsrecovery <-na.omit(res2_rec)
sig_acutevsrecovery <- sig_acutevsrecovery[sig_acutevsrecovery$padj < 0.05,]

# Convert Ensemblb IDs to gene name and place gene_name as first column
# acutevscontrol
sig_acutevscontrol.df <- as.data.frame(sig_acutevscontrol)

sig_acutevscontrol.df <- data.frame(gene_name = mapIds(org.Hs.eg.db, 
                                                       keys = rownames(sig_acutevscontrol.df),
                                                       keytype = "ENSEMBL",
                                                       column = "SYMBOL"),
                                    sig_acutevscontrol.df)

sig_acutevscontrol.df
# acutevsrecovery
sig_acutevsrecovery.df <- as.data.frame(sig_acutevsrecovery)

sig_acutevsrecovery.df <- data.frame(gene_name = mapIds(org.Hs.eg.db, 
                                                        keys = rownames(sig_acutevsrecovery.df),
                                                        keytype = "ENSEMBL",
                                                        column = "SYMBOL"),
                                     sig_acutevsrecovery.df)

# gene names output to file
writeLines(unlist(gene_name), con = "gene_names.txt")

# target gene list
target_genes = c("APOBEC4", "APOBEC3A", "APOBEC3B", "APOBEC3B-AS1", 
                 "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", 
                 "APOBEC3H", "APOBEC2", "ADAR", "ADARB2", "ADARB1")

# Subsetting results of RNA editing genes
subset_res <- sig_acutevscontrol.df[sig_acutevscontrol.df$gene_name %in% target_genes, ]
subset_res_recovery <- sig_acutevsrecovery.df[sig_acutevsrecovery.df$gene_name %in% target_genes, ]

# Bar plots
ggplot(subset_res, aes(x = gene_name, y = log2FoldChange, fill = gene_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "                Significant Log2fold changes
                Acute vs Control",
       x = "Gene Name",
       y = "log2FoldChange") +
  theme_minimal()

ggplot(subset_res_recovery, aes(x = gene_name, y = log2FoldChange, fill = gene_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "                Significant Log2fold changes
                Acute vs Recovery",
       x = "Gene Name",
       y = "log2FoldChange") +
  theme_minimal()

# combining plots
subset_res$Comparison <- "Acute vs Control"
subset_res_recovery$Comparison <- "Acute vs Recovery"

combined_res <- rbind(subset_res, subset_res_recovery)

# Plot the combined data with facet_wrap
ggplot(combined_res, aes(x = gene_name, y = log2FoldChange, fill = gene_name)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_wrap(~ Comparison, scales = "free_x") +  # Side-by-side plots based on Comparison
  labs(title = "Up-regulated RNA editing genes in acute samples",
       x = "Gene Name",
       y = "log2FoldChange") +
  theme_dark() +
  theme(legend.position = "right")

# Extracting gene expression level of ENSG00000100298 (APOBEC3H)
# in the three groups 

normalized_counts <- counts(DEdds2, normalized = TRUE)
apobec3h <- normalized_counts["ENSG00000100298", ]

#---------------------#
# Extract the normalised counts of all the RNA editing genes.
normalized_counts <- t(normalized_counts)
normalized_counts <- as.data.frame(normalized_counts)
rna_edit_genes_counts <- normalized_counts %>% 
  select(ENSG00000173627, ENSG00000128383, ENSG00000179750, 
         ENSG00000244509, ENSG00000100298, ENSG00000124701, 
         ENSG00000160710, ENSG00000185736, ENSG00000197381)

# merging RNA editing gene counts with the rest of the dataframe
rna_edit_genes_counts <- cbind(tdf, rna_edit_genes_counts)

# Select only numeric variables
numeric_rna_edit_genes <- rna_edit_genes_counts[sapply(rna_edit_genes_counts, is.numeric)]

rna_edit_counts <- tdf %>%
  select(starts_with("AtoI_Ch"), starts_with("CtoU_"), Time_point) %>% 
  select(-ends_with("matches"))
rna_edit_counts$Time_point<-factor(rna_edit_counts$Time_point)
str(rna_edit_counts)

# Perform PCA
#counts per chromosome a to i all and c to u all
pca_per_chrom <- prcomp(rna_edit_counts[,1:50], center = TRUE, scale. = TRUE)

#counts per chromosome a to i all
pca_per_chr_atoi <- prcomp(rna_edit_counts[,1:25], center = TRUE, scale. = TRUE)

#counts per chromosome c to u all
pca_per_chr_ctou <- prcomp(rna_edit_counts[,26:50], center = TRUE, scale. = TRUE)

summary(pca_per_chrom)

# Checking the loading scores to determine which
# variables have the largest effect on where the
# samples are plotted in the PCA plot
loadings_pc1 <- round(pca_per_chrom$rotation[, 1], 2)
loadings_pc1
loadings_pc2 <- round(pca_per_chrom$rotation[, 2], 2)
loadings_pc2
# Variables with large negative loading scores will push samples to the
# left of the graph. Variables with large positive scores will push
# samples to the right.


# Create a dataframe for plotting
plot_df <- data.frame(PC1 = pca_per_chrom$x[,1], PC2 = pca_per_chrom$x[,2], Group = tdf$Time_point)

# Create the plot using ggplot2
# Amend titles for each plot

ggplot(plot_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point() +
  theme_minimal() +
  labs(x = "PC1: 71% of variance", y = "PC2: 9% of variance", color = "Group", title = 
      "PCA Plot of first two principal components:
       RNA edits per chromosome")

# biplot
library(devtools)
library(ggbiplot)

loadings <- pca_per_chrom$rotation

# Subset the loadings to include only the vectors of interest
selected_vectors <- c("CtoU_chr22", "CtoU_chr1", "AtoI_Chr22", "AtoI_Chr2", "AtoI_Chr8", "AtoI_Chr21")

# Ensure the PCA object is correctly defined
pca_subset <- pca_per_chrom
pca_subset$rotation <- loadings[selected_vectors,]

g <- ggbiplot(pca_subset, obs.scale = 1, var.scale = 1,
              groups = rna_edit_counts$Time_point, ellipse = TRUE,
              circle = TRUE, var.axes = FALSE)

# Customize the plot
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top') +
  
  theme_dark()
g
# PCA of cytokines, A to I, C to U and RNA editing gene expression
# total of 15 covariates
cyto_edits_rnagenes<-numeric_rna_edit_genes %>% 
  select("LPS_IL_10", "LPS_IL_1b", "LPS_IL_6",
         "LPS_TNF_alpha","AtoI", "CtoU", starts_with("ENSG"))

cyto_edits_rnagenes<-numeric_rna_edit_genes %>% 
  select(-Participant)


# removing records with no data  
cyto_edits_rnagenes <- cyto_edits_rnagenes %>% na.omit()

# adding the group variable
groups_for_pca <- tdf$Time_point[!is.na(rowSums(cyto_edits_rnagenes))]

str(cyto_edits_rnagenes)
pca_cyto_edits_rnagenes <- prcomp(cyto_edits_rnagenes, center = TRUE, scale. = TRUE)
summary(pca_cyto_edits_rnagenes)

# Compute correlation matrix

cor_matrix_rna_edit <- cor(numeric_rna_edit_genes, method = "spearman", use = "pairwise.complete.obs")
str(cor_matrix_rna_edit)
write.table(cor_matrix_rna_edit, file = "cor_matrix_rna_edit_strict_acute.csv", sep = "\t", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)
str(numeric_rna_edit_genes)

#checking significance level of certain correlations:
rcorr(numeric_rna_edit_genes$LPS_IL_10, numeric_rna_edit_genes$AtoI_ChrMT,
      type = "spearman")

ggplot(numeric_rna_edit_genes, aes(x = LPS_IL_10, y = AtoI_ChrMT)) +
  geom_point(aes(colour = factor(tdf$Time_point))) +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  theme_minimal() +
  labs(x = "LPS_IL_10", y = "AtoI_ChrMT", 
       title = "LPS_IL_10 vs AtoI_ChrMT",
       colour = "Group")

#RNA editing genes vs LPS_IL_10

ggplot(numeric_rna_edit_genes, aes(x = ENSG00000160710, y = AtoI_Chr22)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black", alpha = 0.7) +  # Black outline for points with group fill
  geom_smooth(aes(colour = factor(tdf$Time_point)), method = "lm", se = FALSE) +  # Trendlines for each group with confidence interval
  stat_cor(aes(color = factor(tdf$Time_point)), method = "spearman", label.x.npc = "centre", label.y.npc = "bottom") +  # Spearman correlation for each group
  theme_dark() +
  labs(x = "ADAR", y = "A-to-I Chr 22 RNA edit count", 
       title = "ADAR normalised experssion vs A-to-I edit count Chromosome 22",
       fill = "Group", colour = "Group") +
  theme(legend.position = "bottom")

ggplot(numeric_rna_edit_genes, aes(x = ENSG00000244509, y = CtoU_chr3)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black", alpha = 0.7) +  # Black outline for points with group fill
  geom_smooth(aes(colour = factor(tdf$Time_point)), method = "lm", se = FALSE) +  # Trendlines for each group with confidence interval
  stat_cor(aes(color = factor(tdf$Time_point)), method = "spearman", label.x.npc = "left", label.y.npc = "top") +  # Spearman correlation for each group
  theme_dark() +
  labs(x = "APOBEC3C normalised expression", y = "C-to-U Chr 3 RNA edit count", 
       title = "APOBEC3C gene experssion vs C-to-U edit count Chromosome 3",
       fill = "Group", colour = "Group") +
  theme(legend.position = "bottom")

# Il-10 secretion across groups




ggplot(numeric_rna_edit_genes, aes(x = LPS_IL_1b, y = AtoI_Ch16matches)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black", alpha = 0.7) +  # Black outline for points with group fill
  geom_smooth(aes(colour = factor(tdf$Time_point)), method = "lm", se = TRUE) +  # Trendlines for each group with confidence interval
  stat_cor(aes(color = factor(tdf$Time_point)), method = "spearman", label.x.npc = "right", label.y.npc = "top") +  # Spearman correlation for each group
  theme_minimal() +
  labs(x = "LPS_IL_10", y = "ENSG00000197381", 
       title = "LPS_IL_10 vs ENSG00000197381",
       fill = "Group", colour = "Group") +
  theme(legend.position = "bottom")

ggplot(numeric_rna_edit_genes, aes(x = ENSG00000197381, y = LPS_IL_10)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black") +  # Black outline for points with group fill
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = 500, label.y = 4000) +  # Add Spearman's correlation
  theme_dark() +
  labs(x = "ADARB1 normalised counts", y = "Il-10 secretion (ng/L)", 
       title = "ADARB1 expression vs Il-10 secretion",
       fill = "Group") +
  scale_x_continuous(limits = c(0, 1000)) +
  theme(legend.position = "none")

ggplot(numeric_rna_edit_genes, aes(x = ENSG00000128383, y = LPS_IL_10)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black") +  # Black outline for points with group fill
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = 7500, label.y = 2000) +  # Add Spearman's correlation
  theme_dark() +
  labs(x = "APOBEC3A normalised counts", y = "Il-10 secretion (ng/L)", 
       title = "APOBEC3A expression vs Il-10 secretion",
       fill = "Group") +
  scale_x_continuous(limits = c(0, 15000)) +
  theme(legend.position = "none")

ggplot(numeric_rna_edit_genes, aes(x = ENSG00000244509, y = LPS_IL_10)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black") +  # Black outline for points with group fill
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = 3500, label.y = 2000) +  # Add Spearman's correlation
  theme_dark() +
  labs(x = "APOBEC3C normalised counts", y = "Il-10 secretion (ng/L)", 
       title = "APOBEC3C expression vs Il-10 secretion",
       fill = "Group") +
  scale_x_continuous(limits = c(2000, 5000)) +
  theme(legend.position = "none")

ggplot(numeric_rna_edit_genes, aes(x = ENSG00000197381, y = LPS_IL_1b)) +
  geom_point(aes(fill = factor(tdf$Time_point)), shape = 21, size = 3, stroke = 0.5, color = "black") +  # Black outline for points with group fill
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = 600, label.y = 25000) +  # Add Spearman's correlation
  theme_dark() +
  labs(x = "ADARB1 normalised counts", y = "Il-1b secretion (ng/L)", 
       title = "ADARB1 expression vs Il-1b secretion",
       fill = "Group") +
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0, 1000))+
  scale_y_continuous(limits = c(0, 40000))


ggplot(numeric_rna_edit_genes, aes(x = ENSG00000128383, y = LPS_IL_10)) +
  geom_point(aes(colour = factor(tdf$Time_point))) +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  theme_minimal() +
  labs(x = "APOBEC3A normalised counts", y = "Il-10 secretion", 
       title = "APOBEC3A expression vs Il-10 secretion",
       colour = "Group") +
  scale_y_continuous(limits = c(0, 6000))


ggplot(numeric_rna_edit_genes, aes(x = ENSG00000244509, y = LPS_IL_10)) +
  geom_point(aes(colour = factor(tdf$Time_point))) +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  theme_minimal() +
  labs(x = "APOBEC3C normalised counts", y = "Il-10 secretion", 
       title = "APOBEC3C expression vs Il-10 secretion",
       colour = "Group")

cor(numeric_rna_edit_genes$LPS_IL_10, numeric_rna_edit_genes$CtoU, method = "spearman")
cor(numeric_rna_edit_genes$LPS_IL_10, numeric_rna_edit_genes$ENSG00000100298, method = "spearman")
cor(numeric_rna_edit_genes$LPS_IL_10, numeric_rna_edit_genes$ENSG00000128383, method = "spearman")

# Il-1b
ggplot(numeric_rna_edit_genes, aes(x = ENSG00000197381, y = LPS_IL_1b)) +
  geom_point(aes(colour = factor(tdf$Time_point))) +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  theme_minimal() +
  labs(x = "ADARB1 normalised counts", y = "Il-1b secretion", 
       title = "ADARB1 expression vs Il-1b secretion",
       colour = "Group")

# TNF-alpha in acute only

numeric_rna_edit_genes$Time_point <- tdf$Time_point
acute_numeric_rna_edit_genes <- numeric_rna_edit_genes %>% 
  filter(Time_point == "acute")

ggplot(acute_numeric_rna_edit_genes, aes(x = ENSG00000179750, y = LPS_TNF_alpha)) +
  geom_point(aes(colour = factor(acute_numeric_rna_edit_genes$Time_point))) +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  theme_minimal() +
  labs(x = "APOBEC3B normalised counts in acute", y = "TNF-alpha secretion", 
       title = "APOBEC3B expression in acute vs TNF-alpha secretion",
       colour = "Group")


#PCA with cytokine data and all the RNA editing genes
numeric_rna_edit_genes <- numeric_rna_edit_genes %>% 
  select(-Participant)
numeric_rna_edit_genes <- numeric_rna_edit_genes %>% 
  select(-starts_with("AtoI"))
numeric_rna_edit_genes <- numeric_rna_edit_genes %>% 
  select(-starts_with("CtoU"))
numeric_rna_edit_genes <- numeric_rna_edit_genes %>% 
  select(-age, -TotalSnpIndels, -uncommon_snps, -uncommon_snps_qc)

numeric_rna_edit_genes$Time_point <- tdf$Time_point

numeric_rna_edit_genes <- na.omit(numeric_rna_edit_genes)
numeric_rna_edit_genes$Time_point <-as.factor(numeric_rna_edit_genes$Time_point)

pca_cyto_genes <- prcomp(numeric_rna_edit_genes[,1:23], center = TRUE, scale. = TRUE)
summary(pca_cyto_genes)

autoplot(pca_cyto_genes, data = numeric_rna_edit_genes, colour = 'Time_point')

pca_cyto_genes$rotation

# Subset the loadings to include only the vectors I'm interested in

loadings2 <- pca_cyto_genes$rotation
selected_vectors2 <- c("LPS_IL_12p70","IFN_gamma", "ENSG00000128383",
                       "ENSG00000160710", "ENSG00000244509")
pca_cyto_genes$rotation <- loadings2[selected_vectors2, ]

g2 <- ggbiplot(pca_cyto_genes, obs.scale = 1, var.scale = 1,
              groups = numeric_rna_edit_genes$Time_point, ellipse = TRUE,
              circle = TRUE)
g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal',
               legend.position = 'top')
g2

#getting normalised counts of all significant rna editing genes in acute vs control
rna_edit_genes_apobec <- normalized_counts[c("ENSG00000128383", "ENSG00000244509", "ENSG00000100298"), ]

# Transpose data for plotting
rna_edit_genes_apobec <- t(rna_edit_genes_apobec)

rna_edit_genes_apobec <- as.data.frame(rna_edit_genes_apobec)
rna_edit_genes_apobec$condition <-coldata$condition
rna_edit_genes_apobec$id <- thesis_datav4$Participant_ID
rna_edit_genes_apobec$age <- thesis_datav4$age
rna_edit_genes_apobec$ctou <- thesis_datav4$CtoU
str(rna_edit_genes_apobec)

selected_columns_ctou <- c("ENSG00000128383", "ENSG00000244509", "ENSG00000100298", "ctou")
selected_columns_age <- c("ENSG00000128383", "ENSG00000244509", "ENSG00000100298", "age")

# Spearman correlation matrix between counts of rna_editing genes and c to u changes
# and age.
cor_matrix_ctou <- rcorr(as.matrix(rna_edit_genes_apobec[, selected_columns_ctou]), type = "spearman")
cor_matrix_age <- rcorr(as.matrix(rna_edit_genes_apobec[, selected_columns_age]), type = "spearman")

# Extract the correlation matrix
cor_matrix_ctou$r
cor_matrix_age$r

# Filter data for the acute condition
acute_exp_data_apobec <- subset(rna_edit_genes_apobec, condition == "acute")

# Spearman correlation matrix between counts of rna editing genes and 
# c to u and age only in acute group
acute_cor_matrix_ctou <- rcorr(as.matrix(acute_exp_data_apobec[, selected_columns_ctou]), type = "spearman")
acute_cor_matrix_age <- rcorr(as.matrix(acute_exp_data_apobec[, selected_columns_age]), type = "spearman")

# Extract the correlation matrix
acute_cor_matrix_ctou$r
acute_cor_matrix_age$r

# Comparing the gex of ENSG00000100298 (apobec3h) in runs in which a variant in 
# RELL1 (RELT-like protein 1) were detected.

# Define the sample names for acute condition
rell1_var_samples <- c("SRR12926256", "SRR12926262", "SRR12926269", "SRR12926284",
                       "SRR12926292", "SRR12926339", "SRR12926347", "SRR12926351")

rell1_var_samples_paired <- c("SRR12926256", "SRR12926262", "SRR12926269",
                              "SRR12926292", "SRR12926339", "SRR12926347")


rell1_var_samples_recovery <- c("SRR12926257", "SRR12926263", "SRR12926270", "SRR12926293",
                                "SRR12926340", "SRR12926348")

# Extract the normalised counts for ENSG00000100298 in acute condition
rell1_counts <- normalized_counts["ENSG00000100298", rell1_var_samples]
rell1_counts_paired <- normalized_counts["ENSG00000100298", rell1_var_samples_paired]
rell1_counts_rec <-normalized_counts["ENSG00000100298", rell1_var_samples_recovery]


# Extract the normalised counts for ENSG00000100298 in all other samples
no_rell1_counts <- normalized_counts["ENSG00000100298", !colnames(normalized_counts) %in% rell1_var_samples]



boxplot(log2(rell1_counts), log2(no_rell1_counts))
boxplot(log2(rell1_counts), log2(rell1_counts_rec))
boxplot(log2(rell1_counts_paired), log2(rell1_counts_rec))


# Create a data frame for each set of values
df_rell1 <- data.frame(
  Condition = rep("Acute", length(rell1_counts)),
  Log2Counts = log2(rell1_counts)
)

df_rell1_paired <- data.frame(
  Condition = rep("Acute", length(rell1_counts_paired)),
  Log2Counts = rell1_counts_paired
)

df_no_rell1 <- data.frame(
  Condition = rep("Other", length(no_rell1_counts)),
  Log2Counts = log2(no_rell1_counts)
)

df_rell1_counts_rec <- data.frame(
  Condition = rep("recovery", length(rell1_counts_rec)),
  Log2Counts = rell1_counts_rec
)

# Combine the data frames
combined_df <- rbind(df_rell1, df_no_rell1)
combined_df <- rbind(df_rell1, df_rell1_counts_rec)
combined_df <- rbind(df_rell1_paired, df_rell1_counts_rec)

combined_df

ggplot(combined_df, aes(x = Condition, y = Log2Counts, fill = Condition)) +
  geom_boxplot() +
  labs(title = "apobec3h gene expression in individuals with RELL1 variant",
       x = "Condition",
       y = "Normalised counts of Apobec3h")

combined_df$Run <- rownames(combined_df)

# Merge the data frames based on the "Run" column
merged_df <- merge(combined_df, thesis_datav4, by = "Run", all.x = TRUE)

ggplot(merged_df, aes(x = Condition, y = CtoU)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = Participant_ID), position = position_dodge(width = 0.75), vjust = -0.5) +  # Add labels
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.5) +
  labs(title = "Comparison of condition vs count of C to U variants",
       x = "Condition",
       y = "C to U variants")

ctou_in_paired_samples <- merged_df %>%
  select(Condition, CtoU) %>% 
  group_by(Condition) %>% 
  summarise(mean_ctou = mean(CtoU))

# displaying mean counts
ctou_in_paired_samples

ctou_in_paired_samples_age <- merged_df %>%
  select(Condition, CtoU, age) %>% 
  group_by(Condition) %>% 
  summarise(mean_age = mean(age))
ctou_in_paired_samples_age

# mean ages of acute, recovery and controls
ages<-thesis_datav4 %>% 
  select(Time_point, age) %>% 
  group_by(Time_point) %>% 
  summarise(mean_age = mean(age))
ages

# The 2:44175023 G>T variant in PPM1B 
PPM1B_var_samples <- c(
  "SRR12926207",
  "SRR12926212",
  "SRR12926224",
  "SRR12926264",
  "SRR12926266",
  "SRR12926271",
  "SRR12926278",
  "SRR12926349",
  "SRR12926352",
  "SRR12926354",
  "SRR12926366"
)

PPM1B_var_samples_paired<- c(
  "SRR12926212",
  "SRR12926264",
  "SRR12926271",
  "SRR12926278",
  "SRR12926349",
  "SRR12926352",
  "SRR12926366"
)

PPM1B_var_samples_paired_rec<- c(
  "SRR12926213",
  "SRR12926265",
  "SRR12926272",
  "SRR12926279",
  "SRR12926350",
  "SRR12926353",
  "SRR12926367"
)


# Extract the normalised counts for apobec3h in acute and rec condition
PPM1B_counts <- normalized_counts["ENSG00000100298", PPM1B_var_samples]
PPM1B_counts_paired <- normalized_counts["ENSG00000100298", PPM1B_var_samples_paired]
PPM1B_counts_paired_rec <- normalized_counts["ENSG00000100298", PPM1B_var_samples_paired_rec]

# Extract the normalised counts for apobec3h in all other samples
no_PPM1B_counts <- normalized_counts["ENSG00000100298", !colnames(normalized_counts) %in% PPM1B_var_samples]

boxplot(PPM1B_counts, no_PPM1B_counts)
boxplot(PPM1B_counts_paired, PPM1B_counts_paired_rec)

PPM1B_counts_combined <- data.frame(
  Sample = rep(c(names(PPM1B_counts_paired), names(PPM1B_counts_paired_rec)), each = 1),
  Counts = c(PPM1B_counts_paired, PPM1B_counts_paired_rec),
  Condition = rep(c("Acute", "Recovery"), each = length(PPM1B_counts_paired))
)

PPM1B_counts_combined

# Plotting of apobec3h exp and C to U variants counts in PPM1B variant individuals
ggplot(PPM1B_counts_combined, aes(x = Condition, y = Counts, fill = Condition)) +
  geom_boxplot() +
  labs(title = "apobec3h gene expression in individuals with PPM1B variant",
       x = "Condition",
       y = "Normalised counts of Apobec3h ")

PPM1B_counts_combined <- merge(PPM1B_counts_combined, thesis_datav4, by.x = "Sample", by.y = "Run", all.x = TRUE)

ggplot(PPM1B_counts_combined, aes(x = Condition, y = CtoU, fill = Condition)) +
  geom_boxplot() +
  labs(title = "C to U variants in individuals with PPM1B variant",
       x = "Condition",
       y = "Count of C to U variants")

# Finding mean C to U variant count
ctou_in_paired_samples <- PPM1B_counts_combined %>%
  select(Condition, CtoU) %>% 
  group_by(Condition) %>% 
  summarise(mean_ctou = mean(CtoU))

# displaying mean counts
ctou_in_paired_samples

# Removing samples without PPM1B and RELL1 variants

samples_to_remove <- c(
  "SRR12926256", "SRR12926262", "SRR12926269", "SRR12926284",
  "SRR12926292", "SRR12926339", "SRR12926347", "SRR12926351",
  "SRR12926212", "SRR12926264", "SRR12926271", "SRR12926278",
  "SRR12926349", "SRR12926352", "SRR12926366"
)

samples <- c(
  "SRR12926202", "SRR12926204", "SRR12926207", "SRR12926208",
  "SRR12926209", "SRR12926212", "SRR12926214", "SRR12926216",
  "SRR12926218", "SRR12926220", "SRR12926222", "SRR12926224",
  "SRR12926225", "SRR12926227", "SRR12926229", "SRR12926230",
  "SRR12926233", "SRR12926235", "SRR12926237", "SRR12926239",
  "SRR12926240", "SRR12926242", "SRR12926244", "SRR12926246",
  "SRR12926249", "SRR12926251", "SRR12926253", "SRR12926255",
  "SRR12926256", "SRR12926258", "SRR12926260", "SRR12926262",
  "SRR12926264", "SRR12926266", "SRR12926267", "SRR12926269",
  "SRR12926271", "SRR12926273", "SRR12926275", "SRR12926276",
  "SRR12926278", "SRR12926280", "SRR12926281", "SRR12926282",
  "SRR12926283", "SRR12926284", "SRR12926285", "SRR12926287",
  "SRR12926289", "SRR12926291", "SRR12926292", "SRR12926335",
  "SRR12926337", "SRR12926339", "SRR12926341", "SRR12926343",
  "SRR12926345", "SRR12926347", "SRR12926349", "SRR12926351",
  "SRR12926352", "SRR12926354", "SRR12926355", "SRR12926357",
  "SRR12926359", "SRR12926361", "SRR12926363", "SRR12926365",
  "SRR12926366"
)

remaining_samples <- setdiff(samples, samples_to_remove)
# Extract the normalised counts for apobec3h in acute and rec condition
apobec3h_counts_non_variants <- normalized_counts["ENSG00000100298", remaining_samples]
# Convert to data frame
apobec3h_counts_non_variants <- data.frame(
  Sample = names(apobec3h_counts_non_variants),
  Counts = as.numeric(apobec3h_counts_non_variants),
  row.names = NULL
)

# Add rest of data
apobec3h_counts_non_variants <- merge(apobec3h_counts_non_variants, thesis_datav4, by.x = "Sample", by.y = "Run", all.x = TRUE) %>%
  filter(!grepl("single", Paired_single))

id_to_extract <- apobec3h_counts_non_variants %>% 
  select(Participant_ID)

samples_to_extract <- thesis_datav4 %>% 
  select(Run, Participant_ID, Time_point) %>%
  filter(grepl("recovery", Time_point)) %>% 
  filter(Participant_ID %in% id_to_extract$Participant_ID)

samples_to_extract <- samples_to_extract %>% 
  select(Run)

samples_to_extract <- c(samples_to_extract$Run)
samples_to_extract

# Extract the normalised counts for apobec3h in all other samples
no_variants_counts_rec <- normalized_counts["ENSG00000100298", colnames(normalized_counts) %in% samples_to_extract]
str(no_variants_counts_rec)

#convert to dataframe
no_variants_counts_rec <- data.frame(
  Sample = names(no_variants_counts_rec),
  Counts = as.numeric(no_variants_counts_rec),
  row.names = NULL
)

no_variants_counts_rec
apobec3h_counts_non_variants

apobec3h_non_variants <-rbind(no_variants_counts_rec, apobec3h_counts_non_variants)

# Add rest of data
apobec3h_non_variants <- merge(apobec3h_non_variants, thesis_datav4, by.x = "Sample", by.y = "Run", all.x = TRUE) %>%
  filter(!grepl("single", Paired_single))

ggplot(apobec3h_non_variants, aes(x = Time_point, y = Counts, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "apobec3h gene expression in individuals without PPM1B oR RELL1 variant",
       x = "Condition",
       y = "Normalised counts of Apobec3h ")

ggplot(apobec3h_non_variants, aes(x = Time_point, y = CtoU, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "C to U variants in individuals without PPM1B oR RELL1 variant",
       x = "Condition",
       y = "Count of C to U variants")

#------------------------------------------------------------#
# Combining gene expression levels with the rest of the data

thesis_datav4_df <- as.data.frame(thesis_datav4)
normalized_counts_DEdds2 <- t(counts(DEdds2, normalized = TRUE))
str(normalized_counts_DEdds2)

# Extract the normalised counts data from the matrix
matrix_data <- as.matrix(normalized_counts_DEdds2)

# Create a data frame
matrix_dataframe <- as.data.frame(matrix_data)

# Set row names from the first list in dimnames
rownames(matrix_dataframe) <- dimnames(matrix_data)[[1]]

# Set column names from the second list in dimnames
colnames(matrix_dataframe) <- dimnames(matrix_data)[[2]]

# Convert row names of matrix_dataframe to a column named "Run"
matrix_dataframe$Run <- rownames(matrix_dataframe)

# Merged everything!!!
merged_df <- merge(thesis_datav4_df, matrix_dataframe, by = "Run", all.x = TRUE)


# plotting box plots of cytokine release data
ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IFN_gamma), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "LPS_IFN_gamma",
       x = "Condition",
       y = "log2 of delta LPS_IFN_gamma")

ggplot(tdf, aes(x = Time_point, y = (LPS_IL_10/1000), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "       Delta IL-10 secretion with LPS vs 
       control medium across groups",
       x = "Group",
       y = "delta Il-10 (ng/mL)") +
  scale_y_continuous(limits = c(0, 4))+
  theme_dark()+
  theme(legend.position = "bottom")

ggplot(tdf, aes(x = Time_point, y = (LPS_IL_1b/1000), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "       Delta IL-1b secretion with LPS vs 
       control medium across groups",
       x = "Group",
       y = "delta IL-1b (ng/mL)",
       fill = "Groups") +  # Custom label for the legend
  scale_y_continuous(limits = c(0, 30)) +  # Limit y-axis
  theme_dark() +
  theme(legend.position = "bottom")

kruskal.test(LPS_IL_10 ~ Time_point, data = merged_df)
dunn.test(merged_df$LPS_IL_10, merged_df$Time_point, method="bh")

ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IL_27), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log 2 of LPS_IL_27",
       x = "Condition",
       y = "log 2 delta LPS_IL_27")

ggplot(merged_df, aes(x = Time_point, y = LPS_IL_1b, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log 2 of LPS_IL_1b",
       x = "Condition",
       y = "log 2 of delta LPS_IL_1b")

ggplot(merged_df, aes(x = Time_point, y = LPS_IL_6, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "LPS_IL_6",
       x = "Condition",
       y = "delta LPS_IL_6")

ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IL_12p70), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 LPS_IL_12p70",
       x = "Condition",
       y = "log2 delta LPS_IL_12p70")

ggplot(merged_df, aes(x = Time_point, y = (log2(LPS_TNF_alpha)+1), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 LPS_TNF_alpha",
       x = "Condition",
       y = "log2 delta LPS_TNF_alpha")

#plots of acute important cytokine-release data
cytokines_acute <- merged_df %>%
  filter(Time_point == "acute") %>%
  select(Participant_ID, LPS_TNF_alpha, LPS_IL_1b, LPS_IL_6, LPS_IL_10)

cytokines_acute
plot(cytokines_acute$LPS_TNF_alpha)
plot(merged_df$LPS_TNF_alpha, col = merged_df$Time_point)
str(merged_df)
plot(cytokines_acute$LPS_IL_1b)
plot(cytokines_acute$LPS_IL_6)
plot(cytokines_acute$LPS_IL_10)
summary(cytokines_acute[,-1])


ggplot(merged_df, aes(x = Time_point, y = log2(IFN_gamma), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IFN_gamma",
       x = "Condition",
       y = "log2 delta Kp IFN_gamma")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_10), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_10",
       x = "Condition",
       y = "log2 delta Kp IL_10")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_27), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_27",
       x = "Condition",
       y = "log2 delta Kp IL_27")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_1b), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_1b",
       x = "Condition",
       y = "log2 delta Kp IL_1b")

ggplot(merged_df, aes(x = Time_point, y = IL_6, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "Kp IL_6",
       x = "Condition",
       y = "delta Kp IL_6")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_12p70), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_12p70",
       x = "Condition",
       y = "log2 delta Kp IL_12p70")

ggplot(merged_df, aes(x = Time_point, y = log2(TNF_alpha), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp TNF_alpha",
       x = "Condition",
       y = "log2 delta Kp TNF_alpha")

#Plots of the distribution of cytokines in the different groups:
vars_of_int <- merged_df %>% 
  select(LPS_TNF_alpha, ENSG00000100298 )
# Checking for correlation between TNF alpha release and CtoU variant count
# and Apobec3h normalised gene counts.
TNFvsApobec3h<-merged_df %>% 
  select(Time_point, Participant_ID, LPS_TNF_alpha, ENSG00000100298, CtoU) %>% 
  filter(Time_point == "acute") %>% 
  filter(LPS_TNF_alpha < 40000) # remove outliers

scatter.smooth(TNFvsApobec3h$CtoU, TNFvsApobec3h$LPS_TNF_alpha) # no relationship
scatter.smooth(TNFvsApobec3h$CtoU, TNFvsApobec3h$LPS_TNF_alpha)
rcorr(TNFvsApobec3h$CtoU, TNFvsApobec3h$LPS_TNF_alpha, type = "spearman")

merged_df %>% 
  select(Time_point, ENSG00000100298, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000, ENSG00000100298 < 50) %>% 
  select(LPS_TNF_alpha, ENSG00000100298) %>% 
  plot()

merged_df %>% 
  select(Time_point, CtoU, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000) %>% 
  select(LPS_TNF_alpha, CtoU) %>% 
  plot()

plot(merged_df$LPS_TNF_alpha)
plot(y=log2(merged_df$LPS_TNF_alpha), x=log2(merged_df$ENSG00000100298))

merged_df %>% 
  select(Time_point, ENSG00000100298, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000, ENSG00000100298 < 50) %>% 
  select(LPS_TNF_alpha, ENSG00000100298) %>% 
  plot()

merged_df %>% 
  select(Time_point, CtoU, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000) %>% 
  select(LPS_TNF_alpha, CtoU) %>% 
  plot()

# Is there a relationship between any of the APOBEC genes of interest 
# and the cytokine release data?
# ALL SAMPLES
covariates <- merged_df %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000128383, ENSG00000179750, ENSG00000244509,
         ENSG00000100298, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)
covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

#Acute stage only:
covariates <- merged_df %>% 
  select(Time_point, CtoU, AtoI_TOTALmatches, 
         ENSG00000128383, ENSG00000179750, ENSG00000244509,
         ENSG00000100298, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>% 
  filter(Time_point == "acute") %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000128383, ENSG00000179750, ENSG00000244509,
         ENSG00000100298, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)

covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

plot(covariates$CtoU, covariates$AtoI_TOTALmatches)

# Is there a relationship between any of the ADAR genes of interest 
# and the cytokine release data?
# ALL SAMPLES

covariates <- merged_df %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000160710,
         ENSG00000197381,
         LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)
covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

#Acute stage only:
covariates <- merged_df %>% 
  select(Time_point, CtoU, AtoI_TOTALmatches, 
         ENSG00000160710,
         ENSG00000197381, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>% 
  filter(Time_point == "acute") %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000160710,
         ENSG00000197381, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)
covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

# Box plots of cytokine release 
# working with LPS data
# Acute Individuals with and without variants of interest
# rell1

with_rell1_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
with_rell1_var <- with_rell1_var[with_rell1_var$Run %in% rell1_var_samples, ]
without_rell1_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
without_rell1_var <- without_rell1_var[!(without_rell1_var$Run %in% rell1_var_samples), ]
without_rell1_var <- na.omit(without_rell1_var)
without_rell1_var

# Getting cytokine release data of the other groups as well

cytokines_controls <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>%
  filter(Time_point == "control")
cytokines_controls <- na.omit(cytokines_controls)

cytokines_recovery <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>%
  filter(Time_point == "recovery")
cytokines_recovery <- na.omit(cytokines_recovery)

boxplot(with_rell1_var$LPS_IL_10, without_rell1_var$LPS_IL_10, 
        cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10,  
        main = "Comparison of IL10 release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_10", 
        names = c("Var", "!Var","Rec","Cont")) # IL-10

t_test_result <- t.test(with_rell1_var$LPS_IL_10, without_rell1_var$LPS_IL_10)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_rell1_var$LPS_IL_10, without_rell1_var$LPS_IL_10)
print(wilcox_test_result)

boxplot(log2(with_rell1_var$LPS_IL_1b), log2(without_rell1_var$LPS_IL_1b), 
        log2(cytokines_recovery$LPS_IL_1b), log2(cytokines_controls$LPS_IL_1b),  
        main = "Comparison of IL-1b release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "log2(delta LPS_IL_1b)", 
        names = c("Var", "!Var","Rec","Cont")) # IL-1b

boxplot(with_rell1_var$LPS_IL_6, without_rell1_var$LPS_IL_6, 
        cytokines_recovery$LPS_IL_6, cytokines_controls$LPS_IL_6,  
        main = "Comparison of IL6 release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_6", 
        names = c("Var", "!Var","Rec","Cont")) # IL-6

t_test_result <- t.test(with_rell1_var$LPS_IL_6, without_rell1_var$LPS_IL_6)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_rell1_var$LPS_IL_6, without_rell1_var$LPS_IL_6)
print(wilcox_test_result)

boxplot(log2(with_rell1_var$LPS_TNF_alpha), log2(without_rell1_var$LPS_TNF_alpha), 
        log2(cytokines_recovery$LPS_TNF_alpha), log2(cytokines_controls$LPS_TNF_alpha),  
        main = "Comparison of TNF-a release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "log 2 delta TNF-a", 
        names = c("Var", "!Var","Rec","Cont")) # TNF-a

t_test_result <- t.test(with_rell1_var$LPS_TNF_alpha, without_rell1_var$LPS_TNF_alpha)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_rell1_var$LPS_TNF_alpha, without_rell1_var$LPS_TNF_alpha)
print(wilcox_test_result)

# PPM1B (PPM1B_var_samples)
with_PPM1B_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
with_PPM1B_var <- with_PPM1B_var[with_PPM1B_var$Run %in% PPM1B_var_samples, ]
with_PPM1B_var <- na.omit(with_PPM1B_var)

without_PPM1B_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
without_PPM1B_var <- without_PPM1B_var[!(without_PPM1B_var$Run %in% PPM1B_var_samples), ]
without_PPM1B_var <- na.omit(without_PPM1B_var)
without_PPM1B_var

#removing a extremely high Il-10 outlier from the acute group
without_PPM1B_var <- without_PPM1B_var[without_PPM1B_var$Run != "SRR12926276", ]

boxplot(with_PPM1B_var$LPS_IL_10, without_PPM1B_var$LPS_IL_10, 
        cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10,  
        main = "Comparison of IL10 release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_10", 
        names = c("Var", "!Var","Rec","Cont")) # IL-10

boxplot(with_PPM1B_var$LPS_IL_1b, without_PPM1B_var$LPS_IL_1b, 
        cytokines_recovery$LPS_IL_1b, cytokines_controls$LPS_IL_1b,  
        main = "Comparison of IL1b release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_1b", 
        names = c("Var", "!Var","Rec","Cont")) # IL-1b

boxplot(with_PPM1B_var$LPS_IL_6, without_PPM1B_var$LPS_IL_6, 
        cytokines_recovery$LPS_IL_6, cytokines_controls$LPS_IL_6,  
        main = "Comparison of IL6 release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_6", 
        names = c("Var", "!Var","Rec","Cont")) # IL-6

t_test_result <- t.test(with_PPM1B_var$LPS_IL_6, without_PPM1B_var$LPS_IL_6)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_PPM1B_var$LPS_IL_6, without_PPM1B_var$LPS_IL_6)
print(wilcox_test_result)

boxplot(log2(with_PPM1B_var$LPS_TNF_alpha), log2(without_PPM1B_var$LPS_TNF_alpha), 
        log2(cytokines_recovery$LPS_TNF_alpha), log2(cytokines_controls$LPS_TNF_alpha),  
        main = "Comparison of TNF-a release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "log2 delta TNF-a", 
        names = c("Var", "!Var","Rec","Cont")) # TNF-a

t_test_result <- t.test(with_PPM1B_var$LPS_TNF_alpha, without_PPM1B_var$LPS_TNF_alpha)
print(t_test_result)

# Age

plot(merged_df$age, merged_df$AtoIandCtoU, main = "No of C to U variants vs Age", 
     xlab = "Age", ylab = "CtoU")

# No. of C to U vs age

plot(merged_df$age, merged_df$CtoU, main = "No of C to U variants vs Age", 
     xlab = "Age", ylab = "CtoU")

# Add a linear trend line
abline(lm(merged_df$CtoU ~ merged_df$age), col = "red")
rcorr(merged_df$CtoU, merged_df$age, type = "spearman")


merged_df %>%
  filter(Time_point == "acute") %>%
  ggplot(aes(x = age, y = CtoU)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  labs(title = "Scatter Plot with Trend Line",
       x = "Age",
       y = "CtoU")

#------------------ Responders and non-responders --------------#
cytokine_column_names <- colnames(merged_df[,10:23])

cytokine_responses

cytokine_responses <- tdf %>% 
  select(Participant, Time_point, Paired_single, all_of(cytokine_column_names)) %>% 
  filter(Paired_single == "paired")

cytokine_responders <- cytokine_responses %>%
  group_by(Participant) %>%
  mutate(across(all_of(cytokine_column_names), 
                ~ ifelse(.[Time_point == "acute"] == .[Time_point == "recovery"], NA, 
                         ifelse(.[Time_point == "acute"] < .[Time_point == "recovery"],
                                "lower", "higher"))))


cytokine_responders <- cytokine_responders %>% 
  filter(Time_point == "acute")

cytokine_responders <- as.data.frame(cytokine_responders)


#counting the number of "lower" and "higher" cytokine releases:
cytokine_responders$lower_count <- apply(cytokine_responders == "lower", 1, function(x) sum(x, na.rm = TRUE))
cytokine_responders$higher_count <- apply(cytokine_responders == "higher", 1, function(x) sum(x, na.rm = TRUE))

# Getting runs SRR numbers
acute_paired_runs <- tdf %>%
  select(Run, Time_point, Paired_single) %>% 
  filter(Time_point == "acute") %>%
  filter(Paired_single == "paired") %>% 
  select(Run)

cytokine_responders$Run <-acute_paired_runs$Run

#Move Run as the first column
cytokine_responders <- cytokine_responders %>% relocate(Run)

cytokine_responders$state <- ifelse(cytokine_responders$lower_count > cytokine_responders$higher_count,
                                    "tolerant", "normal")

cytokine_responders$state <- ifelse(cytokine_responders$lower_count == 0 & cytokine_responders$higher_count == 0,
                                    "no_data", cytokine_responders$state)


#moving state as first column
cytokine_responders <- cytokine_responders %>% 
  select(state, everything())

#Adding tolerant label to tdf

tolerant_runs <- cytokine_responders$Run[cytokine_responders$state == "tolerant"]

cytokine_state_runs <- cytokine_responders %>%
  select(Run, state)

# Join the new dataframe with 'tdf' and create 'cytokine_state' column
tdf <- tdf %>%
  left_join(cytokine_state_runs, by = "Run") %>%
  mutate(cytokine_state = ifelse(is.na(state), Time_point, state)) %>% 
  select(-state)

tdf$cytokine_state <- as.factor(tdf$cytokine_state)

write.table(cytokine_responders, file = "cytokine_responders.csv", sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

#--------- now that I added some more group information, I'll re-plot some
# of the important things:
rna_edit_genes_counts$cytokine_state <- tdf$cytokine_state 

#By RNA-editing genes:
ggplot(rna_edit_genes_counts, aes(x = cytokine_state,
                y = rna_edit_genes_counts$ENSG00000128383,
                color = cytokine_state)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "APOBEC3A expression",
       x = "Group",
       y = "APOBEC3A expression",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rna_edit_genes_counts, aes(x = cytokine_state,
                                  y = rna_edit_genes_counts$ENSG00000197381,
                                  color = cytokine_state)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "ADARB1 expression",
       x = "Group",
       y = "ADARB1 expression",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rna_edit_genes_counts, aes(x = cytokine_state,
                                  y = rna_edit_genes_counts$ENSG00000244509,
                                  color = cytokine_state)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "APOBEC3C expression",
       x = "Group",
       y = "APOBEC3C expression",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rna_edit_genes_counts, aes(x = cytokine_state,
                                  y = rna_edit_genes_counts$ENSG00000100298,
                                  color = cytokine_state)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "APOBEC3H expression",
       x = "Group",
       y = "APOBEC3H expression",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(rna_edit_genes_counts, aes(x = cytokine_state,
                                  y = rna_edit_genes_counts$ENSG00000160710,
                                  color = cytokine_state)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75),
               alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "ADAR expression",
       x = "Group",
       y = "ADAR expression",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# removing acute samples with no paired cytokine release
# data at recovery stage

paired_samples <- rna_edit_genes_counts %>% 
  filter(cytokine_state != "acute") %>% 
  filter(cytokine_state != "no_data") %>%
  filter(cytokine_state != "recovery") %>% 
  filter(cytokine_state != "control") 
  

# removing non-numeric variables:
numeric_paired_samples <- paired_samples %>%
  select_if(is.numeric)

numeric_paired_samples$cytokine_state <- paired_samples$cytokine_state

numeric_paired_samples <- numeric_paired_samples
numeric_paired_samples$
#performing wilcoxson tests on normal vs tolerant:

perform_wilcox_test <- function(var) {
  formula <- as.formula(paste(var, "~ cytokine_state"))
  wilcox_result <- wilcox.test(formula, data = numeric_paired_samples)
  tidy_result <- broom::tidy(wilcox_result)
  tidy_result$variable <- var
  return(tidy_result)
}

# Get numeric column names
numeric_vars <- numeric_paired_samples %>%
  select_if(is.numeric) %>%
  colnames()

# Apply the function to each numeric variable
results <- map_df(numeric_vars, perform_wilcox_test)

# Sort the results by p-value
sorted_results <- results %>%
  arrange(p.value)
sorted_results

#------------- checking some other correlations---------#

plot(tdf$LPS_IL_6, tdf$AtoI_Ch7matches)

ggplot(tdf, aes(x=LPS_IL_10, y=CtoU_chr3, color = Time_point)) +
  geom_point() +
  theme_minimal()

plot(tdf$LPS_IL_6, tdf$CtoU_chr12)

# Perform PCA
ctou_cytokines <- tdf %>% 
  select(3:17, starts_with("CtoU_chr"))
str(ctou_cytokines)
ctou_cytokines <- na.omit(ctou_cytokines)
pca_ctou_cytokines <- prcomp(ctou_cytokines[,2:37], center = TRUE, scale. = TRUE)
summary(pca_ctou_cytokines)
pca_ctou_cytokines
autoplot(pca_ctou_cytokines, data = ctou_cytokines, colour = 'Time_point',
         alpha = 0.7, size = 2)

# correlation between AtoI22 and ENSG00000160710 (ADAR gene)
numeric_rna_edit_genes$AtoI_Ch22matches
ggplot(numeric_rna_edit_genes, aes(x = ENSG00000160710, y = AtoI_Ch22matches)) +
  geom_point(aes(color = Time_point)) +  # Color points by Time_point
  geom_smooth(method = lm, se = FALSE, color = "black") +  # Add linear trend line
  labs(title = "Correlation between A to I Chr 22 and ADAR gene expression",
       x = "ADAR normalised counts",
       y = "A to I Chr 22 counts")

rcorr(rna_edit_genes_counts$ENSG00000160710,
      rna_edit_genes_counts$AtoI_Ch22, type = "spearman")


ggplot(rna_edit_genes_counts, aes(x = ENSG00000197381, y = AtoI_Ch12matches)) +
  geom_point(aes(color = Time_point)) +  # Color points by Time_point
  geom_smooth(method = lm, se = FALSE, color = "black") +  # Add linear trend line
  labs(title = "Correlation between A to I Chr 12 and ADARB1 gene expression",
       x = "ADARB1 normalised counts",
       y = "A to I Ch12 counts")

rcorr(rna_edit_genes_counts$ENSG00000197381,
      rna_edit_genes_counts$AtoI_Ch12matches, type = "spearman")


ggplot(rna_edit_genes_counts, aes(x = ENSG00000244509, y = CtoU_chr3)) +
  geom_point(aes(color = Time_point)) +  # Color points by Time_point
  geom_smooth(method = lm, se = FALSE, color = "black") +  # Add linear trend line
  labs(title = "Correlation between CtoU_chr3 and APOBEC3C gene expression",
       x = "APOBEC3C normalised counts",
       y = "CtoU_chr3 counts")

rcorr(rna_edit_genes_counts$ENSG00000244509,
      rna_edit_genes_counts$CtoU_chr3, type = "spearman")

rcorr(rna_edit_genes_counts$LPS_IL_10,
      rna_edit_genes_counts$CtoU_chr3, type = "spearman")

#------------ merging all the normalised counts with tdf----#

numeric_tdf <- tdf[sapply(tdf, is.numeric)]
numeric_tdf <- numeric_counts_tdf[, 2:73]
merged_norm_tdf <- cbind(numeric_tdf, normalized_counts)
merged_norm_tdf$Time_point <- tdf$Time_point

rna_edit_genes_counts
# Compute correlation
global_cor_matrix <- cor(numeric_tdf, normalized_counts,
                         method = "spearman", use = "pairwise.complete.obs")
transposed_cor_matrix <- t(global_cor_matrix)

write.table(transposed_cor_matrix, file = "transposed_cor_matrix.csv", sep = "\t", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)

#--------visualising some relationships:

ggplot(merged_norm_tdf, aes(x = ENSG00000145016, y = AtoI_Ch22matches)) +
  geom_point(aes(color = Time_point)) +  # Color points by Time_point
  geom_smooth(method = lm, se = FALSE, color = "black") +  # Add linear trend line
  labs(title = "Correlation between AtoI_Ch22 and RUBCN gene expression",
       x = "RUBCN normalised counts",
       y = "AtoI_Ch22matches counts")

merged_norm_tdf$CtoU_chr4

ggplot(merged_norm_tdf, aes(x = ENSG00000122566, y = CtoU)) +
  geom_point(aes(color = Time_point)) +  # Color points by Time_point
  geom_smooth(method = lm, se = FALSE, color = "black") +  # Add linear trend line
  labs(title = "Correlation between LPS_IL_10 and FPR2 gene expression",
       x = "FPR2 normalised counts",
       y = "LPS_IL_10")

cor(merged_norm_tdf$ENSG00000122566, merged_norm_tdf$CtoU, method = 'spearman')

#---------------------variant analysis
# Normal filtration
ctou_ann_stats_norm <- read_delim("annovar_output/ctou_ann_stats_norm2.csv", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_types = cols(Start = col_number(), End = col_number(), 
                                                   perc_total = col_number(), perc_acute = col_number(), 
                                                   perc_con = col_number(), perc_rec = col_number()), 
                                  trim_ws = TRUE)

atoi_ann_stats_norm <- read_delim("annovar_output/atoi_ann_stats_norm2.csv", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_types = cols(Start = col_number(), End = col_number(), 
                                                   perc_total = col_number(), perc_acute = col_number(), 
                                                   perc_con = col_number(), perc_rec = col_number()), 
                                  trim_ws = TRUE)

ctou_ann_stats_norm <- as.data.frame(ctou_ann_stats_norm)
atoi_ann_stats_norm <- as.data.frame(atoi_ann_stats_norm)

#strict filtration
ctou_ann_stats_strict <- read_delim("annovar_output/ctou_ann_stats_strict.csv", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_types = cols(Start = col_number(), End = col_number(), 
                                                   perc_total = col_number(), perc_acute = col_number(), 
                                                   perc_con = col_number(), perc_rec = col_number()), 
                                  trim_ws = TRUE)

atoi_ann_stats_strict <- read_delim("annovar_output/atoi_ann_stats_strict.csv", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_types = cols(Start = col_number(), End = col_number(), 
                                                   perc_total = col_number(), perc_acute = col_number(), 
                                                   perc_con = col_number(), perc_rec = col_number()), 
                                  trim_ws = TRUE)

ctou_ann_stats_strict <- as.data.frame(ctou_ann_stats_strict)
atoi_ann_stats_strict <- as.data.frame(atoi_ann_stats_strict)

ctou_ann_stats_strict

# creating unique_var_ids

ctou_ann_stats_norm$unique_var_id <- paste(ctou_ann_stats_norm$Chr,
                                           ctou_ann_stats_norm$Start, sep=":")

atoi_ann_stats_norm$unique_var_id <- paste(atoi_ann_stats_norm$Chr,
                                           atoi_ann_stats_norm$Start, sep=":")

# moving it as the first column
ctou_ann_stats_norm <- ctou_ann_stats_norm[, c("unique_var_id",
                                               setdiff(names(ctou_ann_stats_norm),
                                                       "unique_var_id"))]

atoi_ann_stats_norm <- atoi_ann_stats_norm[, c("unique_var_id",
                                               setdiff(names(atoi_ann_stats_norm),
                                                       "unique_var_id"))]

#-------------Chi square tests on strict edit list of C-to-U
# c to u or a to i...just change name of datasets and results variable
p_values_fisher <- numeric(nrow(ctou_ann_stats_strict))

# Loop through each RNA edit and perform chi-square test
for (i in 1:nrow(ctou_ann_stats_strict)) {
  # Create a 2x2 contingency table for the current RNA edit
  contingency_table <- matrix(c(
    ctou_ann_stats_strict$a_count[i], ctou_ann_stats_strict$a_noncount[i],
    ctou_ann_stats_strict$c_count[i], ctou_ann_stats_strict$c_noncount[i]
  ), nrow = 2)
  
  # Perform chi-square test
  test_result <- fisher.test(contingency_table)
  
  # Store the p-value
  p_values_fisher[i] <- test_result$p.value
}

# Add p-values as the 7th column
ctou_ann_stats_strict$p_value_fisher <- p_values_fisher

# Adjust the p-values using Benjamini-Hochberg (BH) correction
ctou_ann_stats_strict$BH_adjusted_p_fisher <- p.adjust(p_values_fisher, method = "BH")

write.table(ctou_ann_stats_strict, file = "ctou_ann_stats_strict_p.csv", sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(atoi_ann_stats_strict, file = "atoi_ann_stats_strict_p.csv", sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

ctou_strict_volcano<-ctou_ann_stats_strict %>% 
  select(Gene.refGene, log_adj_fisher, log_odds_ratio)

atoi_strict_volcano<-atoi_ann_stats_strict %>% 
  select(Gene.refGene, log_adj_fisher, log_odds_ratio)

write.table(ctou_strict_volcano, file = "ctou_strict_volcano.csv", sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(atoi_strict_volcano, file = "atoi_strict_volcano.csv", sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

# calculating the odds ratio:

total_disease <- 69
total_control <- 56

# Initialize a vector to store log odds ratios with continuity correction
log_odds_ratios <- numeric(nrow(atoi_ann_stats_strict))

# Loop through each RNA edit and compute the log odds ratio with continuity correction
for (i in 1:nrow(atoi_ann_stats_strict)) {
  # Number of disease and control samples with and without the edit
  X1 <- atoi_ann_stats_strict$a_count[i]
  X2 <- atoi_ann_stats_strict$c_count[i]
  
  # Apply continuity correction by adding 0.5 to all values
  X1 <- X1 + 0.5
  X2 <- X2 + 0.5
  no_edit_disease <- (total_disease - atoi_ann_stats_strict$a_count[i]) + 0.5
  no_edit_control <- (total_control - atoi_ann_stats_strict$c_count[i]) + 0.5
  
  # Compute the odds in disease and control groups
  odds_disease <- (X1 / no_edit_disease)
  odds_control <- (X2 / no_edit_control)
  
  # Compute the log odds ratio
  log_odds_ratios[i] <- log(odds_disease / odds_control)
}

# Add the log odds ratio as a new column in the dataframe
atoi_ann_stats_strict$log_odds_ratio <- log_odds_ratios

# add the -log10 of the p adjusted values
atoi_ann_stats_strict$log_adj_fisher <- -log10(atoi_ann_stats_strict$BH_adjusted_p_fisher)
ctou_ann_stats_strict$log_adj_fisher <- -log10(ctou_ann_stats_strict$BH_adjusted_p_fisher)


#-------------Fisher's exact tests
# c to u or a to i...just change name of datasets and results variable
# Function to perform Fisher's exact test
perform_test <- function(df, var1, var2) {
  mat <- matrix(c(df[[paste0(var1, "_count")]], df[[paste0(var2, "_count")]], df[[paste0(var1, "_noncount")]], df[[paste0(var2, "_noncount")]]), ncol = 2)
  test <- fisher.test(mat)
  return(c(test$p.value, test$estimate))
}

# Define the groups
groups <- c("a", "c", "r")

# Perform Fisher's exact test for each pairwise comparison
fishers_atoi <- combn(groups, 2, simplify = FALSE, FUN = function(x) {
  atoi_ann_stats_norm %>%
    rowwise() %>%
    mutate(
      p_value = perform_test(cur_data(), x[1], x[2])[1],
      odds_ratio = perform_test(cur_data(), x[1], x[2])[2]
    ) %>%
    ungroup() %>%
    mutate(
      p_adjust = p.adjust(p_value, method = "BH"),
      comparison = paste(x, collapse = "_vs_")
    ) %>%
    select(unique_var_id, comparison, p_value, p_adjust, odds_ratio)
}) %>% bind_rows()

# attempting a volcano plot

fishers_atoi <- fishers_atoi %>%
  mutate(sig = ifelse(p_adjust < 0.05 & abs(log2(odds_ratio)) > 1, "yes", "no"))

fishers_ctou <- fishers_ctou %>%
  mutate(sig = ifelse(p_adjust < 0.05 & abs(log2(odds_ratio)) > 1, "yes", "no"))

fishers_results <- rbind(fishers_atoi, fishers_ctou)
fishers_results <- fishers_results %>% 
  filter(p_value < 0.5)
# Create the volcano plot
ggplot(fishers_results, aes(x = log2(odds_ratio), y = -log10(p_adjust))) +
  geom_point(aes(color = sig), size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  theme_classic() +
  labs(x = "Log2 Odds Ratio", y = "-Log10 Adjusted P-Value", title = "Volcano Plot")

#--- Checking relationship of those with this variant

EIF2AK2_variant <- c('SRR12926202', 'SRR12926233', 'SRR12926244', 'SRR12926251', 'SRR12926337', 'SRR12926339', 'SRR12926343', 'SRR12926351')
MAVS_variant <-c('SRR12926202', 'SRR12926239', 'SRR12926251', 'SRR12926273', 'SRR12926282', 'SRR12926283', 'SRR12926337', 'SRR12926349', 'SRR12926354', 'SRR12926361')
GNB4_variant<-c('SRR12926202', 'SRR12926207', 'SRR12926212', 'SRR12926218', 'SRR12926230', 'SRR12926235', 'SRR12926244', 'SRR12926246', 'SRR12926253', 'SRR12926282', 'SRR12926287', 'SRR12926289', 'SRR12926339', 'SRR12926345', 'SRR12926347', 'SRR12926361', 'SRR12926363', 'SRR12926365')


# Box plots of cytokine release 
# working with LPS data
# Acute Individuals with and without variants of interest
# EIF2AK2

with_EIF2AK2_var <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

with_EIF2AK2_var <- with_EIF2AK2_var[with_EIF2AK2_var$Run %in% EIF2AK2_variant, ]

without_EIF2AK2_var <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

without_EIF2AK2_var <- without_EIF2AK2_var[!(without_EIF2AK2_var$Run %in% EIF2AK2_variant), ]

without_EIF2AK2_var <- na.omit(without_EIF2AK2_var)
with_EIF2AK2_var <- na.omit(with_EIF2AK2_var)


# Getting cytokine release data of the other groups as well

cytokines_controls <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>%
  filter(Time_point == "control")
cytokines_controls <- na.omit(cytokines_controls)

cytokines_recovery <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>%
  filter(Time_point == "recovery")
cytokines_recovery <- na.omit(cytokines_recovery)

boxplot(with_EIF2AK2_var$LPS_IL_10, without_EIF2AK2_var$LPS_IL_10, 
        cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10,  
        main = "Comparison of IL10 release
        with EIF2AK2 variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_10", 
        names = c("Var", "!Var","Rec","Cont")) # IL-10

# converting to a ggplot:

EIF2AK2_data <- data.frame(
  LPS_IL_10 = c(with_EIF2AK2_var$LPS_IL_10, without_EIF2AK2_var$LPS_IL_10, 
                cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10),
  Group = factor(rep(c("EIF2AK2 edit", "acute","recovery","control"), 
                     times = c(length(with_EIF2AK2_var$LPS_IL_10), length(without_EIF2AK2_var$LPS_IL_10), 
                               length(cytokines_recovery$LPS_IL_10), length(cytokines_controls$LPS_IL_10))))
)

ggplot(EIF2AK2_data, aes(x = Group, y = LPS_IL_10/1000, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "IL-10 release with EIF2AK2 edit",
       x = "Group", y = "delta IL-10(ng/mL)") +
  theme_dark() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 4))

# Perform Kruskal-Wallis test
kruskal.test(LPS_IL_10 ~ Group, data = EIF2AK2_data)

EIF2AK2_dunn_result <- dunn.test(EIF2AK2_data$LPS_IL_10, EIF2AK2_data$Group, method="bonferroni")

# Acute Individuals with and without variants of interest
# MAVS

with_MAVS_variant <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

with_MAVS_variant <- with_MAVS_variant[with_MAVS_variant$Run %in% MAVS_variant, ]

without_MAVS_var <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

without_MAVS_var <- without_MAVS_var[!(without_MAVS_var$Run %in% MAVS_variant), ]

without_MAVS_var <- na.omit(without_MAVS_var)
with_MAVS_variant <- na.omit(with_MAVS_variant)

# converting to a ggplot:

MAVS_data <- data.frame(
  LPS_IL_10 = c(with_MAVS_variant$LPS_IL_10, without_MAVS_var$LPS_IL_10, 
                cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10),
  Group = factor(rep(c("MAVS edit", "acute","recovery","control"), 
                     times = c(length(with_MAVS_variant$LPS_IL_10), length(without_MAVS_var$LPS_IL_10), 
                               length(cytokines_recovery$LPS_IL_10), length(cytokines_controls$LPS_IL_10))))
)

ggplot(MAVS_data, aes(x = Group, y = LPS_IL_10/1000, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "IL-10 release with MAVS edit",
       x = "Group", y = "delta IL-10 (ng/mL)") +
  theme_dark() +
  scale_y_continuous(limits = c(0, 3.5)) +
  theme(legend.position = "none")

# Acute Individuals with and without variants of interest
# GNB4

with_GNB4_variant <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

with_GNB4_variant <- with_GNB4_variant[with_GNB4_variant$Run %in% GNB4_variant, ]

without_GNB4_var <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

without_GNB4_var <- without_GNB4_var[!(without_GNB4_var$Run %in% GNB4_variant), ]

without_GNB4_var <- na.omit(without_GNB4_var)
with_GNB4_variant <- na.omit(with_GNB4_variant)

# converting to a ggplot:

GNB4_data <- data.frame(
  LPS_IL_1b = c(with_GNB4_variant$LPS_IL_1b, without_GNB4_var$LPS_IL_1b, 
                cytokines_recovery$LPS_IL_1b, cytokines_controls$LPS_IL_1b),
  Group = factor(rep(c("GNB4 edit", "acute","recovery","control"), 
                     times = c(length(with_GNB4_variant$LPS_IL_1b), length(without_GNB4_var$LPS_IL_1b), 
                               length(cytokines_recovery$LPS_IL_1b), length(cytokines_controls$LPS_IL_1b))))
)

ggplot(GNB4_data, aes(x = Group, y = (LPS_IL_1b/1000), fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "IL-1b release with GNB4 edit",
       x = "Group", y = "delta IL-1b (ng/mL)") +
  theme_dark() +
  scale_y_continuous(limits = c(-3, 30)) +
  theme(legend.position = 'none')

# GPR65

GPR65<-c('SRR12926202', 'SRR12926209', 'SRR12926220', 'SRR12926224', 'SRR12926225', 'SRR12926230', 'SRR12926235', 'SRR12926244', 'SRR12926246', 'SRR12926249', 'SRR12926251', 'SRR12926253', 'SRR12926264', 'SRR12926269', 'SRR12926271', 'SRR12926273', 'SRR12926275', 'SRR12926278', 'SRR12926283', 'SRR12926284', 'SRR12926289', 'SRR12926337', 'SRR12926339', 'SRR12926341', 'SRR12926343', 'SRR12926349', 'SRR12926351', 'SRR12926352', 'SRR12926354', 'SRR12926355', 'SRR12926361', 'SRR12926365')

with_GPR65_variant <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

with_GPR65_variant <- with_GPR65_variant[with_GPR65_variant$Run %in% GPR65, ]

without_GPR65_var <- tdf %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")

without_GPR65_var <- without_GPR65_var[!(without_GPR65_var$Run %in% GPR65), ]

without_GPR65_var <- na.omit(without_GPR65_var)
with_GPR65_variant <- na.omit(with_GPR65_variant)

# converting to a ggplot:

GPR65_data <- data.frame(
  LPS_IL_10 = c(with_GPR65_variant$LPS_IL_10, without_GPR65_var$LPS_IL_10, 
                cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10),
  Group = factor(rep(c("GPR65 edit", "acute","recovery","control"), 
                     times = c(length(with_GPR65_variant$LPS_IL_10), length(without_GPR65_var$LPS_IL_10), 
                               length(cytokines_recovery$LPS_IL_10), length(cytokines_controls$LPS_IL_10))))
)

ggplot(GPR65_data, aes(x = Group, y = (LPS_IL_10)/1000, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "IL-10 release with GPR65 edit",
       x = "Group", y = "delta Il-10 (ng/mL)") +
  theme_dark() +
  scale_y_continuous(limits = c(0, 4)) +
  theme(legend.position = 'none')

# Perform Kruskal-Wallis test
kruskal.test(LPS_IL_10 ~ Group, data = EIF2AK2_data)

EIF2AK2_dunn_result <- dunn.test(EIF2AK2_data$LPS_IL_10, EIF2AK2_data$Group, method="bonferroni")

# ['SRR12926229', 'SRR12926246', 'SRR12926249', 'SRR12926276', 'SRR12926349', 'SRR12926359'] 

# Trying to find the mean of cytokines per common acute variant

# Read the .txt file
lines <- readLines("atoi_acute_srr.txt")

# Initialize a list to store the means
means <- list()
medians <- list()
# Loop over each line in the file
for (i in seq_along(lines)) {
  # Extract the SRR numbers
  srr_numbers <- gsub("\\[|\\]|'| ", "", lines[i])  # Remove brackets, apostrophes, and spaces
  srr_numbers <- strsplit(srr_numbers, ",")[[1]]  # Split the string into individual SRR numbers
  
  # Subset the dataframe to include only these SRR numbers
  subset_tdf <- subset(tdf, Run %in% srr_numbers)
  
  # Calculate the mean/median of the LPS_IL_10 values for these SRR numbers
  #mean_cytokine <- mean(subset_tdf$LPS_TNF_alpha, na.rm = TRUE)
  median_cytokine <- median(subset_tdf$LPS_TNF_alpha, na.rm = TRUE)

  # Store the mean or median in the list
  #means[[i]] <- mean_cytokine
  medians[[i]] <- median_cytokine
}

# Convert the list of means to a character vector
means_char <- sapply(means, as.character)
medians_char <- sapply(medians, as.character)

# Write the means to a file, with each mean on a new line
writeLines(means_char, "LPS_TNF_alpha_atoi.txt")
writeLines(medians_char, "median_LPS_TNF_alpha_atoi.txt")

# Checking if there is a difference in 
tdf_il_10_sexes_acute <- tdf %>% 
  select(LPS_IL_10, sex, Time_point) %>% 
  filter(Time_point == "acute")

ggplot(tdf_il_10_sexes_acute, aes(x = sex, y = LPS_IL_10)) +
  geom_boxplot() +
  labs(title = "Boxplot of LPS_IL_10 by Sex",
       x = "Sex",
       y = "LPS_IL_10")

# mean and median of cytokines by group:
mean_median_cytokines <- tdf %>%
  group_by(Time_point) %>%
  summarise(
    mean_LPS_IL_10 = mean(LPS_IL_10, na.rm = TRUE),
    median_LPS_IL_10 = median(LPS_IL_10, na.rm = TRUE),
    mean_LPS_IL_6 = mean(LPS_IL_6, na.rm = TRUE),
    median_LPS_IL_6 = median(LPS_IL_6, na.rm = TRUE),
    mean_LPS_IL_1b = mean(LPS_IL_1b, na.rm = TRUE),
    median_LPS_IL_1b = median(LPS_IL_1b, na.rm = TRUE),
    mean_LPS_TNF_alpha = mean(LPS_TNF_alpha, na.rm = TRUE),
    median_LPS_TNF_alpha = median(LPS_TNF_alpha, na.rm = TRUE)
  )

# Print the summary dataframe
print(mean_median_cytokines)

# Plotting the last box plots of cytokine release with coloured points for 
# variant-carrying samples
interesting_variants <- read_csv("interesting_variants.csv", col_names = FALSE)

interesting_variants <-as.data.frame(interesting_variants)

# get gene names with variants:
gene_names <- interesting_variants %>% pull(X1)

get_srr_numbers <- function(df, gene_name) {
  srr_string <- df %>% 
    filter(X1 == gene_name) %>% 
    mutate(X2 = gsub("\\['", "", X2),
           X2 = gsub("'\\]", "", X2),
           X2 = gsub("'", "", X2),
           X2 = gsub(" ", "", X2)) %>%
    pull(X2)
  srr_numbers <- strsplit(srr_string, ",")[[1]]
  return(srr_numbers)
}

add_gene_column <- function(df, gene_name, srr_numbers) {
  df[[gene_name]] <- ifelse(df$Run %in% srr_numbers, gene_name, "no")
  return(df)
}

for (gene in gene_names) {
  srr_numbers <- get_srr_numbers(interesting_variants, gene)
  tdf <- add_gene_column(tdf, gene, srr_numbers)
}

#adding ADAM gene variant samples:
ADAMTS2 <- c('SRR12926208', 'SRR12926222', 'SRR12926224', 'SRR12926229', 
             'SRR12926230', 'SRR12926235', 'SRR12926246', 'SRR12926253', 
             'SRR12926283', 'SRR12926284', 'SRR12926337', 'SRR12926339', 
             'SRR12926345', 'SRR12926346', 'SRR12926349', 'SRR12926352',
             'SRR12926361')

tdf$ADAMTS2 <- ifelse(tdf$Run %in% ADAMTS2, "ADAMTS2", "no")

#adding LILRA2 gene variant samples
LILRA2 <- c('SRR12926202', 'SRR12926209', 'SRR12926230', 'SRR12926244',
            'SRR12926273', 'SRR12926283', 'SRR12926285', 'SRR12926339',
            'SRR12926341', 'SRR12926351')

tdf$LILRA2 <- ifelse(tdf$Run %in% LILRA2, "LILRA2", "no")

#adding ADGRE3 gene variant samples
ADGRE3 <- c('SRR12926204', 'SRR12926249', 'SRR12926283', 'SRR12926351')

tdf$ADGRE3 <- ifelse(tdf$Run %in% ADGRE3, "ADGRE3", "no")

#adding HPSE gene variant samples
HPSE <- c('SRR12926202', 'SRR12926218', 'SRR12926233', 'SRR12926235',
          'SRR12926262', 'SRR12926272', 'SRR12926282', 'SRR12926339',
          'SRR12926347', 'SRR12926351', 'SRR12926354', 'SRR12926361',
          'SRR12926367')

tdf$HPSE <- ifelse(tdf$Run %in% HPSE, "HPSE", "no")

#adding TRAF3IP2-AS1 gene variant samples
TRAF3IP2_AS1 <- c('SRR12926216', 'SRR12926230', 'SRR12926244', 'SRR12926256',
                  'SRR12926260', 'SRR12926271', 'SRR12926278', 'SRR12926282',
                  'SRR12926337', 'SRR12926361')
tdf$TRAF3IP2_AS1 <- ifelse(tdf$Run %in% TRAF3IP2_AS1, "TRAF3IP2_AS1", "no")

#adding DHCR24 gene variant samples
DHCR24 <- c('SRR12926202', 'SRR12926204', 'SRR12926207', 'SRR12926208',
            'SRR12926209', 'SRR12926212', 'SRR12926213', 'SRR12926214',
            'SRR12926215', 'SRR12926216', 'SRR12926218', 'SRR12926220',
            'SRR12926223', 'SRR12926224', 'SRR12926225', 'SRR12926226',
            'SRR12926227', 'SRR12926229', 'SRR12926230', 'SRR12926232',
            'SRR12926233', 'SRR12926234', 'SRR12926235', 'SRR12926236',
            'SRR12926237', 'SRR12926239', 'SRR12926240', 'SRR12926242',
            'SRR12926243', 'SRR12926244', 'SRR12926246', 'SRR12926247',
            'SRR12926249', 'SRR12926250', 'SRR12926251', 'SRR12926253',
            'SRR12926255', 'SRR12926256', 'SRR12926258', 'SRR12926262',
            'SRR12926264', 'SRR12926266', 'SRR12926267', 'SRR12926268',
            'SRR12926269', 'SRR12926271', 'SRR12926273', 'SRR12926275',
            'SRR12926276', 'SRR12926277', 'SRR12926278', 'SRR12926280',
            'SRR12926281', 'SRR12926282', 'SRR12926283', 'SRR12926284',
            'SRR12926285', 'SRR12926287', 'SRR12926288', 'SRR12926289',
            'SRR12926290', 'SRR12926291', 'SRR12926292', 'SRR12926294',
            'SRR12926296', 'SRR12926298', 'SRR12926303', 'SRR12926305',
            'SRR12926307', 'SRR12926312', 'SRR12926324', 'SRR12926325',
            'SRR12926328', 'SRR12926335', 'SRR12926337', 'SRR12926338',
            'SRR12926339', 'SRR12926340', 'SRR12926343', 'SRR12926347',
            'SRR12926349', 'SRR12926351', 'SRR12926352', 'SRR12926354',
            'SRR12926355', 'SRR12926357', 'SRR12926359', 'SRR12926360',
            'SRR12926361', 'SRR12926362', 'SRR12926363', 'SRR12926365',
            'SRR12926366')

tdf$DHCR24 <- ifelse(tdf$Run %in% DHCR24, "DHCR24", "no")

#adding CST7 gene variant samples
CST7 <- c('SRR12926202', 'SRR12926207', 'SRR12926208', 'SRR12926209', 'SRR12926212', 'SRR12926214', 'SRR12926218', 'SRR12926220', 'SRR12926225', 'SRR12926229', 'SRR12926230', 'SRR12926234', 'SRR12926235', 'SRR12926239', 'SRR12926240', 'SRR12926242', 'SRR12926244', 'SRR12926249', 'SRR12926253', 'SRR12926256', 'SRR12926257', 'SRR12926258', 'SRR12926260', 'SRR12926262', 'SRR12926264', 'SRR12926266', 'SRR12926267', 'SRR12926269', 'SRR12926271', 'SRR12926273', 'SRR12926275', 'SRR12926276', 'SRR12926278', 'SRR12926281', 'SRR12926283', 'SRR12926285', 'SRR12926287', 'SRR12926291', 'SRR12926294', 'SRR12926316', 'SRR12926323', 'SRR12926335', 'SRR12926337', 'SRR12926339', 'SRR12926340', 'SRR12926341', 'SRR12926342', 'SRR12926343', 'SRR12926345', 'SRR12926347', 'SRR12926349', 'SRR12926351', 'SRR12926352', 'SRR12926353', 'SRR12926354', 'SRR12926355', 'SRR12926356', 'SRR12926357', 'SRR12926358', 'SRR12926359', 'SRR12926361', 'SRR12926363', 'SRR12926364', 'SRR12926365', 'SRR12926366')
tdf$CST7 <- ifelse(tdf$Run %in% CST7, "CST7", "no")

HLA_DOA <- c('SRR12926206', 'SRR12926213', 'SRR12926219', 'SRR12926221', 'SRR12926226', 'SRR12926232', 'SRR12926236', 'SRR12926238', 'SRR12926245', 'SRR12926247', 'SRR12926252', 'SRR12926257', 'SRR12926265', 'SRR12926268', 'SRR12926270', 'SRR12926272', 'SRR12926274', 'SRR12926279', 'SRR12926286', 'SRR12926298', 'SRR12926299', 'SRR12926301', 'SRR12926307', 'SRR12926309', 'SRR12926312', 'SRR12926314', 'SRR12926315', 'SRR12926316', 'SRR12926319', 'SRR12926320', 'SRR12926321', 'SRR12926323', 'SRR12926324', 'SRR12926327', 'SRR12926330', 'SRR12926333', 'SRR12926340', 'SRR12926344', 'SRR12926348', 'SRR12926353', 'SRR12926364')
tdf$HLA_DOA <- ifelse(tdf$Run %in% HLA_DOA, "HLA_DOA", "no")

PMF1<-c('SRR12926204', 'SRR12926206', 'SRR12926207', 'SRR12926208', 'SRR12926209', 'SRR12926210', 'SRR12926211', 'SRR12926212', 'SRR12926213', 'SRR12926214', 'SRR12926215', 'SRR12926216', 'SRR12926217', 'SRR12926218', 'SRR12926219', 'SRR12926222', 'SRR12926223', 'SRR12926225', 'SRR12926226', 'SRR12926227', 'SRR12926228', 'SRR12926229', 'SRR12926230', 'SRR12926231', 'SRR12926232', 'SRR12926233', 'SRR12926234', 'SRR12926238', 'SRR12926239', 'SRR12926240', 'SRR12926241', 'SRR12926242', 'SRR12926243', 'SRR12926244', 'SRR12926245', 'SRR12926246', 'SRR12926247', 'SRR12926248', 'SRR12926249', 'SRR12926250', 'SRR12926256', 'SRR12926257', 'SRR12926258', 'SRR12926259', 'SRR12926262', 'SRR12926263', 'SRR12926267', 'SRR12926268', 'SRR12926269', 'SRR12926270', 'SRR12926271', 'SRR12926272', 'SRR12926273', 'SRR12926274', 'SRR12926276', 'SRR12926277', 'SRR12926280', 'SRR12926283', 'SRR12926284', 'SRR12926287', 'SRR12926289', 'SRR12926290', 'SRR12926291', 'SRR12926292', 'SRR12926293', 'SRR12926294', 'SRR12926296', 'SRR12926298', 'SRR12926299', 'SRR12926302', 'SRR12926303', 'SRR12926307', 'SRR12926308', 'SRR12926309', 'SRR12926310', 'SRR12926311', 'SRR12926313', 'SRR12926314', 'SRR12926316', 'SRR12926317', 'SRR12926318', 'SRR12926320', 'SRR12926322', 'SRR12926323', 'SRR12926324', 'SRR12926325', 'SRR12926327', 'SRR12926328', 'SRR12926330', 'SRR12926333', 'SRR12926334', 'SRR12926335', 'SRR12926336', 'SRR12926343', 'SRR12926344', 'SRR12926345', 'SRR12926346', 'SRR12926352', 'SRR12926353', 'SRR12926355', 'SRR12926356', 'SRR12926359', 'SRR12926363', 'SRR12926364', 'SRR12926365', 'SRR12926367')
tdf$PMF1 <- ifelse(tdf$Run %in% PMF1, "PMF1", "no")

HLA_DQB1<-c('SRR12926206', 'SRR12926211', 'SRR12926215', 'SRR12926220', 'SRR12926221', 'SRR12926226', 'SRR12926233', 'SRR12926234', 'SRR12926244', 'SRR12926245', 'SRR12926248', 'SRR12926250', 'SRR12926269', 'SRR12926270', 'SRR12926271', 'SRR12926272', 'SRR12926277', 'SRR12926282', 'SRR12926284', 'SRR12926311', 'SRR12926313', 'SRR12926317', 'SRR12926361', 'SRR12926362', 'SRR12926367')
tdf$HLA_DQB1 <- ifelse(tdf$Run %in% HLA_DQB1, "HLA_DQB1", "no")

SLC26A2<-c('SRR12926203', 'SRR12926204', 'SRR12926207', 'SRR12926208', 'SRR12926211', 'SRR12926212', 'SRR12926214', 'SRR12926215', 'SRR12926216', 'SRR12926217', 'SRR12926218', 'SRR12926221', 'SRR12926224', 'SRR12926226', 'SRR12926227', 'SRR12926232', 'SRR12926233', 'SRR12926236', 'SRR12926237', 'SRR12926240', 'SRR12926244', 'SRR12926247', 'SRR12926248', 'SRR12926249', 'SRR12926253', 'SRR12926254', 'SRR12926255', 'SRR12926258', 'SRR12926266', 'SRR12926267', 'SRR12926273', 'SRR12926274', 'SRR12926275', 'SRR12926278', 'SRR12926281', 'SRR12926283', 'SRR12926285', 'SRR12926288', 'SRR12926291', 'SRR12926294', 'SRR12926298', 'SRR12926299', 'SRR12926303', 'SRR12926307', 'SRR12926308', 'SRR12926311', 'SRR12926312', 'SRR12926316', 'SRR12926321', 'SRR12926322', 'SRR12926326', 'SRR12926333', 'SRR12926338', 'SRR12926341', 'SRR12926342', 'SRR12926344', 'SRR12926346', 'SRR12926349', 'SRR12926356', 'SRR12926357', 'SRR12926358', 'SRR12926359', 'SRR12926360', 'SRR12926362', 'SRR12926364', 'SRR12926365')
tdf$SLC26A2<- ifelse(tdf$Run %in% SLC26A2, "SLC26A2", "no")

MAVS<-c('SRR12926202', 'SRR12926206', 'SRR12926238', 'SRR12926239', 'SRR12926251', 'SRR12926272', 'SRR12926273', 'SRR12926282', 'SRR12926283', 'SRR12926297', 'SRR12926337', 'SRR12926349', 'SRR12926354', 'SRR12926361')
tdf$MAVS<- ifelse(tdf$Run %in% MAVS, "MAVS", "no")

NSUN3<-c('SRR12926214', 'SRR12926215', 'SRR12926229', 'SRR12926255', 'SRR12926264', 'SRR12926290', 'SRR12926291', 'SRR12926352', 'SRR12926359')
tdf$NSUN3<- ifelse(tdf$Run %in% NSUN3, "NSUN3", "no")

SFT2D2<-c('SRR12926202', 'SRR12926204', 'SRR12926206', 'SRR12926214', 'SRR12926215', 'SRR12926220', 'SRR12926223', 'SRR12926224', 'SRR12926226', 'SRR12926227', 'SRR12926230', 'SRR12926234', 'SRR12926235', 'SRR12926238', 'SRR12926239', 'SRR12926240', 'SRR12926242', 'SRR12926243', 'SRR12926244', 'SRR12926251', 'SRR12926255', 'SRR12926256', 'SRR12926259', 'SRR12926260', 'SRR12926261', 'SRR12926262', 'SRR12926267', 'SRR12926268', 'SRR12926269', 'SRR12926271', 'SRR12926272', 'SRR12926273', 'SRR12926274', 'SRR12926275', 'SRR12926278', 'SRR12926279', 'SRR12926281', 'SRR12926284', 'SRR12926286', 'SRR12926287', 'SRR12926289', 'SRR12926291', 'SRR12926297', 'SRR12926299', 'SRR12926300', 'SRR12926304', 'SRR12926309', 'SRR12926311', 'SRR12926312', 'SRR12926315', 'SRR12926320', 'SRR12926325', 'SRR12926326', 'SRR12926332', 'SRR12926335', 'SRR12926337', 'SRR12926339', 'SRR12926343', 'SRR12926349', 'SRR12926351', 'SRR12926352', 'SRR12926355', 'SRR12926356', 'SRR12926357', 'SRR12926360', 'SRR12926364', 'SRR12926366', 'SRR12926367')
tdf$SFT2D2<- ifelse(tdf$Run %in% SFT2D2, "SFT2D2", "no")

GPR65<-c('SRR12926202', 'SRR12926209', 'SRR12926215', 'SRR12926220', 'SRR12926224', 'SRR12926225', 'SRR12926226', 'SRR12926230', 'SRR12926232', 'SRR12926235', 'SRR12926243', 'SRR12926244', 'SRR12926246', 'SRR12926247', 'SRR12926248', 'SRR12926249', 'SRR12926251', 'SRR12926253', 'SRR12926264', 'SRR12926269', 'SRR12926271', 'SRR12926272', 'SRR12926273', 'SRR12926275', 'SRR12926278', 'SRR12926283', 'SRR12926284', 'SRR12926286', 'SRR12926288', 'SRR12926289', 'SRR12926297', 'SRR12926300', 'SRR12926301', 'SRR12926304', 'SRR12926308', 'SRR12926310', 'SRR12926313', 'SRR12926337', 'SRR12926339', 'SRR12926340', 'SRR12926341', 'SRR12926343', 'SRR12926349', 'SRR12926350', 'SRR12926351', 'SRR12926352', 'SRR12926354', 'SRR12926355', 'SRR12926358', 'SRR12926361', 'SRR12926362', 'SRR12926364', 'SRR12926365')
tdf$GPR65<- ifelse(tdf$Run %in% GPR65, "GPR65", "no")

EIF2AK2<-c('SRR12926202', 'SRR12926233', 'SRR12926244', 'SRR12926251', 'SRR12926337', 'SRR12926339', 'SRR12926343', 'SRR12926351')
tdf$EIF2AK2<- ifelse(tdf$Run %in% EIF2AK2, "EIF2AK2", "no")

TTC9C<-c('SRR12926214', 'SRR12926215', 'SRR12926223', 'SRR12926226', 'SRR12926230', 'SRR12926235', 'SRR12926246', 'SRR12926248', 'SRR12926249', 'SRR12926255', 'SRR12926256', 'SRR12926258', 'SRR12926265', 'SRR12926266', 'SRR12926269', 'SRR12926270', 'SRR12926277', 'SRR12926278', 'SRR12926280', 'SRR12926284', 'SRR12926285', 'SRR12926294', 'SRR12926298', 'SRR12926321', 'SRR12926339', 'SRR12926349', 'SRR12926350', 'SRR12926351', 'SRR12926354', 'SRR12926358', 'SRR12926359', 'SRR12926361')
tdf$TTC9C<- ifelse(tdf$Run %in% TTC9C, "TTC9C", "no")

SYNJ1<-c('SRR12926209', 'SRR12926215', 'SRR12926220', 'SRR12926229', 'SRR12926230', 'SRR12926244', 'SRR12926246', 'SRR12926249', 'SRR12926260', 'SRR12926261', 'SRR12926264', 'SRR12926267', 'SRR12926269', 'SRR12926273', 'SRR12926275', 'SRR12926281', 'SRR12926282', 'SRR12926285', 'SRR12926290', 'SRR12926291', 'SRR12926309', 'SRR12926312', 'SRR12926323', 'SRR12926341', 'SRR12926342', 'SRR12926343', 'SRR12926349', 'SRR12926350', 'SRR12926351', 'SRR12926352', 'SRR12926354', 'SRR12926359', 'SRR12926361')
tdf$SYNJ1<- ifelse(tdf$Run %in% SYNJ1, "SYNJ1", "no")

ADAMTS2<-c('SRR12926202', 'SRR12926204', 'SRR12926208', 'SRR12926218', 'SRR12926224', 'SRR12926226', 'SRR12926230', 'SRR12926234', 'SRR12926235', 'SRR12926236', 'SRR12926237', 'SRR12926239', 'SRR12926242', 'SRR12926249', 'SRR12926256', 'SRR12926281', 'SRR12926283', 'SRR12926289', 'SRR12926337', 'SRR12926339', 'SRR12926343', 'SRR12926345', 'SRR12926349', 'SRR12926355', 'SRR12926361')
tdf$ADAMTS2<- ifelse(tdf$Run %in% ADAMTS2, "ADAMTS2", "no")

MAVS<-c('SRR12926202', 'SRR12926206', 'SRR12926238', 'SRR12926239', 'SRR12926251', 'SRR12926272', 'SRR12926273', 'SRR12926282', 'SRR12926283', 'SRR12926297', 'SRR12926337', 'SRR12926349', 'SRR12926354', 'SRR12926361')

# Box plots of cytokines by group with coloured points of variant samples:
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = MAVS), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("MAVS" = "red", "no" = "black"))

# TRAF3IP2_AS1
# Box plots of cytokines by group with coloured points of variant samples:
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = TRAF3IP2_AS1), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("TRAF3IP2_AS1" = "red", "no" = "black"))

# ATP6AP1
# OK
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ATP6AP1), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("ATP6AP1" = "red", "no" = "black"))

ggplot(tdf, aes(x = Time_point, y = log2(LPS_IL_1b))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ADAMTS2), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_1b by Group",
       x = "Group",
       y = "LPS_IL_1b") +
  scale_color_manual(values = c("ADAMTS2" = "red", "no" = "black"))

ggplot(tdf, aes(x = Time_point, y = log2(LPS_IL_1b))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = SYNJ1), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_1b by Group",
       x = "Group",
       y = "LPS_IL_1b") +
  scale_color_manual(values = c("SYNJ1" = "red", "no" = "black"))

ggplot(tdf, aes(x = Time_point, y = log2(LPS_TNF_alpha))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ATP6AP1), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_TNF_alpha by Group",
       x = "Group",
       y = "LPS_IL_1b") +
  scale_color_manual(values = c("ATP6AP1" = "red", "no" = "black"))

#ADAMTS2
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ADAMTS2), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("ADAMTS2" = "red", "no" = "black"))

#HSPA1L
ggplot(tdf, aes(x = Time_point, y = log2(LPS_TNF_alpha))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = HSPA1L), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_TNF_alpha by Group",
       x = "Group",
       y = "LPS_TNF_alpha") +
  scale_color_manual(values = c("HSPA1L" = "red", "no" = "black"))

#FAM160B1
ggplot(tdf, aes(x = Time_point, y = log2(LPS_TNF_alpha))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = FAM160B1), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_TNF_alpha by Group",
       x = "Group",
       y = "LPS_TNF_alpha") +
  scale_color_manual(values = c("FAM160B1" = "red", "no" = "black"))

ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = FAM160B1), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("FAM160B1" = "red", "no" = "black"))

#RICTOR
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = RICTOR), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("RICTOR" = "red", "no" = "black"))

ggplot(tdf, aes(x = Time_point, y = log2(LPS_IL_1b))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = RICTOR), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_1b by Group",
       x = "Group",
       y = "LPS_IL_1b") +
  scale_color_manual(values = c("RICTOR" = "red", "no" = "black"))

# LILRA2
# OK
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = LILRA2), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("LILRA2" = "red", "no" = "black"))

# ADGRE3
# OK but too few values
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ADGRE3), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS_IL_10 by Group",
       x = "Group",
       y = "LPS_IL_10") +
  scale_color_manual(values = c("ADGRE3" = "red", "no" = "black"))

# HPSE
# OK
ggplot(tdf, aes(x = Time_point, y = log2(LPS_TNF_alpha))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = HPSE), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of TNF-alpha by Group",
       x = "Group",
       y = "log 2 TNF-alpha") +
  scale_color_manual(name = "RNA-edit status",
                     labels = c("HPSE RNA-edit present", "no edit"),
                     values = c("HPSE" = "red", "no" = "grey")) +
  scale_fill_discrete(name = "Time point")


# LILRA2
# OK
ggplot(tdf, aes(x = Time_point, y = LPS_IL_10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = LILRA2), position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot of LPS IL-10 by Group",
       x = "Group",
       y = "LPS IL-10") +
  scale_color_manual(name = "RNA-edit status",
                     labels = c("LILRA2 RNA-edit present", "no edit"),
                     values = c("LILRA2" = "red", "no" = "grey")) +
  scale_fill_discrete(name = "Time point")

# plot IL-10 levels only in acute samples and stratified by variant.
# annotate variant presence in all samples
# List of variants:
# HPSE
# LILRA2
# TRAF3IP2-AS1 / TRAF3IP2_AS1
# VHL 
# MARCHF1
# RICTOR
# ADGRE3

tdf$HPSE_group <- ifelse(tdf$Time_point == "acute" & tdf$HPSE == "HPSE", "Acute with HPSE edit",
                    ifelse(tdf$Time_point == "acute" & tdf$HPSE == "no", "Acute without edit",
                           ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$LILRA2_group <- ifelse(tdf$Time_point == "acute" & tdf$LILRA2 == "LILRA2", "Acute with LILRA2 edit",
                         ifelse(tdf$Time_point == "acute" & tdf$LILRA2 == "no", "Acute without edit",
                                ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$TRAF3IP2_group <- ifelse(tdf$Time_point == "acute" & tdf$TRAF3IP2_AS1 == "TRAF3IP2_AS1", "Acute with TRAF3IP2_AS1 edit",
                           ifelse(tdf$Time_point == "acute" & tdf$TRAF3IP2_AS1 == "no", "Acute without edit",
                                  ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$VHL_group <- ifelse(tdf$Time_point == "acute" & tdf$VHL == "VHL", "Acute with VHL edit",
                             ifelse(tdf$Time_point == "acute" & tdf$VHL == "no", "Acute without edit",
                                    ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$MARCHF1_group <- ifelse(tdf$Time_point == "acute" & tdf$MARCHF1 == "MARCHF1", "Acute with MARCHF1 edit",
                        ifelse(tdf$Time_point == "acute" & tdf$MARCHF1 == "no", "Acute without edit",
                               ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$RICTOR_group <- ifelse(tdf$Time_point == "acute" & tdf$RICTOR == "RICTOR", "Acute with RICTOR edit",
                            ifelse(tdf$Time_point == "acute" & tdf$RICTOR == "no", "Acute without edit",
                                   ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$ADGRE3_group <- ifelse(tdf$Time_point == "acute" & tdf$ADGRE3 == "ADGRE3", "Acute with ADGRE3 edit",
                           ifelse(tdf$Time_point == "acute" & tdf$ADGRE3 == "no", "Acute without edit",
                                  ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf$ADAMTS2_group <- ifelse(tdf$Time_point == "acute" & tdf$ADAMTS2 == "ADAMTS2", "Acute with ADAMTS2 edit",
                           ifelse(tdf$Time_point == "acute" & tdf$ADAMTS2 == "no", "Acute without edit",
                                  ifelse(tdf$Time_point == "control", "Controls", "Recovery")))

tdf_lastvariants_acute <- tdf %>% 
  filter(Time_point == "acute")

# HPSE plots and statistical analysis
ggplot(tdf_lastvariants_acute, aes(x = HPSE_group, y = LPS_IL_10, color = HPSE_group)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5) +
  labs(title = "Il-10 release across groups", y = "Il-10", color = "HPSE_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 2000)) +
  guides(color = FALSE)

acute_HPSE_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), HPSE) %>%
  filter(HPSE == "HPSE")

acute_no_HPSE_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), HPSE) %>%
  filter(HPSE == "no")

wilcox.test(acute_HPSE_cytokines$LPS_IL_10, acute_no_HPSE_cytokines$LPS_IL_10)
            
ggplot(tdf_lastvariants_acute, aes(x = HPSE_group, y = LPS_TNF_alpha, color = HPSE_group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "TNF-alpha release across groups", y = "TNF-alpha", color = "HPSE_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 20000)) + 
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

tdf$log2LPS_TNF_alpha <- log2(tdf$LPS_TNF_alpha)

ggplot(tdf_lastvariants_acute, aes(x = HPSE_group, y = LPS_IL_1b, color = HPSE_group)) +
  geom_boxplot() +
  labs(title = "Il-1b release across groups", y = "Il-1b", color = "HPSE_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 20000)) +
  guides(color = FALSE)+
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

wilcox.test(acute_HPSE_cytokines$LPS_IL_1b, acute_no_HPSE_cytokines$LPS_IL_1b)
dunn.test(tdf$LPS_IL_1b, tdf$HPSE_group, method = "bonferroni")


# LILRA2 plots and statistical analysis
ggplot(tdf_lastvariants_acute, aes(x = LILRA2_group, y = LPS_IL_10, color = LILRA2_group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Il-10 release across groups", y = "Il-10", color = "LILRA2_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 2000)) +
  guides(color = FALSE)+
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

acute_LILRA2_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), LILRA2) %>%
  filter(LILRA2 == "LILRA2")

acute_no_LILRA2_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), LILRA2) %>%
  filter(LILRA2 == "no")

wilcox.test(acute_LILRA2_cytokines$LPS_IL_10,
            acute_no_LILRA2_cytokines$LPS_IL_10)


# VHL plots and statistical analysis
ggplot(tdf_lastvariants_acute, aes(x = VHL_group, y = LPS_IL_10, color = VHL_group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Il-10 release across groups", y = "Il-10", color = "VHL_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 2000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

ggplot(tdf_lastvariants_acute, aes(x = VHL_group, y = LPS_IL_6, color = VHL_group)) +
  geom_boxplot() +
  labs(title = "Il-6 release across groups", y = "Il-6", color = "VHL_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(-10000, 100000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

ggplot(tdf_lastvariants_acute, aes(x = VHL_group, y = LPS_TNF_alpha, color = VHL_group)) +
  geom_boxplot() +
  labs(title = "TNF-alpha release across groups", y = "TNF-alpha", color = "VHL_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 20000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

# last one I hope:
ggplot(tdf_lastvariants_acute, aes(x = ADAMTS2_group, y = LPS_IL_1b, color = ADAMTS2_group)) +
  geom_boxplot() +
  labs(title = "Il-1b release across groups", y = "Il-1b", color = "ADAMTS2_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 20000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

acute_VHL_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), VHL) %>%
  filter(VHL == "VHL")

acute_no_VHL_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), VHL) %>%
  filter(VHL == "no")

wilcox.test(acute_VHL_cytokines$LPS_IL_10,
            acute_no_VHL_cytokines$LPS_IL_10)

# TRAF3IP2 plots and statistical analysis

ggplot(tdf_lastvariants_acute, aes(x = TRAF3IP2_group, y = LPS_IL_10, color = TRAF3IP2_group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Il-10 release across groups", y = "Il-10", color = "TRAF3IP2_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 2000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

ggplot(tdf_lastvariants_acute, aes(x = TRAF3IP2_group, y = LPS_IL_1b, color = TRAF3IP2_group)) +
  geom_boxplot() +
  labs(title = "Il-1b release across groups", y = "Il-1b", color = "TRAF3IP2_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 20000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

ggplot(tdf_lastvariants_acute, aes(x = TRAF3IP2_group, y = LPS_TNF_alpha, color = TRAF3IP2_group)) +
  geom_boxplot() +
  labs(title = "TNF-alpha release across groups", y = "TNF-alpha", color = "TRAF3IP2_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(-1000, 20000)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

acute_TRAF3IP2_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), TRAF3IP2_AS1) %>%
  filter(TRAF3IP2_AS1 == "TRAF3IP2_AS1")

acute_no_TRAF3IP2_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), TRAF3IP2_AS1) %>%
  filter(TRAF3IP2_AS1 == "no")

wilcox.test(acute_TRAF3IP2_cytokines$LPS_IL_1b,
            acute_no_TRAF3IP2_cytokines$LPS_IL_1b)

# MARCHF1 plots and statistical analysis
ggplot(tdf_lastvariants_acute, aes(x = MARCHF1_group, y = LPS_IL_10, color = MARCHF1_group)) +
  geom_boxplot() +
  labs(title = "Il-10 release across groups", y = "Il-10", color = "MARCHF1_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)

acute_MARCHF1_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), MARCHF1) %>%
  filter(MARCHF1 == "MARCHF1")

acute_no_MARCHF1_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), MARCHF1) %>%
  filter(MARCHF1 == "no")

wilcox.test(acute_MARCHF1_cytokines$LPS_IL_10,
            acute_no_MARCHF1_cytokines$LPS_IL_10)

# RICTOR plots and statistical analysis
ggplot(tdf_lastvariants_acute, aes(x = RICTOR_group, y = LPS_TNF_alpha, color = RICTOR_group)) +
  geom_boxplot() +
  labs(title = "TNF_alpharelease across groups", y = "TNF_alpha", color = "RICTOR_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5) +
  coord_cartesian(ylim = c(-10000, 75000))

acute_RICTOR_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), RICTOR) %>%
  filter(RICTOR == "RICTOR")

acute_no_RICTOR_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), RICTOR) %>%
  filter(RICTOR == "no")

wilcox.test(acute_RICTOR_cytokines$LPS_IL_10,
            acute_no_RICTOR_cytokines$LPS_IL_10)

acute_no_RICTOR_cytokines$LPS_IL_1b

# Trying a permutation test
# Remove NA values
acute_RICTOR_cytokines <- na.omit(acute_RICTOR_cytokines)
acute_no_RICTOR_cytokines <- na.omit(acute_no_RICTOR_cytokines)

# Combine the two groups into one vector
combined <- c(acute_RICTOR_cytokines$LPS_IL_6, acute_no_RICTOR_cytokines$LPS_IL_6)

# Create a grouping variable
group <- c(rep("RICTOR", length(acute_RICTOR_cytokines$LPS_IL_6)), 
           rep("no_RICTOR", length(acute_no_RICTOR_cytokines$LPS_IL_6)))
group <- as.factor(group)

# Ensure the data is numeric
combined <- as.numeric(combined)

# Perform the permutation test

test_result <- oneway_test(combined ~ group, distribution = approximate(nresample = 9999))

# Print the test result
print(test_result)

# ADGRE3 plots and statistical analysis
ggplot(tdf_lastvariants_acute, aes(x = ADGRE3_group, y = LPS_IL_10, color = ADGRE3_group)) +
  geom_boxplot() +
  labs(title = "LPS_IL_10 release across groups", y = "LPS_IL_10", color = "ADGRE3_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5)
  #coord_cartesian(ylim = c(-10000, 75000))

ggplot(tdf_lastvariants_acute, aes(x = ADGRE3_group, y = LPS_IL_1b, color = ADGRE3_group)) +
  geom_boxplot() +
  labs(title = "LPS_IL_1b release across groups", y = "LPS_IL_1b", color = "ADGRE3_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 150000))

ggplot(tdf_lastvariants_acute, aes(x = ADGRE3_group, y = LPS_IL_6, color = ADGRE3_group)) +
  geom_boxplot() +
  labs(title = "LPS_IL_6 release across groups", y = "LPS_IL_6", color = "ADGRE3_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 150000))

ggplot(tdf_lastvariants_acute, aes(x = ADGRE3_group, y = LPS_TNF_alpha, color = ADGRE3_group)) +
  geom_boxplot() +
  labs(title = "LPS_TNF_alpha release across groups", y = "LPS_TNF_alpha", color = "ADGRE3_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE) +
  stat_summary(fun = median, geom = "text", aes(label = ..y..), vjust = -0.5) +
  coord_cartesian(ylim = c(-1000, 70000))

acute_ADGRE3_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), ADGRE3) %>%
  filter(ADGRE3 == "ADGRE3")

acute_no_ADGRE3_cytokines <- tdf %>%
  filter(Time_point == "acute") %>% 
  select(starts_with("LPS"), ADGRE3) %>%
  filter(ADGRE3 == "no")

wilcox.test(acute_ADGRE3_cytokines$LPS_TNF_alpha,
            acute_no_ADGRE3_cytokines$LPS_TNF_alpha)

# Trying a permutation test
# Remove NA values
acute_ADGRE3_cytokines <- na.omit(acute_ADGRE3_cytokines)
acute_no_ADGRE3_cytokines <- na.omit(acute_no_ADGRE3_cytokines)

# Combine the two groups into one vector
combined <- c(acute_ADGRE3_cytokines$LPS_IL_10, acute_no_ADGRE3_cytokines$LPS_IL_10)

# Create a grouping variable
group <- c(rep("ADGRE3", length(acute_ADGRE3_cytokines$LPS_IL_10)), 
           rep("no_ADGRE3", length(acute_no_ADGRE3_cytokines$LPS_IL_10)))
group <- as.factor(group)

# Ensure the data is numeric
combined <- as.numeric(combined)

# Perform the permutation test

test_result <- oneway_test(combined ~ group, distribution = approximate(nresample = 9999))

# Print the test result
print(test_result)

# adjusting p-values
exonic_a_to_i_for_volcano_plot$Adjusted_P_Value <- p.adjust(exonic_a_to_i_for_volcano_plot$p_value_fisher,
                                                            method = "BH")
exonic_a_to_i_for_volcano_plot$Adjusted_P_Value_neglog <- -log10(exonic_a_to_i_for_volcano_plot$Adjusted_P_Value)
output_volcanoatoiexonic <- "adjusted_p_atoi_exonic.csv"

# Export tables as tab-delimited text files.
write.table(exonic_a_to_i_for_volcano_plot, file = output_volcanoatoiexonic, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)



