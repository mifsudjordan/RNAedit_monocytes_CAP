# New script from the moment I created the large variant dataframe.

ctou[ctou == "none"] <- 0
ctou[ctou == "CtoU"] <- 1
ctou <- t(ctou)

new_ctou <- ctou[-1, ]
colnames(new_ctou) <- ctou[1, ]
rownames(new_ctou) <- thesis_datav4$Run
ctou <- new_ctou
rm(new_ctou)

# checking normalised counts of PPM1B
PPM1b_gene_counts <- normalized_counts["ENSG00000138032", ]
PPM1b_gene_counts <- t(t(PPM1b_gene_counts))

# loading variant counts after more stringent filtration
new_var_counts <- as.data.frame(read_csv("stats_table_all_new.csv"))

# box plots of variant counts by Time_point
ggplot(new_var_counts, aes(x = Time_point, y = CtoU, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot and Scatter Plot of C to U by Time_point",
       x = "Time_point",
       y = "count of C to U variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(new_var_counts, aes(x = Time_point, y = AtoI_TOTALmatches, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot and Scatter Plot of A to I matched to ADAR db by Time_point",
       x = "Time_point",
       y = "count of A to I variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(new_var_counts, aes(x = Time_point, y = AtoI, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot and Scatter Plot of A to I by Time_point",
       x = "Time_point",
       y = "count of A to I variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hist(new_var_counts$AtoI_TOTALmatches)

#Negative binomial test results
nb_test_result_CtoU_new <- glm.nb(CtoU ~ Time_point, data = new_var_counts)
nb_test_result_AtoI_new <- glm.nb(AtoI ~ Time_point, data = new_var_counts)
summary(nb_test_result_CtoU_new)
summary(nb_test_result_AtoI_new)
