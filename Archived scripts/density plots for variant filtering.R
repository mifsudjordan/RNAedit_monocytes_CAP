# density plots for variant filtering:

library(ggplot2)
library(BiocGenerics)
library(VariantAnnotation)

# Load the VCF file
vcf <- readVcf("sample_haplo.vcf", "hg19")

# Extract the QD field
QD <- info(vcf)$QD

# Create a data frame
qd_df <- data.frame(QD)

# Create the density plot
ggplot(qd_df, aes(x = QD)) + 
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(x = "QD", y = "Density", title = "Density plot of QD field")

# Extract the FS field
FS <- info(vcf)$FS

# Create a data frame
fs_df <- data.frame(FS)
summary(fs_df$FS)
# Create the density plot
ggplot(fs_df, aes(x = FS)) + 
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(x = "FS", y = "Density", title = "Density plot of FS field")

# Extract the SOR field
SOR <- info(vcf)$SOR

# Create a data frame
sor_df <- data.frame(SOR)

# Create the density plot
ggplot(sor_df, aes(x = SOR)) + 
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(x = "SOR", y = "Density", title = "Density plot of SOR field")

# Extract the SOR field
MQ <- info(vcf)$MQ

# Create a data frame
mq_df <- data.frame(MQ)
summary(mq_df)
# Create the density plot
ggplot(mq_df, aes(x = MQ)) + 
  geom_density(fill = "lightblue") +
  labs(x = "MQ", y = "Density", title = "Density plot of MQ field")

# Extract the SOR field
MQRankSum <- info(vcf)$ReadPosRankSum

# Create a data frame
mq_ranksum_df <- data.frame(MQRankSum)
summary(mq_ranksum_df)
# Create the density plot
ggplot(mq_ranksum_df, aes(x = MQRankSum)) + 
  geom_density(fill = "lightblue") +
  labs(x = "MQ", y = "Density", title = "Density plot of MQ field")
