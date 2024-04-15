library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

# Read in data
sp <- readRDS("./Data/Species numbers for plotting.RDS")
sp <- sp[!sp %in% c(6, 9, 26, 23)]
length(unique(sp))

sp.info <- readRDS("./data/Species info for plotting.RDS")

## Remove 4 species with too few observations 
sp.info <- subset(sp.info, !sp.info$sp.name %in% c("TSUCAN", "QUERUB", "CORALT", "AMELAE"))
n.sp <- length(unique(sp.info$sp.name))

load("./data/Average posterior comparison/common exponents/ForestGEO_Wabikon_Legacy_CDD_LogisticSurvival_cdd_matrices_20231110.RData")

# Calculations
realLivingCDD.spAPC <- apply(realLivingCDD.sp, 2, mean)
realLegacyCDD.spAPC <- apply(realLegacyCDD.sp, 2, mean)

sp.w <- tapply(sp, sp, length)

realLivingCDD.plot <- apply(realLivingCDD.sp, 1, weighted.mean, w = sp.w)
realLegacyCDD.plot <- apply(realLegacyCDD.sp, 1, weighted.mean, w = sp.w)

# Data frame creation
df1 <- data.frame(
  group = c("Legacy CDD", "Living CDD"),
  lower = c(quantile(realLegacyCDD.plot, 0.025), quantile(realLivingCDD.plot, 0.025)),
  upper = c(quantile(realLegacyCDD.plot, 0.975), quantile(realLivingCDD.plot, 0.975)),
  lower2 = c(quantile(realLegacyCDD.plot, 0.25), quantile(realLivingCDD.plot, 0.25)),
  upper2 = c(quantile(realLegacyCDD.plot, 0.75), quantile(realLivingCDD.plot, 0.75)),
  middle = c(quantile(realLegacyCDD.plot, 0.5), quantile(realLivingCDD.plot, 0.5)),
  color = c("red", "blue")
)

# Plot using ggplot
fig2 <- ggplot(df1, aes(x = middle, y = group)) +
  geom_point(aes(color = color), shape = 15, size = 4) +
  geom_segment(aes(xend = lower, x = upper, yend = group, y = group, color = color), linewidth = 1.5) +
  geom_segment(aes(xend = lower2, x = upper2, yend = group, y = group, color = color), linewidth = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  labs(x = "Change in annual survival probability", y = NULL, title = "") +
  scale_color_manual(values = c("#FF5733", "#3399FF")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ggsave(fig2, file = "./Figures/figure2.jpg", width = 6, height = 4)

# ... Previous code blocks ...

# Species level estimates
i = 1
df <- data.frame(
  group = c("Legacy CDD", "Living CDD"),
  lower = c(quantile(realLegacyCDD.sp[,i], 0.025), quantile(realLivingCDD.sp[,i], 0.025)),
  upper = c(quantile(realLegacyCDD.sp[,i], 0.975), quantile(realLivingCDD.sp[,i], 0.975)),
  lower2 = c(quantile(realLegacyCDD.sp[,i], 0.25), quantile(realLivingCDD.sp[,i], 0.25)),
  upper2 = c(quantile(realLegacyCDD.sp[,i], 0.75), quantile(realLivingCDD.sp[,i], 0.75)),
  middle = c(quantile(realLegacyCDD.sp[,i], 0.5), quantile(realLivingCDD.sp[,i], 0.5)),
  sp = sp.info$sp.name[i]
)

for (i in 2:n.sp) {
  df <- rbind(df,
              data.frame(
                group = c("Legacy CDD", "Living CDD"),
                lower = c(quantile(realLegacyCDD.sp[,i], 0.025), quantile(realLivingCDD.sp[,i], 0.025)),
                upper = c(quantile(realLegacyCDD.sp[,i], 0.975), quantile(realLivingCDD.sp[,i], 0.975)),
                lower2 = c(quantile(realLegacyCDD.sp[,i], 0.25), quantile(realLivingCDD.sp[,i], 0.25)),
                upper2 = c(quantile(realLegacyCDD.sp[,i], 0.75), quantile(realLivingCDD.sp[,i], 0.75)),
                middle = c(quantile(realLegacyCDD.sp[,i], 0.5), quantile(realLivingCDD.sp[,i], 0.5)),
                sp = sp.info$sp.name[i])
  )
}

spt <- fread("./data/sps_wabikon.txt", header = T)

df <- data.table(merge(df, spt, by = 'sp', all = F))
df$latin <- as.factor(df$latin)
df$latin <- reorder(df$latin, df$middle)

count_species_below_zero <- df %>%
  filter(upper < 0) %>%
  group_by(group) %>%
  summarise(count_species = n_distinct(sp))

ran_eff_plot_sp <- ggplot(df, aes(x = middle, y = latin, group = group, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1, color = "gray25") +
  geom_point(shape = 16, size = 4) +
  geom_segment(aes(xend = lower, x = upper, yend = latin, y = latin), linewidth = 1.5) +
  geom_segment(aes(xend = lower2, x = upper2, yend = latin, y = latin), linewidth = 2) +
  scale_color_manual(values = c("#FF5733", "#3399FF")) +
  labs(x = "Change in annual survival probability", y = NULL) +
  facet_grid(group ~ .) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )



# ggsave(ran_eff_plot_sp, file = "./Figures/Random slopes CDD.jpg", width = 173, height = 173, units = 'mm')

# beta.m <- readRDS("./data/beta_dot_m.RDS")
beta.m <- read.csv("./data/betas.csv", header = T)

beta.m <- data.table(beta.m)

colnames(beta.m) <- c("Intercept", "DBH", "Cdens", "Hdens", "LCdens", "LHdens")

beta.long <- data.table(beta.m %>%
                          pivot_longer(everything(), names_to = 'beta', values_to = 'est'))
beta.comp <- beta.long[, .(med = median(est), lower = quantile(est, 0.025), lower2 = quantile(est, 0.25),
                           upper = quantile(est, 0.975), upper2 = quantile(est, 0.75)), by = beta]
beta.comp$beta <- factor(beta.comp$beta, levels =  c("LHdens", "LCdens", "Hdens","Cdens",  "DBH","Intercept" ))

fixed_effects_plot <- ggplot(beta.comp, aes(x = med, y = beta)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = lower, x = upper, yend = beta, y = beta), linewidth = 1.5) +
  geom_segment(aes(xend = lower2, x = upper2, yend = beta, y = beta), linewidth = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  labs(x = "Standardized parameter estimates", y = NULL) +
  theme_bw()

# ggsave(filename = "./Figures/fixed effects all.jpg", width = 82, height = 82, units = 'mm')

