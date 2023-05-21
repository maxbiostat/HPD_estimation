library(ggplot2)
library(tidyverse)

results.df <-
  read_csv(file = "../saved_data/MCMC_lognormal_m=0_s=1_results.csv")

target <- "high"

results.df <-
  results.df %>% mutate(rae = abs((point - true) / true))
results.df <- results.df %>% mutate(sqe = (point - true) ^ 2)

###############
MRAE <- aggregate(rae ~ quantity + target_ESS,
                  results.df, mean)
RMSE <-
  aggregate(sqe ~ quantity + target_ESS, results.df, function(x)
    sqrt(mean(x)))
RMSE <- RMSE %>% rename(rmse = sqe)

Var <- aggregate(point ~ quantity + target_ESS, results.df, var)
Var <- Var %>% rename(var = point)

True <- aggregate(true ~ quantity + target_ESS, results.df, mean)

Bias <- aggregate(point ~ quantity + target_ESS, results.df, mean)
Bias <- merge(Bias, True, by = c("quantity", "target_ESS"))
Bias <-  Bias %>% mutate(bias = point - true)
Bias <-  Bias %>% mutate(bias_sq = bias^2)
Bias$point <- NULL
Bias$true <- NULL

RMSE.final <- Reduce(function(x, y)
  merge(x, y, by = c("quantity", "target_ESS")),
  list(RMSE, Bias, Var))

RMSE.final <- RMSE.final %>% mutate(rmse_b = sqrt(bias_sq + var))
RMSE.final

write.csv(MRAE,
          file = paste0("../saved_data/MRAE_MCMC_",
                        target, ".csv"),
          row.names = FALSE)

write.csv(RMSE.final,
          file = paste0("../saved_data/RMSE_MCMC_",
                        target, ".csv"),
          row.names = FALSE)


library(ggplot2)

sampling_dists <-
  ggplot(data = results.df,
         aes(x = point, fill = quantity)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous("", expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  geom_vline(data = results.df,
             aes(xintercept = true), linetype = "longdash") +
  facet_grid(target_ESS ~ quantity, scales = "free") +
  # stat_theodensity(distri = "norm",
  #                  size = .7,  aes(color = "blue")) +
  # stat_theodensity(distri = "lnorm",
  #                  size = .7,  aes(color = "gold")) +
  # scale_color_identity(
  #   name = "Distribution",
  #   breaks = c("blue", "gold"),
  #   labels = c("normal", "log-normal"),
  #   guide = "legend"
  # ) +
  theme_bw(base_size = 16)

sampling_dists

ggsave(
  plot = sampling_dists,
  filename = paste0("../figures/MCMC_simu_estimator_distributions_",
                    target, ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)

mrae <- ggplot(data = MRAE,
       aes(
         x = target_ESS,
         y = rae,
         colour = quantity,
         fill = quantity
       )) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous("Target effective sample size", expand = c(0, 0 )) + 
  scale_y_continuous("Mean absolute relative error (MRAE)",
                     expand = c(0, 0)) +
  facet_grid(quantity ~ ., scales = "free_y") +
  geom_hline(yintercept = .05, linetype = "longdash") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

mrae

ggsave(
  plot = mrae,
  filename = paste0("../figures/MCMC_simu_MRAE_",
                    target, ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
#######
