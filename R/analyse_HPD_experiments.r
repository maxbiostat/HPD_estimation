expfs <- system("ls HPD_experiments/*.csv", intern = TRUE)

results <- tibble::tibble(do.call(rbind, lapply(expfs, read.csv)))

library(ggplot2)

cov_p <- ggplot(results,
                aes(
                  x = method,
                  y = est_coverage,
                  colour = method,
                  shape = quantity
                )) +
  geom_point(size = 3.5) +
  facet_grid(asymmetry ~ efficiency, scales = "free") +
  geom_hline(yintercept = 0.95,
             linetype = "longdash") +
  theme_bw()
cov_p


cov_pL <- ggplot(subset(results, quantity == "HPD_L"),
                aes(
                  x = method,
                  y = est_coverage,
                  colour = method
                )) +
  geom_point(size = 3.5) +
  facet_grid(asymmetry ~ efficiency, scales = "free") +
  geom_hline(yintercept = 0.95,
             linetype = "longdash") +
  ggtitle("HPD lower end") + 
  theme_bw()
cov_pL

cov_pU <- ggplot(subset(results, quantity == "HPD_U"),
                aes(
                  x = method,
                  y = est_coverage,
                  colour = method,
                  fill = method
                )) +
  geom_point(size = 3.5, shape = 24) +
  facet_grid(asymmetry ~ efficiency, scales = "free") +
  geom_hline(yintercept = 0.95,
             linetype = "longdash") +
  ggtitle("HPD upper end") + 
  theme_bw()
cov_pU


bias_p <- ggplot(results,
                 aes(
                   x = method,
                   y = abs(bias),
                   colour = method,
                   shape = quantity
                 )) +
  geom_point() +
  facet_grid(asymmetry ~ efficiency, scales = "free") +
  theme_bw()
bias_p



# ggsave(
#   plot = pp,
#   filename =
#     paste0("../figures/BCI_experiment_",
#            exper,
#            ".pdf"),
#   scale = 1,
#   width = 297,
#   height = 210,
#   units = "mm",
#   dpi = 300
# )
