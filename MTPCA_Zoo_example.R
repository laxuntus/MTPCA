# Richard Forsyth 1990. “Zoo.” UCI Machine Learning Repository,
# DOI: https://doi.org/10.24432/C5R59V.
zoo <- read.csv("zoo.csv")

animals = zoo[,-c(1, 14, 18)] # Removing non-binary variables

mod = MTPCA(animals, rep("B", 15), 3, 5)
lambda <- mod$lambda
U <- mod$U
rownames(U) = colnames(animals)
U # The loadings

mu = mod$mu

Z <- gradient_descent_all(animals, rep("B", 15), mu, lambda, U, 1)

library(ggplot2)

# Component scores
Z_scatter = data.frame(PC2 = Z[,2], PC3 = Z[,3], name = zoo$animal_name)

ggplot(Z_scatter, aes(PC2, PC3, label=name)) +
  geom_point() +
  geom_text_repel(aes(label = name), max.overlaps = Inf, size = 9 / .pt) +
  theme_bw()

# Scree plot
mod2 = MTPCA(animals, rep("B", 15), 6, 5)
dat_scree = data.frame(component=seq(1,6), lambda=mod2$lambda)

ggplot(dat_scree, aes(component, lambda)) +
  geom_point(size=4) +
  geom_line() +
  xlab("Principal component") +
  ylab("Variance explained") +
  scale_x_continuous(n.breaks=6) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("Scree plot") +
  theme(plot.title = element_text(hjust = 0.5))