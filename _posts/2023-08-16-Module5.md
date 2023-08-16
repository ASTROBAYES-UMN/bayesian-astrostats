---
layout: post
title: "Lognormal Distribution"
date: 2023-07-10
categories: modules
tags: test demo pleaseWork
pin: false
---

# Our Data

```

library(tidyverse)
library(modelsummary)
library(FITSio)

pr_dat <- read_csv("Planet_radii.csv")

summary(pr_dat)

ggplot(pr_dat, aes(x = Radius)) + 
  geom_histogram(aes(y =..density..),
                   colour = "black", fill = "white") + 
  geom_density(colour = "blue")

```
![Desktop View](/assets/img/Module-5/fig_1.png){: w="700" h="400" }


## Lognormal Distribution
Distribution of a continuous random variable whose logarithm is normally doistributed.

::: {#def-LN}
Say $X \sim LogNormal (\mu, \lambda)$ if $Y = log X \sim N(\mu, \lambda)$. Then the pdf of $X$ is $$p(x \mid \mu, \lambda) = \frac{1}{x \sqrt{2 \pi \lambda}} e^{- \frac{1}{2 \lambda} (\log x - \mu)^{2}} .$$
:::

Then the mean and varaince of $X$ are $$E[X \mid \mu, \lambda] = e^{\mu + \frac{\lambda}{2}}$$ and $$Var(X \mid \mu, \lambda) = \left[ e^{\lambda} - 1\right] e^{2\mu + \lambda}$$

## What does LogNormal data look like?

```

set.seed(5731)

msim <- 1e3

sim_dat1 <- rlnorm(msim, meanlog = 0, sdlog = sqrt(1))

sim_dat2 <- rlnorm(msim, meanlog = -1, sdlog = sqrt(1))

sim_dat3 <- rlnorm(msim, meanlog = 0, sdlog = sqrt(100))
 
sim_dat4 <- rlnorm(msim, meanlog = 5, sdlog = sqrt(4))

summary(sim_dat1)

signif(t.test(sim_dat1, conf.level = 0.95)conf.int, digits = 4)

summary(sim_dat2)

signif(t.test(sim_dat2, conf.level = 0.95)conf.int, digits = 4)

summary(sim_dat3)

signif(t.test(sim_dat3, conf.level = 0.95)conf.int, digits = 4)

summary(sim_dat4)

signif(t.test(sim_dat4, conf.level = 0.95)conf.int, digits = 4)

```

```
ggplot(as_tibble(sim_dat1), aes(x=value)) + 
  geom_histogram() +
  labs(x = "Simulated values")

ggplot(as_tibble(sim_dat2), aes(x=value)) + 
  geom_histogram() +
  labs(x = "Simulated values")

ggplot(as_tibble(sim_dat3), aes(x=value)) + 
  geom_histogram() + 
  labs( x = "Simulated values")

ggplot(as_tibble(sim_dat4), aes(x=value)) + 
  geom_histogram() +
  labs(x = "Simulated values")
```

![Desktop View](/assets/img/Module-3/fig_2.png){: w="700" h="400" }
_Figure 1: Histograms of 1000 Log Normal Simulations_

