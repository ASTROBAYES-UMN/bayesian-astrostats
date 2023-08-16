---
layout: post
title: "Introduction to MCMC"
date: 2023-07-10
categories: modules
tags: test demo pleaseWork
pin: false
---

# Planet Radii

```


library(tidyverse)
library(modelsummary)
library(FITSio)

pr_dat <- read_csv("Planet_radii.csv")

summary(pr_dat)

ggplot(pr_dat, aes(x = Radius)) + theme_classic() +
  geom_histogram(aes(y =..density..),
                   colour = "black", fill = "white") + 
  geom_density(colour = "blue") + labs(x="Earth Radii")

```
![Desktop View](/assets/img/Module-6/fig_1.png){: w="700" h="400" }


# Statistical Model, Part 1

1.  $\theta$, mean planet radius

2.  $\Theta = [0, \infty)$

3.  $x$ observed data, described above

4.  $\mathcal{X}$, all possible vectors of `r length(pr_dat$Radius)` non-negative scalars

5.  $p(x \mid \theta)$, a description of how the data arise in nature

6.  $\nu(\theta)$, assumptions about the mean palnet radius

## Lognormal Distribution

::: {#def-LN}
Say $X \sim LogNormal (\mu, \lambda)$ if $Y = log X \sim N(\mu, \lambda)$. Then the pdf of $X$ is $$p(x \mid \mu, \lambda) = \frac{1}{x \sqrt{2 \pi \lambda}} e^{- \frac{1}{2 \lambda} (\log x - \mu)^{2}} .$$
:::

Then the mean and varaince of $X$ are $$E[X \mid \mu, \lambda] = e^{\mu + \frac{\lambda}{2}}$$ and $$Var(X \mid \mu, \lambda) = \left[ e^{\lambda} - 1\right] e^{2\mu + \lambda}$$

# Statistical Model, Part 2

1.  $\theta= e^{\mu + \frac{\lambda}{2}}$ or $(\mu, \lambda)$, parameters from LogNormal

2.  $\Theta = [0, \infty)$ 

3.  $x$ observed data, described above

4.  $\mathcal{X}$, all possible vectors of $n=$ `r length(pr_dat$Radius)` non-negative scalars

5.  $p(x \mid \mu, \lambda)$, a description of how the data arise in nature

$$p(x \mid \mu, \lambda) = \prod_{i=1}^{n}\frac{1}{x_i \sqrt{2 \pi \lambda}} e^{- \frac{1}{2 \lambda} (\log x_i - \mu)^{2}}$$

6.  $\nu(\mu, \lambda)$, assumptions about $(\mu, \lambda)$

Consider two priors:

1. $\nu(\mu, \lambda) = 1/\lambda$

2. $\nu(\mu, \lambda) = \nu(\mu) \nu(\lambda)$ with
$$\mu \sim N(a,b) \quad\quad\lambda \sim IG(c, d).$$
# MCMC for Default (Improper) Prior

```

library(invgamma)
library(SimTools)
library(mcmcse)

set.seed(5731)

y <- log(pr_datRadius)
n <-length(pr_datRadius)
ybar<-mean(y)
sv_y <- var(y)
```

```

# Minimum Monte Carlo sample size
minESS(2)
```


```
msim <- 1e4

markov <- matrix(NA_real_, nrow = msim, ncol = 2)

markov[1,] <- c(ybar, sv_y)

for (iter in 2:msim){
  markov[iter, 1] <- rnorm(1, mean = ybar, sd=sqrt(markov[iter-1,2]/n))
  markov[iter, 2] <- rinvgamma(1, shape = n/2, 
                               rate = (n*(markov[iter,1] - ybar)^2 +
                                         (n-1)*sv_y)/2)
}


plot(ts(markov[,1:2]), 
     main="Time series (trace) plots of each component")

plot(ts(markov[(0.9*msim):msim,1:2]), 
     main="Time series of last 10% for each component")


# Autocorrelation plots
acf(markov[,1])
acf(markov[,2])

# observed effective sample size
multiESS(markov)

# marginal posteriors
S_post<-Smcmc(markov)
plot(S_post, Q=c(0.05,0.95))

```


```
apply(markov, 2, mean)

apply(markov, 2, quantile, probs = c(0.025,0.975))

exp(mean(markov[,1]) + mean(markov[,2])/2)
```
![Desktop View](/assets/img/Module-6/fig_2.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_3.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_4.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_5.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_6.png){: w="700" h="400" }


# Inference
```

theta <- exp(markov[,1] + markov[,2]/2)

lb <- quantile(theta, 0.025)
ub <- quantile(theta, 0.975)
```

The mean radii is estimated to be `r signif(mean(theta), digits=4)` and the 0.95-credible interval is (`r signif(lb, digits=4)`, `r signif(ub, digits=4)`).

```
den <- density(theta)

plot(den, main="Estimated Posterior Density of Theta", xlab="Mean Radii")

# Lower and higher indices on the X-axis
l <- min(which(den$x >= lb))
h <- max(which(den$x < ub))

polygon(c(den$x[c(l, l:h, h)]),
        c(0, den$y[l:h], 0),
        col = "azure2", border=NA)

```
![Desktop View](/assets/img/Module-6/fig_7.png){: w="700" h="400" }


# MCMC for Proper Prior

```

library(invgamma)
library(SimTools)

a<-0
b<-1
c<-3
d<-2

markov <- matrix(NA_real_, nrow = msim, ncol = 2)

markov[1,] <- c(ybar, sv_y)

for (iter in 2:msim){
  markov[iter, 1] <- rnorm(1, mean = (n*ybar*b + 
                      a*markov[iter-1,2])/(n*b+markov[iter-1,2]), 
                      sd = sqrt(b*markov[iter-1,2] / (n*b + markov[iter-1,2])))
  markov[iter, 2] <- rinvgamma(1, shape = c+n/2, 
                               rate = (n*(markov[iter,1] - ybar)^2 +
                                         (n-1)*sv_y  + 2*d)/2)
}


plot(ts(markov[,1:2]), 
     main="Time series (trace) plots of each component")

plot(ts(markov[(0.9*msim):msim,1:2]), 
     main="Time series of last 10% for each component")


# Autocorrelation plots
acf(markov[,1])
acf(markov[,2])

# observed effective sample size
multiESS(markov)

# marginal posteriors
S_post<-Smcmc(markov)
plot(S_post, Q=c(0.05,0.95))

```
![Desktop View](/assets/img/Module-6/fig_8.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_9.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_10.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_11.png){: w="700" h="400" }
![Desktop View](/assets/img/Module-6/fig_12.png){: w="700" h="400" }


# Inference
```

theta <- exp(markov[,1] + markov[,2]/2)

lb <- quantile(theta, 0.025)
ub <- quantile(theta, 0.975)
```

The mean radii is estimated to be `r signif(mean(theta), digits=4)` and the 0.95-credible interval is (`r signif(lb, digits=4)`, `r signif(ub, digits=4)`).

```

den <- density(theta)

plot(den, main="Estimated Posterior Density of Theta", xlab="Mean Radii")

# Lower and higher indices on the X-axis
l <- min(which(denx >= lb))
h <- max(which(denx < ub))

polygon(c(denx[c(l, l:h, h)]),
        c(0, deny[l:h], 0),
        col = "azure2", border=NA)

```
![Desktop View](/assets/img/Module-6/fig_13.png){: w="700" h="400" }
