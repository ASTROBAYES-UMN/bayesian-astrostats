---
layout: post
title: "Module 4"
date: 2023-08-10
categories: modules
tags: Module 
pin: false
math: true
---

# Galaxy Counts

The WISP survey observed hundreds of regions in the sky, each covering an area of 3.5 square arcminutes, and recorded the number of galaxies. 

*Question* If we were to choose a region of the same size how many galaxies should we expect to observe?

```

library(tidyverse)
library(modelsummary)

# Read the data in
GC_dat <- read_table("GalaxyCounts_MT.txt", col_names = "GCount")

head(GC_dat)

datasummary_skim(GC_dat)

nobs <- length(GC_dat GCount)
nobs
```


```

GC_box <- ggplot(data =  GC_dat, aes(x = GCount))  +
  geom_boxplot(outlier.color = "blue") +
  labs(x = "Galaxy Counts",  
       title = "Boxplot of Galaxy Counts")

GC_box

```

```

GC_hist <- ggplot(GC_dat, aes(x = GCount)) + 
  geom_histogram(aes(y =..density..),
                   colour = "black", fill = "white") + 
  geom_density(colour = "blue") +
  labs(x = "Galaxy Counts", y = "Density", 
       title = "Histogram of Galaxy Counts")

GC_hist

```

# Statistical Model, Part 1

1.  $\theta$, mean surface density of galaxies

2.  $\Theta = [0, \infty)$

3.  $x$ observed data, described above

4.  $\mathcal{X}$, all possible vectors of `r length(GC_dat&GCount)` non-negative integers

5.  $p(x \mid \theta)$, a description of how the data arise

6.  $\nu(\theta)$, assumptions about the mean surface density of galaxies

# Likelihood

A natural starting point is to assume $X_i \mid \theta \stackrel{ind}{\sim} \text{Poisson}(\theta)$ for $i=1,\ldots,n$  As we saw previously, this implies
$$E[X \mid \theta] = Var[X \mid \theta] = \theta.$$

Is this a reasonable assumption for the galaxy count data?
```


GC_mean <- mean(GC_datGCount) 

GC_mean / var(GC_datGCount)

```

The sample variance is larger than the sample mean, so the Poisson may not be the best choice.  This phenomenon is called *overdispersion*, that is, more variability than would be expected under the Poisson model. 

The Negative Binomial distribution provides one model which may work better than the Poisson in the presence of overdispersion.

## Negative Binomial

::: {#def-neg_bin}
If $X\mid a, \theta \sim \text{NegBin}(a,\theta)$ where $0 < \theta < 1$, then

$$p(x \mid a, \theta) = \frac{\Gamma(a+x)}{x! \Gamma(a)} \theta^{a}(1-\theta)^x \quad x =0,1,2,3,\ldots$$
:::

$$E[X\mid a,\theta]=a(1-\theta)/\theta$$ 

$$Var[X\mid a,\theta] = a(1-\theta)/\theta^2$$ 

Notice that the variance is larger than the mean since $0 < \theta < 1$ and
$$Var[X\mid a, \theta] = \frac{E[X\mid a, \theta]}{\theta}$$



### Poisson-Gamma Mixture
There are many interpretations of the Negative Binomial.  For example, the Negative Binomial mass function ca be viewed in terms of the number of failures which occur in a sequence of Bernoulli trials before a target number of successes is reached.

A more modern interpretation is that it is the result of a Poisson-Gamma mixture. Suppose $X \mid \lambda \sim \text{Poisson}(\lambda)$ and $\lambda \mid a, \theta \sim \text{Gamma}(a,\theta/(1-\theta))$. Then

$$p(x\mid \lambda)p(\lambda \mid a, \theta)$$

will define a joint probability function for $(X, \lambda) \mid \theta$ so if we integrate with respect to $\lambda$ we will obtain a marginal pmf for $X \mid \theta$.

\begin{align*}
  p(x|a, \theta) & = \int p(x|\lambda)p(\lambda \mid a, \theta) d\lambda \\
                 & = \int  \frac{e^{-\lambda} \lambda^x}{x!} \frac{\theta^a}{\Gamma(a)
                   (1-\theta)^a}  \lambda^{a-1} e^{-\lambda \theta/(1-\theta)} d\lambda\\
                 & = \frac{\theta^a}{x! \Gamma(a)
                   (1-\theta)^a} \int  \lambda^{a+x-1} e^{-\lambda/(1-\theta)} d\lambda \\
                 &= \frac{\Gamma(a+x)}{x! \Gamma(a)} \theta^{a}
                   \left(1-\theta\right)^x ,
\end{align*} which is the Negative Binomial pmf.

### Negative Binomial Data

```{r}
#| echo: true
#| eval: true
#| warning: false
#| message: false

set.seed(5731)

msim <-nobs

sim_dat1 <- rnbinom(msim, size = 1, prob = 1/2)

sim_dat2 <- rnbinom(msim, size = 5, prob = 1/2)

sim_dat3 <- rnbinom(msim, size = 1, prob = 1/4)

sim_dat4 <- rnbinom(msim, size = 1, prob = 3/4)

```


::: {#fig-pois_dat layout="[[2,2]]"}
```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

sim_dat1 <- as_tibble(sim_dat1)

nb_plot1 <- ggplot(sim_dat1, aes(x=value)) + geom_histogram() +
  labs(title = "a=1, theta=1/2", x = "Simulated values")

nb_plot1
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

sim_dat2 <- as_tibble(sim_dat2)

nb_plot2 <- ggplot(sim_dat2, aes(x=value)) + geom_histogram() +
  labs(title = "a=5, theta=1/2", x = "Simulated values")

nb_plot2
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

sim_dat3 <- as_tibble(sim_dat3)

nb_plot3 <- ggplot(sim_dat3, aes(x=value)) + geom_histogram() +
  labs(title = "a=1, theta=1/4", x = "Simulated values")

nb_plot3
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

sim_dat4 <- as_tibble(sim_dat4)

nb_plot4 <- ggplot(sim_dat4, aes(x=value)) + geom_histogram() +
  labs(title = "a =1, theta=3/4", x = "Simulated values")

nb_plot4
```

Histograms of `r msim` Negative Binomial simulations
:::

# Statistical Model, Part 2

1.  $E[X]$, mean surface density of galaxies

2.  $\Theta = [0, \infty)$

3.  $x$ observed data, described above

4.  $\mathcal{X}$, all possible vectors of `r length(GC_dat&GCount)` non-negative integers

5. Description of how the data arise, either

- $X_{i} \mid \theta \stackrel{ind}{\sim} \text{Negative Binomial}(a,\theta)$ so $E[X] = a\theta / (1-\theta)$


or 

- $X_{i} \mid \theta \stackrel{ind}{\sim} \text{Poisson}(\theta)$, so $E[X] = \theta$.

6.  $\nu(\theta)$, assumptions about $\theta$ or mean surface density of galaxies

::: {.callout-important}
Parameters only have meaning in the context of a given model.
:::

1. If $X_{i} \mid \theta \stackrel{ind}{\sim} \text{Negative Binomial}(a,\theta)$, then the conjugate prior is $\theta \sim \text{Beta}(b,c)$:

$$q(\theta \mid x_1,\ldots, x_n) \propto \left[ \prod_{i=1}^{n} \theta^a (1-\theta)^{x_i}\right] \theta^{b-1} (1-\theta)^{c-1} = \theta^{na + b -1} (1-\theta)^{n\bar{x} + c - 1}$$
Hence the posterior is $\text{Beta}(na + b, n\bar{x} + c)$.

2. If $X_{i} \mid \theta \stackrel{ind}{\sim} \text{Poisson}(\theta)$, then the conjugate prior is $\theta \sim \text{Gamma}(a,b)$:

$$q(\theta \mid x_1,\ldots, x_n) \propto \left[ \prod_{i=1}^{n} \theta^{x_{i}} e^{-\theta}\right]  \theta^{a-1} e^{-b\theta} = \theta^{a+n\bar{x} - 1} e^{-(n+b)\theta}$$
Hence the posterior is $\text{Gamma}(a + n \bar{x}, n+b)$.

# Prior Predictive Check

## Poisson-Gamma

```{r}
#| echo: true
#| eval: true
#| warning: false
#| message: false

msim <- nobs

theta <- rgamma(msim, shape = 1, rate = 1)
priorpred_dat1 <- rpois(msim, theta)

theta <- rgamma(msim, shape = 2, rate = 0.5)
priorpred_dat2 <- rpois(msim, theta)
```

::: {#fig-pg_ppd layout="[[1,1]]"}

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false


pg_ppd_plot1 <- ggplot(as_tibble(priorpred_dat1), aes(x=value)) + geom_histogram() +
  labs(title = "Prior shape=1, Prior rate=1", x = "Simulated values")

pg_ppd_plot1
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

pg_ppd_plot2 <- ggplot(as_tibble(priorpred_dat2), aes(x=value)) + geom_histogram() +
  labs(title = "Prior shape=2, Prior rate=0.5", x = "Simulated values")

pg_ppd_plot2
```

Prior predictive simulations of length `r msim` for Poisson-Gamma hierarchy.
:::

## Negative Binomial-Beta

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

theta <- rbeta(msim, shape1 = 1, shape2 = 1)
ppnb_dat1 <- rnbinom(msim, size = 1, prob = theta)


theta <- rbeta(msim, shape1 = 10, shape2 = 10)
ppnb_dat2 <- rnbinom(msim, size = 5, prob = theta)

```

::: {#fig-pg_ppd layout="[[1,1]]"}

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

nb_ppd_plot1 <- ggplot(as_tibble(ppnb_dat1), aes(x=value)) + geom_histogram() +
  labs(title = "Prior shape1=1, Prior shape2=1", x = "Simulated values")

nb_ppd_plot1
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

nb_ppd_plot2 <- ggplot(as_tibble(ppnb_dat2), aes(x=value)) + geom_histogram() +
  labs(title = "Prior shape1=10, Prior shape2=10", x = "Simulated values")

nb_ppd_plot2
```

Prior predictive simulations of length `r msim` for Negative Binomial-Beta hierarchy.
:::

# Galaxy Counts

## Poisson-Gamma

Taking $a=2.5$ and $b=0.5$ with $n=$ `r nobs` and $\bar{x}=$ `r GC_mean`, the posterior is Gamma(`r 2.5+nobs*GC_mean`, `r 0.5 + nobs`).

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

prior_shape <- 2.5
prior_rate <- 0.5

post_shape <- prior_shape + nobs*GC_mean
post_rate <- prior_rate + nobs

post_mean <- post_shape / post_rate

post_lb <- qgamma(.055, shape = post_shape, rate = post_rate)
post_ub <- qgamma(.945, shape = post_shape, rate = post_rate)
```


```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

pg_prior_plot <- ggplot() + xlim(0,10) + theme_classic() + 
  geom_function(fun = dgamma, args = 
                  list(shape=prior_shape, rate=prior_rate)) + 
  labs(x = "Theta=Mean surface density", title = "Prior (shape=2.5, rate=0.5)")

pg_post_plot <- ggplot()  + 
  xlim(5,7)  + theme_classic() +
 stat_function(fun = dgamma, 
               args = list(shape = post_shape, rate = post_rate)) + 
 labs(x = "Theta=Mean surface density", 
      y = "Posterior density", title = "Posterior")  +
  stat_function(fun = dgamma, 
                args = list(shape = post_shape, rate = post_rate),
                xlim = c(post_lb,post_ub),
                geom = "area", fill = "azure2") +
  geom_vline(xintercept = post_mean, colour = "blue", linetype = "dashed")

```

::: {#fig-pg_ppd layout="[[1,1]]"}
```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

pg_prior_plot
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

pg_post_plot
```

Prior and posterior for Poisson-Gamma model with posterior mean `r signif(post_mean, digits = 4)` and 89\% credible interval (`r signif(post_lb, digits = 4)`, `r signif(post_ub, digits = 4)`).
:::

## Negative Binomial-Beta

Taking $a=5$ and $b=c=10$ with $n=$ `r length(GC_dat&GCount)` and $\bar{x}=$ `r mean(GC_dat$GCount)`, the posterior is Beta(`r 10+5*length(GC_dat$GCount)`, `r 10 + length(GC_dat$GCount)*mean(GC_dat$GCount)`).

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

a<-5

prior_shape1 <- 10
prior_shape2 <- 10

post_shape1 <- a + prior_shape1*nobs
post_shape2 <- prior_shape2 + nobs*GC_mean

post_mean <- post_shape1 / (post_shape1 + post_shape2)

post_lb <- qbeta(.055, shape1 =  post_shape1, shape2 = post_shape2)
post_ub <- qbeta(.945, shape1 = post_shape1, shape2 = post_shape2)
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

nb_prior_plot <- ggplot() + xlim(0,1) + theme_classic() + 
  geom_function(fun = dbeta, args = 
                  list(shape1=prior_shape1, shape2=prior_shape2)) + 
  labs(x = "Theta=probability", y= "Prior density", title = "Prior (shape1=10, shape2=10)")

nb_post_plot <- ggplot()  + 
  xlim(0.55,.70)  + theme_classic() +
 stat_function(fun = dbeta, 
               args = list(shape1 = post_shape1, shape2 = post_shape2)) + 
 labs(x = "Theta=probability", 
      y = "Posterior density", title = "Posterior")  +
  stat_function(fun = dbeta, 
                args = list(shape1 = post_shape1, shape2=post_shape2),
                xlim = c(post_lb,post_ub),
                geom = "area", fill = "azure2") +
  geom_vline(xintercept = post_mean, colour = "blue", linetype = "dashed")

```

::: {#fig-pg_ppd layout="[[1,1]]"}
```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

nb_prior_plot
```

```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

nb_post_plot
```

Prior and posterior for Negative Binomial-Beta model with posterior mean `r signif(post_mean, digits = 4)` and 89\% credible interval (`r signif(post_lb, digits = 4)`, `r signif(post_ub, digits = 4)`).
:::

# Posterior Predictive Distribution

Recall our original question.

*Question* If we were to choose a region of the same size how many galaxies should we expect to observe?

Another way of restating the questions is as follows.

*Question* Given observed data $x_{1}, \ldots, x_{n}$ what can we say about a new (as yet unobserved)  rv  $\tilde{X}$?

::: {#def-post_pred}
The *posterior predictive density* is 
$$p(\tilde{x} | x_{1}, \ldots, x_{n})  = \int_{\Theta}
p(\tilde{x} | \theta) q(\theta | x_{1}, \ldots, x_{n}) d\theta .$$
:::

It is worth emphasizing that the posterior predictive distribution
depends only on the observed data and hyperparameters.  That is, it
does not depend on any unknown quantities and thus it is appropriate
to use it for prediction of future observations.

::: {#exm-pg_ppd}
Suppose $X_{1}, \ldots, X_{n} \stackrel{iid}{\sim}
\text{Poisson}(\theta)$ and $\theta \sim \text{Gamma}(a,b)$.
Then the posterior is
$$q(\theta | x_{1}, \ldots, x_{n}) \propto \left[ \prod_{i=1}^{n} e^{-\theta} \theta^{x_{i}} \right] \theta^{a-1} e^{-b\theta} = \theta^{a + n\bar{x} -1} e^{-(b+n)\theta}$$
and hence
$\theta | x_{1},\ldots, x_{n} \sim \text{Gamma}(a + n\bar{x},
b+n)$. 

The posterior predictive distribution is
$\tilde{X} | x_{1}, \ldots, x_n \sim \text{NegBin}(a + n \bar{x},
(b+n)/(b+n+1))$ since
\begin{align*}
p(\tilde{x}|x_1 \ldots, x_n) & = \int_{\Theta}
p(\tilde{x} | \theta) q(\theta | x_{1}, \ldots, x_{n}) d \theta \\
& = \frac{1}{\tilde{x}!} \frac{(b+n)^{a + n \bar{x}}}{\Gamma(a+n\bar{x})}
  \int_{0}^{\infty} \theta^{a + n\bar{x} + \tilde{x} -1}
  e^{-(b+n+1)\theta}   d \theta\\
& = \frac{1}{\tilde{x}!} \frac{(b+n)^{a + n
  \bar{x}}}{\Gamma(a+n\bar{x})} \frac{\Gamma(a + n\bar{x} +
  \tilde{x})}{(b+n+1)^{a + n\bar{x} + \tilde{x}}} \\
& = \frac{\Gamma(a + n\bar{x} +
  \tilde{x})}{\Gamma(\tilde{x}+1) \Gamma(a+n\bar{x})} \left(
  \frac{b+n}{b+n+1}\right)^{a + n \bar{x}} \left(
  \frac{1}{b+n+1}\right)^{\tilde{x}}  .
\end{align*} 
:::

::: {#exr-negbin_ppd}
$X_{1}, \ldots, X_{n} \stackrel{iid}{\sim} \text{NegBin}(a, \theta)$
where $a>0$ is known. If $\theta \sim \text{Beta}(b, c)$ find the posterior predictive distribution for each of the
  priors. Hint: Consider the
    [Beta-Negative
      Binomial distribution](https://en.wikipedia.org/wiki/Beta_negative_binomial_distribution).
:::

::: {#exr-calc}
 Suppose $X|\theta \sim p(x|\theta)$ such that $E[X|\theta]=\theta$
  and the prior $\pi(\theta)$ yields a posterior $q(\theta|x)$.

1.  Compare the posterior predictive mean with the posterior mean.

2. Show that the posterior variance is bounded above by the
  posterior predictive variance.
  
3. Suppose $X|\theta \sim \text{N}(\theta, 1)$ and
  $\theta \sim \text{N}(0,1)$.  Compare the posterior mean with the
  posterior predictive mean and the posterior variance with the
  posterior predictive variance.
  

Hint: The easy way to do this problem is via the [Law of Iterated Expectation](https://en.wikipedia.org/wiki/Law_of_total_expectation) and the [Law of Total Variance](https://en.wikipedia.org/wiki/Law_of_total_variance), but it can be solved in a direct, if somewhat messy ,fashion.


## Sampling

If the posterior predictive distribution is a name-brand distribution,
then we can often simulate directly from it.  Recall that the poisterior predictive distribution for the Poisson-Gamma hierarchy is a Negtaive Binomial.

```{r}
#| echo: true
#| eval: true
#| warning: false
#| message: false

mc_sims <-1e4

postpred_pg <- rnbinom(mc_sims, prior_shape + nobs*GC_mean, (prior_rate + nobs)/(prior_rate + nobs +1))

postpred_pg_mean <- mean(postpred_pg)

pg_postpred_plot <- ggplot(as_tibble(postpred_pg), aes(x=value))  + 
  xlim(0,20)  + theme_classic() +
  geom_histogram(aes(y=..density..), 
                 colour = "black", fill = "white") +
  geom_vline(xintercept = postpred_pg_mean, colour = "blue", linetype = "dashed") + 
  labs(x = "Simulated values", y = "Density") 

```


::: {#fig-postpred_pg}
```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

pg_postpred_plot
```

Histogram of `r mc_sims` realizations from posterior predictive distribution based on Poisson-Gamma hierarchy.  The posterior predicitve mean is `r signif(mean(postpred_pg), digits=3)` and an 89\% prediction interval is (`r quantile(postpred_pg, 0.045)`, `r quantile(postpred_pg, 0.945)`).
:::

If the posterior predictive distribution is non-standard, then it is often easy to simulate from
with the following algorithm.  Notice the similarity with the
algorithm for sampling the prior predictive distribution.

::: {#def-post_pred_samp}

The posterior predictive sampling algorithm:

For $j \in \{1,\ldots, m\}$

1. Draw $\theta_j$ from $q(\cdot | x)$.  

2. Draw $X_j$ from $p(\cdot|\theta_j)$ and call the observed
     value $x_j$.  
     
3. Save $x_j$ and discard $\theta_j$.

4. Set $j=j+1$.
:::

```{r}
#| echo: true
#| eval: true
#| warning: false
#| message: false

theta_sims <- rbeta(mc_sims, shape1=nobs*prior_shape1 + a, shape2=nobs*GC_mean + prior_shape2)
postpred_nb <- rnbinom(mc_sims, a, theta_sims)


```

::: {#fig-postpred_pnb
```{r}
#| echo: false
#| eval: true
#| warning: false
#| message: false

nb_postpred_plot <- ggplot(as_tibble(postpred_nb), aes(x=value))  + 
  xlim(0,20)  + theme_classic() +
  geom_histogram(aes(y=..density..), 
                 colour = "black", fill = "white") +
  geom_vline(xintercept = mean(postpred_nb), colour = "blue", linetype = "dashed") + 
  labs(x = "Simulated values", y = "Density") 

nb_postpred_plot
```

Histogram of `r mc_sims` realizations from posterior predictive distribution based on Negative Binomial-Beta hierarchy.  The posterior predictive mean is `r signif(mean(postpred_nb), digits=3)` and an 89\% prediction interval is (`r quantile(postpred_nb, 0.045)`, `r quantile(postpred_nb, 0.945)`).
:::

## Monte Carlo Sample Size

In the prior section the mean of the posterior predictive distribution for the Poisson-Gamma hierarchy was estimated to be `r signif(mean(postpred_pg), digits=3)` based on `r mc_sims` simulations.  Is this enough?

```{r}
#| echo: true
#| eval: true
#| warning: false
#| message: false

ci <- signif(t.test(postpred_pg, conf.level = 0.95)&conf.int, digits = 4)
ci

# mc_sims <- 1e3
postpred_pg <- rnbinom(1e3, prior_shape + nobs*GC_mean, (prior_rate + nobs)/(prior_rate + nobs +1))

mean(postpred_pg)

ci <- signif(t.test(postpred_pg, conf.level = 0.95)&conf.int, digits = 4)
ci

# mc_sims <- 1e4
postpred_pg <- rnbinom(1e4, prior_shape + nobs*GC_mean, (prior_rate + nobs)/(prior_rate + nobs +1))

mean(postpred_pg)

ci <- signif(t.test(postpred_pg, conf.level = 0.95)&conf.int, digits = 4)
ci

# mc_sims <- 1e5
postpred_pg <- rnbinom(1e5, prior_shape + nobs*GC_mean, (prior_rate + nobs)/(prior_rate + nobs +1))

mean(postpred_pg)

ci <- signif(t.test(postpred_pg, conf.level = 0.95)$conf.int, digits = 4)
ci


# mc_sims <- 1e6
postpred_pg <- rnbinom(1e6, prior_shape + nobs*GC_mean, (prior_rate + nobs)/(prior_rate + nobs +1))

mean(postpred_pg)

ci <- signif(t.test(postpred_pg, conf.level = 0.95)&conf.int, digits = 4)
ci
```


::: {#exr-nb_mc_size}
How many simulations should be used to estimate the posterior predictive mean for the Negative Binomial-Beta hierarchy?
:::