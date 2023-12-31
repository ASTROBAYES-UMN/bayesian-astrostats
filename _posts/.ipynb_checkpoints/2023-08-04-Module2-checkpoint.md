---
layout: post
title: "Introductory Probability and Statistics"
date: 2023-07-10
categories: modules
tags: test demo pleaseWork
pin: false
---

```{python}
#| warning: false
#| message: false
#| eval: false
#| echo: false
import pandas as pd
import scipy as sp
import pymc as pm
```

## Introduction to Quarto

This document was produced using [Quarto](https://quarto.org/) in [Rstudio](https://www.rstudio.com/). Quarto is fairly new to me so I make no claims about the optimality of my approach. I would greatly appreciate being made aware of any typos or outright errors in what follows so that I can correct them.

The advantage of Quarto is that it is straightforward to incorporate a variety of code in a single document that makes all of the embedded computations fully reproducible. For example, we can write mathematics with [LaTeX](https://www.latex-project.org/) and include both `R` and `Python`, among other languages in a single document. The document can be written in [Markdown](https://en.wikipedia.org/wiki/Markdown) or [Jupyter](https://jupyter.org/) notebooks. I am using Markdown.

The document produced can be a pdf, html, or docx file. And this can be changed at any point.

We can easily include images, tables, and figures. See Figure 1, for example.


We can write displayed equations as 
$$\frac{dN}{d M_{V}} \propto \frac{1}{\sigma \sqrt{2 \pi}} \exp \left\{ - \frac{1}{2 \sigma^2} \left( M_{V} - M_{V,0} \right)^2 \right\} .$$ 
or we can write in-line mathematics as $y = m \cos (x) + b$. For displayed equations simply begin and end with \$\$ while in-line math begins and ends with a single \$. We can also include cross-referenced theorems, figures, definitions, lemmas, and so on.

::: {#thm-pi_over_e}
$$\int_{\mathbb{R}} \frac{\cos(x)}{1 + x^2}dx = \frac{\pi}{e}.$$
:::

See @thm-pi_over_e for an interesting integral.

::: callout-note
So far I haven't been able to figure out how to change the numbering so that everything doesn't start with 0. Help would be appreciated.
:::

### Languages Supported in Quarto

[Bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) is fully supported.

```{bash}
ls *.qmd
```

Which shows the Quarto documents in my current directory.

Quarto also supports [Julia](https://quarto.org/docs/computations/julia.html) and [Observable JS](https://quarto.org/docs/computations/ojs.html), but I will focus on [Python](https://quarto.org/docs/computations/python.html) and [R](https://quarto.org/docs/computations/r.html).

## Getting Started with R and Python

Let's begin with some simple `R`

```{r}
# This is a comment
x <- 100
x <- x + 1
```

and follow it up with some simple `Python`

```{python}
# This is a comment
x = 100
x = x + 1
```

The code is quite similar, but the assignment operator in `R` is `<-` while in `Python` it is `=`. Technically, you can also use `=` in `R`, but this is not advised. A nice discussion of the implications can be found on [stackoverflow](https://stackoverflow.com/questions/1741820/what-are-the-differences-between-and-assignment-operators).

In both languages the left and right sides of the command `x = y` or `x <- y` are different: the left side of an assignment must be a name that can be assigned to the result of evaluating the right side.

Both languages are [object-oriented](https://en.wikipedia.org/wiki/Object-oriented_programming) programming languages. Everything is an object. An object is a combination of data and functions. When processing an assignment statement they attach (bind) a name to an object. These names are called variables.

If we want to evaluate whether a variable has a given value, we use the operator `==`.

In `R`

```{r}
y <- (x==100)
y
```

In `Python`

```{python}
y = (x==100)
y
```

In both cases, the variable y is assigned a logical value (True or False) after comparing the numerical value of x with 100.

In `Python` you can assign multiple variables with a single `=`.

```{python}
a , b , c = 1 , 2 , 3
print(a,b,c)
```

but, to my knowledge, this is not a thing in `R`. Instead you would need to create them separately.

```{r}
a <- 1
b <- 2
c <- 3

cat(a,b,c)
```

## Getting Help

The command `help(...)` in `R` will give the information it has about the expression in parenthesis. But the following will *not* print the results of `help` in this document.

```{r}
help(sqrt)
```

Instead it will appear in the Help tab in R studio. Figure 2 shows what you would see.


The command `help(...)` in `Python` will also give the information it has about the expression in parenthesis.

Try `help(sqrt)` to see what happens when Python does not know a name.

```{python}
#| eval: false
help(sqrt)
```

If you try the above code, `Python` will barf. I used the `#| eval: false` option to prevent evaluation of it. In the next section, we will see why and how to fix this.

## Libraries and Packages

A substantial difference between `R` and `Python` is that in `R` the `base` package is loaded automagically. This package contains many basic functions and built-in constants such as `sqrt` and `pi`, respectively. In `Python` this is not the case.

In `Python` you will need to use an external library that does not come standard with `Python`. These libraries are imported with the `import` command.

An often-used module is [NumPy](https://numpy.org/doc/stable/reference/), which is a collection of tools for numerical calculations.

```{python}
import numpy as np
```

Now the following will recognize `np.sqrt` and compute $\sqrt{9}$.

```{python}
np.sqrt(9)
```

In both `R` and `Python` trigonometric functions use angles expressed in radians. First, `R`

```{r}
pi

exp(1)

sin(pi/2)
```

then `Python`

```{python}
np.pi

np.e

np.sin(np.pi/2.)
```

While `R` comes with more available on startup, we almost always want to do things beyond what is available in `base R` so we need to load packages in `R` too. Here is the syntax for that.

```{r}
#| message: false
library(tidyverse)
library(ggplot2)
library(modelsummary)
```

::: callout-note
There are many more similarities between `R` and `Python`.

-   Neither offers protection against changing the value of a symbol. (I.e., you can set $\pi$ to 5.)

-   Name collision can occur if you inadvertently reuse variable names for different quantities. It is something to look out for, especially in longer programs.

-   Variable names are case sensitive.

-   Best practice is to use long, meaningful names for variables. A few variable names are forbidden (e.g., if, and, TRUE,...). You can find all the forbidden variable names with a Web search for `R` or `Python` reserved words.
:::

## Printing

Let's begin with `R`.

```{r}
a <- 2e7
b <- 123
print(a)
cat(a,b)
cat("Hello", "World!")
my_string <- "Which is larger"
cat(my_string, a, "or", b, "?")
```

Now for the same thing in `Python`.

```{python}
a=2e7
b=123
print(a)
print(a,b)
print("Hello" + "World!")
my_str1 = "Which is larger"
print(my_str1, a, "or", b, "?")
```

## Plots

For a demonstration of a line plot on a polar axis, see @fig-polar.

```{python}
#| label: fig-polar
#| fig-cap: "A line plot on a polar axis"

import matplotlib.pyplot as plt

r = np.arange(0, 2, 0.01)
theta = 2 * np.pi * r
fig, ax = plt.subplots(
  subplot_kw = {'projection': 'polar'} 
)
ax.plot(theta, r)
ax.set_rticks([0.5, 1, 1.5, 2])
ax.grid(True)
plt.show()
```
