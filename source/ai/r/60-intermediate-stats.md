---
tags: r tutorial statistics survey regression
---

# Intermediate Statistics: From Surveys to Statistics

**Author:** Kyle Monahan
**Date:** April 2020

---

**Learning objectives:** By the end of this tutorial you will understand how to import and clean survey data, recode and construct new variables, specify a complex survey design object, compute survey-weighted descriptive statistics, fit simple and multiple survey-weighted regression models (linear and logistic), and conduct basic regression diagnostics.

**Dataset:** This tutorial uses the [DHS Model Dataset](https://dhsprogram.com/data/Download-Model-Datasets.cfm?flag=1) : specifically the "Children's Data - Children's Recode (KR)" file in Stata format (`.DTA`). Download this file from the DHS website and place it in the same directory as this tutorial. Contact DataLab at datalab-support@elist.tufts.edu for assistance accessing the data.

---

## Getting Started

R is a statistical platform similar to Stata, SAS, and SPSS. It allows you to manipulate data, perform descriptive statistics, recode variables, and bring in your own data.

For additional online statistics resources, see [go.tufts.edu/stats](https://go.tufts.edu/stats).

## Set Working Directory

Go to **Session > Set Working Directory > To Source File Location** before running any code. All file paths are relative to the working directory.

## R Basics

```r
n <- 200
n
```

Variables are stored in RAM and appear in the **Environment** pane. You can combine values into a vector using `c()`:

```r
grades <- c(94, 96, 72, 92)
```

Access a specific element:

```r
grades[2]
```

Compute descriptive statistics:

```r
mean(grades)
median(grades)
sd(grades)
```

---

## Background: The DHS Dataset

We will explore the relationship between **breastfeeding** and **weight-for-height z-score (WHZ)** among children aged 12–23 months, using the DHS Children's Recode (KR) dataset.

**Unit of analysis:** Children born to interviewed women in the five years preceding the survey (0–59 months).

The research question: *Is breastfeeding associated with wasting (low weight-for-height) status among children aged 12–23 months?*

---

## Install and Load Packages

```r
# Uncomment the line below to install if needed:
# install.packages(c("tidyverse", "haven", "survey"), repos = "https://cloud.r-project.org")
library("tidyverse")
library("haven")
library("survey")
```

Packages used:

- `tidyverse`: `dplyr` for data manipulation, `ggplot2` for visualization
- `haven`: importing Stata/SAS/SPSS data files
- `survey`: analyzing complex survey data with proper weights and design

---

## Data Science Pipelines

The general workflow for a data analysis pipeline is:

**Import → Tidy → Transform → Visualize → Model → Communicate**

We will follow this workflow through the analysis.

## Import Data

Download the DHS model dataset from [https://dhsprogram.com/data/Download-Model-Datasets.cfm?flag=1](https://dhsprogram.com/data/Download-Model-Datasets.cfm?flag=1) and place the `ZZKR62FL.DTA` file in your working directory. Then import it:

```r
dhs_full <- haven::read_dta("ZZKR62FL.DTA")
```

---

## Data Management

Limit the dataset to children aged 12–23 months whose anthropometric measurements were taken, and select the variables relevant to our analysis.

**Key variables:**

| Variable               | Description                                                                  |
| ---------------------- | ---------------------------------------------------------------------------- |
| `caseid`               | Identification number                                                        |
| `m4`                   | Breastfeeding status (93=ever/not current; 94=never; 95=still breastfeeding) |
| `hw11`                 | Weight-for-Height Z-score × 100 (WHZ); 9998 = inconsistent                   |
| `hw1`                  | Child's age in months                                                        |
| `hw13`                 | Result of anthropometric measurement (0 = measured)                          |
| `v411`–`v414w`         | Complementary food consumption indicators                                    |
| `v190`                 | Household wealth index (quintiles 1–5)                                       |
| `v149`                 | Mother's educational attainment (0–5)                                        |
| `h31`, `h22`, `h11`    | Cough, fever, diarrhea in past 2 weeks                                       |
| `v005`                 | Sample weight (divide by 1,000,000)                                          |
| `v021`, `v002`, `v023` | PSU, SSU, and stratification variable                                        |

Use the pipe operator `%>%` to chain multiple data manipulation steps:

```r
dhs_example <- dhs_full %>%
  filter(between(hw1, 12, 23), hw13 == 0,
         hw11 != 9998 & !is.na(hw11),
         m4 %in% c(93, 94, 95)) %>%
  mutate(whz           = hw11 / 100,
         sample_weight = v005 / 1000000) %>%
  select(caseid, age = hw1, whz, breastfeed = m4,
         v411, v411a, v412, v412a,
         starts_with("v414"),
         wealth_index = v190, mother_edu = v149,
         h31, h22, h11,
         sample_weight, psu = v021, ssu = v002, strat = v023)
```

```r
head(dhs_full)
```

---

## Recoding Variables and Creating New Variables

### Factor Categorical Variables

```r
cfvarname <- names(select(dhs_full, starts_with("v414")))
catvar <- c("breastfeed", cfvarname, "v411", "v412", "v411a", "v412a",
            "wealth_index", "h31", "h22", "h11", "psu", "ssu", "strat")
dhs_example[catvar] <- lapply(dhs_example[catvar], factor)
summary(dhs_example)
```

### Recode Breastfeeding as Binary

The `m4 = 94` (never breastfed) category has very few observations in the model dataset. Recode breastfeeding as a binary variable (1 = currently breastfeeding, 0 = not):

```r
dhs_example <- dhs_example %>%
  mutate(breastfeed = factor(case_when(
    breastfeed == 93 | breastfeed == 94 ~ 0,
    TRUE ~ 1)))
```

### Construct Minimum Dietary Diversity (MDD)

MDD is defined as consuming foods from ≥ 4 food groups in the 24 hours preceding the interview:

```r
dhs_example <- dhs_example %>%
  mutate(
    grain      = case_when(v412a == 1 | v414e == 1 | v414f == 1 ~ 1, TRUE ~ 0),
    legume     = case_when(v414o == 1 ~ 1, TRUE ~ 0),
    dairy      = case_when(v411 == 1 | v411a == 1 | v412 == 1 | v414v == 1 | v414p == 1 ~ 1, TRUE ~ 0),
    flesh      = case_when(v414h == 1 | v414m == 1 | v414n == 1 ~ 1, TRUE ~ 0),
    egg        = case_when(v414g == 1 ~ 1, TRUE ~ 0),
    vavegfru   = case_when(v414i == 1 | v414j == 1 | v414k == 1 ~ 1, TRUE ~ 0),
    othervegfru = case_when(v414l == 1 ~ 1, TRUE ~ 0),
    mdd = factor(case_when(
      grain + legume + dairy + flesh + egg + vavegfru + othervegfru >= 4 ~ 1,
      TRUE ~ 0))
  ) %>%
  select(-starts_with("v414"), -v412a, -v411, -v411a)
```

```r
head(dhs_example)
```

### Construct Morbidity Variable

Morbidity is defined as presence of cough, diarrhea, or fever in the past two weeks:

```r
dhs_example <- mutate(dhs_example,
  morbidity = factor(case_when(
    h31 == 1 | h31 == 2 | h11 == 1 | h11 == 2 | h22 == 1 ~ 1,
    TRUE ~ 0)))
```

### Construct Wasting Variable

Wasting is defined as WHZ ≤ -2:

```r
dhs_example <- mutate(dhs_example,
  wasting = factor(case_when(
    whz <= -2 ~ 1,
    TRUE ~ 0)))
```

---

## Specifying the Survey Design

Before analysis, specify the complex survey design. This ensures that all subsequent estimates properly account for the DHS sampling units (PSU and SSU), stratification, and sampling weights:

```r
design <- svydesign(ids     = ~psu + ssu,
                    strata  = ~strat,
                    weights = ~sample_weight,
                    data    = dhs_example)
summary(design)
```

---

## Descriptive Statistics

Use `svymean()` for survey-weighted means across multiple variables:

```r
vars <- c("whz", "age", "mother_edu", "wealth_index", "breastfeed",
          "mdd", "grain", "legume", "dairy", "flesh", "egg",
          "vavegfru", "othervegfru", "wasting")

svymean(make.formula(vars), design, na.rm = TRUE)
```

Stratify by breastfeeding status using `svyby()`:

```r
svyby(~whz + age + mother_edu + wealth_index + mdd + wasting,
      ~breastfeed, design, svymean)
```

Some differences in mean values by breastfeeding status are apparent, suggesting we should control for these covariates in the statistical model.

---

## Correlations and Exploratory Analysis

Load `jtools` for survey-weighted correlation:

```r
# install.packages("jtools")
library(jtools)
```

```r
cor_d <- jtools::svycor(~whz + dairy, design = design)
cor_d$cors
```

---

## Statistical Modeling: Simple Linear Regression

Fit a survey-weighted simple linear regression of WHZ on breastfeeding status:

```r
slr <- svyglm(whz ~ breastfeed, design)
summary(slr)
```

## Regression Diagnostics

Examine the distribution of residuals to check the normality assumption:

```r
# Histogram of residuals
svyhist(~residuals(slr), design)

# Shapiro-Wilk test for normality
shapiro.test(svytable(~residuals(slr), design))

# Q-Q plot
survey::svyqqmath(~residuals(slr), design)
```

The right tail may be long by visual inspection, suggesting that the normality assumption is only weakly satisfied in this model dataset (which is intended for practice, not publication).

---

## Statistical Modeling: Multiple Linear Regression

Add covariates to the model:

```r
mlr <- svyglm(whz ~ breastfeed + age + mother_edu + wealth_index + mdd, design)
summary(mlr)
```

## Optional: Interaction Terms

Explore whether the relationship between breastfeeding and WHZ differs by MDD status:

```r
mlr_int <- svyglm(whz ~ breastfeed * mdd + age + mother_edu + wealth_index, design)
summary(mlr_int)
```

---

## Logistic Regression

Use the binary `wasting` variable as the outcome. Specify `family = binomial()` for a logistic model:

```r
logit_bi <- svyglm(wasting ~ breastfeed, design, family = binomial())
summary(logit_bi)

logit <- svyglm(wasting ~ breastfeed * mdd + age + mother_edu + wealth_index,
                design, family = binomial())
summary(logit)
```

---

> **Note:** These models are exploratory demonstrations using a practice dataset. A full analysis for publication would require a more thorough review of model assumptions, variable selection strategy, and the specific DHS sampling design documentation.

For questions, contact DataLab at **datalab-support@elist.tufts.edu**.
