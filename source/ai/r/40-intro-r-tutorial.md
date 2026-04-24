---
tags: r tutorial data-science
---

# Introduction to R

**Authors:** Uku-Kaspar Uustalu & Kyle Monahan
**Date:** 2024-02-16

---

**Learning objectives:** By the end of this tutorial you will be able to assign variables and vectors, perform descriptive statistics, import and explore a CSV dataset, clean and subset data, create scatter plots with base R and ggplot2, use the tidyverse pipe operator, and create basic interactive visualizations with plotly.

**Dataset:** This tutorial uses the [HURDAT2](https://www.nhc.noaa.gov/data/#hurdat) Atlantic Hurricane Database (1851–2022) from NOAA's National Hurricane Center. The data file `atlantic.csv` is included in this tutorial's directory.

---

## Getting Started

R is a statistical platform similar to Stata, SAS, and SPSS. It allows you to manipulate data, perform descriptive statistics, recode variables, and bring in your own data.

If you see code that looks like this: `code`, this refers to code you could write in the R console.

## Set Working Directory

The first thing we need to do is set a working directory. Go to **Session > Set Working Directory > To Source File Location**. This is always good practice!

## Comments

You write a comment by adding a `#` to the start of a line in an R script or code section. Or select the lines you wish to comment and press **Ctrl/Cmd+Shift+C**. Comments appear in green (in the default theme).

## Assigning Variables

The first thing we can do is assign a variable. This looks like an arrow placing a value into a variable: `n <- 200`.

```r
n <- 200
# The variable n is set equal to 200.

# You can write a comment in an R script by writing a "#"
```

This created a variable `n` which appears in the **Environment** tab on the upper right-hand side of the window.

> **Tip:** If you can't see the Console window, just click on the word "Console" below. The `>` tells you that R is waiting for your input!

Now if you call `n`, you will see the console report its value.

```r
n
```

## Console Window

**Try it:** Set `n <- 300` in the Console window, then print `n`.

Try typing `N` in the Console window. What happens? The variable `N` is not found, because R is case-sensitive!

## Functions and Vectors

What if you have multiple numbers to store as a variable? For example, storing exam scores of 94, 96, 72, and 92. We use the combine function `c()`.

```r
scores <- c(94, 96, 72, 92)
```

If we call the variable, the values are printed. This is called a **vector**.

```r
scores
```

Access individual values using square brackets `[]`. To get the second score:

```r
scores[2]
```

Note that R is one-indexed. To get the first score, use `[1]`, not `[0]`.

## Operators

To check if a vector contains an element, use the `%in%` operator.

```r
96 %in% scores
```

It returns `TRUE` if the vector contains the element.

```r
100 %in% scores
```

Returns `FALSE` if it does not. `TRUE` and `FALSE` are boolean (logical) data types in R.

## Descriptive Statistics

Most statistical functions in R easily take a whole vector as input.

```r
mean(scores)    # mean
median(scores)  # median
sd(scores)      # standard deviation
```

## The Environment

Click on the **Environment** tab. You will see the objects you have created. You can also call `objects()` for a list of all object names.

```r
objects()
```

## Modifying a Vector

Suppose we need to replace the value of 72 (the lowest score) with a new score of 85.

```r
scores[3] <- 85
```

## Saving Variables and Creating Tables

Re-calculate stats and create variables to store them. Always save anything you might need going forward.

```r
scores_mean   <- mean(scores)
scores_median <- median(scores)
scores_sd     <- sd(scores)
```

Create a table with the `rbind()` function:

```r
scores_table <- rbind(Mean   = scores_mean,
                      Median = scores_median,
                      SD     = scores_sd)
scores_table
```

Check the class of the object:

```r
class(scores_table)
```

This is a **matrix**, the most primitive form of table in R. Note how our column is unnamed. We can fix this:

```r
colnames(scores_table) <- c("Value")
scores_table
```

You can round values using the `round()` function. Use `help(round)` to see all its arguments:

```r
round(scores_table, digits = 2)
```

## Working with Real Data

We have a CSV file of the Atlantic Hurricane Database (HURDAT2) 1851–2022. It is included as `atlantic.csv` in this tutorial's directory.

## Working Directory

All relative file paths are relative to the current working directory. Use `getwd()` and `dir()` to explore it.

```r
getwd()   # Current working directory
dir()     # Files in the working directory
```

Verify that `atlantic.csv` is present:

```r
"atlantic.csv" %in% dir()
```

If this returns `TRUE`, you are all set. If `FALSE`, set the working directory via **Session > Set Working Directory > To Source File Location**.

## Importing Data

The `atlantic.csv` file is in CSV format. Import it with `read.csv()`.

```r
hurrdata <- read.csv("atlantic.csv")
```

Explore additional arguments with `?read.csv`.

## Exploring Data

A new variable `hurrdata` has been added to the environment.

```r
head(hurrdata)
```

Use `summary()` to get descriptive statistics.

```r
summary(hurrdata)
```

Check the data type:

```r
class(hurrdata)
```

This is a **data frame**, a table with observations as rows and variables as columns.

## Working with Data Frames

Access the maximum wind speed of the sixth observation (Maximum.Wind is the ninth column):

```r
hurrdata[[6, 9]]     # [[row, column]]
hurrdata[[9]][6]     # [[column]][row]
```

Or use the column name:

```r
hurrdata[["Maximum.Wind"]][6]
hurrdata$Maximum.Wind[6]
```

The `$` operator is the preferred way to access data frame columns by name.

## Selecting Data Based on Conditions

Create a frequency table of all unique values in a column:

```r
table(hurrdata$Name)
```

Select rows based on a condition (all hurricanes named Nicole):

```r
hurrdata[hurrdata$Name == "NICOLE", ]
```

Get the maximum wind speed for all hurricanes named Nicole:

```r
max(hurrdata$Maximum.Wind[hurrdata$Name == "NICOLE"])
```

## Data Cleaning: Dates & Strings

The date is stored as a number in YYYYMMDD format. Extract year and month using `substr()`.

```r
date_strings <- as.character(hurrdata$Date)
```

Extract the year from positions 1–4:

```r
hurrdata$Year <- substr(date_strings, start = 1, stop = 4)
```

Extract the month from positions 5–6:

```r
hurrdata$Month <- substr(date_strings, 5, 6)
```

Convert both to numeric:

```r
hurrdata$Month <- as.numeric(hurrdata$Month)
hurrdata$Year  <- as.numeric(hurrdata$Year)
```

## Making a Simple Scatter Plot

Explore how maximum wind speed has changed over time:

```r
plot(hurrdata$Year, hurrdata$Maximum.Wind)
```

## Data Cleaning: Missing Data

Some wind speeds are negative; these are missing value codes. Replace them with `NA`:

```r
hurrdata$Maximum.Wind[hurrdata$Maximum.Wind < 0] <- NA
```

Get the minimum wind speed excluding NAs:

```r
min(hurrdata$Maximum.Wind, na.rm = TRUE)
```

## Sampling the Data

Sample 1000 rows to declutter the scatter plot:

```r
hurrdata2 <- hurrdata[sample(nrow(hurrdata), 1000), ]
```

## Customizing the Scatter Plot

```r
plot(x = hurrdata2$Year,
     y = hurrdata2$Maximum.Wind,
     main = "Selected Annual Hurricane Data, 1851-2022",
     xlab = "Year",
     ylab = "Maximum Wind Speed (knots)",
     pch  = 21,
     col  = "blue",
     bg   = "lightblue")
```

Add a regression line:

```r
abline(reg = lm(formula = Maximum.Wind ~ Year, data = hurrdata2),
       col = "red",
       lty = "dashed",
       lwd = 3)
```

## Base R and the Tidyverse

Thus far we have been working with **base R**, R without any community-developed packages. The most popular collection of R packages for data science is called the **Tidyverse**. Key tidyverse packages used in this tutorial:

- `readr`: reading and writing tabular data
- `lubridate`: working with dates and times
- `dplyr`: easy modification of data frames
- `ggplot2`: streamlined data visualization

## Checking Installed Packages

```r
"tidyverse" %in% installed.packages()
```

Returns `TRUE` if tidyverse is installed.

## Installing Packages

Install a package only if it is not already present:

```r
if ("tidyverse" %in% installed.packages()) {
  message("Tidyverse already installed!")
} else {
  install.packages("tidyverse")
}
```

You only need to install packages once per machine.

## Calling a Library

Before using a package, load it with `library()`:

```r
library(tidyverse)
```

Note any reported conflicts; if two packages have functions with the same name, use the `package::function()` notation to call a specific one.

## Advanced: Package Management Using Librarian

The `librarian` package simplifies package management. The `shelf()` function installs a package if not present and then loads it.

```r
if (! "librarian" %in% installed.packages()) {
  install.packages("librarian")
}
```

```r
librarian::shelf(janitor)
```

## Rewriting the Analysis Using Tidyverse

Use `readr::read_csv()` to import the data (faster than `read.csv()`):

```r
hurrdata3 <- readr::read_csv("atlantic.csv")
```

Clean column names to snake_case with `janitor`:

```r
hurrdata3 <- janitor::clean_names(hurrdata3)
```

Convert the date column with lubridate and add year/month columns using `dplyr::mutate()`:

```r
hurrdata3 <- dplyr::mutate(hurrdata3, date = lubridate::ymd(date))

hurrdata3 <- dplyr::mutate(hurrdata3,
                           year  = lubridate::year(date),
                           month = lubridate::month(date))
```

Replace missing wind speed codes and sample rows:

```r
hurrdata3 <- dplyr::mutate(hurrdata3,
                           maximum_wind = dplyr::na_if(maximum_wind, -99))

hurrdata4 <- dplyr::sample_n(hurrdata3, 1000)
```

## Alternate (Better) Data Visualizations

Use `ggplot2` to create a scatter plot. Build the plot in layers using `+`:

```r
hurrgraph <- ggplot(data = hurrdata4, aes(x = year, y = maximum_wind)) +
  geom_point()
```

```r
hurrgraph
```

Add color, labels, and formatting:

```r
hurrgraph2 <- ggplot(data = hurrdata4,
                     aes(x = year, y = maximum_wind, color = maximum_wind)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme(legend.position = "bottom") +
  ylab("Maximum Wind Speed (knots)") +
  ggtitle("Selected Annual Hurricane Data, 1851-2022") +
  theme(plot.title = element_text(lineheight = 0.8, face = "bold")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
```

```r
hurrgraph2
```

## Advanced: Piping and Grouping

Was sampling the correct approach? Actually each hurricane appears multiple times in the data at different intensities. To see how hurricane intensity has changed over time, we should look at the peak wind speed per hurricane.

Use the pipe operator `%>%` with `dplyr::group_by()` and `dplyr::summarize()`:

```r
hurrdata5 <- hurrdata3 %>%
  dplyr::group_by(name, year) %>%
  dplyr::summarize(maximum_wind = max(maximum_wind), .groups = "drop")
```

Now visualize correctly with a trend line:

```r
hurrgraph3 <- ggplot(data = hurrdata5, aes(x = year, y = maximum_wind)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x) +
  ylab("Maximum Wind Speed (knots)") +
  ggtitle("Hurricane Wind Speeds at Highest Intensity 1851-2021")
```

```r
hurrgraph3
```

## Exercise: Fix the Analysis

Take a look at the year 1950. What causes the sudden change? Is it fair to compare all recorded hurricanes from 1851 to 2022? Should we instead be looking at the most intense hurricane from each year?

Fix the analysis above and produce a graph that answers: *How has the intensity of hurricanes changed over time?*

## Optional: Interactive Graphics

Create interactive visualizations using `plotly`:

```r
librarian::shelf(plotly)
```

```r
hurrdata6 <- hurrdata3 %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(count = dplyr::n_distinct(id), .groups = "drop")
```

```r
hurrgraph4 <- plotly::plot_ly(data = hurrdata6,
                              x = ~year,
                              y = ~count,
                              type = "bar",
                              marker = base::list(
                                color = "lightblue",
                                line  = base::list(color = "blue", width = 0.5))) %>%
  plotly::layout(title = "Number of Hurricanes 1851 - 2022",
                 xaxis = base::list(title = "Year", tickangle = 45),
                 yaxis = base::list(title = "Recorded Number of Hurricanes"))
```

```r
hurrgraph4
```

## Bonus: Defining Functions

Define your own functions in R using the `function()` keyword:

```r
get_x <- function(t) {
  x = 16 * sin(t)^3
  return(x)
}

get_y <- function(t) 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t)
```

Store results in a data frame and plot:

```r
t <- seq(0, 2 * pi, by = 0.1)

parameq   <- data.frame(t = t)
parameq$x <- get_x(parameq$t)
parameq$y <- get_y(parameq$t)

finalplot <- ggplot(data = parameq, aes(x = x, y = y)) +
  geom_polygon(linetype = "solid", linewidth = 1.5,
               color = "black", fill = "red")
```

## FAQ

**What is the difference between `=` and `<-`?**

Both assign values, but they differ in scope. `<-` puts a value in the user's global workspace, making it accessible anywhere. `=` declares it in the scope of a function. In most cases for assigning top-level variables, `<-` is conventional in R.
