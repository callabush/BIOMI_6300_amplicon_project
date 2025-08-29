---
title: "BioDiversity"
author: "Calla Bush St George"
date: "2025-08-26"
output:
  html_document: 
    code_folding: show
    theme: spacelab
    highlight: pygments
    keep_md: yes
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 3
  keep_md: true  
editor_options: 
  chunk_output_type: console
---


``` r
knitr::opts_chunk$set(echo = TRUE,
                      fig.align = "center",
                      fig.path = "../figures/04_Biodiversity/",
                      dev = "png", dpi = 200) 

# send any figure output to this folder 
```

# Setting the Environment 

### Set my seed

``` r
# Any number can be chose
set.seed(567890)
```

## Load Libraries 

``` r
pacman::p_load(tidyverse, devtools, patchwork, iNEXT, phyloseq, rstatix, 
               ggpubr,
               install = FALSE)
```

## Load in Data 

``` r
load("data/02_PreProcessing/raw_preprocessed_physeq.RData")
raw_preprocessed_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1736 taxa and 11 samples ]
## sample_data() Sample Data:       [ 11 samples by 23 sample variables ]
## tax_table()   Taxonomy Table:    [ 1736 taxa by 9 taxonomic ranks ]
```

``` r
# Intuition Check 
min(sample_sums(raw_preprocessed_physeq))
```

```
## [1] 7142
```

``` r
# Setting colors for gut sections 
gutsection_colors <- c(
  "IV" = "dodgerblue4",
  "V" = "#FF5733")

# Make metadata dataframe
metadata_df <-
  raw_preprocessed_physeq %>%
  sample_data() %>%
  data.frame()
```



# Goals

1. Calculate the Hill Diversity of the samples. 
2. Evaluate the rarefaction curves. 
3. Evaluate the Diversity values. 
4. Makes notes of specific samples and their seq depth. 

# Diversity Calculations with iNEXT 


``` r
# prepare input data 
iNEXT_input_df <- 
  raw_preprocessed_physeq %>%
  otu_table() %>%
  data.frame()
# Quick check
dim(iNEXT_input_df)
```

```
## [1] 1736   11
```

``` r
# Run iNEXT: Calculate the Hill Numbers 
# Note that: Species in ROWS, Samples in COLUMNS 
# Remember to set the seed! 
#iNEXT_data <- iNEXT(iNEXT_input_df, 
#                   q = c(0,1,2), datatype = "abundance")

# Save the file
#save(iNEXT_data, file = "data/04_Biodiversity/iNEXT_data.RData")
```

# Evaluate the Diversity! 

``` r
load("data/04_Biodiversity/iNEXT_data.RData")
str(iNEXT_data)
```

```
## List of 3
##  $ DataInfo:'data.frame':	11 obs. of  14 variables:
##   ..$ Assemblage: chr [1:11] "X568_4" "X568_5" "X581_5" "X611_5" ...
##   ..$ n         : num [1:11] 24503 17303 15963 13893 45378 ...
##   ..$ S.obs     : num [1:11] 47 434 441 321 603 52 408 378 520 132 ...
##   ..$ SC        : num [1:11] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ f1        : num [1:11] 0 0 0 0 0 0 0 0 0 0 ...
##   ..$ f2        : num [1:11] 4 10 5 8 5 8 10 15 8 12 ...
##   ..$ f3        : num [1:11] 5 15 15 7 9 3 9 19 11 17 ...
##   ..$ f4        : num [1:11] 0 23 18 13 20 4 10 15 12 17 ...
##   ..$ f5        : num [1:11] 2 18 27 19 18 0 17 18 17 8 ...
##   ..$ f6        : num [1:11] 4 24 21 19 16 2 17 16 14 6 ...
##   ..$ f7        : num [1:11] 1 24 20 13 19 0 15 11 22 7 ...
##   ..$ f8        : num [1:11] 0 17 12 9 15 2 16 16 15 5 ...
##   ..$ f9        : num [1:11] 0 16 19 15 15 0 9 13 11 4 ...
##   ..$ f10       : num [1:11] 2 13 18 7 15 0 14 11 11 4 ...
##  $ iNextEst:List of 2
##   ..$ size_based    :'data.frame':	1320 obs. of  10 variables:
##   .. ..$ Assemblage: chr [1:1320] "X568_4" "X568_4" "X568_4" "X568_4" ...
##   .. ..$ m         : num [1:1320] 1 1362 2723 4084 5445 ...
##   .. ..$ Method    : chr [1:1320] "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##   .. ..$ Order.q   : num [1:1320] 0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ qD        : num [1:1320] 1 29.7 35.1 38.2 40.3 ...
##   .. ..$ qD.LCL    : num [1:1320] 1 28.7 33.8 36.7 38.7 ...
##   .. ..$ qD.UCL    : num [1:1320] 1 30.6 36.3 39.6 41.8 ...
##   .. ..$ SC        : num [1:1320] 0.304 0.994 0.997 0.998 0.999 ...
##   .. ..$ SC.LCL    : num [1:1320] 0.3 0.994 0.997 0.998 0.999 ...
##   .. ..$ SC.UCL    : num [1:1320] 0.309 0.995 0.997 0.998 0.999 ...
##   ..$ coverage_based:'data.frame':	627 obs. of  8 variables:
##   .. ..$ Assemblage: chr [1:627] "X568_4" "X568_4" "X568_4" "X568_4" ...
##   .. ..$ SC        : num [1:627] 0.304 0.994 0.997 0.998 0.999 ...
##   .. ..$ m         : num [1:627] 1 1362 2723 4084 5445 ...
##   .. ..$ Method    : chr [1:627] "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##   .. ..$ Order.q   : num [1:627] 0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ qD        : num [1:627] 1 29.7 35.1 38.2 40.3 ...
##   .. ..$ qD.LCL    : num [1:627] 0.994 28.243 33.305 36.3 38.39 ...
##   .. ..$ qD.UCL    : num [1:627] 1.01 31.13 36.84 40.04 42.17 ...
##  $ AsyEst  :'data.frame':	33 obs. of  7 variables:
##   ..$ Assemblage: chr [1:33] "E05_5" "E05_5" "E05_5" "E1_4" ...
##   ..$ Diversity : chr [1:33] "Species richness" "Shannon diversity" "Simpson diversity" "Species richness" ...
##   ..$ Observed  : num [1:33] 603 163.49 49.5 52 6.21 ...
##   ..$ Estimator : num [1:33] 603 164.59 49.55 52 6.21 ...
##   ..$ s.e.      : num [1:33] 2.5167 1.3929 0.7081 4.5052 0.0437 ...
##   ..$ LCL       : num [1:33] 603 161.86 48.17 52 6.13 ...
##   ..$ UCL       : num [1:33] 607.9 167.3 50.9 60.8 6.3 ...
##  - attr(*, "class")= chr "iNEXT"
```

``` r
typeof(iNEXT_data)
```

```
## [1] "list"
```

# Plot Diversity 

``` r
# Prepare Colors 
color_df <- 
  iNEXT_input_df %>%
  colnames() %>%
  data.frame()
# Check
head(color_df)
```

```
##        .
## 1 X568_4
## 2 X568_5
## 3 X581_5
## 4 X611_5
## 5  E05_5
## 6   E1_4
```

``` r
# Rename the column 
colnames(color_df)[1] <- "names"
# Check
head(color_df)
```

```
##    names
## 1 X568_4
## 2 X568_5
## 3 X581_5
## 4 X611_5
## 5  E05_5
## 6   E1_4
```

``` r
# Make a helper dataframe for plotting with colors 
iNEXT_color_df <- 
  color_df %>%
  # Fix the names for merging
  mutate(names = gsub(names, pattern = "[.]", replace = "-"),
         names = gsub(names, pattern = "X",  replace = "")) %>%
  # Merge with metadata
  left_join(metadata_df, by = "names") %>%
   #Merge with colors for plotting with ggiNEXT
  left_join(y=data.frame(gutsection_colors = gutsection_colors,
            gut_section = names(gutsection_colors)),
            by = "gut_section")
```

# Plot Rarefaction with `ggiNEXT`


``` r
# Plot rarefaction! 
# rarefaction/extrapolation curve, type = 1 

# Order q: 
  # 0 = Richness/ Number of Total taxa
  # 1 = Exponential Shannon / Number of "Common" taxa
  # 2 = Inverse Simpson / Number of "Dominant" taxa 

ggiNEXT(iNEXT_data, type = 1, facet.var = "Order.q") + 
  facet_wrap(~Order.q, scales = "fixed") + 
  scale_color_manual(values = iNEXT_color_df$gutsection_colors, guide = FALSE) + 
  scale_fill_manual(values = iNEXT_color_df$gutsection_colors, guide = FALSE) + 
  scale_shape_manual(values = base::rep(17, nsamples(raw_preprocessed_physeq)),
                     guide = FALSE) + 
  theme(legend.position = "none")
```

```
## Scale for colour is already present.
## Adding another scale for colour, which will replace the existing scale.
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

<img src="../figures/04_Biodiversity/ggiNEXT-1.png" style="display: block; margin: auto;" />


My rarefaction curves and extrapolation looks like I would expect. Species diversity decreases as the hill numbers increase, which makes sense since the higher hill numbers represent taking into account more abundant samples. The range of species diversity is greater between my samples when I don't take into account the abundance of my samples.

# Manually plot Diversity 

## Rarefaction

``` r
iNEXT_manual_df <- 
  iNEXT_data$iNextEst$size_based %>%
  dplyr::rename(names = Assemblage) %>%
  # Fix the samples names 
  mutate(names = gsub(names, pattern = "[.]", replace = "-"),
         names = gsub(names, pattern = "X", replace = "")) %>%
  # join with metadata 
  left_join(., metadata_df, by = "names") %>%
  # Add colors to data frame
  left_join(., data.frame(gutsection_colors = gutsection_colors,
                          gut_section = names(gutsection_colors)),
            by = "gut_section") 

# Inspect 
dim(iNEXT_manual_df)
```

```
## [1] 1320   33
```

``` r
str(iNEXT_manual_df)
```

```
## 'data.frame':	1320 obs. of  33 variables:
##  $ names              : chr  "568_4" "568_4" "568_4" "568_4" ...
##  $ m                  : num  1 1362 2723 4084 5445 ...
##  $ Method             : chr  "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##  $ Order.q            : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ qD                 : num  1 29.7 35.1 38.2 40.3 ...
##  $ qD.LCL             : num  1 28.7 33.8 36.7 38.7 ...
##  $ qD.UCL             : num  1 30.6 36.3 39.6 41.8 ...
##  $ SC                 : num  0.304 0.994 0.997 0.998 0.999 ...
##  $ SC.LCL             : num  0.3 0.994 0.997 0.998 0.999 ...
##  $ SC.UCL             : num  0.309 0.995 0.997 0.998 0.999 ...
##  $ host_species       : chr  "Pomacanthus sexstriatus" "Pomacanthus sexstriatus" "Pomacanthus sexstriatus" "Pomacanthus sexstriatus" ...
##  $ gut_section        : chr  "IV" "IV" "IV" "IV" ...
##  $ region             : chr  "Great Barrier Reef" "Great Barrier Reef" "Great Barrier Reef" "Great Barrier Reef" ...
##  $ location           : chr  "South Island South" "South Island South" "South Island South" "South Island South" ...
##  $ year               : int  2014 2014 2014 2014 2014 2014 2014 2014 2014 2014 ...
##  $ month              : int  12 12 12 12 12 12 12 12 12 12 ...
##  $ day                : int  12 12 12 12 12 12 12 12 12 12 ...
##  $ sample_lab         : chr  "S246" "S246" "S246" "S246" ...
##  $ Time               : chr  "15:30" "15:30" "15:30" "15:30" ...
##  $ SL..mm.            : int  236 236 236 236 236 236 236 236 236 236 ...
##  $ FL..mm.            : int  285 285 285 285 285 285 285 285 285 285 ...
##  $ TW..g.             : num  654 654 654 654 654 ...
##  $ GW..g.             : num  52.5 52.5 52.5 52.5 52.5 52.5 52.5 52.5 52.5 52.5 ...
##  $ Sex                : chr  "M" "M" "M" "M" ...
##  $ Sample_or_Control  : chr  "Sample" "Sample" "Sample" "Sample" ...
##  $ input              : num  52502 52502 52502 52502 52502 ...
##  $ filtered           : num  25983 25983 25983 25983 25983 ...
##  $ denoisedF          : num  25713 25713 25713 25713 25713 ...
##  $ denoisedR          : num  25724 25724 25724 25724 25724 ...
##  $ merged             : num  24787 24787 24787 24787 24787 ...
##  $ nochim             : num  24522 24522 24522 24522 24522 ...
##  $ perc_reads_retained: num  46.7 46.7 46.7 46.7 46.7 ...
##  $ gutsection_colors  : chr  "dodgerblue4" "dodgerblue4" "dodgerblue4" "dodgerblue4" ...
```

``` r
# Plot it - Rarefaction Curve 
iNEXT_manual_df %>%
  # Filter out rows that are calcaulted by rarefaction from iNEXT
  dplyr::filter(Method == "Extrapolation") %>%
  # Make the actual rarefaction plot with 
  # the # of sequences on the x-axis and diversity on the y-axis
  # You can choose to pick one diversity value or plot all three 
  ggplot(aes(x = m, y= qD, color = gut_section, group = names)) + 
  # line 
  geom_line() + 
  #geom_point() + 
  # Challenge: Facet with gut section
  facet_grid(Order.q~gut_section, scales = "fixed") + 
  scale_color_manual(values = gutsection_colors) + 
  theme(legend.position = "bottom")
```

<img src="../figures/04_Biodiversity/iNEXT-manual-1.png" style="display: block; margin: auto;" />


# Diversity vs Gut Section 


``` r
#### 
# Stats
obs_div_df <- 
  iNEXT_manual_df %>%
  dplyr::filter(Method == "Observed") 
# check it
glimpse(obs_div_df)
```

```
## Rows: 33
## Columns: 33
## $ names               <chr> "568_4", "568_4", "568_4", "568_5", "568_5", "568_…
## $ m                   <dbl> 24503, 24503, 24503, 17303, 17303, 17303, 15963, 1…
## $ Method              <chr> "Observed", "Observed", "Observed", "Observed", "O…
## $ Order.q             <dbl> 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1,…
## $ qD                  <dbl> 47.000000, 5.044894, 3.285552, 434.000000, 152.518…
## $ qD.LCL              <dbl> 45.425542, 4.953346, 3.234355, 431.132546, 148.439…
## $ qD.UCL              <dbl> 48.574458, 5.136443, 3.336748, 436.867454, 156.597…
## $ SC                  <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ SC.LCL              <dbl> 0.9999120, 0.9999120, 0.9999120, 0.9997024, 0.9997…
## $ SC.UCL              <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ host_species        <chr> "Pomacanthus sexstriatus", "Pomacanthus sexstriatu…
## $ gut_section         <chr> "IV", "IV", "IV", "V", "V", "V", "V", "V", "V", "V…
## $ region              <chr> "Great Barrier Reef", "Great Barrier Reef", "Great…
## $ location            <chr> "South Island South", "South Island South", "South…
## $ year                <int> 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 20…
## $ month               <int> 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3,…
## $ day                 <int> 12, 12, 12, 12, 12, 12, 14, 14, 14, 18, 18, 18, 22…
## $ sample_lab          <chr> "S246", "S246", "S246", "S346", "S346", "S346", "S…
## $ Time                <chr> "15:30", "15:30", "15:30", "15:30", "15:30", "15:3…
## $ SL..mm.             <int> 236, 236, 236, 236, 236, 236, 190, 190, 190, 210, …
## $ FL..mm.             <int> 285, 285, 285, 285, 285, 285, 235, 235, 235, 256, …
## $ TW..g.              <dbl> 654.3, 654.3, 654.3, 654.3, 654.3, 654.3, 381.1, 3…
## $ GW..g.              <dbl> 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 308.6, 308.6, …
## $ Sex                 <chr> "M", "M", "M", "M", "M", "M", "F", "F", "F", "M", …
## $ Sample_or_Control   <chr> "Sample", "Sample", "Sample", "Sample", "Sample", …
## $ input               <dbl> 52502, 52502, 52502, 38610, 38610, 38610, 35144, 3…
## $ filtered            <dbl> 25983, 25983, 25983, 18739, 18739, 18739, 17699, 1…
## $ denoisedF           <dbl> 25713, 25713, 25713, 18165, 18165, 18165, 17222, 1…
## $ denoisedR           <dbl> 25724, 25724, 25724, 18254, 18254, 18254, 17193, 1…
## $ merged              <dbl> 24787, 24787, 24787, 17389, 17389, 17389, 16071, 1…
## $ nochim              <dbl> 24522, 24522, 24522, 17305, 17305, 17305, 15963, 1…
## $ perc_reads_retained <dbl> 46.70679, 46.70679, 46.70679, 44.81999, 44.81999, …
## $ gutsection_colors   <chr> "dodgerblue4", "dodgerblue4", "dodgerblue4", "#FF5…
```

``` r
# Pull out unique data from the three fractions of samples 
obs_whole_rich_df <- 
  obs_div_df %>%
  dplyr::filter(Order.q == 0)

# Test of the data is normal for the continuous value of richness
shapiro.test(obs_whole_rich_df$qD)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  obs_whole_rich_df$qD
## W = 0.89554, p-value = 0.1629
```

``` r
# 	Shapiro-Wilk normality test
# 
# data:  obs_whole_rich_df$qD
# W = 0.89554, p-value = 0.1629
# Hill number data is normal we can use t-test
wx_station_salinity <- 
  obs_whole_rich_df %>%
  wilcox_test(qD ~ gut_section) 

# Look at it 
wx_station_salinity
```

```
## # A tibble: 1 × 7
##   .y.   group1 group2    n1    n2 statistic      p
## * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl>
## 1 qD    IV     V          4     7         1 0.0121
```

``` r
summary(wx_station_salinity)
```

```
##      .y.               group1             group2                n1   
##  Length:1           Length:1           Length:1           Min.   :4  
##  Class :character   Class :character   Class :character   1st Qu.:4  
##  Mode  :character   Mode  :character   Mode  :character   Median :4  
##                                                           Mean   :4  
##                                                           3rd Qu.:4  
##                                                           Max.   :4  
##        n2      statistic       p         
##  Min.   :7   Min.   :1   Min.   :0.0121  
##  1st Qu.:7   1st Qu.:1   1st Qu.:0.0121  
##  Median :7   Median :1   Median :0.0121  
##  Mean   :7   Mean   :1   Mean   :0.0121  
##  3rd Qu.:7   3rd Qu.:1   3rd Qu.:0.0121  
##  Max.   :7   Max.   :1   Max.   :0.0121
```

``` r
# Pull out unique data from the three fractions of samples 
obs_whole_rich_df <- 
  obs_div_df %>%
  dplyr::filter(Order.q == 1)

# Test of the data is normal for the continuous value of richness
shapiro.test(obs_whole_rich_df$qD)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  obs_whole_rich_df$qD
## W = 0.86171, p-value = 0.06059
```

``` r
# 	Shapiro-Wilk normality test
# 
# data:  obs_whole_rich_df$qD
# W = 0.89554, p-value = 0.1629
# Hill number data is normal we can use t-test
wx_station_salinity <- 
  obs_whole_rich_df %>%
  wilcox_test(qD ~ gut_section) 

# Look at it 
wx_station_salinity
```

```
## # A tibble: 1 × 7
##   .y.   group1 group2    n1    n2 statistic       p
## * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>
## 1 qD    IV     V          4     7         0 0.00606
```

``` r
summary(wx_station_salinity)
```

```
##      .y.               group1             group2                n1   
##  Length:1           Length:1           Length:1           Min.   :4  
##  Class :character   Class :character   Class :character   1st Qu.:4  
##  Mode  :character   Mode  :character   Mode  :character   Median :4  
##                                                           Mean   :4  
##                                                           3rd Qu.:4  
##                                                           Max.   :4  
##        n2      statistic       p          
##  Min.   :7   Min.   :0   Min.   :0.00606  
##  1st Qu.:7   1st Qu.:0   1st Qu.:0.00606  
##  Median :7   Median :0   Median :0.00606  
##  Mean   :7   Mean   :0   Mean   :0.00606  
##  3rd Qu.:7   3rd Qu.:0   3rd Qu.:0.00606  
##  Max.   :7   Max.   :0   Max.   :0.00606
```

``` r
# Pull out unique data from the three fractions of samples 
obs_whole_rich_df <- 
  obs_div_df %>%
  dplyr::filter(Order.q == 2)

# Test of the data is normal for the continuous value of richness
shapiro.test(obs_whole_rich_df$qD)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  obs_whole_rich_df$qD
## W = 0.91519, p-value = 0.2806
```

``` r
# 	Shapiro-Wilk normality test
# 
# data:  obs_whole_rich_df$qD
# W = 0.89554, p-value = 0.1629
# Hill number data is normal we can use t-test
wx_station_salinity <- 
  obs_whole_rich_df %>%
  wilcox_test(qD ~ gut_section) 

# Look at it 
wx_station_salinity
```

```
## # A tibble: 1 × 7
##   .y.   group1 group2    n1    n2 statistic       p
## * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>
## 1 qD    IV     V          4     7         0 0.00606
```

``` r
summary(wx_station_salinity)
```

```
##      .y.               group1             group2                n1   
##  Length:1           Length:1           Length:1           Min.   :4  
##  Class :character   Class :character   Class :character   1st Qu.:4  
##  Mode  :character   Mode  :character   Mode  :character   Median :4  
##                                                           Mean   :4  
##                                                           3rd Qu.:4  
##                                                           Max.   :4  
##        n2      statistic       p          
##  Min.   :7   Min.   :0   Min.   :0.00606  
##  1st Qu.:7   1st Qu.:0   1st Qu.:0.00606  
##  Median :7   Median :0   Median :0.00606  
##  Mean   :7   Mean   :0   Mean   :0.00606  
##  3rd Qu.:7   3rd Qu.:0   3rd Qu.:0.00606  
##  Max.   :7   Max.   :0   Max.   :0.00606
```

``` r
#######
# Plotting

 indexes <- c("Richness", "Shannon Index", "Reverse Simpson Index")
names(indexes) <- c("0", "1", "2")

iNEXT_manual_df %>%
  dplyr::filter(Method == "Observed") %>%
  ggplot(aes(x = gut_section, y = qD)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(.~Order.q, scales = "free",  
             labeller = labeller(Order.q = indexes)) + 
  stat_smooth() + 
  labs(x = "Section", y = "Alpha Diversity Measure") + 
  theme(legend.position = "bottom")
```

```
## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
```

<img src="../figures/04_Biodiversity/div-vs--1.png" style="display: block; margin: auto;" />

This plot compares the relative number of ASVs to gut section (IV or V). Consistantly, gut section V has a higher amount of ASVs. This could be random, however it contradicts what we previously thought about the data. Our previous assumption was that gut section V was less diverse than gut section IV. I'm not sure how to parse this with gut section V having more ASVs. I thought that the number of ASVs was due to the illumina run but all the samples were run at the same time, to my knowledge.

# Session Information 

``` r
# Ensure reproducibility 
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.4.2 (2024-10-31)
##  os       macOS Sequoia 15.6.1
##  system   x86_64, darwin20
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2025-08-26
##  pandoc   3.4 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64/ (via rmarkdown)
##  quarto   1.6.42 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version date (UTC) lib source
##  abind              1.4-8   2024-09-12 [1] CRAN (R 4.4.1)
##  ade4               1.7-23  2025-02-14 [1] CRAN (R 4.4.1)
##  ape                5.8-1   2024-12-16 [1] CRAN (R 4.4.1)
##  backports          1.5.0   2024-05-23 [1] CRAN (R 4.4.0)
##  Biobase            2.66.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  BiocGenerics       0.52.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  biomformat         1.34.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  Biostrings         2.74.1  2024-12-16 [1] Bioconductor 3.20 (R 4.4.2)
##  broom              1.0.9   2025-07-28 [1] CRAN (R 4.4.1)
##  bslib              0.9.0   2025-01-30 [1] CRAN (R 4.4.1)
##  cachem             1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
##  car                3.1-3   2024-09-27 [1] CRAN (R 4.4.1)
##  carData            3.0-5   2022-01-06 [1] CRAN (R 4.4.0)
##  cli                3.6.5   2025-04-23 [1] CRAN (R 4.4.1)
##  cluster            2.1.6   2023-12-01 [2] CRAN (R 4.4.2)
##  codetools          0.2-20  2024-03-31 [2] CRAN (R 4.4.2)
##  colorspace         2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
##  crayon             1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
##  data.table         1.17.8  2025-07-10 [1] CRAN (R 4.4.1)
##  devtools         * 2.4.5   2022-10-11 [1] CRAN (R 4.4.0)
##  digest             0.6.37  2024-08-19 [1] CRAN (R 4.4.1)
##  dplyr            * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
##  ellipsis           0.3.2   2021-04-29 [1] CRAN (R 4.4.0)
##  evaluate           1.0.4   2025-06-18 [1] CRAN (R 4.4.1)
##  farver             2.1.2   2024-05-13 [1] CRAN (R 4.4.0)
##  fastmap            1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
##  forcats          * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
##  foreach            1.5.2   2022-02-02 [1] CRAN (R 4.4.0)
##  Formula            1.2-5   2023-02-24 [1] CRAN (R 4.4.0)
##  fs                 1.6.6   2025-04-12 [1] CRAN (R 4.4.1)
##  generics           0.1.4   2025-05-09 [1] CRAN (R 4.4.1)
##  GenomeInfoDb       1.42.3  2025-01-27 [1] Bioconductor 3.20 (R 4.4.2)
##  GenomeInfoDbData   1.2.13  2025-08-12 [1] Bioconductor
##  ggplot2          * 3.5.2   2025-04-09 [1] CRAN (R 4.4.1)
##  ggpubr           * 0.6.1   2025-06-27 [1] CRAN (R 4.4.1)
##  ggsignif           0.6.4   2022-10-13 [1] CRAN (R 4.4.0)
##  glue               1.8.0   2024-09-30 [1] CRAN (R 4.4.1)
##  gtable             0.3.6   2024-10-25 [1] CRAN (R 4.4.1)
##  hms                1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
##  htmltools          0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
##  htmlwidgets        1.6.4   2023-12-06 [1] CRAN (R 4.4.0)
##  httpuv             1.6.16  2025-04-16 [1] CRAN (R 4.4.1)
##  httr               1.4.7   2023-08-15 [1] CRAN (R 4.4.0)
##  igraph             2.1.4   2025-01-23 [1] CRAN (R 4.4.1)
##  iNEXT            * 3.0.2   2025-07-30 [1] CRAN (R 4.4.1)
##  IRanges            2.40.1  2024-12-05 [1] Bioconductor 3.20 (R 4.4.2)
##  iterators          1.0.14  2022-02-05 [1] CRAN (R 4.4.0)
##  jquerylib          0.1.4   2021-04-26 [1] CRAN (R 4.4.0)
##  jsonlite           2.0.0   2025-03-27 [1] CRAN (R 4.4.1)
##  knitr              1.50    2025-03-16 [1] CRAN (R 4.4.1)
##  labeling           0.4.3   2023-08-29 [1] CRAN (R 4.4.0)
##  later              1.4.2   2025-04-08 [1] CRAN (R 4.4.1)
##  lattice            0.22-6  2024-03-20 [2] CRAN (R 4.4.2)
##  lifecycle          1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
##  lubridate        * 1.9.4   2024-12-08 [1] CRAN (R 4.4.1)
##  magrittr           2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
##  MASS               7.3-61  2024-06-13 [2] CRAN (R 4.4.2)
##  Matrix             1.7-1   2024-10-18 [2] CRAN (R 4.4.2)
##  memoise            2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
##  mgcv               1.9-1   2023-12-21 [2] CRAN (R 4.4.2)
##  mime               0.13    2025-03-17 [1] CRAN (R 4.4.1)
##  miniUI             0.1.2   2025-04-17 [1] CRAN (R 4.4.1)
##  multtest           2.62.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  nlme               3.1-166 2024-08-14 [2] CRAN (R 4.4.2)
##  pacman             0.5.1   2019-03-11 [1] CRAN (R 4.4.0)
##  patchwork        * 1.3.1   2025-06-21 [1] CRAN (R 4.4.1)
##  permute            0.9-8   2025-06-25 [1] CRAN (R 4.4.1)
##  phyloseq         * 1.50.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  pillar             1.11.0  2025-07-04 [1] CRAN (R 4.4.1)
##  pkgbuild           1.4.8   2025-05-26 [1] CRAN (R 4.4.1)
##  pkgconfig          2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
##  pkgload            1.4.0   2024-06-28 [1] CRAN (R 4.4.0)
##  plyr               1.8.9   2023-10-02 [1] CRAN (R 4.4.0)
##  profvis            0.4.0   2024-09-20 [1] CRAN (R 4.4.1)
##  promises           1.3.3   2025-05-29 [1] CRAN (R 4.4.1)
##  purrr            * 1.1.0   2025-07-10 [1] CRAN (R 4.4.1)
##  R6                 2.6.1   2025-02-15 [1] CRAN (R 4.4.1)
##  RColorBrewer       1.1-3   2022-04-03 [1] CRAN (R 4.4.0)
##  Rcpp               1.1.0   2025-07-02 [1] CRAN (R 4.4.1)
##  readr            * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
##  remotes            2.5.0   2024-03-17 [1] CRAN (R 4.4.0)
##  reshape2           1.4.4   2020-04-09 [1] CRAN (R 4.4.0)
##  rhdf5              2.50.2  2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
##  rhdf5filters       1.18.1  2025-03-06 [1] Bioconductor 3.20 (R 4.4.3)
##  Rhdf5lib           1.28.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  rlang              1.1.6   2025-04-11 [1] CRAN (R 4.4.1)
##  rmarkdown          2.29    2024-11-04 [1] CRAN (R 4.4.1)
##  rstatix          * 0.7.2   2023-02-01 [1] CRAN (R 4.4.0)
##  rstudioapi         0.17.1  2024-10-22 [1] CRAN (R 4.4.1)
##  S4Vectors          0.44.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  sass               0.4.10  2025-04-11 [1] CRAN (R 4.4.1)
##  scales             1.4.0   2025-04-24 [1] CRAN (R 4.4.1)
##  sessioninfo        1.2.3   2025-02-05 [1] CRAN (R 4.4.1)
##  shiny              1.11.1  2025-07-03 [1] CRAN (R 4.4.1)
##  stringi            1.8.7   2025-03-27 [1] CRAN (R 4.4.1)
##  stringr          * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
##  survival           3.7-0   2024-06-05 [2] CRAN (R 4.4.2)
##  tibble           * 3.3.0   2025-06-08 [1] CRAN (R 4.4.1)
##  tidyr            * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
##  tidyselect         1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
##  tidyverse        * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
##  timechange         0.3.0   2024-01-18 [1] CRAN (R 4.4.0)
##  tzdb               0.5.0   2025-03-15 [1] CRAN (R 4.4.1)
##  UCSC.utils         1.2.0   2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  urlchecker         1.0.1   2021-11-30 [1] CRAN (R 4.4.0)
##  usethis          * 3.1.0   2024-11-26 [1] CRAN (R 4.4.1)
##  utf8               1.2.6   2025-06-08 [1] CRAN (R 4.4.1)
##  vctrs              0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
##  vegan              2.7-1   2025-06-05 [1] CRAN (R 4.4.1)
##  withr              3.0.2   2024-10-28 [1] CRAN (R 4.4.1)
##  xfun               0.52    2025-04-02 [1] CRAN (R 4.4.1)
##  xtable             1.8-4   2019-04-21 [1] CRAN (R 4.4.0)
##  XVector            0.46.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
##  yaml               2.3.10  2024-07-26 [1] CRAN (R 4.4.0)
##  zlibbioc           1.52.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
## 
##  [1] /Users/cab565/Library/R/x86_64/4.4/library
##  [2] /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/library
##  * ── Packages attached to the search path.
## 
## ──────────────────────────────────────────────────────────────────────────────
```
