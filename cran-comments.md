---
title: "cran-comments"
author: "Hessa Al-Thani"
date: "3/19/2020"
output: md_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Test environments
* local OS X install, (devel and release) R 3.5.2

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
  
  ❯ checking installed package size ... NOTE
installed size is 27.5Mb
sub-directories of 1Mb or more:
  data  27.0Mb

I used xz compression on the larger datasets and 
they show to be around 6.5MB in my directory, I'm not sure
why devtools::check() has interpred it as 27MB. If you can help
me find the issue I will attempt to release the data as 1 or 2 packages 
based on the memory usage but I'd like to know why it's stating 27MB. 
Your insight would be deeply appreciated.


