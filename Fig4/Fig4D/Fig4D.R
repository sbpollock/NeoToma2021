library(tidyverse)
library(googlesheets4)
library(reshape2)
library(gridExtra)
library(wesanderson)

`%notin%` = function(x,y) !(x %in% y)

# Import sheets
rep_100 <- read.csv("./toma1.csv")
rep_200 <- read.csv("./toma2.csv")
rep_300 <- read.csv("./toma3.csv")
