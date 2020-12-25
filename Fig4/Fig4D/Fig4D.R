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

# Combine data and clean
toma <- rbind(rep_100,rep_200,rep_300) %>% 
select(Protein,
         "Control"=MS3Quant2,
         "Dox"=MS3Quant6,
         "Both"=MS3Quant9,
         "dTag"=MS3Quant10) %>%
  drop_na() %>%
  mutate("Signal"=rowSums(.[2:5])) %>%
  filter(Signal>120) %>% 
  #filter_at(vars(Control:dTag),all_vars((.) != 0)) %>% 
  mutate("FCcontrol"=Control/Control) %>%
  mutate("FCdox"=Dox/Control) %>% 
  mutate("FCboth"=Both/Control) %>% 
  mutate("FCdtag"=dTag/Control) %>% 
  select(Protein,FCcontrol,FCdox,FCboth,FCdtag)

toma2 <- toma %>% 
  melt(id.vars=1,measure.vars=2:5)
  
# Get overall correction factors
corfacs <- summarise(group_by(toma2, variable), 
                     MD = as.numeric(format(median(value,na.rm = T),digits=2)))

# Correct data
toma_corr <- toma %>% 
  mutate("Control"=FCcontrol/corfacs$MD[1]) %>%
  mutate("dTag"=FCdtag/corfacs$MD[4]) %>%
  mutate("Dox"=FCdox/corfacs$MD[2]) %>% 
  mutate("Both"=FCboth/corfacs$MD[3]) %>% 
  select(Protein,Control,dTag,Dox,Both)

# Subset for targets of interest
toma_adpgk <- toma_corr %>%
  filter(Protein=="Adpgk") %>% 
  melt(id.vars=1,measure.vars=2:5)

adpgk_dox <- toma_adpgk$value[which(toma_adpgk$variable=="Dox")]
adpgk_both <- toma_adpgk$value[which(toma_adpgk$variable=="Both")]
t.test(adpgk_both,adpgk_dox, paired = T)

toma_b2m <- toma_corr %>%
  filter(Protein=="B2M") %>% 
  melt(id.vars=1,measure.vars=2:5)

toma_tgfp <- toma_corr %>%
  filter(Protein=="tGFP") %>% 
  melt(id.vars=1,measure.vars=2:5)

tgfp_dox <- toma_tgfp$value[which(toma_tgfp$variable=="Dox")]
tgfp_dtag <- toma_tgfp$value[which(toma_tgfp$variable=="dTag")]
tgfp_both <- toma_tgfp$value[which(toma_tgfp$variable=="Both")]
t.test(tgfp_dox,tgfp_dtag, paired = T)

toma_Reps1 <- toma_corr %>%
  filter(Protein=="Reps1") %>% 
  melt(id.vars=1,measure.vars=2:5)

toma_Aatf <- toma_corr %>%
  filter(Protein=="Aatf") %>% 
  melt(id.vars=1,measure.vars=2:5)

toma_Rpl18 <- toma_corr %>%
  filter(Protein=="Rpl18") %>% 
  melt(id.vars=1,measure.vars=2:5)

toma_Gtf2i <- toma_corr %>%
  filter(Protein=="Gtf2i") %>% 
  melt(id.vars=1,measure.vars=2:5)

toma_Cpne1 <- toma_corr %>%
  filter(Protein=="Cpne1") %>% 
  melt(id.vars=1,measure.vars=2:5)
