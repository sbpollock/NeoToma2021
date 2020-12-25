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

# Plot 4plex abundance with bars overlaying points and median value
p_adpgk <- ggplot(toma_adpgk,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..),color=variable),fun = "median", geom="text",position = position_nudge(y=15),size=10)+
  labs(y="Signal",title="Adpgk")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )
  

p_b2m <- ggplot(toma_b2m,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="B2M")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

p_tgfp <- ggplot(toma_tgfp,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="tGFP")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

p_Reps1 <- ggplot(toma_Reps1,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="Reps1")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

p_Aatf <- ggplot(toma_Aatf,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="Aatf")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

p_Rpl18 <- ggplot(toma_Rpl18,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="Rpl18")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

p_Gtf2i <- ggplot(toma_Gtf2i,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="Gtf2i")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

p_Cpne1 <- ggplot(toma_Cpne1,aes(x=variable,y=value,color=variable))+
  stat_summary(fun =median, geom="col",fill=NA,size=3)+
  theme_classic()+
  scale_color_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  scale_fill_manual(values=c("peru","violetred1","turquoise3","darkorchid2"))+
  geom_jitter(aes(y=value,fill=variable),size=3,color="white",pch=21)+
  stat_summary_bin(aes(label = round(..y..,digits=2),color=variable),fun = "median", geom="text",position = position_nudge(y=1),size=10)+
  labs(y="Signal",title="Cpne1")+
  scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=22),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        plot.title = element_text(size=32)
        )

# Save images  
ggsave(filename="Fig4D_Adpgk.png",
       plot=p_adpgk,
       width = 6, height = 4)

ggsave(filename="Fig4D_B2M.png",
       plot=p_b2m,
       width = 6, height = 4)

ggsave(filename="Fig4D_tGFP.png",
       plot=p_tgfp,
       width = 6, height = 4)

ggsave(filename="Fig4D_Reps1.png",
       plot=p_Reps1,
       width = 6, height = 4)

ggsave(filename="Fig4D_Aatf.png",
       plot=p_Aatf,
       width = 6, height = 4)

ggsave(filename="Fig4D_Rpl18.png",
       plot=p_Rpl18,
       width = 6, height = 4)

ggsave(filename="Fig4D_Gtf2i.png",
       plot=p_Gtf2i,
       width = 6, height = 4)

ggsave(filename="Fig4D_Cpne1.png",
       plot=p_Cpne1,
       width = 6, height = 4)
