# Load packages and set working directory
library(tidyverse)
library(rstudioapi)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Import data
sht <- read.csv("data_reuse.csv") %>% 
mutate(Elution=as.factor(Elution))

# Plot data
p <- ggplot(sht,aes(x=Elution,y=Count,fill=Type,color=Type))+
stat_summary(aes(group=Type),fun = "mean", geom = "line",size=3)+
scale_y_continuous(limits = c(0,3000))+
stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2)+
labs(y="Unique Peptide Count",x="Consecutive Uses")+
theme(axis.text = element_text(size=14),
axis.title = element_text(size=14),
legend.text = element_text(size=14),
legend.title = element_text(size=14))

# Save image
ggsave(filename="Fig1E.png",
plot=p,
width = 5, height = 4)