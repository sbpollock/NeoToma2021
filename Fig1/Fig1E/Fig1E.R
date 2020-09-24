library(tidyverse)
library(googlesheets4)
url <- "https://docs.google.com/spreadsheets/d/1ckQf5HFtVQqYj7j0HsYKJ5TcdXlXNoB-0TxZDguELd8/edit#gid=743349605"
sht <- read_sheet(url,sheet="result - reuse") %>% 
  mutate(Elution=as.factor(Elution))

p <- ggplot(sht,aes(x=Elution,y=Count,fill=Type,color=Type))+
stat_summary(aes(group=Type),fun = "mean", geom = "line",size=3)+
scale_y_continuous(limits = c(0,3000))+
stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2)+
labs(y="Unique Peptide Count",x="Consecutive Uses")+
theme(axis.text = element_text(size=14),
axis.title = element_text(size=14),
legend.text = element_text(size=14),
legend.title = element_text(size=14))

ggsave(filename="/Volumes/RLadies/Manuscript_AssayMAP/Figures/Fig1/Fig1E.png",
plot=p,
width = 5, height = 4)
