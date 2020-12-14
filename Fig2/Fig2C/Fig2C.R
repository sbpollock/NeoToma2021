# Load packages
library(tidyverse)
library(googlesheets4)
library(reshape2)

# Import data
sheet <- read.csv("./oxidation.csv")

# Separate TIC by Ox %
Ox_percent <- as.factor(c(0,0.01,0.05,0.1,0.25,0.5,1,1.5))

df0 <- subset(sheet,Ox_percent==0)
df0.01 <- subset(sheet,Ox_percent==0.01)
df0.05 <- subset(sheet,Ox_percent==0.05)
df0.1 <- subset(sheet,Ox_percent==0.1)
df0.25 <- subset(sheet,Ox_percent==0.25)
df0.5 <- subset(sheet,Ox_percent==0.5)
df1 <- subset(sheet,Ox_percent==1)
df1.5 <- subset(sheet,Ox_percent==1.5)

TIC <- c(
sum(df0$Area),
sum(df0.01$Area),
sum(df0.05$Area),
sum(df0.1$Area),
sum(df0.25$Area),
sum(df0.5$Area),
sum(df1$Area),
sum(df1.5$Area)
)

df_ox <- data.frame(Ox_percent,TIC) %>% 
  mutate(log10=log(TIC,base=10))
      
# Filter dataframe for Neoantigens
clean <- function(df) {
  df %>%
    filter(str_detect(Accession,"NEO_")) %>%
    filter(!str_detect(Accession,"#DECOY")) %>%
    filter(!str_detect(Accession,"sp")) %>%
    filter(!str_detect(Accession,"tr")) %>%
    filter(Area > 0)
}

df0_c <- clean(df0)
df0.01_c <- clean(df0.01)
df0.05_c <- clean(df0.05)
df0.1_c <- clean(df0.1)
df0.25_c <- clean(df0.25)
df0.5_c <- clean(df0.5)
df1_c <- clean(df1)
df1.5_c <- clean(df1.5)

neo_count <- c(
length(unique(df0_c$Accession)),
length(unique(df0.01_c$Accession)),
length(unique(df0.05_c$Accession)),
length(unique(df0.1_c$Accession)),
length(unique(df0.25_c$Accession)),
length(unique(df0.5_c$Accession)),
length(unique(df1_c$Accession)),
length(unique(df1.5_c$Accession))
)

df_ox2 <- df_ox %>%
  mutate(neo_count)

# Plot
p2 <- ggplot(df_ox2,
             aes(x=Ox_percent,y=neo_count,fill=Ox_percent))+
      geom_col()+
      geom_text(label=neo_count,nudge_y = 7,size=6)+
      coord_cartesian(ylim=c(0, 223))+
      scale_fill_hue(h=c(90,0),c=75)+
      ylab("Neoantigens detected")+
      xlab("Treatment Condition (% Oxidant)")+
      #ggtitle("Neoantigen count with increasing oxidant","Out of 223 neoantigens")+
      labs(fill="% Oxidant")+
      theme(
axis.text = element_text(size=18),
axis.title.x = element_text(size=18),
axis.title.y = element_text(size=18),
axis.text.x = element_text(angle=45,hjust=1),
legend.position = "none"
)

# Save image
ggsave(filename="Fig2C.png",
plot=p2,
width = 5, height = 5
)
