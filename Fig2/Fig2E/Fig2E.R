# Load packages
library(tidyverse)
library(googlesheets4)
library(reshape2)

# Import data
sheet <- read.csv("./oxidation.csv")

# Define functions
clean <- function(df) {
  df %>%
    filter(str_detect(Accession,"NEO_")) %>%
    filter(!str_detect(Accession,"#DECOY")) %>%
    filter(!str_detect(Accession,"sp")) %>%
    filter(!str_detect(Accession,"tr")) %>%
    filter(Area > 0)
}

met <- function(df) {
  df %>%
    filter(str_detect(Peptide,"M"))
}

`%notin%` = function(x,y) !(x %in% y)

# Find all neoepitopes that have Mets
df_c <- clean(sheet) 
df_cm <- met(df_c)
neos <- unique(df_cm$Accession)

# Simplify dataframe and select max values
df_forms <- df_cm %>% 
  select(Accession,Ox_percent,Peptide,Area) %>%
  mutate(Ox_percent=as.factor(Ox_percent)) %>%
  mutate(Accession=gsub("NEO_","",Accession)) %>% 
  group_by(Accession,Ox_percent) %>% 
  filter(Area == max(Area)) %>%
  arrange(Accession,Ox_percent,Peptide) %>% 
  mutate("Met_count"=str_count(Peptide,"M")) %>% 
  arrange(Met_count,Ox_percent,desc(Area)) 

# Relevel Ox_percent
df_forms$Ox_percent <- factor(df_forms$Ox_percent,levels = rev(levels(factor(df_forms$Ox_percent))))

# Relevel accession names
df_forms$Accession <- factor(df_forms$Accession, levels = unique(df_forms$Accession))

# Make lists of helped by oxidation and hurt by oxidation
# Cast df
df_cast <- df_forms %>% 
  select(Accession,Ox_percent,Area) %>% 
  dcast(Accession ~ Ox_percent, fun.aggregate = sum)

helped <- vector()
for (i in 1:nrow(df_cast)){
  if (df_cast[i,9]<df_cast[i,2]|
      df_cast[i,9]<df_cast[i,3]|
      df_cast[i,9]<df_cast[i,4]|
      df_cast[i,9]<df_cast[i,5]|
      df_cast[i,9]<df_cast[i,6]|
      df_cast[i,9]<df_cast[i,7]|
      df_cast[i,9]<df_cast[i,8]
  ){
    helped <- c(helped,as.character(df_cast[i,1]))
  } else {
        next
      }
}

# Plot
p_m1h <- ggplot(subset(df_forms,Met_count==1 & Accession %in% helped),
aes(x=Accession,y=Ox_percent,fill=Area))+
geom_tile()+
scale_x_discrete(position = "top")+
scale_fill_gradient(low = "white", high = "black",limits=c(0,4E8))+
labs(y="% Oxidant")+
theme_classic()+
theme(
axis.text.x.top = element_text(angle=90,hjust=0, vjust=0.5, size = 12),
axis.title.x = element_blank(),
axis.title.y = element_text(size=12),
axis.text.y = element_text(size = 12),
legend.text = element_text(size = 12),
legend.title = element_text(size = 12))

p_m1nh <- ggplot(subset(df_forms,Met_count==1 & Accession %notin% helped),
aes(x=Accession,y=Ox_percent,fill=Area))+
geom_tile()+
scale_x_discrete(position = "top")+
scale_fill_gradient(low = "white", high = "black",limits=c(0,4E8))+
labs(y="% Oxidant")+
theme_classic()+
theme(
axis.text.x.top = element_text(angle=90,hjust=0, vjust=0.5, size = 12),
axis.title.x = element_blank(),
axis.title.y = element_text(size=12),
axis.text.y = element_text(size = 12),
legend.position = "none")

p_m23h <- ggplot(subset(df_forms,Met_count!=1 & Accession %in% helped),
aes(x=Accession,y=Ox_percent,fill=Area))+
geom_tile()+
scale_x_discrete(position = "top")+
scale_fill_gradient(low = "white", high = "black",limits=c(0,4E8))+
labs(y="% Oxidant")+
theme_classic()+
theme(
axis.text.x.top = element_text(angle=90,hjust=0, vjust=0.5, size = 12),
axis.title.x = element_blank(),
axis.title.y = element_text(size=12),
axis.text.y = element_text(size = 12),
legend.position = "none")

p_m23nh <- ggplot(subset(df_forms,Met_count!=1 & Accession %notin% helped),
aes(x=Accession,y=Ox_percent,fill=Area))+
geom_tile()+
scale_x_discrete(position = "top")+
scale_fill_gradient(low = "white", high = "black",limits=c(0,4E8))+
labs(y="% Oxidant")+
theme_classic()+
theme(
axis.text.x.top = element_text(angle=90,hjust=0, vjust=0.5, size = 12),
axis.title.x = element_blank(),
axis.title.y = element_text(size=12),
axis.text.y = element_text(size = 12),
legend.position = "none")

# Save images
ggsave(filename="Fig2E1_m1h.png",
       plot=p_m1h,
       width = 7.5, height = 2.5)

ggsave(filename="Fig2E2_m1nh.png",
       plot=p_m1nh,
       width = 3.75, height = 2.5)

ggsave(filename="Fig2E3_m23h.png",
       plot=p_m23h,
       width = 1.875, height = 2.5)

ggsave(filename="Fig2E4_m23nh.png",
       plot=p_m23nh,
       width = 1.875, height = 2.5)
