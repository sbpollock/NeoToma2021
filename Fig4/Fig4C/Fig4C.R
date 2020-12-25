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

# Combine data, clean, and smush
toma_smush <- rbind(rep_100,rep_200,rep_300) %>% 
  select(Protein,"UT"=MS3Quant2,"SPS"=MS3SPSIons) %>%
  drop_na() 

toma_smush_hi <- toma_smush %>%
  filter(SPS>1) %>% 
  group_by(Protein) %>% # Same when smushed
  summarize_each(list(sum)) %>% # What to do with values being smushed
  filter(UT>30)

toma_smush_lo <- toma_smush %>% 
  filter(SPS==1) %>% 
  filter(.$Protein %notin% toma_smush_hi$Protein) %>% 
  group_by(Protein) %>% # Same when smushed
  summarize_each(list(sum)) %>% # What to do with values being smushed
  filter(UT>30)
  
# Make type dataframe 
name <- c("B2M","B2M_2","CCS","Adpgk_wt",
          "Adpgk","Adpgk_3","tGFP","tGFP_OX",
          "Reps1","Cpne1","Dpagt1","Aatf","Aatf_OX1","Aatf_OX2",
          "Irgq","Med12","Rpl18","Gtf2i")

type <- c(rep("Control",4),
          rep("Construct",4),
          rep("Previous",10)
          )
df_type <- data.frame(name,type)

# Make complex melting vector (observed Tm >= 40C by DSF)
forms_complex <- c("Dnmt3a","Irgq","Atg9a","Sfi1","Reps1","Adpgk","Huwe1","N4bp2l2","Gpt","Pop1","Flrt2","Ampd2",
                   "Actr1b","Fat1","Tmem135","Itpr3","Ndfip2","Car7","Trrap","Pam","Lats1","Yipf4","Hp1bp3",
                   "Tmub1","Cspg4","Cenpe","Rpl18","Crnkl1","Zbtb40","Slc5a6","Slc4a3","Fam46b","Aatf",
                   "2410127L17Rik","Hcfc2","Lmo7","Spire1","Dpagt1","Smad2","Ivd","Tada3","Phrf1","Cpne1",
                   "Anubl1","Med12","Gtf2i","Cdh19","Nipbl","Atm","Lrrk2","Slc35e4","Prex1","Pdgfrb","D030016E14Rik",
                   "Pkd1","Hnrnpl","Rnpep","Hp1bp3_2","Cenpe_OX","Fat1_OX")


# Assign a type to each protein, otherwise assign "New" for high confidence set
type_hi <- vector()
for (i in 1:nrow(toma_smush_hi)){
  if(sum(grepl(paste0("\\b",toma_smush_hi[i,1],"\\b"),df_type$name))>0){
    temp <- df_type[which(grepl(paste0("\\b",toma_smush_hi[i,1],"\\b"),df_type$name)==T),]
    type_hi <- c(type_hi,as.character(temp[1,2]))
  }else{
    type_hi <- c(type_hi,"New")
  }
}

color_hi <- vector()
for (j in 1:nrow(toma_smush_hi)){
  if(sum(grepl(paste0("\\b",toma_smush_hi[j,1],"\\b"),forms_complex))>0){
    temp <- forms_complex[which(grepl(paste0("\\b",toma_smush_hi[j,1],"\\b"),forms_complex)==T)]
    color_hi <- c(color_hi,"red")
  }else{
    color_hi <- c(color_hi,"black")
  }
   
}

toma_smush_hi2 <- cbind(toma_smush_hi,"type"=type_hi,"color"=color_hi) %>% filter(type=="New"|type=="Previous") %>% 
  select(Protein,UT,type,color)
toma_smush_hi2$Protein <- gsub("OX1","OX",toma_smush_hi2$Protein)
toma_smush_hi2$Protein <- gsub("OX2","OX",toma_smush_hi2$Protein)
toma_smush_hi2 <- toma_smush_hi2 %>% distinct() %>% arrange(desc(UT))
toma_smush_hi2$type <- fct_relevel(toma_smush_hi2$type,c("Previous","New"))
toma_smush_hi2$Protein <- fct_rev(fct_reorder(toma_smush_hi2$Protein,toma_smush_hi2$UT))
color_hi2 <- toma_smush_hi2$color

# Assign a type to each protein, otherwise assign "New" for low confidence set
type_lo <- vector()
for (i in 1:nrow(toma_smush_lo)){
  if(sum(grepl(paste0("\\b",toma_smush_lo[i,1],"\\b"),df_type$name))>0){
    temp <- df_type[which(grepl(paste0("\\b",toma_smush_lo[i,1],"\\b"),df_type$name)==T),]
    type_lo <- c(type_lo,as.character(temp[1,2]))
  }else{
    type_lo <- c(type_lo,"New")
  }
}

color_lo <- vector()
for (j in 1:nrow(toma_smush_lo)){
  if(sum(grepl(paste0("\\b",toma_smush_lo[j,1],"\\b"),forms_complex))>0){
    temp <- forms_complex[which(grepl(paste0("\\b",toma_smush_lo[j,1],"\\b"),forms_complex)==T)]
    color_lo <- c(color_lo,"red")
  }else{
    color_lo <- c(color_lo,"black")
  }
   
}

toma_smush_lo2 <- cbind(toma_smush_lo,"type"=type_lo,"color"=color_lo) %>% filter(type=="New"|type=="Previous") %>% 
  select(Protein,UT,type,color)
toma_smush_lo2$Protein <- gsub("OX1","OX",toma_smush_lo2$Protein)
toma_smush_lo2$Protein <- gsub("OX2","OX",toma_smush_lo2$Protein)
toma_smush_lo2 <- toma_smush_lo2 %>% distinct() %>% arrange(desc(UT))
toma_smush_lo2$type <- fct_relevel(toma_smush_lo2$type,c("Previous","New"))
toma_smush_lo2$Protein <- fct_rev(fct_reorder(toma_smush_lo2$Protein,toma_smush_lo2$UT))
color_lo2 <- toma_smush_lo2$color

# Plot
p_smush_hi <- ggplot(toma_smush_hi2,aes(x=Protein,y=UT,fill=type))+
geom_col()+
labs(y="Signal",fill="Type\n(# Neoantigens)")+
  scale_fill_manual(labels=c("Previous (5)","New (54)"),
                               values=c("cornflowerblue","orange"))+
  scale_y_continuous(limits=c(0,3200))+
theme(axis.text.x=element_text(color=color_hi2,angle = 90,size=18, hjust=1,vjust=0.5),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      legend.position = "none")

p_smush_lo <- ggplot(toma_smush_lo2,aes(x=Protein,y=UT,fill=type))+
geom_col()+
labs(y="Signal",fill="Type\n(# Neoantigens)")+
  scale_fill_manual(values=c("lightgoldenrod2"))+
  scale_y_continuous(limits=c(0,3200))+
theme(axis.text.x=element_text(color=color_lo2, angle = 90,size=18, hjust=1,vjust=0.5),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none")

# Save images
ggsave(filename="Fig4C_hi.png",
       plot=p_smush_hi,
       width = 10, height = 7)

ggsave(filename="Fig4C_lo.png",
       plot=p_smush_lo,
       width = 4, height = 6.8)
