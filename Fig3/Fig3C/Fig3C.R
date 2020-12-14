# Load packages
library(tidyverse)
library(googlesheets4)
library(reshape2)

# Import data
url <- "https://docs.google.com/spreadsheets/d/1Szq9Am0pg4orMjy0Ig23TzGS50emVHk-q7jsJzUSI_4/edit#gid=412290502"

sht <- read_sheet(url,sheet = "data - MS2") 

# Clean data
shtc <- sht %>% 
  select("fmol"=concentration,
         Peptide,
         "Protein"=`Protein Accession`,
         "126"=`Intensity TMT6-126`,
         "127"=`Intensity TMT6-127`,
         "128"=`Intensity TMT6-128`,
         "129"=`Intensity TMT6-129`,
         "130"=`Intensity TMT6-130`,
         "131"=`Intensity TMT6-131`) %>% 
  filter(str_detect(Protein,"NEO_")) %>%
    filter(!str_detect(Protein,"#DECOY")) %>%
    filter(!str_detect(Protein,"sp")) %>%
    filter(!str_detect(Protein,"tr")) %>%
  mutate("sum"=rowSums(.[4:9])) %>% 
group_by(fmol,Protein) %>% # Same when smushed
top_n(1,sum)
  
# Add CV (uncorrected)
SD <- apply(shtc[,4:9],1,sd)
MN <- apply(shtc[,4:9],1,mean)
CV <- (SD/MN)*100
shtc2 <- as.data.frame(cbind(shtc,"CV"=CV))
shtc3 <- shtc2 %>%
  filter_at(vars('126':'131'),all_vars((.) != 0)) %>% 
  filter(CV<100) %>% 
  add_row(fmol=0.1) %>%
  add_row(fmol=0.01) %>% 
add_row(fmol=0.001) 

# Order points from left to right using bins by level, then CV

# Split data by fmol level
lev1 <- shtc3 %>% filter(fmol==0.001) %>% arrange(CV)
lev2 <- shtc3 %>% filter(fmol==0.01) %>% arrange(CV)
lev3 <- shtc3 %>% filter(fmol==0.1) %>% arrange(CV)
lev4 <- shtc3 %>% filter(fmol==1) %>% arrange(CV)
lev5 <- shtc3 %>% filter(fmol==10) %>% arrange(CV)
lev6 <- shtc3 %>% filter(fmol==100) %>% arrange(CV)

# Grab arranged names for each level
l1p <- lev1$Protein
l2p <- lev2$Protein
l3p <- lev3$Protein
l4p <- lev4$Protein
l5p <- lev5$Protein
l6p <- lev6$Protein

# Order proteins by removing those that have already appeared on lower levels
`%notin%` = function(x,y) !(x %in% y)
ni2 <- subset(l2p,l2p %notin% l1p)
ni3 <- subset(l3p,l3p %notin% ni2)
ni4 <- subset(l4p,l4p %notin% ni3)
ni5 <- subset(l5p,l5p %notin% ni4)
ni6 <- subset(l6p,l6p %notin% ni5)

# Create protein order
prot_order <- unique(c(l1p,ni2,ni3,ni4,ni5,ni6))

#Order proteins in graph
shtc3$Protein <- as.factor(shtc3$Protein)
shtc3$Protein <- fct_relevel(shtc3$Protein,prot_order)
  
# Name variable for later plotting
sht_MS2 <- shtc3





# Work up TOMAHAQ data
library(tidyverse)
library(googlesheets4)
library(reshape2)

url <- "https://docs.google.com/spreadsheets/d/1Szq9Am0pg4orMjy0Ig23TzGS50emVHk-q7jsJzUSI_4/edit#gid=412290502"

sht1 <- read_sheet(url,sheet = "data - 100") 
sht2 <- read_sheet(url,sheet = "data - 200") 
sht3 <- read_sheet(url,sheet = "data - 300") 

# Combine data
sht <- rbind(sht1,sht2,sht3) %>% 
  select("fmol"=concentration,
         Peptide,
         Protein,
         "126"=MS3Quant1,
         "127"=MS3Quant2,
         "128"=MS3Quant5,
         "129"=MS3Quant6,
         "130"=MS3Quant9,
         "131"=MS3Quant10) %>% 
  mutate("Signal"=rowSums(.[4:9])) %>% 
  group_by(fmol,Peptide,Protein) %>% # Same when smushed
  summarize_each(list(sum)) %>% 
  filter(Signal>90)
  #filter_at(vars('126':'131'),all_vars((.) >= 10))

# Count targets
test <- sht %>% filter(fmol==0.001)
namez <- vector()
for (a in 1:nrow(test)){
  namez <- c(namez,strsplit(test$Protein[a],"_")[[1]][1])
}
length(unique(namez))

# Add CV
SD <- apply(sht[,4:9],1,sd)
MN <- apply(sht[,4:9],1,mean)
CV <- (SD/MN)*100
sht2 <- cbind(sht,"CV"=CV) %>% 
filter(CV<100) 


# Order points from left to right using bins by level, then CV

# Split data by fmol level
lev1 <- sht2 %>% filter(fmol==0.001) %>% arrange(CV)
lev2 <- sht2 %>% filter(fmol==0.01) %>% arrange(CV)
lev3 <- sht2 %>% filter(fmol==0.1) %>% arrange(CV)
lev4 <- sht2 %>% filter(fmol==1) %>% arrange(CV)
lev5 <- sht2 %>% filter(fmol==10) %>% arrange(CV)
lev6 <- sht2 %>% filter(fmol==100) %>% arrange(CV)

# Grab arranged names for each level
l1p <- lev1$Protein
l2p <- lev2$Protein
l3p <- lev3$Protein
l4p <- lev4$Protein
l5p <- lev5$Protein
l6p <- lev6$Protein

# Order proteins by removing those that have already appeared on lower levels
`%notin%` = function(x,y) !(x %in% y)
ni2 <- subset(l2p,l2p %notin% l1p)
ni3 <- subset(l3p,l3p %notin% ni2)
ni4 <- subset(l4p,l4p %notin% ni3)
ni5 <- subset(l5p,l5p %notin% ni4)
ni6 <- subset(l6p,l6p %notin% ni5)

# Create protein order
prot_order <- unique(c(l1p,ni2,ni3,ni4,ni5,ni6))

#Order proteins in graph
sht2$Protein <- as.factor(sht2$Protein)
sht2$Protein <- fct_relevel(sht2$Protein,prot_order)

# Plot
p_toma <- ggplot(sht2,aes(x=Protein,y=as.factor(fmol),size=Signal,color=CV))+
  geom_point()+
  scale_size_area(max_size=10,
                  breaks=c(1e2,1e3,1e4,1e5,1e6))+
  scale_color_gradient2(low="blue",
                        mid="grey",
                        high="red",
                        midpoint = 20)+
  ylab("fmol peptide")+
  xlab("Proteins")+
  labs(color="CV %",title="TOMAHAQ")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))

ggsave(filename="/Volumes/RLadies/Manuscript_AssayMAP/Figures/Fig3/FigS3_TOMA_SN90.png",
       plot=p_toma,
       width = 12, height = 4)

# Split out names and filter table so that it's comparable to MS2
sht_TOMA <- sht2 %>% 
  mutate("simple"=str_split(Protein,"_")[[1]][1]) %>% 
  filter(simple!="B2M") %>% 
  filter(simple!="CCS") %>% 
  filter(simple!="tGFP") %>% 
  filter(Protein!="ITFHFTPV") %>% 
  filter(Protein!="ITFHFTPVL")







# Make detected df
Concentration <- rep(c("100","10","1","0.1","0.01","0.001"),2)
Data <- c(rep("untargeted",6),rep("TOMAHAQ",6))

Detected <- c(
length(unique((subset(sht_MS2,sht_MS2$fmol==100)$Protein))),
length(unique((subset(sht_MS2,sht_MS2$fmol==10)$Protein))),
length(unique((subset(sht_MS2,sht_MS2$fmol==1)$Protein))),
0,
0,
0,
length(unique((subset(sht_TOMA,sht_TOMA$fmol==100)$simple))),
length(unique((subset(sht_TOMA,sht_TOMA$fmol==10)$simple))),
length(unique((subset(sht_TOMA,sht_TOMA$fmol==1)$simple))),
length(unique((subset(sht_TOMA,sht_TOMA$fmol==0.1)$simple))),
length(unique((subset(sht_TOMA,sht_TOMA$fmol==0.01)$simple))),
length(unique((subset(sht_TOMA,sht_TOMA$fmol==0.001)$simple)))
)

df_bars <- data.frame("Concentration"=as.factor(Concentration),Data,Detected)
df_bars$Concentration <- fct_relevel(df_bars$Concentration,c("100","10","1","0.1","0.01","0.001"))
df_bars$Data <- fct_relevel(df_bars$Data,c("untargeted","TOMAHAQ"))

# Make neoantigens detected graph
p_count <- ggplot(df_bars,aes(x=Concentration,y=Detected,fill=Data))+
geom_col(position = position_dodge())+
geom_text(label=Detected,position = position_dodge(width = 1),vjust=-0.5,size=6)+
labs(y="Neoantigens Detected",x="Concentration (fmol)") +
theme(
axis.text = element_text(size=16),
axis.title = element_text(size=14),
legend.title = element_blank(),
legend.text = element_text(size=14)
)+
scale_y_continuous(limits=c(0,165))

ggsave(filename="/Volumes/RLadies/Manuscripts/Manuscript_2020_MCP/Figures/Fig3/Fig3C.png",
       plot=p_count,
       width = 6, height = 3, units = "in")
