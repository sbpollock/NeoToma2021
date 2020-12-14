# Load packages
library(tidyverse)
library(googlesheets4)
library(reshape2)

# Import data
url <- "https://docs.google.com/spreadsheets/d/1CUggPwOFysCPVGKmF4AdluJ2BYQpSLjMSf-Ou-QGI8A/edit#gid=1191321520"
sht_desalt <- read_sheet(url,sheet="data - desalt") %>% mutate("type"="desalt")
sht_TMT <- read_sheet(url,sheet="data - TMT") %>% mutate("type"="TMT")
sht_RAT <- read_sheet(url,sheet="data - RAT") %>% mutate("type"="RAT")

# Rename Area columns to make consistent
colnames(sht_desalt)[8] <- "Area"
colnames(sht_TMT)[8] <- "Area"
colnames(sht_RAT)[8] <- "Area"

# Remove Spec column to allow rbind
sht_desalt <- sht_desalt[,-13]
sht_TMT <- sht_TMT[,-13]
sht_RAT <- sht_RAT[,-13]

df0 <- rbind(sht_desalt,sht_TMT,sht_RAT)

# Clean data
clean <- function(df) {
  df %>%
    filter(str_detect(Accession,"NEO_")) %>%
    filter(!str_detect(Accession,"#DECOY")) %>%
    filter(!str_detect(Accession,"sp")) %>%
    filter(!str_detect(Accession,"tr")) %>%
    filter(Area > 0)
}

df1 <- clean(df0)

# Count neoantigens detected
treatment <- c("Desalt","TMT","RAT")

df1_ds <- df1 %>% filter(type=="desalt")
df1_tmt <- df1 %>% filter(type=="TMT")
df1_rat <- df1 %>% filter(type=="RAT")

ds_n <- df1_ds %>% select(Accession) %>% distinct()
tmt_n <- df1_tmt %>% select(Accession) %>% distinct()
rat_n <- df1_rat %>% select(Accession) %>% distinct()

lost_in_TMT_accession <- setdiff(ds_n,tmt_n)
lost_then_recovered <- setdiff(intersect(ds_n,rat_n),tmt_n)

ds_n2 <- df1_ds %>% select(Peptide) %>% distinct()
tmt_n2 <- df1_tmt %>% select(Peptide) %>% distinct()
rat_n2 <- df1_rat %>% select(Peptide) %>% distinct()

gained_with_RAT <- setdiff(rat_n2,tmt_n2)

lst <- list(ds_n,tmt_n,rat_n)

neo_count <- unlist(map(lst,nrow))
df_e <- data.frame("treatment"=as.factor(treatment),neo_count)
df_e$treatment <- fct_relevel(df_e$treatment,c("Desalt","TMT","RAT"))
```

# Plot
p_count <- ggplot(df_e,aes(x=treatment,y=neo_count,fill=treatment)) +
geom_col()+
geom_text(label=neo_count, nudge_y = 7, size=6)+
coord_cartesian(ylim=c(0, 223))+
scale_fill_brewer(palette = "Oranges")+
ylab("Neoantigens detected")+
xlab("Treatment")+
labs(fill="% Oxidant")+
theme(
axis.text = element_text(size=16),
axis.text.x = element_text(angle=45,hjust = 1),
axis.title = element_text(size=18),
axis.title.x = element_blank(),
legend.position = "none"
)

# Save image
ggsave(plot=p_count,
filename="/Volumes/RLadies/Manuscripts/2020_MCP/Figures/Fig2/Fig2B1.png",
width = 2, height = 5)
