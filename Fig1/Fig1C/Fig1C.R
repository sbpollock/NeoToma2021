# Load packages
library(tidyverse)

# Import data
sht <- read.csv("./data_peptide_counts.csv")
sht2 <- sht %>% select("Peptide_Count"=MHC_Peptides,Cell_Line,Species,"Cell_Count"=Cell_Equivalents)

sht2$Peptide_Count <- as.numeric(sht2$Peptide_Count)
sht2$Cell_Count <-as.numeric(sht2$Cell_Count)

# Plot
p <- ggplot(sht2,aes(x=Species,y=Peptide_Count,color=Cell_Line)) +
geom_point(aes(size=as.numeric(Cell_Count),fill=Cell_Line),color="white",pch=21) +
scale_y_continuous(limits=c(0,6500))+
scale_size_area()+
scale_fill_brewer(palette="Dark2")+
labs(y = "Unique Peptide Count",
fill = "Cell Line",
size = "Effective \nCell Count")+
theme(
axis.title.x = element_blank(),
axis.text = element_text(size=18),
axis.text.x = element_text(angle=45,hjust=1),
axis.title.y = element_text(size=14),
legend.title = element_text(size=14),
legend.text = element_text(size=12))+
guides(size = guide_legend(override.aes = list(fill="black")))+
guides(fill = guide_legend(override.aes = list(size=5)))

# Save image
ggsave(filename="Fig1C.png",width = 3.5, height = 5)
