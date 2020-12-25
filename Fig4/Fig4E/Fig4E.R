library(ggrepel)

df_pep <- read.csv(file = "/Volumes/RLadies/Manuscripts/2020_MCP/Figures/Fig4/Saundra1 - id-Adpgk(m)-g in MC38 treatment to MS - data - MS2.csv", stringsAsFactors = F)

# Break apart df into a list where each element is a unique protein accession
df_pep_l <- split(df_pep,df_pep$Protein.Accession)

# Filter dataframes within list for rows that do not have signal > 10,000 in 126 and 128 channels
dfl2 <- list()
for (b in 1:length(df_pep_l)){
  dfl2[[b]] <- df_pep_l[[b]] %>% 
    filter(Intensity.TMT6.126<10000 & Intensity.TMT6.128<(0.1*Intensity.TMT6.127))
  names(dfl2)[b] <- names(df_pep_l)[b]
}

# Filter list for dataframes that have more than one peptide
dfl3 <- list()
n <- 1
for (a in 1:length(dfl2)){
  if (sum(dfl2[[a]]$X.Spec)>1){
    dfl3[[n]] <- dfl2[[a]]
    names(dfl3)[n] <- names(dfl2)[a]
    n <- n+1
  } else {
    next
  }
}

# (sum(dfl2[[a]]$X.Spec)
# (nrow(dfl2[[a]])>1)

# Create clean dataframes within a new list
abundance_list <- list()
for (i in 1:length(dfl3)){
df_temp <- dfl3[[i]] %>%
select(
"Protein"=Protein.Accession,
"ut"=Intensity.TMT6.127,
"dox"=Intensity.TMT6.129,
"both"=Intensity.TMT6.130,
"dtag"=Intensity.TMT6.131) #%>% 
#melt(id.vars=1,measure.vars=2:5)
  
abundance_list[[i]] <- df_temp
names(abundance_list)[i] <- names(dfl3)[i] 
}

# Bind into a single large dataframe
big_df <- bind_rows(abundance_list) %>% 
  group_by(Protein) %>% # Same when smushed
  summarize_each(list(sum)) # Different when smushed

# Simplify names
big_df$Protein <- gsub("CONSTRUCT_FKBP12F36V_AdpgkR304M_tGFP","AdpgkR304M",big_df$Protein)

for (c in 2:nrow(big_df)){
big_df[c,1] <- str_split(big_df[c,1],"\\|")[[1]][3]
big_df[c,1] <- str_split(big_df[c,1],"_")[[1]][1]
}

# UT vs. Dox abundance
p_ut_dox <- ggplot(big_df)+
  geom_abline(slope=1)+
  geom_point(aes(x=ut,y=dox))+
  labs(x="Control (signal intensity)",y="Dox (signal intensity)")+
  theme_minimal()+
    scale_x_log10(limits=c(1E4,1E7))+
    scale_y_log10(limits=c(1E4,1E7))+
  geom_label_repel(data=subset(big_df,Protein=="RBP2"|Protein=="PSMD3"),
      nudge_y=0.4,aes(x=ut,y=dox,label=Protein),size=6)+
    geom_point(data=subset(big_df,dox>=3*ut),
      aes(x=ut,y=dox,color="red"))+
    theme(legend.position = "none",
          axis.title = element_text(size=18),
          axis.text = element_text(size=14),
  axis.text.x = element_text(angle=45))+
    geom_label_repel(data=subset(big_df,Protein=="AdpgkR304M"),fill="green",
      nudge_y = 0.5, size=6,
      aes(x=ut,y=dox,label=Protein))+
    geom_point(data=subset(big_df,Protein=="AdpgkR304M"),color="green",size=6,
      aes(x=ut,y=dox))
  
p_ut_dtag <- ggplot(big_df)+
  geom_abline(slope=1)+
  geom_point(aes(x=ut,y=dtag))+
  labs(x="Control (signal intensity)",y="dTag (signal intensity)")+
  theme_minimal()+
    scale_x_log10(limits=c(1E4,1E7))+
    scale_y_log10(limits=c(1E4,1E7))+
  geom_label_repel(data=subset(big_df,Protein=="RBP2"|Protein=="PSMD2"|Protein=="PSMD3"),
      nudge_y=0.65,aes(x=ut,y=dtag,label=Protein),size=6)+
    geom_point(data=subset(big_df,dtag>=3*ut),
      aes(x=ut,y=dtag,color="red"))+
    theme(legend.position = "none",
          axis.title = element_text(size=18),
          axis.text = element_text(size=14),
          axis.text.x = element_text(angle=45))

p_dox_both <- ggplot(big_df)+
  geom_abline(slope=1)+
  geom_point(aes(x=dox,y=both))+
  labs(x="Dox (signal intensity)",y="Dox+dTag (signal intensity)")+
  theme_minimal()+
    scale_x_log10(limits=c(1E4,1E7))+
    scale_y_log10(limits=c(1E4,1E7))+
  geom_label_repel(data=subset(big_df,both>=2*dox),
      nudge_y=0.4,aes(x=dox,y=both,label=Protein),size=6)+
    geom_point(data=subset(big_df,both>=2*dox),
      aes(x=dox,y=both,color="red"))+
    theme(legend.position = "none",
          axis.title = element_text(size=18),
          axis.text = element_text(size=14),
  axis.text.x = element_text(angle=45))
    

ggsave(plot = p_ut_dox,
filename = "/Volumes/RLadies/Manuscripts/2020_MCP/Figures/Fig4/Fig4E/Fig4E_dox_2.png",
width = 4.5, height = 6)

ggsave(plot = p_ut_dtag,
filename = "/Volumes/RLadies/Manuscripts/2020_MCP/Figures/Fig4/Fig4E/Fig4E_dtag_2.png",
width = 4.5, height = 6)

ggsave(plot = p_dox_both,
filename = "Fig4E_dox_both.png",
width = 4.5, height = 6)
