############################################ 
## % unoxidized vs. sulfoxide vs. sulfone ##
############################################

# Make a function to subset for Met peptides
met <- function(df) {
  df %>%
    filter(str_detect(Peptide,"M"))
}

# Subset for Met peptides
df0_cm <- met(df0_c)
df0.01_cm <- met(df0.01_c)
df0.05_cm <- met(df0.05_c)
df0.1_cm <- met(df0.1_c)
df0.25_cm <- met(df0.25_c)
df0.5_cm <- met(df0.5_c)
df1_cm <- met(df1_c)
df1.5_cm <- met(df1.5_c)

# Get total area for all df_cm's
total_areas <- c(
sum(df0_cm$Area),
sum(df0.01_cm$Area),
sum(df0.05_cm$Area),
sum(df0.1_cm$Area),
sum(df0.25_cm$Area),
sum(df0.5_cm$Area),
sum(df1_cm$Area),
sum(df1.5_cm$Area)
)

# Create functions for filtering oxidation states
filter_unox <- function(df) {
  df %>%
    filter(!str_detect(PTM,"Oxidation")) %>% 
    filter(!str_detect(PTM,"Dioxidation"))
}

filter_sulfox <- function(df) {
  df %>%
    filter(str_detect(PTM,"Oxidation")) 
}

filter_sulfone <- function(df) {
  df %>%
    filter(str_detect(PTM,"Dioxidation"))
}

# Find fractions of each oxidation state
unoxidized <- c(
  sum(filter_unox(df0_cm)$Area/total_areas[1]),
  sum(filter_unox(df0.01_cm)$Area/total_areas[2]),
  sum(filter_unox(df0.05_cm)$Area/total_areas[3]),
  sum(filter_unox(df0.1_cm)$Area/total_areas[4]),
  sum(filter_unox(df0.25_cm)$Area/total_areas[5]),
  sum(filter_unox(df0.5_cm)$Area/total_areas[6]),
  sum(filter_unox(df1_cm)$Area)/total_areas[7],
  sum(filter_unox(df1.5_cm)$Area)/total_areas[8]
)

sulfoxide <- c(
  sum(filter_sulfox(df0_cm)$Area/total_areas[1]),
  sum(filter_sulfox(df0.01_cm)$Area/total_areas[2]),
  sum(filter_sulfox(df0.05_cm)$Area/total_areas[3]),
  sum(filter_sulfox(df0.1_cm)$Area/total_areas[4]),
  sum(filter_sulfox(df0.25_cm)$Area/total_areas[5]),
  sum(filter_sulfox(df0.5_cm)$Area/total_areas[6]),
  sum(filter_sulfox(df1_cm)$Area)/total_areas[7],
  sum(filter_sulfox(df1.5_cm)$Area)/total_areas[8]
)

sulfone <- c(
  sum(filter_sulfone(df0_cm)$Area/total_areas[1]),
  sum(filter_sulfone(df0.01_cm)$Area/total_areas[2]),
  sum(filter_sulfone(df0.05_cm)$Area/total_areas[3]),
  sum(filter_sulfone(df0.1_cm)$Area/total_areas[4]),
  sum(filter_sulfone(df0.25_cm)$Area/total_areas[5]),
  sum(filter_sulfone(df0.5_cm)$Area/total_areas[6]),
  sum(filter_sulfone(df1_cm)$Area)/total_areas[7],
  sum(filter_sulfone(df1.5_cm)$Area)/total_areas[8]
)

df_unox <- data.frame(Ox_percent,"state"="unoxidized","fraction"=unoxidized)
df_sulfox <- data.frame(Ox_percent,"state"="sulfoxide","fraction"=sulfoxide)
df_sulfone <- data.frame(Ox_percent,"state"="sulfone","fraction"=sulfone)

df_uss <- rbind(df_unox,df_sulfox,df_sulfone)

df_0 <- df_uss %>% filter(Ox_percent==0) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_0.01 <- df_uss %>% filter(Ox_percent==0.01) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_0.05 <- df_uss %>% filter(Ox_percent==0.05) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_0.1 <- df_uss %>% filter(Ox_percent==0.1) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_0.25 <- df_uss %>% filter(Ox_percent==0.25) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_0.5 <- df_uss %>% filter(Ox_percent==0.5) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_1 <- df_uss %>% filter(Ox_percent==1) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)
df_1.5 <- df_uss %>% filter(Ox_percent==1.5) %>% mutate(sum=sum(fraction)) %>% 
  mutate("fractionC"=fraction/sum)

df_uss2 <- rbind(df_0,df_0.01,df_0.05,df_0.1,df_0.25,df_0.5,df_1,df_1.5)

p_states <- ggplot(df_uss2,aes(x=Ox_percent,y=fractionC,fill=state))+
geom_col()+
ylab("Fraction")+
xlab("Treatment Condition (% Oxidant)")+
#ggtitle("Peptide oxidation state breakdown with increasing oxidant")+
labs(fill="Oxidation state")+
theme(axis.text = element_text(size=14),
axis.title = element_text(size=14),
legend.title = element_text(size=14),
legend.text = element_text(size=14),
axis.text.x = element_text(angle=45,hjust=1))

# Save image
ggsave(plot=p_states,
filename="Fig2D.png",
width = 4.5, height = 5)
