# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_pdf.pdf
library(knitr)
library(kableExtra)
Small <- c("130ug","65ug")
Large <- c("450ug","125ug")
Characteristic <- c("Capacity, Antibody (Small)",
                    "Capacity, Antibody (Large)",
                    "Capacity, Complex (Large)",
                    "MHC-I in 250M GRANTA",
                    "MHC-I/mL = MHC-I/cartridge")
Value <- c("130",
           "450",
           "73",
           "36",
           "5")
#rownames <- c("Antibody","Recombinant MHC*")
colnames <- c("Characteristic","Value (ug)")
#df <- data.frame(Small,Large)
df <- data.frame(Characteristic,Value)
#rownames(df) <- rownames
colnames(df) <- colnames
kable(df,"latex", booktabs = T, escape = F) %>%
  kable_styling(full_width = F) %>%
#footnote(symbol = "Against antibody cartridges at full capacity.") %>% 
  as_image()
