abu <- read_delim(file = "data/talk06/relative_abundance_for_RUN_ERR1072629_taxonlevel_species.txt",
                  delim = "\t", quote = "", comment = "#");

library(tidytidbits);
abu.dat <- 
  abu %>% arrange( desc( relative_abundance ) ) %>% 
  lump_rows( scientific_name, relative_abundance, n = 10, other_level = "Others" );

## 难点：先确定排序，再计算cumsum --
abu.dat$scientific_name <- reorder(abu.dat$scientific_name, abu.dat$relative_abundance, sort);
lvls <- levels( abu.dat$scientific_name );
lvls <- c( lvls[ ! lvls %in% c("Others", "Unknown") ], c("Others", "Unknown") );

abu.dat$scientific_name <- factor( abu.dat$scientific_name, levels = lvls );

abu.dat <- abu.dat %>% 
  arrange( scientific_name ) %>% 
  mutate(cumsum = cumsum( relative_abundance ) - 0.5 * relative_abundance);

ggplot(abu.dat, aes(x = 2, y = relative_abundance, 
                    fill = scientific_name)) +
  geom_bar( stat = "identity", color = "white" ) +
  coord_polar(theta = "y", start = 30)+
  geom_text(aes(y = 100 - cumsum, label = str_c( round(relative_abundance, 2), "%")), color = "white")+
  theme_void()+
  xlim(0.5, 2.5)

