library(tidyverse);
mt <- mtcars %>% rownames_to_column() %>% arrange( -cyl, -mpg  );

mt$rowname <- factor(mt$rowname, levels = as.character(mt$rowname));

mt %>% ggplot( aes( x = rowname, y = mpg, fill = factor(cyl) ) ) + geom_bar( stat = "identity" );

