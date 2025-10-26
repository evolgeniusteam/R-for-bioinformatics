library(tidyverse);
library(modeldata);

## show me how to use dplyr group_by function
data("mtcars");
mtcars %>% group_by(cyl) %>% summarise( mean_mpg = mean(mpg), n = n() );

