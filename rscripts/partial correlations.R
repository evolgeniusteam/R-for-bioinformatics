
if(!require(ppcor)){
  install.packages("ppcor");
}
library(ppcor);


with( mtcars, cor.test(mpg, cyl) );
with( mtcars, cor.test(mpg, wt) );
with( mtcars, pcor.test(mpg, cyl, wt, method = "spearman") );
with( mtcars, pcor.test(mpg, cyl, wt, method = "pearson") );

## get data ready --
mpg_resid <- resid(lm(data=mtcars, mpg ~ wt))
cyl_resid <- resid(lm(data=mtcars, cyl ~ wt))
cor.test(mpg_resid, cyl_resid, method = "spearman");

## plot
library(ggplot2);
ggplot(data=NULL) + geom_point(aes(x=mpg_resid, y=cyl_resid))

