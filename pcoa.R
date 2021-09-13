require(vegan);
dat.dist <- vegdist( dat, method = "bray" );

require(ape);
PCOA <- pcoa(dat.dist);
rel.eig <- PCOA$values$Relative_eig[1:10];
val <- data.frame( PCOA$vectors[, 1:2] );

plot_pcoa <- 
  ggscatter(val, x = "Axis.1", y = "Axis.2", 
            ellipse = TRUE, mean.point = TRUE, star.plot = TRUE, 
            alpha = 0.6, ggtheme = theme_bw()) + 
  xlab( paste0( "Axis.1 (", round( rel.eig[1] * 100, digits = 2), "%)") ) +
  ylab( paste0( "Axis.2 (", round( rel.eig[2] * 100, digits = 2), "%)") );
