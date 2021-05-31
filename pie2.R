library(ggplot2)
library(dplyr)
library(tidytidbits)
library(readr)
zy4t<- read_delim(file = "data/talk06/relative_abundance_for_RUN_ERR1072629_taxonlevel_species.txt",delim="\t",quote="",comment="#")
zy4<-zy4t%>%arrange(desc(relative_abundance))%>%lump_rows(scientific_name,relative_abundance,n=10,other_level = "others")
zy4<-zy4[c(1,3:10,2,11),]
zy4$fr<-zy4$relative_abundance/sum(zy4$relative_abundance)
zy4$ymax<-cumsum(zy4$fr)
zy4$ymin<-c(0,head(zy4$ymax,n=-1))
zy4$lp<-(zy4$ymax+zy4$ymin)/2
zy4$l<-paste(round(as.vector(zy4$relative_abundance),3),"%",sep = " ")
ggplot(zy4,aes(ymax=ymax,ymin=ymin,xmax=4,xmin=3))+geom_rect(aes(fill=scientific_name))+geom_text(x=3.5,aes(y=lp,label=l),size=3)+coord_polar(theta = "y")+xlim(2,4)+theme_void()
