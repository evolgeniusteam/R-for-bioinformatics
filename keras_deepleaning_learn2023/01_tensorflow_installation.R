## install mini conda; if not done
reticulate::install_miniconda(); ## /Users/wchen/Library/r-miniconda-arm64

## install tensorflow
Sys.setenv("RETICULATE_PYTHON" = "/Users/wchen/Library/r-miniconda-arm64/envs/r-reticulate/bin/python");
library(reticulate);

## install tensorflow --
library(tensorflow); ## install this package if not done!!
install_tensorflow(envname = "r-reticulate", version = "cpu");

## done!!