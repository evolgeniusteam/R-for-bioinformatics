Sys.setenv("RETICULATE_PYTHON" = "/Users/wchen/Library/r-miniconda-arm64/envs/r-reticulate/bin/python");
library(reticulate);

## -- 
library(tensorflow);
library(keras);
mnist <- dataset_mnist();

train_images <- mnist$train$x;
train_labels <- mnist$train$y;


test_images <- mnist$test$x;
test_labels <- mnist$test$y;

## -- contents of the variables 
str(train_images);
dim(train_images);

str(train_labels);

network <- keras_model_sequential() %>% 
  layer_dense(units = 512, activation = "relu", input_shape = c(28 * 28)) %>% 
  layer_dense(units = 10, activation = "softmax");

network %>% compile(
  optimizer = "rmsprop", 
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
);

## 2.4 -
train_images <- array_reshape(train_images, c(60000, 28 * 28));
train_images <- train_images / 255;

test_images <- array_reshape(test_images, c(10000, 28 * 28));
test_images <- test_images / 255;

## 2.5 --
train_labels <- to_categorical(train_labels);
test_labels <- to_categorical(test_labels);

network %>% fit(train_images, train_labels, epochs = 5, batch_size = 128);
