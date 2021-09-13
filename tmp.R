library(caret);
cv_lm.swiss <- train(
  Fertility ~ ., 
  data = swiss, 
  method = "lm",
  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 10)
);

cv_mars.swiss <- train(
  Fertility ~.,
  data = swiss,
  metric = "RMSE", 
  method = "earth",
  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 10)
)

data.frame(
  RMSE = c( RMSE( swiss$Fertility, predict( cv_lm.swiss, swiss ) ), RMSE( swiss$Fertility, predict( cv_mars.swiss, swiss ) ) ),
  R2 = c( R2( swiss$Fertility, predict( cv_lm.swiss, swiss ) ), R2( swiss$Fertility, predict( cv_mars.swiss, swiss ) ) ),
  row.names = c("LM", "MARS")
)

library(vip);
p1 <- vip(cv_lm.swiss, num_features = 5, geom = "point", value = "gcv") + ggtitle("LM:GCV")
p2 <- vip(cv_mars.swiss, num_features = 5, geom = "point", value = "gcv") + ggtitle("MARS:GCV")

gridExtra::grid.arrange(p1, p2, ncol = 2)
