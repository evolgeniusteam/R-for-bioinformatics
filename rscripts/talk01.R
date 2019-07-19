### --  装入必要的libraries --
library(RMySQL);
library(dplyr);

### -- 指定用户名，密码，数据库名
mysql.dbname = "r4ds_test";
dbCon <- dbConnect(MySQL(), user="r4ds", password="r4ds", dbname=mysql.dbname );

# Run query to get results as dataframe
dat <- dbGetQuery(dbCon, "SELECT * FROM grades");

## -- 任务： 为每个人计算：
## -- 平均成绩、上课总数、及格门数、不及格门数

stats <- dat %>% group_by(name) %>% 
  summarise( avg_grade = mean(grade), course_count = n(), 
             course_pass = sum(  grade >= 60  ), course_fail = sum( grade  < 60  ) ) %>% 
  arrange( -avg_grade );

## -- 显示计算结果 --
stats;
