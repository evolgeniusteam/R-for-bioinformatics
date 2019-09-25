Input = ("
Time    Student  Score
Before  a         65
Before  b         75
Before  c         86
Before  d         69
Before  e         60
Before  f         81
Before  g         88
Before  h         53
Before  i         75
Before  j         73
After   a         77
After   b         98
After   c         92
After   d         77
After   e         65
After   f         77
After   g        100
After   h         73
After   i         93
After   j         75
")

scores = read.table(textConnection(Input),header=TRUE)