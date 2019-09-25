Input = ("
Student  Sex     Teacher  Steps  Rating
a        female  Catbus    8000   7
b        female  Catbus    9000  10
c        female  Catbus   10000   9
d        female  Catbus    7000   5
e        female  Catbus    6000   4
f        female  Catbus    8000   8
g        male    Catbus    7000   6
h        male    Catbus    5000   5
i        male    Catbus    9000  10
j        male    Catbus    7000   8
k        female  Satsuki   8000   7
l        female  Satsuki   9000   8
m        female  Satsuki   9000   8
n        female  Satsuki   8000   9
o        male    Satsuki   6000   5
p        male    Satsuki   8000   9
q        male    Satsuki   7000   6
r        female  Totoro   10000  10
s        female  Totoro    9000  10
t        female  Totoro    8000   8
u        female  Totoro    8000   7
v        female  Totoro    6000   7
w        male    Totoro    6000   8
x        male    Totoro    8000  10
y        male    Totoro    7000   7
z        male    Totoro    7000   7
")

Data = read.table(textConnection(Input),header=TRUE);

