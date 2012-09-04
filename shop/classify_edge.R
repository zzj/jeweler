library('MASS')
library('glmnet')
library('e1071')
library('class')
library('LiblineaR')


data = read.table("../training", stringsAsFactors = F)
print(dim(data))
Y = data[,1]
X = data[,2:dim(data)[2]]
a = LiblineaR(scale(X), (Y),cross = 10, type = 4)

a = LiblineaR(scale(X), (Y),cross = 0, type = 4)
b = predict(a,scale(X))
