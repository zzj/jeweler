from __future__ import print_function
from sklearn import svm
from sklearn import tree
from sklearn import cross_validation
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from sklearn.metrics import zero_one_score
from sklearn.neighbors.nearest_centroid import NearestCentroid
from sklearn.ensemble import RandomForestClassifier

import os
import sys
import glob
import pickle

X = []
Y = []
for infile in glob.glob(os.path.join("simulation_data/", "*")):
    (x, y) = pickle.load(open(infile, "rb"))
    X.extend(x)
    Y.extend(y)
half = len(X) * 9 / 10
clf = svm.LinearSVC()
clf = tree.DecisionTreeClassifier()
clf = RandomForestClassifier(n_estimators=20)
fit = clf.fit(X[0:half], Y[0:half])
X = X[half:]
Y = Y[half:]
predict = fit.predict(X)
print(zero_one_score(Y, predict))
print(precision_score(Y, predict, pos_label=0))
print(recall_score(Y, predict, pos_label=0))
print(f1_score(Y, predict, pos_label=0))
pickle.dump(fit, open("learning_model", "wb"))

# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=precision_score)
# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=recall_score)
# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=f1_score)
# print(sum(scores)/len(scores))
