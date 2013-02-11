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
import random


def get_data(files):
    X = []
    Y = []
    for infile in files:
        (x, y) = pickle.load(open(infile, "rb"))
        X.extend(x)
        Y.extend(y)
    return (X, Y)
data = glob.glob(os.path.join("simulation_data/", "*"))
accuracy = 0
precision = 0
recall = 0
f1 = 0
area = 0
print(len(data))
for i in range(5):
    learning = data[:-55]
    training = data[-55:]
    data = data[-55:] + data[:-55]
    X, Y = get_data(learning)
    half = len(X)
    clf = RandomForestClassifier(n_estimators=20)
    clf = tree.DecisionTreeClassifier()
    clf = svm.LinearSVC()
    clf = svm.SVC(kernel='linear', probability=True)
    print(len(Y))
    fit = clf.fit(X, Y)
    X, Y = get_data(training)
    predict = fit.predict(X)
    accuracy += zero_one_score(Y, predict)
    precision += precision_score(Y, predict, pos_label=0)
    recall += recall_score(Y, predict, pos_label=0)
    f1 += f1_score(Y, predict, pos_label=0)
    print(zero_one_score(Y, predict))
    print(precision_score(Y, predict, pos_label=0))
    print(recall_score(Y, predict, pos_label=0))
    print(f1_score(Y, predict, pos_label=0))
    pickle.dump(fit, open("learning_model_cv" + str(i), "wb"))
    predict_prob = fit.predict_proba(X)
    precision_curve, recall_curve, thresholds = precision_recall_curve(Y, predict_prob[:,1])
    area += auc(recall_curve, precision_curve)
    print("area is %lf" % area)

print(accuracy / 5)
print(precision / 5)
print(recall / 5)
print(f1 / 5)
print(area / 5)
#pickle.dump(fit, open("learning_model", "wb"))

# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=precision_score)
# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=recall_score)
# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=f1_score)
# print(sum(scores)/len(scores))
