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
from sklearn.ensemble import ExtraTreesClassifier

import os
import sys
import glob
import pickle
import random
import numpy

def get_data(files):
    X = []
    Y = []
    for infile in files:
        (x, y) = pickle.load(open(infile, "rb"))
        X.extend(x)
        Y.extend(y)
    return (X, Y)

def get_features(X, selected):
    if selected:
        new_X = numpy.matrix(X)
        return new_X[:,selected]
    return X

def get_clf():
    clf = tree.DecisionTreeClassifier()
    clf = svm.SVC(kernel='linear', probability=True)
    clf = svm.LinearSVC()
    clf = RandomForestClassifier(n_estimators=40)
    return clf

data = glob.glob(os.path.join("simulation_data/", "*"))
accuracy = 0
precision = 0
recall = 0
f1 = 0
area = 0
total_files = len(data)
N = 11

selected = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
selected = None

num_folds = total_files / N
for i in range(num_folds):
    train = data[N:]
    test = data[:N]
    data = data[N:] + data[:N]
    X, Y = get_data(train)
    X = get_features(X, selected)
    print(sum(Y))
    print(len(Y))
    predict = [y for y in Y]
    random.shuffle(predict)

    print(precision_score(Y, predict, pos_label=0))
    print(recall_score(Y, predict, pos_label=0))
    print(f1_score(Y, predict, pos_label=0))
    clf = get_clf()
    fit = clf.fit(X, Y)
    X, Y = get_data(test)
    X = get_features(X, selected)
    predict = fit.predict(X)
    accuracy += zero_one_score(Y, predict)
    precision += precision_score(Y, predict, pos_label=0)
    recall += recall_score(Y, predict, pos_label=0)
    f1 += f1_score(Y, predict, pos_label=0)
    pickle.dump(fit, open("learning_model_cv" + str(i), "wb"))
    predict_prob = fit.predict_proba(X)
    precision_curve, recall_curve, thresholds = precision_recall_curve(Y, predict_prob[:,1])
    area += auc(recall_curve, precision_curve)

print(accuracy / num_folds)
print(precision / num_folds)
print(recall / num_folds)
print(f1 / num_folds)
print(area / num_folds)
#pickle.dump(fit, open("learning_model", "wb"))

# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=precision_score)
# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=recall_score)
# scores = cross_validation.cross_val_score(clf, X, Y, cv=10, score_func=f1_score)
# print(sum(scores)/len(scores))
