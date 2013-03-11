from __future__ import print_function
from sklearn import svm
from sklearn import tree
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

def main():
    data = glob.glob(os.path.join("simulation_data/", "*"))
    if (len(sys.argv) > 0):
        for i in range(1, len(sys.argv)):
            data.remove(sys.argv[i])
    X, Y = get_data(data)
    clf = RandomForestClassifier(n_estimators=20)
    fit = clf.fit(X, Y)
    pickle.dump(fit, open("learning_model_all", "wb"))

if __name__ == "__main__":
    main()