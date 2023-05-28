import re

import numpy as np
import pandas as pd

regex = re.compile(r"\[|\]|<", re.IGNORECASE)

from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
)
from sklearn.metrics import *
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_score,
    cross_validate,
    train_test_split,
)

from sklearn.preprocessing import MinMaxScaler
from skopt import BayesSearchCV

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
from warnings import filterwarnings
warnings.filterwarnings('ignore')

filterwarnings("ignore")

seed = 0


import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras.wrappers.scikit_learn import KerasClassifier
from keras.optimizers import *



data = pd.read_csv("training_cleaned.csv", header=0, sep=",")

data['BPlabel_encoded'] = data['BPlabel'].map( {'most likely':0,'probable':1, 'least likely':2})
Y = data["BPlabel_encoded"] 
data = data.drop(["BPlabel"],1)
data.shape 

X = pd.read_csv("selected_features_training_data.csv", header=0)
X.columns = [
    regex.sub("_", col) if any(x in str(col) for x in set(("[", "]", "<"))) else col
    for col in X.columns.values
]

X = MinMaxScaler().fit_transform(X)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=seed)

def baseline_model(optimizer='SGD', learning_rate=0.01, momentum=0, activation='relu', dropout_rate=0.0, weight_constraint=0, neurons=1):
    model = Sequential()
    model.add(Dense(neurons, input_dim=X.shape[1], activation=activation)) #dense layers perform: output = activation(dot(input, kernel) + bias).
    model.add(Dropout(dropout_rate))
    model.add(Dense(neurons, activation=activation)) #8 is the dim/ the number of hidden units (units are the kernel)
    model.add(Dense(3, activation='softmax'))
    model.compile(loss='sparse_categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    return model

keras = KerasClassifier(build_fn=baseline_model, batch_size=32, epochs=10, verbose=0)

optimizer = ['SGD','RMSprop','adam']
activation = ['softmax', 'relu', 'tanh'] 
lr_rate = (0.01, 0.1)
momentum = (0.3, 0.6, 0.9)
dropout_rate = (0.0, 0.1, 0.2)
neurons = (50, 160, 200)
keras_params = dict(optimizer=optimizer, learning_rate=lr_rate, momentum=momentum, activation=activation, 
                    dropout_rate=dropout_rate, neurons=neurons)


inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)

clf = BayesSearchCV(keras, keras_params, cv=inner_cv,iid=False, n_jobs=1)

scoring = ['accuracy', 'balanced_accuracy', 'f1_weighted', 
          'precision_weighted','recall_weighted']

nested_cv_results = cross_validate(clf, X , Y, cv=outer_cv, scoring=scoring, error_score="raise")
nested_cv_results2 = cross_val_score(clf, X , Y, cv=outer_cv, scoring='balanced_accuracy', error_score="raise")

print( 'Nested CV results for all scores:', '\n', nested_cv_results, '\n')
print( 'Accuracy Nested CV Average', np.median(nested_cv_results['test_accuracy']))
print( 'Balanced Accuracy Nested CV Average', np.median(nested_cv_results['test_balanced_accuracy'] ))
print( 'F1 Nested CV Average', np.median(nested_cv_results['test_f1_weighted'] ))
print( 'Precision Nested CV Average', np.median(nested_cv_results['test_precision_weighted'] ))
print( 'Recall Nested CV Average', np.median(nested_cv_results['test_recall_weighted'] ))
clf.fit(X_train, Y_train)
print("Best Parameters: \n{}\n".format(clf.best_params_))
print('Non-nested CV Results:')
y_pred_train = clf.predict(X_train)
y_pred = clf.predict(X_test)
print( 'Train accuracy:', accuracy_score(Y_train, y_pred_train), 'Test accuracy:', accuracy_score(Y_test, y_pred))
print( 'Train balanced accuracy:', balanced_accuracy_score(Y_train, y_pred_train), 'Test balanced accuracy:', balanced_accuracy_score(Y_test, y_pred))
print( 'Train F1', f1_score(Y_train, y_pred_train, average='weighted'), 'Test F1:', f1_score(Y_test, y_pred, average='weighted'))
print( 'Train recall:', recall_score(Y_train, y_pred_train, average='weighted'),'Test recall:', recall_score(Y_test, y_pred,average='weighted'))
print( 'Train precision:', precision_score(Y_train, y_pred_train,average='weighted'), 'Test precision:', precision_score(Y_test, y_pred,average='weighted'))
