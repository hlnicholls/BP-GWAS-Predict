import sys
import pandas as pd
import numpy as np
from numpy import sort
from scipy.stats import spearmanr
from scipy.cluster import hierarchy
import scipy.cluster
from numpy import absolute, mean, sort, std
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

import re
regex = re.compile(r"\[|\]|<", re.IGNORECASE)

from sklearn import datasets, metrics, preprocessing, model_selection
import sklearn.neighbors._base
sys.modules['sklearn.neighbors.base'] = sklearn.neighbors._base
from sklearn.preprocessing import MinMaxScaler,StandardScaler
from sklearn.model_selection import train_test_split, KFold,RepeatedKFold, StratifiedKFold, cross_val_score, cross_validate, cross_val_predict, GridSearchCV, RandomizedSearchCV, validation_curve, learning_curve
from sklearn.metrics import mean_squared_error, r2_score, explained_variance_score, mean_absolute_error, max_error

import skopt
from skopt import BayesSearchCV 

from missingpy import MissForest

import shap
from BorutaShap import BorutaShap

import xgboost
import lightgbm
from catboost import CatBoostClassifier
from lightgbm import LGBMClassifier
from sklearn.linear_model import LinearRegression, Lasso, ElasticNet
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, VotingClassifier, StackingClassifier, BaggingClassifier, ExtraTreesClassifier

from sklearn.svm import SVR


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from warnings import filterwarnings
filterwarnings("ignore")

import matplotlib.pyplot as plt
import missingno as msno

seed = 0

dataset = pd.read_csv("training_data.txt", sep="\t")
dataset = dataset.rename({'IPA_BP':'IPA_BP_annotation'}, axis=1)

data = dataset.drop(
    ["Gene", ], 1
)  

data["BPlabel_encoded"] = data["BPlabel"].map(
    {"most likely": 1, "probable": 2, "least likely": 3}
)
Y = data["BPlabel_encoded"]


data = pd.read_csv("training_cleaned.csv", header=0, sep=",")

data["BPlabel_encoded"] = data["BPlabel"].map(
     {"most likely": 1, "probable": 2, "least likely": 3}
)
Y = data["BPlabel_encoded"]
data = data.drop(["BPlabel"], 1)

X = pd.read_csv("cleaned_imputed_training_data.csv", header=0)
X.columns = [
    regex.sub("_", col) if any(x in str(col) for x in set(("[", "]", "<"))) else col
    for col in X.columns.values
]

X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=0.2, random_state=seed
)

xgbr = xgboost.XGBClassifier(random_state=seed, objective='reg:squarederror', verbosity = 0, eval_metric='mlogloss') 
xgbr_params = {
    'max_depth':  (1, 4), 
    'learning_rate': (0.01, 0.2, 'log-uniform'),  
    'n_estimators':  (10, 50), 
    'reg_alpha':  (1, 10, 'log-uniform'), 
    'reg_lambda':  (1, 10, 'log-uniform')} 

lgbm = LGBMClassifier(random_state=seed)
lgbm_params = {
    "max_depth": (1, 4),
    "learning_rate": (0.01, 0.2, "log-uniform"),
    "n_estimators": (10, 50),
    "reg_alpha": (1, 10, "log-uniform"),
    "reg_lambda": (1, 10, "log-uniform"),
}

catboost = CatBoostClassifier(random_seed=seed, verbose=False)
cat_params = {
     "iterations": (10, 50),
     'learning_rate': (0.01, 0.2, 'log-uniform'), 
     'depth':  (1, 4), 
}


gbr = GradientBoostingClassifier(random_state=seed)
gbr_params = {
    'learning_rate': (0.01, 0.2),
    'max_depth': (1, 4),
    "max_features":["log2","sqrt", "auto"],
    "criterion": ["friedman_mse", "squared_error"],
    'n_estimators': (10, 50)
    }

rfr = RandomForestClassifier(random_state=seed)
rfr_params={'n_estimators': (10, 50), 
             'max_features': ['sqrt', 'log2'],
             'max_depth' : (1, 4),
             'criterion' :['gini', 'entropy']} 

dt = DecisionTreeClassifier(random_state=seed)
dt_params= {"criterion": ['gini', 'entropy'],
            'max_features': ['sqrt', 'log2'],
            'max_depth' : (1, 4)}

extra = ExtraTreesClassifier(random_state=seed)
extra_params ={'n_estimators': (10, 50), 
             'max_features': ['sqrt', 'log2'],
             'max_depth' : (1, 4),
             'criterion' :['gini', 'entropy']}


inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)

models = []

models.append(('XGBR', BayesSearchCV(xgbr, xgbr_params, cv=inner_cv,iid=False,n_jobs=1, random_state=seed))) 
models.append(("LGBM", BayesSearchCV(lgbm, lgbm_params, cv=inner_cv, iid=False, n_jobs=1, random_state=seed)))
models.append(("CB", BayesSearchCV(catboost, cat_params, cv=inner_cv, iid=False, n_jobs=1, random_state=seed)))
models.append(('GBR', BayesSearchCV(gbr, gbr_params, cv=inner_cv,iid=False, n_jobs=1, random_state=seed)))
models.append(('RFR', BayesSearchCV(rfr, rfr_params, cv=inner_cv,iid=False, n_jobs=1, random_state=seed)))
models.append(('DT', BayesSearchCV(dt, dt_params, cv=inner_cv, iid=False, n_jobs=1, random_state=seed)))
models.append(('ExtraTrees', BayesSearchCV(extra, extra_params, cv=inner_cv, iid=False, n_jobs=1, random_state=seed)))


results = []
names = []
medians =[]
scoring = ['accuracy', 'balanced_accuracy', 'f1_weighted', 
          'precision_weighted','recall_weighted']

models_list_balancedac = []


for name, model in models:
    nested_cv_results = model_selection.cross_validate(model, X , Y, cv=outer_cv, scoring=scoring, error_score="raise")
    nested_cv_results2 = model_selection.cross_val_score(model, X , Y, cv=outer_cv, scoring='balanced_accuracy', error_score="raise")
    results.append(nested_cv_results2)
    names.append(name)
    print(name, 'Nested CV results for all scores:', '\n', nested_cv_results, '\n')
    print(name, 'Accuracy Nested CV Average', np.mean(nested_cv_results['test_accuracy']))
    print(name, 'Balanced Accuracy Nested CV Average', np.mean(nested_cv_results['test_balanced_accuracy'] ))
    print(name, 'F1 Nested CV Average', np.mean(nested_cv_results['test_f1_weighted'] ))
    print(name, 'Precision Nested CV Average', np.mean(nested_cv_results['test_precision_weighted'] ))
    print(name, 'Recall Nested CV Average', np.mean(nested_cv_results['test_recall_weighted'] ))
    model.fit(X, Y)
    print('\n')
    print("Best Parameters: \n{}\n".format(model.best_params_))
    print("Best Estimator:", model.best_estimator_)
    best_model = model.best_estimator_
    print('\n')
    print('Non-nested CV Results:')
    best_model.fit(X_train, Y_train)
    y_pred_train = best_model.predict(X_train)
    y_pred = best_model.predict(X_test)
    best_model.fit(X, Y)
    median_balancedac = np.median(nested_cv_results['test_balanced_accuracy'])
    models_list_balancedac.append((best_model, median_balancedac))


print('All r2 results:', results)         

best_model1, best_balancedac = sorted(models_list_balancedac, key = lambda x: x[1], reverse=True)[0]
print('Best model by median balanced accuracy:',best_model1)

Feature_Selector = BorutaShap(model=best_model1, importance_measure="shap", classification=False)

Feature_Selector.fit(X=X, y=Y, n_trials=200, random_state=seed)

subset = Feature_Selector.Subset()
X_boruta_sel = subset
X_boruta_sel.to_csv(r"selected_features_training_data.csv", index=False)

Feature_Selector.plot(which_features="accepted")