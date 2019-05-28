# Load libraries
import numpy
from numpy import arange

#import matplotlib
#matplotlib.use('Qt5Agg')
#import matplotlib.pyplot as pyplot

import pandas as pd
from pandas import read_csv
from pandas import set_option
from pandas.plotting import scatter_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.metrics import mean_squared_error

from sklearn.externals.joblib import dump 
from sklearn.externals.joblib import load




def get_prediction(x, y):

  inputArray = numpy.array([[x, y]])
  
  filename_scaler = 'velt_scaler.sav'
  loaded_scaler = load(filename_scaler)
  rescaledInputArray = loaded_scaler.transform(inputArray)

  filename_model = 'velt_model.sav'
  loaded_model = load(filename_model) 

  # Predict
  prediction = loaded_model.predict(rescaledInputArray)
  
  return prediction
  
