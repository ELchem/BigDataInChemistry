# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 21:11:14 2021

@author: Dell
"""

import matplotlib.pyplot as plt
import seaborn as sns; sns.set() #set the style
import scipy.stats as stats
import numpy as np
import pandas as pd

from sklearn.metrics import mean_squared_error
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

from scipy import stats
from scipy.stats import ttest_ind

import statsmodels.api as sm
from statsmodels.graphics.gofplots import ProbPlot
Filename = "delaney-processed.csv"
ESOL = pd.read_csv(Filename)

#pairplot
#sns.pairplot(ESOL)
#plt.savefig("pairplot.png")
#plt.show()

#OLS
#regression
Y = ESOL['ESOL predicted log solubility in mols per litre']
X = ESOL['measured log solubility in mols per litre']
model = sm.OLS(Y, sm.add_constant(X))
model_fit = model.fit()
print(model_fit.summary())

#scatter + line
#fig,ax = plt.subplots()
#ax.scatter(X,Y)
#ax.plot(X,model_fit.fittedvalues)
#plt.xlabel('Measured Sol')
#plt.ylabel('Predicted Sol')
#plt.title("Measured vs Predicted Sol")
#plt.savefig("regression1.png")
#plt.show()
#residuals, stu residuals and lev
# model values
model_fitted_y = model_fit.fittedvalues
# model residuals
model_residuals = model_fit.resid
# normalized residuals
model_norm_residuals = model_fit.get_influence().resid_studentized_internal
# absolute squared normalized residuals
model_norm_residuals_abs_sqrt = np.sqrt(np.abs(model_norm_residuals))
# absolute residuals
model_abs_resid = np.abs(model_residuals)
# leverage, from statsmodels internals
model_leverage = model_fit.get_influence().hat_matrix_diag

#fig2,axes = plt.subplots(1,2, figsize = (11,8))
#axes[0].scatter(model_fitted_y,model_residuals)
#axes[0].set_xlabel("fitted values")
#axes[0].set_ylabel("residuals")
#axes[0].set_title("fitted vs residuals")
#
#
#axes[1].scatter(model_norm_residuals,model_leverage)
#axes[1].set_xlabel("normalized residuals")
#axes[1].set_ylabel("leverage")
#axes[1].set_title("normalized residuals vs levereage")
#plt.savefig("fig4.png")
#plt.show()
#Multiple linear regression
parameters = ['Molecular Weight','Number of H-Bond Donors','Number of Rings',
              'Number of Rotatable Bonds','Polar Surface Area']
X2 = ESOL[parameters]
Y2 = ESOL['measured log solubility in mols per litre']
model2 = sm.OLS(Y2, sm.add_constant(X2))
model_fit2 = model2.fit()
print(model_fit2.summary())

training_set_x, test_set_x, training_set_y, test_set_y = train_test_split(ESOL[parameters],
                                                                          ESOL['measured log solubility in mols per litre'], 
                                                                          test_size=0.2,
                                                                          random_state=0)


model_fit2 = LinearRegression().fit(training_set_x, training_set_y)
fit2_vals_test = model_fit2.predict(test_set_x)
print("R^2 value for training data:",model_fit2.score(training_set_x,training_set_y))
print("RSE for training data:",mean_squared_error(training_set_y, model_fit2.predict(training_set_x)))
print("R^2 value for testing data:",model_fit2.score(test_set_x,test_set_y))
print("RSE for test data:",mean_squared_error(test_set_y, fit2_vals_test))

for i in range(len(parameters)):
    temp_params = parameters.copy()
    del temp_params[i]
    X3 = ESOL[temp_params]
    Y3 = ESOL['measured log solubility in mols per litre']
    model3 = sm.OLS(Y3, sm.add_constant(X3))
    model_fit3 = model3.fit()
    print(model_fit3.summary())
    
    training_set_x2, test_set_x2, training_set_y2, test_set_y2 = train_test_split(ESOL[temp_params],
                                                                              ESOL['measured log solubility in mols per litre'], 
                                                                              test_size=0.8,
                                                                              random_state=0)
    model_fit2 = LinearRegression().fit(training_set_x2, training_set_y2)
    fit2_vals_test = model_fit2.predict(test_set_x2)
    print(parameters[i],"removed")
    print("R^2 value for training data:",model_fit2.score(training_set_x2,training_set_y2))
    print("RSE for training data:",mean_squared_error(training_set_y2, model_fit2.predict(training_set_x2)))
    print("R^2 value for testing data:",model_fit2.score(test_set_x2,test_set_y2))
    print("RSE for test data:",mean_squared_error(test_set_y2, fit2_vals_test))
    





