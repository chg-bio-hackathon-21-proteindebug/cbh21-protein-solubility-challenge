import os
import sklearn
import joblib
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn import svm

from sklearn.ensemble import RandomForestRegressor as rf


model = svm.SVR()
model = rf()
train_data = pd.read_csv("pdb2csv_v3.csv")
sol_val = train_data["Solubility Score"]
print(train_data.head())
train_data.drop(["Solubility Score", "PDB File", "Sequence"], inplace=True, axis = 1)
fit_mod = model.fit(train_data,sol_val)
cv_score = cross_val_score(model,train_data,sol_val, cv = 5)
print(cv_score)
joblib.dump(fit_mod,"first_model.bin")
