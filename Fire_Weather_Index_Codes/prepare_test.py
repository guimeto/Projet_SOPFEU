# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 10:06:07 2019

@author: guillaume
"""
import pandas as pd      # aliasing as pd
import numpy as np
file1 = "K:/PROJETS/PROJET_FIRE_INDEX/scripts/testBatch.csv"
col_names = ['Date', 'Temp.', 'RH', 'Wind', 'Rain', 'FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI']
df = pd.read_csv(file1, skiprows=1,  names=col_names)

df["Month"] = df.iloc[:,0].str.split('/', expand=True)[0]
df["Day"] = df.iloc[:,0].str.split('/', expand=True)[1]
df["Year"] = df.iloc[:,0].str.split('/', expand=True)[2]

new_dataset = df[['Month','Day','Temp.','RH','Wind','Rain']]

new_dataset.to_csv('data.txt', sep='\t', index=False)
