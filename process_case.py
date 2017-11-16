import pandas as pd
import numpy as np


data = pd.read_csv('dataset/transposed_case_na.csv')
data.replace('NA', np.nan)
data.replace(to_replace=np.nan, value=0.0, inplace=True)
data.replace(to_replace='NA.54E-18', value=0.0, inplace=True)
data.to_csv('dataset/transposed_case_0.csv', index=None)
