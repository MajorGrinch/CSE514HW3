import pandas as pd
import numpy as np

data = pd.read_csv('dataset/temp_ctrl.csv')

# convert = pd.DataFrame(0, index=np.arange(len(data)), columns=data.columns)
for x in data.columns:
    if data[x].dtype == np.object:
        print(x)
data.replace('NaN', np.nan)
data.replace(to_replace='NA.54E-18', value=0.0, inplace=True)
# data['GI_16904380-S'].astype('float64').tolist()
data['WGACON89'] = data['WGACON89'].astype('float64')
data.to_csv('dataset/new_ctrl.csv', index=None)
# eighteenth = [data[x].nlargest(18).tolist()[-1] for x in data.columns]
# # print(data.iloc[:, 0])
# convert = pd.DataFrame(0, index=np.arange(len(data)), columns=data.columns)
# for i in range(1800):
#     convert.iloc[:, i] = data.iloc[:, i] >= eighteenth[i]
# result = convert[:].astype('int')
# # print result
# for i in range(0, 188):
#     for j in range(0, 1800):
#         if result.iloc[i, j] == 1:
#             result.iloc[i, j] = j + 1
#         else:
#             result.iloc[i, j] = 0
# result.to_csv('processed_ctrl.csv', index=None)
# result.iloc[:, :10].to_csv('first10_genes.csv', index=None)
