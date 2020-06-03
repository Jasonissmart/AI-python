from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

dff = pd.read_excel('sw結果.xlsx',sheet_name='日期')

# 要40年SW




lable = dff.columns

df = pd.read_excel('台灣公債資料.xlsx',sheet_name='1-10')
#df = df[:60]
df = pd.DataFrame(df)
name = df.columns

x = [1,2,3,4,5,6,7,8,9,10,30,40] # For T=40 Case
xnew = np.linspace(1, 40, num= 40, endpoint=True)
result = []
for i in range(1, 244, 1):
    y = df[name[i]].tolist()
    y.append(0.0379) #UFR-10bp
    y.append(0.038) # UFR
    f2 = interp1d(x, y, kind='cubic')
    result.append(f2(xnew))


result_array = np.array(result).reshape(243,40)
df2 = pd.DataFrame(result_array)
df2.to_csv('cubic_spline_normal_forward_3.8_40.csv')

df6 = pd.read_csv('cubic_spline_normal_forward_3.8_40.csv')
df6.set_index('Unnamed: 0',inplace=True)
df6 = df6.T
name = df6.columns

data = []

for i in range(0,235,6):
    plt.style.use('ggplot')
    plt.figure(num = 3,figsize= (20,13))
    plt.title('UFR = 3.8% ; Normal Interpolation ; Cubic Spline')
    plt.xlabel('year', fontsize=30)
    plt.ylabel('spot rate', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    data = df6[name[i]].values
    plt.plot(data,label = lable[i])
    #plt.legend()



plt.savefig('UFR = 3.8% ; Normal Interpolation ; Cubic Spline.png')