import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dirName = 'output'
k = 0
rateType = 'forward'
df = pd.read_excel(dirName+'/output'+str(k)+'_'+rateType+'.xlsx')
df.set_index('Year',inplace=True)

name = df.columns

x = '4.4% ; 1510 Interpolation_60' #儲存圖片的檔名＆圖片的大標名稱
for i in range(0,len(name),6):  #間隔6個月畫一條線
    plt.style.use('ggplot')
    plt.figure(num = 3,figsize= (20,13))
    plt.title(x)

    plt.xlabel('year')
    plt.ylabel('forward rate')
    plt.plot(df[name[i]][0:60].values,label = name[i]) #60年代表畫圖畫到60, 如果今天是收斂到40年可以改成40~
    plt.legend(bbox_to_anchor=(1, 0),loc = 3,borderaxespad=0) #看要不要加圖標
plt.savefig(x+'.png')
plt.show()