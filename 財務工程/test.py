import math
import numpy as np
import time
import pandas as pd
import os
import configparser
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

def matmult(a,b): #input兩個矩陣返回乘法結果
    zip_b = zip(*b)
    # uncomment next line if python 3 : 
    # zip_b = list(zip_b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
             for col_b in zip_b] for row_a in a]

def hh(z):
    return (z + math.exp(z*-1))/2

def Hmat(u, v):
    return hh(u+v) - hh(abs(u-v))

def Galfa(alfa, Q, mm, umax, nrofcoup, T2, Tau):
    h = []
    temp1 = []
    # for i in range(umax*nrofcoup):
    #     h.append([])
    #     for j in range(umax*nrofcoup):
    #         h[i].append([])
    # for i in range(mm):
    #     h.append([])
    #     for j in range(1):
    #         h[i].append([])

    for i in range(1, umax*nrofcoup+1):
        h.append([])
        for j in range(1, umax*nrofcoup+1):
            h[i-1].append(Hmat(alfa * i / nrofcoup, alfa * j / nrofcoup))

    for i in range(mm):
        temp1.append([])
        temp1[i].append(1 - sum(Q[i]))

    Q = np.matrix(Q)
    h = np.matrix(h)
    b = np.linalg.inv((Q*h)*Q.transpose())*temp1
    
    gamma = Q.transpose() * b
    temp2 = 0
    temp3 = 0
    for i in range(1, umax*nrofcoup+1):
        # print('temp2=='+str(temp2)+' gamma=='+str(gamma[i-1, 0]))
        temp2 = temp2 + gamma[i-1, 0] * i / nrofcoup
        temp3 = temp3 + gamma[i-1, 0] * np.sinh(alfa * i / nrofcoup)
    kappa = (1 + alfa * temp2) / temp3
    output = []
    output.append(alfa / abs(1 - kappa * math.exp(T2 * alfa)) - Tau)
    output.append(gamma)
    return output

def AlfaScan(lastalfa, stepsize, Q, mm, umax, nrofcoup, T2, Tau):
    '''
    這邊可能有錯 要檢查迴圈
    '''
    alfaReturn = 0
    for alfa in np.arange(lastalfa + stepsize / 10 - stepsize, lastalfa, stepsize / 10):
        alfa = float(round(alfa, 7))
        galfa_output = Galfa(alfa, Q, mm, umax, nrofcoup, T2, Tau)
        alfaReturn = alfa
        if galfa_output[0] <= 0:
            break
    return [alfaReturn, galfa_output[1]]

def SmithWilsonBruteForce(Instrument, DataIn, nrofcoup, CRA, UFRac, alfamin, Tau, T2):
    precision = 6
    Tau = Tau / 10000
    nrofrates = 0
    for i in range(len(DataIn)):
        nrofrates += DataIn[i][0]
    # print('nrofrates=='+str(nrofrates))

    u = []  #len(nrofrates)
    r = []  #len(nrofrates)
    j = 0
    for i in range(len(DataIn)):
        if DataIn[i][0] == 1:
            j = j + 1
            u.append(DataIn[i][1])
            r.append(DataIn[i][2] - CRA / 10000)
    umax = max(u)
    # print('umax==='+str(umax))
    lnUFR = math.log(1+UFRac)
    # print('lnUFR==='+str(lnUFR))

    if Instrument == 'Zero':
        nrofcoup = 1
    
    q = []
    for i in range(nrofrates):
        q.append([])
        for j in range(umax*nrofcoup):
            q[i].append(0)
    # print('q是'+str(nrofrates)+'X'+str(umax*nrofcoup)+'陣列')
    if Instrument == 'Zero':
        for i in range(nrofrates):
            for j in range(nrofrates):
                if i != j:
                    q[i][j] = 0
                else:
                    q[i][u[i]-1] = math.exp(lnUFR*-1 * u[i]) * (math.pow((1 + r[i]),u[i]))
    elif Instrument == 'Swap' or Instrument == 'Bond':
        for i in range(nrofrates):
            for j in range(u[i]*nrofcoup-1):
                q[i][j] = math.exp(lnUFR*-1 * j / nrofcoup) * r[i] / nrofcoup
            q[i][j] = math.exp(lnUFR*-1 * j / nrofcoup) * (1 + r[i] / nrofcoup)
    galfa_output = Galfa(alfamin, q, nrofrates, umax, nrofcoup, T2, Tau)
    # print(galfa_output[0])
    # time.sleep(50)
    if galfa_output[0] <= 0:
        alfaReturn = alfamin
        gamma = galfa_output[1]
    else:
        stepsize = 0.1
        alfaReturn = 0
        for alfa in np.arange(alfamin+stepsize, 20, stepsize):
            alfa = float(round(alfa, 4))
            alfaReturn = alfa
            if Galfa(alfa, q, nrofrates, umax, nrofcoup, T2, Tau)[0] <= 0:
                break
        for i in range(1, precision):
            alfascanoutput = AlfaScan(alfaReturn, stepsize, q, nrofrates, umax, nrofcoup, T2, Tau)
            alfaReturn = alfascanoutput[0]
            stepsize = stepsize / 10
        gamma = alfascanoutput[1]

    h = []
    for i in range(121+1):
        h.append([])
        for j in range(umax*nrofcoup):
            h[i].append([])

    g = []
    for i in range(121+1):
        g.append([])
        for j in range(umax*nrofcoup):
            g[i].append([])

    for i in range(121+1):
        for j in range(1, umax*nrofcoup+1):
            h[i][j-1] = Hmat(alfaReturn * i, alfaReturn * j / nrofcoup)
            if j / nrofcoup > i:
                g[i][j-1] = (alfaReturn * (1 - math.exp(alfaReturn*-1 * j / nrofcoup) * np.cosh(alfaReturn * i)))
            else:
                g[i][j-1] = (alfaReturn * math.exp(alfaReturn*-1 * i) * np.sinh(alfaReturn * j / nrofcoup))
 
    tempDiscount = []
    tempintensity = []
    discount = []
    fwintensity = []
    yldintensity = []
    forwardac = []
    zeroac = []
    forwardcc = []
    zerocc = []

    temptempdiscount = (np.matrix(h)*gamma).transpose()
    temptempintensity = (np.matrix(g)*gamma).transpose()
    for i in range(121+1):
        tempDiscount.append(temptempdiscount[0,i])
        tempintensity.append(temptempintensity[0,i])
    temp = 0
    for i in range(1, umax*nrofcoup+1):
        temp = temp + (1 - math.exp(alfaReturn * -1 * i / nrofcoup)) * gamma[i-1, 0]

    yldintensity.append(lnUFR - alfaReturn * temp)
    fwintensity.append(yldintensity[0])
    discount.append(1)
    yldintensity.append(lnUFR - math.log(1 + tempDiscount[1]))
    fwintensity.append(lnUFR - tempintensity[1] / (1 + tempDiscount[1]))
    discount.append(math.exp(lnUFR*-1) * (1 + tempDiscount[1]))
    zeroac.append(0)
    zeroac.append(1 / discount[1] - 1)
    forwardac.append(0)
    forwardac.append(zeroac[1])
    for i in range(2, 121):
        yldintensity.append(lnUFR - math.log(1 + tempDiscount[i]) / i)
        fwintensity.append(lnUFR - tempintensity[i] / (1 + tempDiscount[i]))
        discount.append(math.exp(lnUFR * -1 * i) * (1 + tempDiscount[i]))
        zeroac.append(math.pow(1 / discount[i], (1 / i)) - 1)            
        forwardac.append(discount[i - 1] / discount[i] - 1)
  
    yldintensity.append(0)
    fwintensity.append(0)
    zeroac.append(0)
    forwardac.append(0)
    discount.append(alfaReturn)

    forwardcc.append(0)
    zerocc.append(0)
    for i in range(1, 121):
        forwardcc.append(math.log(1 + forwardac[i]))
        zerocc.append(math.log(1 + zeroac[i]))
    forwardcc.append(0)
    zerocc.append(0)
    # print(len(discount))
    # print(len(yldintensity))
    # print(len(zeroac))
    # print(len(fwintensity))
    # print(len(forwardcc))
    # print(len(forwardac))
    # print(forwardcc)
    #zeroac forwardcc forwardac(少了第一項 index=1)有誤
    # print(yldintensity)
    return [discount, yldintensity, zeroac, fwintensity, forwardcc, forwardac]

def genDataIn(dataFileName, dataDate, needOneList):
    dataIn = []
    dataFile = open(dataFileName, 'r')
    dataLines = dataFile.readlines()
    for data in dataLines:
        if dataDate == data.split(',')[0]:
            for i in range(10):
                if i+1 in needOneList:
                    dataIn.append([1, i+1, float(data.split(',')[i+1].replace('\n', ''))])
                else:
                    dataIn.append([0, i+1, float(data.split(',')[i+1].replace('\n', ''))])
    return dataIn

def getDateList(dataFileName):
    dateList = []
    dataFile = open(dataFileName, 'r')
    dataLines = dataFile.readlines()
    for data in dataLines:
        dateList.append(data.split(',')[0])
    return dateList


#main
config = configparser.ConfigParser()
config.read('config.ini', encoding='utf-8')
needOneBigList = []
dirName = 'output'
dataFileName = config.get('Parameter', 'DataFileName')
Instrument = config.get('Parameter', 'Instrument')
nrofcoup = float(config.get('Parameter', 'Nrofcoup'))
CRA = float(config.get('Parameter', 'CRA'))
UFRac = float(config.get('Parameter', 'UFRac'))
alfamin = float(config.get('Parameter', 'Alfamin'))
Tau = float(config.get('Parameter', 'Tau'))
T2 = float(config.get('Parameter', 'T2'))
rateType = config.get('Parameter', 'RateType')

for needOne in config.get('Parameter', 'LiquidMaturity').split('|'):
    needOneBigList.append([float(i) for i in needOne.split(',')])

dateList = getDateList(dataFileName)

if os.path.isdir(dirName) == False:
    os.mkdir(dirName)

for k in range(len(needOneBigList)):
    row = 0
    col = 0
    workbook = xlsxwriter.Workbook(dirName+'/output'+str(k)+'_'+rateType+'.xlsx')
    worksheet = workbook.add_worksheet()
    # for i in range(len(dateList)):
    #     if i == 0:
    #         worksheet.write(row, col, 'Year')
    #     else:
    #         worksheet.write(row, col, str(i))
    #     row = row + 1
    worksheet.write(0, col, 'Year')
    for date in dateList:
        row = 0
        col = col + 1
        print('正在運行 日期=='+date)
        worksheet.write_string(row, col, str(date))
        DataIn = genDataIn(dataFileName, date, needOneBigList[k])
        if rateType == 'forward':
            rateList = SmithWilsonBruteForce(Instrument, DataIn, nrofcoup, CRA, UFRac, alfamin, Tau, T2)[5]
        else:
            rateList = SmithWilsonBruteForce(Instrument, DataIn, nrofcoup, CRA, UFRac, alfamin, Tau, T2)[2]
        row = row + 1
        worksheet.write_string(row, col, str(date))
        for i in range(len(rateList)-1):
            # print(rateList[i])
            worksheet.write_string(row, 0, str(i))
            if i == 0:
                continue
            else:
                worksheet.write_string(row, col, str(rateList[i]))
            row = row + 1 
    workbook.close()

    
    # df = pd.read_csv(dirName+'/output'+str(k)+'_'+rateType+'.csv', encoding='utf-8')
    # df2 = df.T
    # df2.to_csv(dirName+'/output'+str(k)+'_'+rateType+'.csv', header=False)

    df = pd.read_excel(dirName+'/output'+str(k)+'_'+rateType+'.xlsx')
    df.set_index('Year',inplace=True)

    name = df.columns

    x = rateType+'_'+str(UFRac)+'_Interpolation_'+str(T2) #儲存圖片的檔名＆圖片的大標名稱
    for i in range(0,len(name),6):  #間隔6個月畫一條線
        plt.style.use('ggplot')
        plt.figure(num = 3,figsize= (20,13))
        plt.title(x)

        plt.xlabel('year',fontsize=30)
        plt.ylabel('forward rate', fontsize=30)
        plt.plot(df[name[i]][0:60].values,label = name[i]) #60年代表畫圖畫到60, 如果今天是收斂到40年可以改成40~
 #       plt.legend(bbox_to_anchor=(1, 0),loc = 3,borderaxespad=0) #看要不要加圖標

    plt.savefig(dirName+'/output'+str(k)+'_'+x+'.png')
