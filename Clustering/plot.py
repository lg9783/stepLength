import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plotdata(filename, ax, tag):
    data_list=[] 
    # # Read in the csv to a pandas dataframe
    for i in range(1,6):
        temp=pd.read_csv('{}_{}.csv'.format(filename,i), index_col=None, header=0)
        data_list.append(temp)

    data = pd.concat(data_list, axis=0, ignore_index=True) # ignore event No column


    numGBP = 10584e-9

    BPdensity = 10584/5.37E+05 #BP/volume in nm3
    refBPdensity = 0.012 #Cancers 2021, 13, 4940
    densityFactor = refBPdensity/BPdensity

    # data.DoseGy = 200*data.LET*10*1.60218e-19/(1000*5.37212e-22)

    # Normailise to per Gy per GBP
    data.TotalSBdirect = densityFactor*data.TotalSBdirect/data.DoseGy/numGBP
    data.TotalSBindirect = densityFactor*data.TotalSBindirect/data.DoseGy/numGBP
    data.TotalSBtotal = densityFactor*data.TotalSBtotal/data.DoseGy/numGBP

    data.DSBdirect = densityFactor*data.DSBdirect/data.DoseGy/numGBP
    data.DSBindirect = densityFactor*data.DSBindirect/data.DoseGy/numGBP
    data.DSBtotal = densityFactor*data.DSBtotal/data.DoseGy/numGBP

    # # LET vs number of strand breaks
    totalSBtotalMean=[]
    totalSBdirectMean=[]
    totalSBindirectMean=[]
    totalSBtotalSTD=[]
    totalSBdirectSTD=[]
    totalSBindirectSTD=[]
    LET = []
    Energy = list(set(data.EnergyMeV))

    for i in Energy:
        totalSBtotalMean.append(data.TotalSBtotal[data.EnergyMeV==i].mean())
        totalSBdirectMean.append(data.TotalSBdirect[data.EnergyMeV==i].mean())
        totalSBindirectMean.append(data.TotalSBindirect[data.EnergyMeV==i].mean())

        totalSBtotalSTD.append(data.TotalSBtotal[data.EnergyMeV==i].std())
        totalSBdirectSTD.append(data.TotalSBdirect[data.EnergyMeV==i].std())
        totalSBindirectSTD.append(data.TotalSBindirect[data.EnergyMeV==i].std())
        LET.append(data.LET[data.EnergyMeV==i].mean())

    # LET = [102,87,76]

    ax[0].errorbar(LET, totalSBtotalMean,yerr=totalSBtotalSTD, fmt='.', label='Total SB - {}'.format(tag))
    ax[0].errorbar(LET, totalSBindirectMean,yerr=totalSBindirectSTD, fmt='.', label='Indirect SB - {}'.format(tag))
    ax[0].errorbar(LET, totalSBdirectMean,yerr=totalSBdirectSTD, fmt='.', label='Direct SB - {}'.format(tag))

    ax[0].set_xlim([0,110])

    ax[0].legend()

    ax[0].set_xlabel('LET keV/$\mu$m')
    ax[0].set_ylabel('Number of strand breaks ($Gy^{-1} Gbp^{-1}$)')


    # # LET vs number of strand breaks
    DSBtotalMean=[]
    DSBtotalSTD=[]


    for i in Energy:

        DSBtotalMean.append(data.DSBtotal[data.EnergyMeV==i].mean())

        DSBtotalSTD.append(data.DSBtotal[data.EnergyMeV==i].std())

    ax[1].errorbar(LET, DSBtotalMean,yerr=DSBtotalSTD, fmt='.', label='DSB - {}'.format(tag))

    ax[1].set_xlim([0,110])

    ax[1].legend()

    ax[1].set_xlabel('LET keV/$\mu$m')
    ax[1].set_ylabel('Number of DSB ($Gy^{-1} Gbp^{-1}$)')

fig, ax = plt.subplots(1,2)

filename = "resultsAlpha"
plotdata(filename, ax, "alpha")

filename = "resultsProton"
plotdata(filename, ax, "proton")

plt.show()