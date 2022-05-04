import numpy as np
from scipy.spatial import cKDTree
from ROOT import TFile

from clustering import clustering

def IsEdepSufficient(pEdep, fEMinDamage, fEMaxDamage):
    if (pEdep<fEMinDamage):
        return False
    if(pEdep>fEMaxDamage):
        return True
    else:
        proba = (pEdep - fEMinDamage)/(fEMaxDamage-fEMinDamage)
        if test:
            return True
        else:
            return (proba>np.random.rand())


def checkPoint(point, T0, T1):
    result = [-1,-1]
    Rdirect = 0.35 #nm
    d0, idx0 = T0.query(point, k=1)
    d1, idx1 = T1.query(point, k=1)

    if d0 < d1:
        currentDist = d0
        currentStrand = 0
        currentCopy = idx0
    else:
        currentDist = d1
        currentStrand = 1
        currentCopy = idx1

    if currentDist < Rdirect:
        result[0]=currentStrand
        result[1]=currentCopy

    return result



filename = "test.csv"
outputFilename = "results.csv"
Rdirect = True
fEMinDamage =  5
fEMaxDamage = 37.5
test = False

with open("sugarPos.csv", "r") as f:
    data = f.readlines()

data = [[float(b)*10**6 for b in a.split("\t")] for a in data]

sugar0=[]
sugar1 = []

for line in data:
    sugar0.append(line[0:3])
    sugar1.append(line[3:6])

sugar0 = np.array(sugar0)
sugar1 = np.array(sugar1)

T0 = cKDTree(sugar0)
T1 = cKDTree(sugar1)

# Load data to analyse

with open(filename, "r") as f:
    input = f.readlines()

input = [a.split(",") for a in input]
dose=[]
sugarEdep=[]
clusteringResults=[]

for line in input:
    fFile = TFile.Open(line[0]) 
    fDirectory = fFile.Get("ntuple")


    tEdep = fDirectory.Get("EventEdep")

    EnergyDeposit=0
    TotalEnergyDeposit=0
    EventIDe=0
    chromatinVolume=5.37212e-22 # in m3
    totalEnergyDepSugar = 0


    entryEdepNumber=tEdep.GetEntries()

    tEdep.GetEntry(entryEdepNumber - 1)

    numEvt = tEdep.EventNo+1 # all events are recorded in EventEdep but only those with edep or reactions are included in direct and indirect trees

    for i in range(0,entryEdepNumber):
    
        tEdep.GetEntry(i)
        TotalEnergyDeposit += tEdep.Edep_J
    
    dose.append(TotalEnergyDeposit / (1000 * chromatinVolume))
    

    # // Read out direct damage
    input_tree = fDirectory.Get("Direct")

    nentries = input_tree.GetEntries()

    cumulatedEnergyDep = {}

    for irow in range(0,nentries):
    
        input_tree.GetEntry(irow)
        if (Rdirect==False):
            if ((input_tree.Strand == 1) or (input_tree.Strand == 2)):
                key = (input_tree.EventNo,input_tree.Strand - 1,input_tree.copyNum)
                if key in cumulatedEnergyDep:
                    cumulatedEnergyDep[key] += input_tree.eDep_eV
                else:
                    cumulatedEnergyDep[key] = input_tree.eDep_eV

                
        if (Rdirect==True):
            result = checkPoint([input_tree.x,input_tree.y,input_tree.z], T0,T1)

            if (result[0] != -1):
                key = (input_tree.EventNo,result[0],result[1])
                if key in cumulatedEnergyDep:
                    cumulatedEnergyDep[key] += input_tree.eDep_eV
                else:
                    cumulatedEnergyDep[key] = input_tree.eDep_eV


    eventsListDirect=[]
    copyListDirect=[]
    strandListDirect=[]


    for key in cumulatedEnergyDep:
        totalEnergyDepSugar+=cumulatedEnergyDep[key]
        if IsEdepSufficient(cumulatedEnergyDep[key],fEMinDamage, fEMaxDamage):
            eventsListDirect.append(key[0])
            copyListDirect.append(key[2])
            strandListDirect.append(key[1])

    # // Read out indirect damage
    input_tree_indirect = fDirectory.Get("Indirect")

    eventsListIndirect=[]
    copyListIndirect=[]
    strandListIndirect=[]

    nentries_indirect = input_tree_indirect.GetEntries()
    input_tree_indirect.GetEntry(nentries_indirect - 1)

    numEvt_indirect = input_tree_indirect.EventNo+1

    indirectDamage = np.zeros((numEvt_indirect,2,10584))

    for irow in range(0,nentries_indirect):
        input_tree_indirect.GetEntry(irow)
        if test:
            eventsListIndirect.append(input_tree_indirect.EventNo)
            copyListIndirect.append(input_tree_indirect.copyNum)
            strandListIndirect.append(input_tree_indirect.Strand-1)
        else:    
            if np.random.rand() <= 0.4:
                eventsListIndirect.append(input_tree_indirect.EventNo)
                copyListIndirect.append(input_tree_indirect.copyNum)
                strandListIndirect.append(input_tree_indirect.Strand-1)

    # clustering
    clusteringResults.append(clustering(numEvt,eventsListDirect,copyListDirect,strandListDirect,eventsListIndirect,copyListIndirect, strandListIndirect))

    sugarEdep.append(totalEnergyDepSugar/ ((numEvt) * 2 * 10584))
    print("Finished: {}".format(line[0]))

# write out results to txt file
with open(outputFilename, "w") as f:
    f.write("Filename,EnergyMeV,LET,DoseGy,MeanEdepPerSugar,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")

    for i, line in enumerate(clusteringResults):
        text = str(line).strip("[").strip("]")
        f.write("{},{},{},{},{},{}".format(input[i][0],input[i][1],input[i][1],dose[i],sugarEdep[i],text))
        f.write('\n')

        if test and Rdirect:
            assert text=="14, 10, 0, 2, 6, 4, 0, 1, 20, 10, 1, 4"
            assert dose[i] == 1606.114532065553
            assert sugarEdep[i] == 0.008914399092970522 
        
        if test and (not Rdirect):
            assert text=="13, 9, 0, 2, 6, 4, 0, 1, 19, 9, 1, 4"
            assert dose[i] == 1606.114532065553
            assert sugarEdep[i] == 0.008127047115142354 
