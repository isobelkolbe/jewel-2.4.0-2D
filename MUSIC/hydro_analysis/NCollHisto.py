# Short python script to bin the data in NCollList.dat
# Outputs a file in the working directory called "NCollHistogram.dat"
# Written by Isobel Kolbe, Feb. 2023    

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib
import seaborn as sns



print("Will now bin NcollList.dat")
# To estimate edges for given A
# nuclA=208                                   
# R=1.12*nuclA**(1/3)-0.86*nuclA**(-1/3)

# Or always do max edges (-10 fm, 10 fm)
R=10

# SigmaNN (get from MUSIC file: parameters_dict_user_IPGlasma_wevo.py)
sigmaNNmb=72  # in mb
sigmaNN=sigmaNNmb*0.1  # in fm^2
nuclA1=208
nuclA2=208
AB=nuclA1*nuclA2

myBinWidth=0.2  # Standard in JEWEL: 0.2 fm                     
nBins=round(2*R/myBinWidth)+1
print("There are ", nBins," bins with width ", myBinWidth,".")
myRange=[[-R,R],[-R,R]]

# Import data into a pandas frame
colNames=['x','y']
df = pd.read_csv("NcollList.dat",delim_whitespace=True,skiprows=[0], names=colNames)

# Create histogram object to write out later
H ,xedges, yedges= np.histogram2d(
    df["x"].values.tolist(),df["y"].values.tolist(),
    bins=nBins,
    range=myRange)


#--------------------------------------------------
#   Uncomment this block to see a plot
#--------------------------------------------------
# myCmap = matplotlib.cm.get_cmap('viridis')
# myColor=myCmap(0.2)
# myPlot = sns.jointplot(data=df,
#     x='x',y='y',
#     # kind='hex',
#     kind='hist',
#     binwidth=myBinWidth,
#     cmap=myCmap,
#     marginal_kws=dict(color=myColor,
#     binwidth=myBinWidth))
# plt.show()

#--------------------------------------------------



# Output histogram data and integrate Ncoll
outputFile= open("NCollHistogram.dat","w")

outputFile.write("Collision system: " + str(nuclA1) + " X " + str(nuclA2) + "\n")
outputFile.write("Sigma_NN used: " + str(sigmaNN) +  "fm^2\n")
# outputFile.write("TAB normalization constant: " + str(normC) + "\n")
outputFile.write("----------------------------------- \n")
outputFile.write("xmid \t ymid \t Ncoll \n")

testInt=0
for i in range((len(xedges)-1)):
    xval=abs(xedges[i+1]-xedges[i])/2+xedges[i]
    for j in range((len(yedges)-1)):
        yval=abs(yedges[j+1]-yedges[j])/2+yedges[j]
        
        # outputFile.write("{:8.3f}".format(xval))
        # outputFile.write("{:8.3f}".format(yval))
        # outputFile.write("{:8.3f}".format(H[i,j]))
        # outputFile.write("{:8.3f}".format(H[i,j]/sigmaNN))
        # outputFile.write("\n")
        # testInt=testInt+normC*H[i,j]/sigmaNN*abs(xedges[i+1]-xedges[i])*abs(yedges[j+1]-yedges[j])
        
        outputFile.write(f"{xval:<10.3f}{yval:<10.3f}{H[i,j]:<10.3f}")

outputFile.close()

# normC=AB/TABint
# print("Normalization constant: ", normC)
# print("integrated TAB: ", testInt, ".  Compare with A1 * A2: ",AB)
