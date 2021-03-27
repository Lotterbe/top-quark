import ROOT
from collections import OrderedDict

import  sys

#from Plotting.Database import config, scaleDB, getScaleFactor
from Plotting import Database
from Plotting import infofile

if len(sys.argv) != 4:
	print("USAGE: %s <input  file > <input  file > <output  file >"%(sys.argv [0]))
	sys.exit (1)

inFileNameone   = sys.argv [1]
inFileNametwo   = sys.argv [2]
outFileName = sys.argv [3]
print("Reading  from", inFileNameone , "and from", inFileNametwo, "and  writing  to", outFileName)

inFileone = ROOT.TFile.Open(inFileNameone ,"READ")
histDataone = inFileone.Get("top_mass")
histDataone.SetDirectory(0)
print(type(histDataone))
#histDataone = Dataone.
inFileone.Close()

inFiletwo = ROOT.TFile.Open(inFileNametwo ,"READ")
histDatatwo = inFiletwo.Get("top_mass")
histDatatwo.SetDirectory(0)
print(type(histDatatwo))
inFiletwo.Close()

check = histDataone.Add(histDatatwo)

print(check)

file = open("backgroundfilelist.txt","r")
filelist = file.readlines()
histFile = ROOT.TFile.Open(str(filelist[0][:-1]), "READ")
background = histFile.Get("top_mass")
background.SetDirectory(0)
histFile.Close()

for filename in filelist[1:]:
	histFile = ROOT.TFile.Open(filename[:-1], "READ")
	backgroundpart = histFile.Get("top_mass")
	backgroundpart.SetDirectory(0)
	pathname = "/home/student/atlas-outreach-PyROOT-framework-13tev/results/"
	dataname = filename[:-1][len(pathname):-5]
	print("dataname", dataname)
	Database.config["Luminosity"] = 1000
	scaling = Database.getScaleFactor(dataname)
	backgroundpart.Scale(scaling)
	background.Add(backgroundpart)
	histFile.Close()
	#with ROOT.TFile.Open(filename,"READ") as backgrounddata:
	#	print(tzpe(backgrounddata))
    #print(line.rstrip())
file.close()

checktwo = histDataone.Add(background, -1)

print(checktwo)

outHistFile = ROOT.TFile.Open(outFileName ,"RECREATE")

outHistFile.cd()
if check == True:
	histDataone.Write()
outHistFile.Close()
