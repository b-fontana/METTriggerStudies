import os
import sys

import ROOT

# Enable multi-threading                                                                                                                                                                                    
# ROOT.ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

core_dir = str(os.getcwd()).split('producer')

sys.path.insert(1, core_dir[0]+"/python")

import argparse
parser = argparse.ArgumentParser(description='Ntuplizer options')
parser.add_argument('-v','--tauid_version', choices=['2p1', '2p5'], dest="tauid_version", default='2p5')
parser.add_argument('-o','--outfile', dest='outfile')
parser.add_argument('-i','--inputFiles', dest='inputFiles', default = False)
parser.add_argument('-mc', '--isMC', action='store_true', default = False)
options = parser.parse_args()

# select the version
TauID_ver = options.tauid_version
outFile   = options.outfile
useFiles  = options.inputFiles 
isMC  = options.isMC

def create_rdataframe(folders, inputFiles=None):
    if not inputFiles:
        inputFiles = []
        for folder in folders:
            files = os.listdir(folder)
            inputFiles += [folder + f for f in files]

    return ROOT.RDataFrame("Events", tuple(inputFiles))

def obtain_picontuple(df):

    branches = []

    # Tag And Probe selection {Obtaining high pure Z -> mu tau Events}
    ## select muon (Tag) candidate

    df = df.Filter("HLT_IsoMu24 > 0")

    df = df.Define("pass_quadJet", "return HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65 > 0;")
    df = df.Define("pass_deepTau", "return HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1 > 0;")
    df = df.Define("pass_both", "return HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65 > 0 && HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1 > 0;" )
    branches = ["dau1_pt", "dau1_eta", "pass_quadJet", "pass_deepTau", "pass_both"]
    branch_list = ROOT.vector('string')()
    for branch_name in branches:
        branch_list.push_back(branch_name)
    df.Snapshot("Events", "outfile.root", branch_list)



if __name__ == '__main__':

    useFiles = str(useFiles)

    if ".txt" in useFiles:
      print("Using files in {}".format(useFiles))
      folders = []
      inputFiles_run3 = []

      with open(useFiles) as f:
        for line in f:
          line = line.replace('\n',"") # trim newline character
          if "root://" not in line and not line.startswith("/eos") and not line.startswith("/pnfs"):
            inputFiles_run3.append('root://cms-xrd-global.cern.ch//'+line) # prepend with redirector
          else:
            inputFiles_run3.append(line)
      print(inputFiles_run3)
      df = ROOT.RDataFrame("Events", tuple(inputFiles_run3))

    # else:
    #   print("Not a valid inputFiles argument")
    #   print("Use 8102, 8136, or a text file of nanoaodfile locations")

    obtain_picontuple(df)
   