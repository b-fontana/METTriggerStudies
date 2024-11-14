#! /usr/bin/env python
# Example of how to read the a correctionlib JSON file

# This is a minimal version of a tau POG example:
# For more information, see the README in
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/TAU
#import sys; sys.path.insert(0,"correctionlib") # add correctionlib to path
from correctionlib import _core

# This corrections file is documented in https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/
# qhere the general correctionlib documentation is kept
fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22/btagging.json.gz"

# Load CorrectionSet
if fname.endswith(".json.gz"):
  import gzip
  with gzip.open(fname,'rt') as file:
    #data = json.load(file)
    data = file.read().strip()
  cset = _core.CorrectionSet.from_string(data)
else:
  cset = _core.CorrectionSet.from_file(fname)

# Load Correction objects that can be evaluated
corr = cset["deepJet_wp_values"]

btagWPl = corr.evaluate("L")
btagWPm = corr.evaluate("M")
btagWPt = corr.evaluate("T")
btagWPxt = corr.evaluate("XT")
btagWPxxt = corr.evaluate("XXT")

print("L WP =", btagWPl)
print("M WP =", btagWPm)
print("T WP =", btagWPt)
print("XT WP =", btagWPxt)
print("XXT WP =", btagWPxxt)
