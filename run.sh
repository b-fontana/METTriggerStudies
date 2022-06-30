
TAG="MultipleDatasets"
BASE="jobs/${TAG}/submission/"

# condor_submit "${BASE}"Histos_SKIM_MET/Histos_SKIM_MET.condor
# condor_submit "${BASE}"Histos_SKIM_EGamma/Histos_SKIM_MET.condor
# condor_submit "${BASE}"Histos_SKIM_TT_fullyLep/Histos_SKIM_TT_fullyLep.condor
# condor_submit "${BASE}"Histos_SKIM_TT_semiLep/Histos_SKIM_TT_semiLep.condor
# condor_submit "${BASE}"Histos_SKIM_TT_fullyHad/Histos_SKIM_TT_fullyHad.condor

# condor_submit "${BASE}"Counts_SKIM_MET/Counts_SKIM_MET.condor
# condor_submit "${BASE}"Counts_SKIM_EGamma/Counts_SKIM_EGamma.condor
# condor_submit "${BASE}"Counts_SKIM_TT_fullyLep/Counts_SKIM_TT_fullyLep.condor
# condor_submit "${BASE}"Counts_SKIM_TT_semiLep/Counts_SKIM_TT_semiLep.condor
# condor_submit "${BASE}"Counts_SKIM_TT_fullyHad/Counts_SKIM_TT_fullyHad.condor

#condor_submit "${BASE}"HaddHistoMET/HaddHistoMET.condor
#condor_submit "${BASE}"HaddHistoTT/HaddHistoTT.condor

# condor_submit "${BASE}"HaddCountsMET/HaddCountsMET.condor
# condor_submit "${BASE}"HaddCountsTT/HaddCountsTT.condor

# condor_submit "${BASE}"HaddHistoAggMET/HaddHistoAggMET.condor
# condor_submit "${BASE}"HaddHistoAggTT/HaddHistoAggTT.condor

# condor_submit "${BASE}"HaddCountsAggMET/HaddCountsAggMET.condor
# condor_submit "${BASE}"HaddCountsAggTT/HaddCountsAggTT.condor

# condor_submit "${BASE}"EffAndScaleFactors/EffAndScaleFactors.condor

# condor_submit "${BASE}"EffAndSFAgg/EffAndSFAgg.condor

# for sub in "${BASE}"Discriminator_*; do
# 	condor_submit "${sub}/${sub}.condor"
# done

# for sub in "${BASE}"UnionWeightsCalculator_*;do
# 	condor_submit "${sub}.condor"
# done

# condor_submit "${BASE}"Closure/Closure.condor
