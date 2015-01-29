#!/bin/bash
# arrays on parameters
declare -a Wt=(0.6 0.8 1.01 1.20 1.40 1.608 1.8 2 2.2 2.4 2.6 2.8 3)
declare -a tSM=(0 1 2 3 4 5 6)
###############################################
# turn off the automatic htmal open
sed -i -e 's/# automatic_html_opening = True/automatic_html_opening = False/g' Cards/me5_configuration.txt
#change the number of events to 3
#sed -i -E -e 's/[0-9]+ = nevents/20000 = nevents/g' Cards/run_card.dat
#change the minimum distance between jets
#sed -i -E -e 's/ [0-9]+?\.?[0-9]+?.* \= drjj    \! min distance between jets/0.0001 \= drjj    \! min distance between jets/g' Cards/run_card.dat
#sed -i -E -e 's/ [0-9]+?\.?[0-9]+?.* \= draa    \! min distance between gammas/ 0.0001 \= draa    \! min distance between gammas/g' Cards/run_card.dat
#sed -i -E -e 's/ [0-9]+?\.?[0-9]+?.* \= draj    \! min distance between gamma and jet/ 0.0001 \= draj    \! min distance between gamma and jet/g' Cards/run_card.dat
#change the CM energy to 2Tev
#sed -i -E -e 's/[0-9]+.*=.*ebeam1/4000     = ebeam1/g' Cards/run_card.dat
#sed -i -E -e 's/[0-9]+.*=.*ebeam2/4000     = ebeam2/g' Cards/run_card.dat
# change the Higgs mass to 125 Gev
#sed -i -E -e 's/25 [0-9]+?\.?[0-9]+?.* .*MH/25 1.25000e+02 \# MH/g' Cards/param_card.dat
# change LR
#sed -i -E -e 's/25 [0-9]+?\.?[0-9]+?.* .*MH/25 1.25000e+02 \# MH/g' Cards/param_card.dat
################################################
mkdir top_Wvary
  #vary on the masses
  for (( i = 0 ; i <${#tSM[@]} ; i++ )); do
# for (( i = 6 ; i < 7 ; i++ )); do
    echo "${Wt[$i]}" " " "${tSM[$i]}"
    #change the parameters
    # sed -i -E -e "s/ 35 [0-9]+?\.?[0-9]+?.* \# MH02/ 35 ${MR[$i]} \# MH02/g" Cards/param_card.dat
    sed -i -E -e "s/.* \# WT/DECAY   6 ${Wt[$i]} \# WT/g" Cards/param_card.dat
    #generate the events
    ./bin/generate_events 0 Wt_${tSM[$i]} run1
    gunzip Events/Wt_${tSM[$i]}/unweighted_events.lhe.gz
    mv Events/Wt_${tSM[$i]}/unweighted_events.lhe top_Wvary/Wt_${tSM[$i]}.lhe
    #take the CX
    echo "Wt" "${tSM[$i]}" "$(grep -E '\#  Integrated weight \(pb\)  \:  ' top_Wvary/Wt_${tSM[$i]}.lhe)" >> top_Wvary/CX_LHC13.txt
  done
