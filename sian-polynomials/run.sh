#!/bin/bash

# running on theta
directories=(new-akt new-cd8-cells new-goodwin-osc new_qwwc new_saird new_seiajrc new_seiddrf new_seiqrdc new_seir new_seir2 new_seirahd new_seirc new_seirp new_seisehr new_siar new_sidarthe new_sijru new_sir_denom new_siraqj new_sircd new_sliqr new_ssaair original_biohydrogenation original_chem_reac_network original_cholera original_daisy_ex3 original_daisy_mamil3 original_daisy_mamil4 original_hiv original_hiv2 original_hpv_group5 original_hpv_mf_group1_subs original_hpv_mf_group4 original_lipolisys original_lv original_nfkb original_oral_glucose original_pharm original_seir original_seir2 original_sir_r0 original_sirsforced original_slowfast original_treatment original_tumor)
for each in "${directories[@]}"; do
    mkdir -p $each/runtimes
    mkdir -p $each/outputs
    for d in $each/new_weights_with_trb_pos_char.mpl; do
        outname=(${d//"/"/ })
        outname="${outname[1]::-4}.out"
        timename="time_${outname::-4}.out"
        echo $d
        (/usr/bin/time -f "elapsed,user+system,memory\n%e,=%U+%S,%M000" /local/maple2022/bin/maple $d) >$each/outputs/$outname 2>$each/runtimes/$timename
        break # run only random weights
    done
    for d in $each/original_with_trb_pos_char.mpl; do
        outname=(${d//"/"/ })
        outname="${outname[1]::-4}.out"
        timename="time_${outname::-4}.out"
        echo $d
        (/usr/bin/time -f "elapsed,user+system,memory\n%e,=%U+%S,%M000" /local/maple2022/bin/maple $d) >$each/outputs/$outname 2>$each/runtimes/$timename
        break # run only random weights
    done
    # break
done
