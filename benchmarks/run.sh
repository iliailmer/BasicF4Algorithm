#!/bin/bash

# running on theta
directories=(alea6 bayes148 cyclic7 cyclic8 cyclic9 eco12 gametwo7 jason katsura7 katsura9 katsura12 mayr42 noon9 reimer7 reimer8 shwarz yang1)
for each in "${directories[@]}"; do
    mkdir -p $each/runtimes
    mkdir -p $each/outputs
    for d in $each/$each.mpl; do
        outname=(${d//"/"/ })
        outname="${outname[1]::-4}.out"
        timename="time_${outname::-4}.out"
        echo $d
        (/usr/bin/time -f "elapsed,user+system,memory\n%e,=%U+%S,%M000" /local/maple2022/bin/maple $d) >$each/outputs/$outname 2>$each/runtimes/$timename
        break # run only random weights
    done
    # for d in $each/original_with_trb_pos_char.mpl; do
    #     outname=(${d//"/"/ })
    #     outname="${outname[1]::-4}.out"
    #     timename="time_${outname::-4}.out"
    #     echo $d
    #     (/usr/bin/time -f "elapsed,user+system,memory\n%e,=%U+%S,%M000" /local/maple2022/bin/maple $d) >$each/outputs/$outname 2>$each/runtimes/$timename
    #     break # run only random weights
    # done
    # break
done
