for each in */original*.mpl; do
    dir_name=$(dirname $each)
    file_name=$(basename $each)
    cp $each $dir_name/new_weights_with_trb_pos_char.mpl
done
