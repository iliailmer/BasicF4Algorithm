for x in ./*.mpl; do
    mkdir "${x%.*}" && mv "$x" "${x%.*}"
done
