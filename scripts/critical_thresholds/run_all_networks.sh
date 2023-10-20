g++ -O3 -I src_eigenvalues/ mpa_perco_threshold.cpp -o mpa_perco_threshold

IFS=$'\n'
FILES="$(find ../../data/networks/ -type d)"
for network in $FILES
do 
    NAME=$(basename "$network")
    echo $NAME
    if [ "$NAME" = "networks" ]; then 
        continue 
    fi
    if [ "$NAME" = "threshold.dat" ]; then 
        continue 
    fi
    #echo "../../data/edgelists/$NAME/.txt"
    #echo "../../data/empirical_networks/$NAME/threshold.dat"
    ./mpa_perco_threshold "../../data/networks/$NAME/edgelist.txt" > "../../data/networks/$NAME/threshold.dat"
done