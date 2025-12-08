#!/bin/bash

arr=("$@")
second_last=("${arr[@]:3}")
process=""
for i in "${second_last[@]}"; do
    case "$i" in
        21) process+="g " ;;
        1)  process+="d " ;;
        2)  process+="u " ;;
	3)  process+="s " ;;
	4)  process+="c " ;;
	5)  process+="b " ;;
	6)  process+="t " ;;
        -1)  process+="d~ " ;;
        -2)  process+="u~ " ;;
        -3)  process+="s~ " ;;
        -4)  process+="c~ " ;;
        -5)  process+="b~ " ;;
        -6)  process+="t~ " ;;
	22) process+="a " ;;
	23) process+="z " ;;
	24) process+="w+ " ;;
	-24) process+="w- " ;;
	25) process+="h " ;;
	11) process+="e- " ;;
	12) process+="ve " ;;
	13) process+="mu- " ;;
	14) process+="vm " ;;
	15) process+="ta- " ;;
	16) process+="vt " ;;
	-11) process+="e+ " ;;
        -12) process+="ve~ " ;;
        -13) process+="mu+ " ;;
        -14) process+="vm~ " ;;
        -15) process+="ta+ " ;;
        -16) process+="vt~ " ;;
    esac
done

arr=($process)
arr=("${arr[@]:0:2}" ">" "${arr[@]:2}")

next=${#arr[@]}
next=$((next - 1)) 

cat > "script.txt" << EOF
generate ${arr[@]}
output standalone temp_dir
launch temp_dir
EOF

chmod +x "script.txt"

> output.txt
./mg5_aMC script.txt > output.txt

filename="momenta_$2_$3.txt"
> $filename

# Print the momenta
awk -v out=$filename -v num_lines=$next '
/px/ {
for (i=0; i<num_lines; i++) {
        if (getline line > 0) {
            n = split(line, words, /[[:space:]]+/)

            # Print words 2 to last
            extracted = ""
	    for (j=3; j<=6 && j<=n; j++) {
                extracted = extracted words[j] " "
            }
            gsub(/[[:space:]]+$/, "", extracted)
            if (extracted != "")
                print extracted >> out
        }
    }
}' "output.txt"

# Print the matrix element
awk '/ANSWER for/ { count++ } END { print count }' output.txt >> $filename
awk '/ANSWER for/ { print $5 }' output.txt  >> $filename

echo "...done for this mg folder"

mv $filename $1/ME_checks






