counter=0
for pos1 in {G,A}; do
	for pos2 in {G,A}; do
		for pos3 in {G,A}; do
			for pos4 in {G,A}; do
				for pos5 in {G,A}; do
					echo ">rbs_${counter}"
					echo "${pos1}${pos2}${pos3}${pos4}${pos5}"
					counter=$(( counter+1 ))
				done
			done
		done
	done
done