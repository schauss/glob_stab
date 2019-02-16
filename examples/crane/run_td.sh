#!/bin/bash
glob_stab_bin_dir="../../bin/"
ch_eq=crane_td.chEq
glob_stab_opts="-jf"
postfix=""

for p in crane_td1_2d2.p  crane_td1_4d.p  crane_td2_2d2.p  crane_td2_4d.p  crane_td3_2d2.p  crane_td3_4d.p  crane_td4_2d2.p  crane_td4_4d.p; do
	echo Processing ${p}...;

	# extract name (without path or extension)
	name=$(basename $p .p)${postfix};

	# call glob stab
	if [ ! -f ${name}.asy -o $p -nt ${name}.asy -o $ch_eq -nt ${name}.asy ]; then
		${glob_stab_bin_dir}/glob_stab -c $ch_eq -p $p $glob_stab_opts -ao ${name}.asy > ${name}.log;
		echo "${name}.asy created!"
	else
		echo "${name}.asy exists => skipping!"
	fi

	# create pdf
	if [ ! -f ${name}.pdf -o ${name}.asy -nt ${name}.pdf ]; then
		asy -f pdf ${name}.asy;
		echo "${name}.pdf created!"
	else
		echo "${name}.pdf exists => skipping!"
	fi
done
