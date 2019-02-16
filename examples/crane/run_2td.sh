#!/bin/bash
glob_stab_bin_dir="../../bin/"
ch_eq=crane_2td.chEq
glob_stab_opts="-jf --tm_order=5"
postfix=""

for p in crane_2td*.p; do
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
