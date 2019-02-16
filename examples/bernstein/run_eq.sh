#!/bin/bash
glob_stab_bin_dir="../../bin/"
glob_stab_opts="-s0"
postfix="_inequality"

for p in *.p; do
	echo Processing ${p}...;

	# extract name (without path or extension)
	name=$(basename $p .p)${postfix};
	ch_eq=$(basename $p .p).eq

	if [ ! -f $ch_eq ]; then
		echo "${ch_eq} doesn't exist => skipping!"
		continue
	fi

	# call glob stab
	if [ ! -f ${name}.asy -o $p -nt ${name}.asy -o $ch_eq -nt ${name}.asy ]; then
		${glob_stab_bin_dir}/glob_stab -c $ch_eq -p $p $glob_stab_opts -ao ${name}.asy > ${name}.log;
	fi

	# create pdf
	if [ ! -f ${name}.pdf -o ${name}.asy -nt ${name}.pdf ]; then
		asy -f pdf ${name}.asy;
		echo "${name}.pdf created!"
	else
		echo "${name}.pdf exists => skipping!"
	fi
done
