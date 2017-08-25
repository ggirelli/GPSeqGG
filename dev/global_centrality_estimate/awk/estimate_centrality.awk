
# Estimate centrality, requires two input variables:
#
# Args:
#	calc (string): 'ratio', 'diff' or 'logratio'
#	type (string): '2p' (two point), 'f' (fixed) or 'g' (global)
#	cumrat (bool): pre-processing for Cumulative of Ratio metric
#	ratcum (bool): pre-processing for Ratio of Cumulative metric

BEGIN {
	OFS = FS = "\t";

}

function estimate(calc, a, b) {
	switch (calc) {
		case "ratio":
			return a / b;
			break;
		case "diff":
			return a - b;
			break;
		case "logratio":
			return log(a / b);
			break;
	}
}

{
	# Cumulative ratio
	if ( cumrat == 1) {
		# Sum probabilities
		for ( i = 2; i <= NF; i++ ) {
			$i = $i + $(i-1);
		}
	}

	# Ratio of cumulatives
	if ( ratcum == 1 ) {
		# Build table
		for ( i = 1; i <= NF; i++ ) {
			nf=split($i, ff, ",");
			for ( j = 1; j <= nf; j++ ) {
				a[i, j] = ff[j];
			}
		}
		# Calculate ratio of cumulatives
		for ( i = 1; i <= NF; i++ ) {
			a[i, 2] += a[i - 1, 2];
			a[i, 1] += a[i - 1, 1];
			$i = a[i, 2] / (a[i, 1] * a[i, 3]);
		}
	}

	# Calculate depending of 
	switch (type) {
		case "2p":
			print estimate(calc, $NF, $1);
			break;
		case "f":
			output = 0;
			for ( i = 2; i <= NF; i++ ) {
				output = output + estimate(calc, $i, $1);
			}
			print output;
			break;
		case "g":
			output = 0;
			for ( i = 2; i <= NF; i++ ) {
				output = output + estimate(calc, $i, $(i-1));
			}
			print output;
			break;
	}
}
