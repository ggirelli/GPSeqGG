
BEGIN {
	OFS = FS = "\t";
	convert["X"] = 23;
	convert["Y"] = 24;
}

{
	chrid = substr($1, 4);
	if ( chrid in convert ) {
		chrid = convert[chrid];
	}
	print chrid OFS $0;
}
