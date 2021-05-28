BEGIN {
	# A record ends on )); so the first record will include the boilerplate
	# at the beginning of the file.
	# No field separator, we split the entire record by hand.
	RS = "));"
	ORS = ""
}
/MOCK_[A-Z_]*METHOD/ {
	# Start by closing the previous record.  We don't do that at the end of
	# individual records to properly handle the end of file.
	if (FNR > 1) {
	  print "));"
	}

	# The prefix is the text before the macro in the current record.
	prefix = $0
	gsub(/MOCK_[A-Z_]*METHOD.*/, "", prefix)
	if (prefix != "") {
	  print prefix
	}

	# Extract the macro and its arguments.
	mock = $0
	gsub(/.* +MOCK_/, "MOCK_", mock)
	macro = mock
	gsub(/\(.*/, "", macro)

	# See if the macro is const.
	if (macro ~ /CONST/) {
	  spec = "const, override"
	} else {
	  spec = "override"
	}

	name = mock
	gsub(/.*METHOD[0-9]+(_T)?\(/, "", name)
	gsub(/,.*/, "", name)
	gsub(/[ \n]/, "", name)

	ret = mock
	gsub(".*\\(.*" name ",", "", ret)
	gsub("\\(.*", "", ret)
	gsub(/\n/, "", ret)
	gsub(/^ */, "", ret)
	gsub(/ *$/, "", ret)

	# The params are everything after the return type.
	params = $0
	# Cannot use gsub here because ret may contain the character *
	retpos = index(params, ret "(")
	params = substr(params, retpos + length(ret) + 1)
	split(params, p, ",\n")
	pparams = ""
	for (i in p) {
	  if (i > 1) {
	    pparams = pparams ",\n"
	  }
	  if (p[i] ~ /<.+,.+>/) {
	    pparams = pparams "(" p[i] ")"
	  } else {
	    pparams = pparams p[i]
	  }
	}
	if (ret ~ /,/) {
	  ret = "(" ret ")"
	}
	print "MOCK_METHOD(" ret ", " name ", (" pparams "), (" spec
}
!/MOCK_[A-Z_]*METHOD/{
	if (FNR > 1) {
	  print "));"
	}
	print $0
}
