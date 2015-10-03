if [[ $# < 2 ]]
then
  echo "Usage: generate_initial_state.sh infile juliandate > outfile"
  exit 1
fi
infile=$1
juliandate=$2
sed 's/= /= +/g;s/=-/= -/g' "${infile}" | \
  awk -f generate_initial_state.awk juliandate="${juliandate}"
