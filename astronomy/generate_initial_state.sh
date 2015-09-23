sed 's/= /= +/g;s/=-/= -/g' HORIZONS*.mbox | \
  awk -f generate_initial_state.awk > \
  initial_state.proto.txt
