sed 's/= /= +/g;s/=-/= -/g' HORIZONS.mbox | awk -f horizons.awk > horizons.cfg
