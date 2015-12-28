import numpy as np
from astropy.time import Time

# Finds the times (in TDB) of events in two sets of lunar eclipses, given
# the more easily found UTC times.

# First Lunar Eclipse was 1950-04-02
# P1 = 18:10:49 UT
# U1 = 19:09:19
# U2 = 20:30:38
# U3 = 20:57:33
# U4 = 22:18:54
# P4 = 23:17:21

# Last Lunar Eclipse will be 2048-01-01 
# P1 = 03:52:39 UT
# U1 = 05:05:17
# U2 = 06:24:27
# U3 = 07:20:23
# U4 = 08:39:33
# P4 = 09:52:05

P1U1234P4 = ['1950-04-02T18:10:49', '1950-04-02T19:09:19', \
             '1950-04-02T20:30:38', '1950-04-02T20:57:33', \
             '1950-04-02T22:18:54', '1950-04-02T23:17:21']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv = open('eclipsetimes.csv', 'w')
csv.write(str(et2.jd))
csv.write('\n')

P1U1234P4 = ['2048-01-01T03:52:39', '2048-01-01T05:05:17', \
             '2048-01-01T06:24:27', '2048-01-01T07:20:23', \
             '2048-01-01T08:39:33', '2048-01-01T09:52:05']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.close()
