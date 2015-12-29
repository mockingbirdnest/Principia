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

# Additional total and 2x penumbral eclipses in 1950/51
# 1950-09-26
# P1 = 01:21:43 UT
# U1 = 02:31:48
# U2 = 03:54:33
# U3 = 04:38:49
# U4 = 06:01:33
# P4 = 07:11:47
# 1951-03-23
# P1 = 08:50:00
# P4 = 12:24:19
# 1951-09-15
# P1 = 10:29:16
# P4 = 14:23:52

# 2x partial eclipses in 1952
# 1952-02-11
# P1 = 22:08:20 UT
# U1 = 00:04:17
# U4 = 01:14:24
# P4 = 03:10:15
# 1952-08-05
# P1 = 17:28:13 UT
# U1 = 18:33:49
# U4 = 21:01:00
# P4 = 22:06:35

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

P1U1234P4 = ['1950-09-26T01:21:43', '1950-09-26T02:31:48', \
             '1950-09-26T03:54:33', '1950-09-26T04:38:49', \
             '1950-09-26T06:01:33', '1950-09-26T07:11:47']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1P4 = ['1951-03-23T08:50:00', '1951-03-23T12:24:19']
eclipseTimes = Time(P1P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1P4 = ['1951-09-15T10:29:16', '1951-09-15T14:23:52']
eclipseTimes = Time(P1P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1U14P4 = ['1952-02-11T22:08:20', '1952-02-12T00:04:17', \
           '1952-02-12T01:14:24', '1952-02-12T03:10:15']
eclipseTimes = Time(P1U14P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1U14P4 = ['1952-08-05T17:28:13', '1952-08-05T18:33:49', \
           '1952-08-05T21:01:00', '1952-08-05T22:06:35']
eclipseTimes = Time(P1U14P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1U1234P4 = ['2048-01-01T03:52:39', '2048-01-01T05:05:17', \
             '2048-01-01T06:24:27', '2048-01-01T07:20:23', \
             '2048-01-01T08:39:33', '2048-01-01T09:52:05']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.close()
