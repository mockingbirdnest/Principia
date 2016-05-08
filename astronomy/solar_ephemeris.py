import numpy as np
from astropy.time import Time

# Finds the times (in TDB) of events in two sets of lunar eclipses, given
# the more easily found UTC times.

# First Solar Eclipse was 1950-03-18 (Annular)
# P1 = 13:10:46.0 UT
# U1 = 15:08:31.9
# U4 = 15:55:12.3
# P4 = 17:52:56.2

# 1950-09-12 (Total)
# P1 = 01:23:12.0
# U1 = 02:49:29.0
# U2 = 02:50:46.7
# U3 = 04:26:15.5
# U4 = 04:27:38.0
# P4 = 05:55:36.3

# 1951-03-07 (Annular)
# P1 = 18:03:56.9
# P2 = 20:12:04.8
# P3 = 21:34:33.7
# P4 = 23:42:34.8
# U1 = 19:05:25.1
# U2 = 19:07:03.1
# U3 = 22:39:23.4
# U4 = 22:41:07.0

#1951-09-01 (Annular)
# P1 = 09:54:27.3 UT
# P2 = 12:04:19.2
# P3 = 13:38:35.0
# P4 = 15:48:10.5
# U1 = 10:57:20.2
# U2 = 11:00:03.8
# U3 = 14:42:44.1
# U4 = 14:45:22.1

csv = open('solareclipsetimes.csv', 'w')

P1U14P4 = ['1950-04-02T13:10:46.0', '1950-04-02T15:08:31.9', \
           '1950-04-02T15:55:12.3', '1950-04-02T17:52:56.2']
eclipseTimes = Time(P1U14P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1U1234P4 = ['1950-09-12T01:23:12.0', '1950-09-12T02:49:29.0', \
             '1950-09-12T02:50:46.7', '1950-09-12T04:26:15.5', \
             '1950-09-12T04:27:38.0', '1950-09-12T05:55:36.3']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1U1234P4 = ['1951-03-07T18:03:56.9', '1951-03-07T19:05:25.1', \
             '1951-03-07T19:07:03.1', '1951-03-07T22:39:23.4', \
             '1951-03-07T22:41:07.0', '1951-03-07T23:42:34.8']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

P1U1234P4 = ['1951-09-1T09:54:27.3', '1951-09-1T10:57:20.2', \
             '1951-09-1T11:00:03.8', '1951-09-1T14:42:44.1', \
             '1951-09-1T14:45:22.1', '1951-09-1T15:48:10.5']
eclipseTimes = Time(P1U1234P4, format='isot', scale='utc')
et2 = eclipseTimes.tdb
csv.write(str(et2.jd))
csv.write('\n')

csv.close()
