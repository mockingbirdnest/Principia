# Satellite orbits in SP3 format

This directory contains satellite orbit files in the standard product 3 orbit
format, used by Principia as test data.

These files come from various agencies; we document their provenance here.

## International DORIS Service (IDS) products

See https://ids-doris.org/ids/data-products/tables-of-data-products.html,
https://ids-doris.org/documents/Data-Structure-Formats.pdf,
https://ids-doris.org/ids/organization/analysis-centers.html,
https://cddis.nasa.gov/docs/2016/cddis_IDSreport_2016.pdf.

- [grgtop03.b97344.e97348.D_S.sp3](grgtop03.b97344.e97348.D_S.sp3):<br>
  ftp://doris.ensg.eu/pub/doris/products/orbits/grg/top/ <br>
  ftp://doris.ign.fr/pub/doris/products/orbits/grg/top/ <br>
  ftp://cddis.gsfc.nasa.gov/doris/products/orbits/grg/top/
- [grgja203.b08243.e08247.D_S.sp3](grgja203.b08243.e08247.D_S.sp3):<br>
  ftp://doris.ensg.eu/pub/doris/products/orbits/grg/ja2/ <br>
  ftp://doris.ign.fr/pub/doris/products/orbits/grg/ja2/ <br>
  ftp://cddis.gsfc.nasa.gov/doris/products/orbits/grg/ja2/
- [ssaja102.b03007.e03017.DGS.sp3](ssaja102.b03007.e03017.DGS.sp3):<br>
  ftp://doris.ensg.eu/pub/doris/products/orbits/ssa/ja2/ <br>
  ftp://doris.ign.fr/pub/doris/products/orbits/ssa/ja2/ <br>
  ftp://cddis.gsfc.nasa.gov/doris/products/orbits/ssa/ja2/

## International GNSS Service (IGS) products

See http://www.igs.org/products,
https://kb.igs.org/hc/en-us/articles/115003935351,
https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/orbit_products.html.

- [esa11802.eph](esa11802.eph):<br>
  ftp://gssc.esa.int/gnss/products/1180/ <br>
  ftp://igs.ensg.ign.fr/pub/igs/products/1180/ <br>
  ftp://lox.ucsd.edu/pub/products/1180/ <br>
  ftp://cddis.gsfc.nasa.gov/gnss/products/1180/

### Multi-GNSS Experiment (MGEX)

See http://mgex.igs.org/IGS_MGEX_Products.php.

- [COD0MGXFIN_20181260000_01D_05M_ORB.SP3](COD0MGXFIN_20181260000_01D_05M_ORB.SP3):<br>
  ftp://igs.ensg.eu/pub/igs/products/mgex/2000/ <br>
  ftp://igs.ign.fr/pub/igs/products/mgex/2000/ <br>
  ftp://cddis.gsfc.nasa.gov/gnss/products/mgex/2000/
- [COD0MGXFIN_20183640000_01D_05M_ORB.SP3](COD0MGXFIN_20183640000_01D_05M_ORB.SP3):<br>
  ftp://igs.ensg.eu/pub/igs/products/mgex/2034/ <br>
  ftp://igs.ign.fr/pub/igs/products/mgex/2034/ <br>
  ftp://cddis.gsfc.nasa.gov/gnss/products/mgex/2034/
- [WUM0MGXFIN_20190270000_01D_15M_ORB.SP3](WUM0MGXFIN_20190270000_01D_15M_ORB.SP3):<br>
  ftp://igs.ensg.eu/pub/igs/products/mgex/2038/ <br>
  ftp://igs.ign.fr/pub/igs/products/mgex/2038/ <br>
  **not found** at ftp://cddis.gsfc.nasa.gov/gnss/products/mgex/2038/


### International GLONASS Experiment (IGEX) / International GLONASS Service Pilot Project (IGLOS-PP)

As far as I can tell, the `glonass` directory of the CDDIS data centre
originates with IGEX and later IGLOS-PP.

- [mcc14000.sp3](mcc14000.sp3):<br>
  ftp://cddis.nasa.gov/glonass/products/1400/

## International Laser Ranging Service (ILRS) products

See https://ilrs.cddis.eosdis.nasa.gov/data_and_products/products/index.html.

- [asi.orb.etalon2.171209.v70.sp3](asi.orb.etalon2.171209.v70.sp3): <br>
  ftp://edc.dgfi.tum.de/pub/slr/products/orbits/etalon2/171209/ <br>
  ftp://cddis.gsfc.nasa.gov/pub/slr/products/orbits/etalon2/171209/
- [ilrsa.orb.lageos2.160319.v35.sp3](ilrsa.orb.lageos2.160319.v35.sp3),<br>
  [ilrsb.orb.lageos2.160319.v35.sp3](ilrsb.orb.lageos2.160319.v35.sp3):<br>
  ftp://edc.dgfi.tum.de/pub/slr/products/orbits/lageos2/160319/ <br>
  ftp://cddis.gsfc.nasa.gov/pub/slr/products/orbits/lageos2/160319/
- [ilrsa.orb.lageos2.180804.v70.sp3](ilrsa.orb.lageos2.180804.v70.sp3):<br>
  ftp://edc.dgfi.tum.de/pub/slr/products/orbits/lageos2/180804/ <br>
  ftp://cddis.gsfc.nasa.gov/pub/slr/products/orbits/lageos2/180804/

## Miscellaneous

### National Geospatial-intelligence Agency (NGA)

See http://earth-info.nga.mil/GandG/update/index.php?dir=gnss&action=gnss#tab_ephemeris.
Note that we use satellite centre of mass ephemerides, rather than antenna
phase centre ephemerides.
The former are found under ftp://ftp.nga.mil/pub2/gps/pedata/ in directories
*yyyy*pe, have the extension `.eph`, and mention `SATCOM` in the comment
records; the latter under ftp://ftp.nga.mil/pub2/gps/apcpe/ in directories
*yyyy*apc, with the extension `.apc`, and `SATAPC` in the comment records.

- [nga20342.eph](nga20342.eph):<br>
  ftp://ftp.nga.mil/pub2/gps/pedata/2019pe/
