# NCPTrans 
*Quick and easy horizontal coordinate system and vertical reference system transformations on the Dutch Continental Shelf.*

# WORK IN PROGRESS

## Functionality
This software can be used to:

 - Transform coordinates between:
	 - RD-xy
	 - ETRS89-latlon
	 - UTM Zone 31N-eastingnorthing
	 - UTM Zone 32N-eastingnorthing
 - Transform vertical references between:
	 - NAP (NLGEO2018)
	 - LAT (NLLAT2018)
	 - GRS80 ellipsoidal height
 - Sample:
	 - NAP (NLGEO2018) - LAT (NLLAT2018) separations
	 - NLGEO2018 - GRS80 separations
	 - NLLAT2018 - GRS80 separations

## RDNAPTRANS
Transformations related to RD or NAP are based on the document '*RDNAPTRANS2018™: Coordinate transformation to and from Stelsel van de Rijksdriehoeksmeting and Normaal Amsterdams Peil*' version 27 June 2022 by NSGI. This document is available for download at [https://nsgi.nl](https://nsgi.nl). Only variant 1 of the transformation is implemented.

Although this software meets the accuracy requirements set by prior document and the NSGI validation service, **THIS SOFTWARE IS NOT APPROVED TO USE THE NAME RDNAPTRANS2018**. No guarantees are provided with regard to the correctness of this software or its results. One should be extremely cautious using the results from this software.

## Installation
This software does not rely on external libraries. Simply download this repository and import `NCPTrans` in the script of your choice.

The grids `rdcorr2018.txt`, `nlgeo2018.txt` and `nllat2018.txt` are not shipped with this software and needs to be downloaded from either [https://nsgi.nl](https://nsgi.nl) or [https://hydro.nl](https://hydro.nl). After downloading, put these three grids in the same folder as `constants.py` and run `helpers.prepare_grid_db()`. This will convert all necessary grids into a sqlite3 database file.

After setting up, you may want to validate the RDNAPTRANS transformation procedure. More information at [NSGI Validatieservice](https://www.nsgi.nl/coordinatenstelsels-en-transformaties/tools/validatieservice). The script `validator.py` may be used to support the validation. Please keep in mind the usage of these transformations is at your own risk. No guarantees are given with respect to to correctness of the results. Even if the validation passes.

## Basic usage
```python
import NCPTrans

# RDNAP -> ETRS89
x = 108360.7859
y = 415757.2741
h = 3.29
lat, lon, h = NCPTrans.rdnap_to_etrsh(x, y, h)

# ETRS89 -> RDNAP
lat = 52.011555
lon = 3.976206
h = 45.3
x, y, nap = NCPTrans.etrsh_to_rdnap(lat, lon, h)
```

## Extended transformations
```python
from NCPTrans import transform

result = transform(x, y, h, src_xy, trg_xy, src_h, trg_h)
# src_xy, trg_xy -> ['RD', 'ETRS89']
# src_h -> ['ETRS89', 'NAP', 'LAT']
# trg_h -> ['ETRS89', 'NAP', 'LAT', 'NAP_LAT_SEPARATION']

# Get NAP-LAT difference at point:
lat = 52.011555
lon = 3.976206
(lon, lat, NAP-LAT-difference) = transform(lon, lat, src_xy='ETRS89', trg_h='NAP_LAT_SEPARATION')
```
## Notes
This software is not intended as a computationally efficient transformation. It should only be used for the examination of a single point or very small dataset.

---  

willem van opstal | 2024