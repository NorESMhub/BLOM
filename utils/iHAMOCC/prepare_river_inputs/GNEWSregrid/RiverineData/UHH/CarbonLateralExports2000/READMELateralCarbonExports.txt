LATERAL CARBON EXPORTS FILES ("NEWS2000")
Emilio Mayorga. mayorga@apl.washington.edu. 3/6/2011


* NEWSOutputExports_2000Carbon*.csv: Basin model output from 
  Global NEWS (DOC, POC, TSS, etc) and Hartmann et al (2009) DIC models.
  basinid is the STN-30p basin id (see below).

* loadsreg_glb_ocea-sea_NEWSOutputExports_2000Carbon_*.csv
  ocean+sea regional aggregations.

* DIC versions:
  H2009 - Hartmann et al, 2009
  Fer25 - H2009 adjusted downward in weathered, tropical ferralsols soils

* STN30v6ngBASINID.shp (and associated files), polygon GIS shape file,
  with every basin represented as a single polygon (a GIS "multi-polygon"
  representation; some basins may have separate polygon components all
  making up a single basin). The only attribute is STN-30p basin IDs
  as used in NEWS 2 (including the NEWS scenarios work).

* STN30v6ngNEWS.csv, a CSV text file showing all main general basin
  attributes as used in NEWS 2. See GlobalNEWS2README.xls for description
  of attributes.

* Endorheic basins: Join the attribute table with the shape file. In the
  attribute "ocean", endorheic basins are ones coded with the value "Land".
