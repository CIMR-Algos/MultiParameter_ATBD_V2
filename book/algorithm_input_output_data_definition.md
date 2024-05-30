# Algorithm Input and Output Data Definition (IODD)

# Input data

The input data for the multi-parameter retrieval includes all channels at
horizontal and vertical polarization and their uncertainties coming from the
L1b processor. In addition, time and location is required for calculation of
field of resampled TBs. The evaluation of the algorithm is planned on the
C-band footprints, with all other channels being resampled to it. The
uncertainties are assumed Gaussian in the retrieval. 

| Field | Description | Shape/Amount |
| ---   | ----------- | ------------ |
| L1B TB | L1B Brightness Temperature of all channels (both H and V polarization) | full swath or section of it (Nscans, Npos) |
| L1B NeΔT | Random radiometric uncertainty of all channels | full swath or section of it (Nscans, Npos) |


The complete list of all output variables output in the EASE2 grid for each hemisphere is given in the table below.

| Abbreviation | CF variable name | Description | Shape/Amount |
| --- | --- | --- | --- |
| WSP | wind_speed | Wind speed | (1440,1440) |
| WSPerr | wind_speed_error | Wind speed error | (1440,1440) |
| TWV | total_water_vapour | Total water vapour | (1440,1440) |
| TWVerr | total_water_vapour_error | Total water vapour error | (1440,1440) |
| CLW | cloud_liquid_water | Cloud liquid water | (1440,1440) |
| CLWerr | cloud_liquid_water_error | Cloud liquid water error | (1440,1440) |
| SST | sea_surface_temperature | Sea surface temperature | (1440,1440) |
| SSTerr | sea_surface_temperature_error | Sea surface temperature error | (1440,1440) |
| IST | ice_surface_temperature | Ice surface temperature | (1440,1440) |
| ISTerr | ice_surface_temperature_error | Ice surface temperature error | (1440,1440) |
| SIC | sea_ice_concentration | Sea ice concentration | (1440,1440) |
| SICerr | sea_ice_concentration_error | Sea ice concentration error | (1440,1440) |
| MYIF | multi_year_ice_fraction | Multi-year ice fraction | (1440,1440) |
| MYIFerr | multi_year_ice_fraction_error | Multi-year ice fraction error | (1440,1440) |
| SIT | sea_ice_thickness | Sea ice thickness | (1440,1440) |
| SITerr | sea_ice_thickness_error | Sea ice thickness error | (1440,1440) |
| SSS | sea_surface_salinity | Sea surface salinity | (1440,1440) |
| SSSerr | sea_surface_salinity_error | Sea surface salinity error | (1440,1440) |









## Output data

The output data consist of the parameters, namely {term}`WSP`, {term}`TWV`,
{term}`CLW`, {term}`SST`, {term}`IST`, {term}`SIC`, {term}`MYIF` and
{term}`SIT`, {term}`SSS`, as well as time and location for the resampling. The uncertainties are assumed
Gaussian in the retrieval and are provided along with the parameters. 

## Auxiliary data
As auxiliary data {term}`ECMWF`
surface analysis data is used for background values in the retrieval. The
variables used are {term}`WSP`, {term}`TWV`, {term}`CLW`, {term}`T2M`,
{term}`TSK`. 


## Ancillary data
No ancillary data is used in the retrieval.
