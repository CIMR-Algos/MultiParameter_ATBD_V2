# Algorithm Input and Output Data Definition (IODD)

## Input data

The input to the retrieval is using the L1B data product from the CIMR instrument, which includes the brightness temperatures of the channels 1.4, 6.9, 10.7, 18.7, and 36.5&nbsp;GHz and their uncertainties.
Technically, missing values are allowed, in case of malfunctioning channels, but the retrieval
uncertainties will be larger. In addition, for a better constrain of the
solution space, {term}`ECMWF` analysis data are highly recommended as additional input (see {ref}`sec:auxiliary_data` below).
The $\mathbf{y}$ and $\mathbf{S}_e$ from {eq}`eq:chi2` are set with the brightness temperatures and their uncertainties from the L1B data product.
The input data for the multi-parameter retrieval includes all channels at
horizontal and vertical polarization and their uncertainties coming from the
L1b processor. In addition, time and location is required for the resampling of the variables. The algorithm is then applied on the C-band footprints, with all other channels being resampled to it. The
uncertainties are assumed Gaussian in the retrieval. 

| Field | Description | Shape/Amount |
| ---   | ----------- | ------------ |
| L1B TB | L1B Brightness Temperature of all channels (both H and V polarization) | full swath or section of it (Nscans, Npos) |
| L1B NeΔT | Random radiometric uncertainty of all channels | full swath or section of it (Nscans, Npos) |


## Output data

The output data will include the retrieved geophysical parameters listed in equation
{eq}`eqxy` and their posterior uncertainties. In addition, quality flags
are derived for each quantity and the retrieval procedure in general.
The parameters are are: {term}`WS`, {term}`TWV`,
{term}`CLW`, {term}`SST`, {term}`IST`, {term}`SIC`, {term}`MYI` and
{term}`SIT`, {term}`SSS`, as well as time and location for the resampling. The uncertainties are assumed
Gaussian in the retrieval and are provided along with the parameters. 
The complete list of all output variables output in the EASE2 grid for each hemisphere is given in the table below.

| Abbreviation | CF variable name | Description | Shape/Amount |
| --- | --- | --- | --- |
| WS | wind_speed | Wind speed | (1440,1440) |
| WSerr | wind_speed_error | Wind speed error | (1440,1440) |
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
| MYI | multi_year_ice_fraction | Multi-year ice fraction | (1440,1440) |
| MYIerr | multi_year_ice_fraction_error | Multi-year ice fraction error | (1440,1440) |
| SIT | sea_ice_thickness | Sea ice thickness | (1440,1440) |
| SITerr | sea_ice_thickness_error | Sea ice thickness error | (1440,1440) |
| SSS | sea_surface_salinity | Sea surface salinity | (1440,1440) |
| SSSerr | sea_surface_salinity_error | Sea surface salinity error | (1440,1440) |


(sec:auxiliary_data)=
## Auxiliary data
{term}`ECMWF` surface analysis data is used as background values for the retrieval. The
variables used are {term}`WS`, {term}`TWV`, {term}`CLW`, {term}`T2M`,
{term}`TSK`. They are used to fill the $\mathbf{S}_a$ matrix and the
$\mathbf{x}_a$ in equation {eq}`eq:chi2`. For near real-time retrieval, the
$\mathbf{S}_a$ and $\mathbf{x}_a$ can use monthly or seasonal values, as
the retrieval is not sensitive to the exact values of the background variables with the rather wide covariance matrix.
```{note}
In this document a retrieval with fixed background values and TB error is used (execpt in {ref}`first-assessment`). This may alter the results of the retrieval on the test cards compared to actual data where the background values are from ECMWF analysis data.
```