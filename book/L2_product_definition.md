# Level-2 product definition

The Level-2b product is the result of the processing from Level-1b brightness
temperatures to Level-2R gridded geophysical parameters. The product is provided in the resampled
grid of the CIMR instrument at resolution of the C-band channel, which means it contains 547
samples per scan for each horn. Multiplied by the number of scans for each orbit gives the
number of retrieved parameter values per orbit for each geophysical parameter. The product is provided in netCDF format and contains the
following variables following the [CF conventions](http://cfconventions.org/):

(product_variables)=
| Variable name | Unit | dimensions |
| --- | ---- | ---| 
| `time` | days since 2000-01-01 00:00:00 | `time` |
| `wind_speed` | m s$^{-1}$ | n_scans * n_samples_earth * n_horns|
| `total_water_vapor` | kg m$^{-2}$ | n_scans * n_samples_earth * n_horns|
| `cloud_liq_water` | kg m$^{-2}$ | n_scans * n_samples_earth * n_horns|
| `sea_surface_temperature` | K | n_scans * n_samples_earth * n_horns|
| `ice_surface_temperature` | K | n_scans * n_samples_earth * n_horns|
| `sea_ice_fraction` | 1 | n_scans * n_samples_earth * n_horns|
| `multi_year_ice_fraction` | 1 | n_scans * n_samples_earth * n_horns|
| `sea_ice_thickness` | m | n_scans * n_samples_earth * n_horns|
| `wind_speed standard_error` | m s$^{-1}$ | n_scans * n_samples_earth * n_horns|
| `total_water_vapor standard_error` | kg m$^{-2}$ | n_scans * n_samples_earth * n_horns|
| `cloud_liq_water standard_error` | kg m$^{-2}$ | n_scans * n_samples_earth * n_horns|
| `sea_surface_temperature standard_error` | K | n_scans * n_samples_earth * n_horns|
| `ice_surface_temperature standard_error` | K | n_scans * n_samples_earth * n_horns|
| `sea_ice_fraction standard_error` | 1 | n_scans * n_samples_earth * n_horns|
| `multi_year_ice_fraction standard_error` | 1 | n_scans * n_samples_earth * n_horns|
| `sea_ice_thickness standard_error` | m | n_scans * n_samples_earth * n_horns|
| `quality_flag` | 1 | n_scans * n_samples_earth * n_horns| 
| `sea_surface_salinity` | g kg$^{-1}$ | n_scans * n_samples_earth * n_horns|
| `sea_surface_salinity standard_error` | g kg$^{-1}$ | n_scans * n_samples_earth * n_horns|


The dimension follow the definition of the original CIMR L1b product.

The quality flag is a 64-bit mask with the following bits:
(product_flags)=
| Bit | Description | state |
| --- | ---- | --- |
| 0 | valid solution | set if valid |
| 1 | default solver convergence flag (<50 iterations) | set if true |
| 2 | fallback solver used | set if true |
| 3 | fallback solver converged | set if fallback solver converged | set if true | 
| 4 | no convergence | set if true (invalid solution) |
| 5 | anomaly detected | set if true (invalid solution) |
| 6 | invalid wind_speed | set if true |
| 7 | invalid total_water_vapor | set if true |
| 8 | invalid cloud_liq_water | set if true |
| 9 | invalid sea_surface_temperature | set if true |
| 10 | invalid ice_surface_temperature | set if true |
| 11 | invalid sea_ice_fraction | set if true |
| 12 | invalid multi_year_ice_fraction | set if true |
| 13 | invalid sea_ice_thickness | set if true |
| 14 | anomaly in residual | set if true |
| 15-23 | Reserved |  |
| 24 | anomaly in L-band | set if true |
| 25 | anomaly in C-band | set if true |
| 26 | anomaly in X-band | set if true |
| 27 | anomaly in Ku-band | set if true |
| 28 | anomaly in K-band | set if true |
| 29-49 | Reserved |  |
| 50 | Landmask | set for land (fixed) |
| 51 | Ice shelf mask | set for ice shelf |
| 52-64 | Reserved |  |

