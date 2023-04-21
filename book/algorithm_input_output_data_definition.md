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
