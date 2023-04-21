# Abstract

This is the ATBD of the multi parameter retrieval. It consists of several
empirical and physical models that are used to retrieve the parameters of
{term}`SIC`, {term}`SIT`, {term}`MYIF`, {term}`WSP`, {term}`TWV`, {term}`CLW`,
{term}`SST`, {term}`IST` and {term}`SSS` in a self-consistent way. The method
respects individual contributions to the brightness temperatures at each
frequency and their uncertainties provided from the {term}`CIMR` L1b data. The
individual error covariances of the different parameters and an a priori state
can be used as input. The a priori is partly taken from {term}`ECMWF` analysis
for compatible quantities. The output of this multi
parameter retrieval is physically consistent and can be used in turn as a priori
for the other retrieval parameters. The method is based on the works of
{cite}`Pedersen1991,Scarlat2017,Scarlat2018,Scarlat2020`.

