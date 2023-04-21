# Introduction, purpose and scope

This is the {term}`ATBD` of the multi parameter retrieval. The method is based
on the works of {cite}`Pedersen1991,Scarlat2017,Scarlat2018,Scarlat2020` and
will be described in this document. It is applied to {term}`CIMR` L1b data which
is resampled to a common footprint for all frequencies (we call L1R). It is
using all {term}`CIMR` frequency channels, namely 1.4, 6.9, 10.7, 18.7, and 36.5
GHz and {term}`ECMWF` Analysis as input. The output of this multi parameter
retrieval is in the same resolution format, i.e., L2R. It is physically
consistent and can be used in turn as a priori for the other retrieval
parameters. 

This document is describing the algorithm and processing steps of the L2R multi
parameter retrieval product. The document is intended for the {term}`CIMR` users and
interested parties. It is not intended to replace a product user guide. The
algorithm is implemented in the *draft* multi parameter retrieval algorithm which
is distributed in the "algorithm" directory of this repository.

To this date, there are no instruments in operation which measure all
frequencies at the same time as the {term}`CIMR`. Therefore, for the development
of the algorithm, data from {term}`AMSR2` and {term}`SMOS` were used as input.
This has a direct impact on the retrieval performance as the collocation of the
measurements is not fixed. For surface parameters, this might be less of an
issue as the surface is not varying in so short time scales. For atmospheric
parameters, however this is a major issue. While the retrieval relies on the
idea of a fixed state observation, the current data for the development of the
algorithm is coming from a potential divergent state. This is when
{term}`AMSR2` and {term}`SMOS` have different overflight times. This may lead
to inconsistent results with higher residuals and may even prevent convergence
of the retrieval. At the moment daily means are used as input, which is only reliable for a stable atmospheric and surface state.
