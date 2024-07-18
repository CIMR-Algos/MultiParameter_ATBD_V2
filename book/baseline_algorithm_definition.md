# Baseline Algorithm Definition
## Introduction

In the following the main algorithm for the multi parameter retrieval is
described in detail. The algorithm is based on the works of
{cite}`Pedersen1991,Scarlat2017,Scarlat2018,Scarlat2020`. The algorithm is
divided into several steps, which are described in the following.


## CIMR Level-1b re-sampling approach


The {term}`CIMR` Level-1b data is first resampled to the C-band footprints using
nearest neighbor resampling. The target footprint with C-band is a compromise between
resolution and accuracy given sensitivities of the inversion of the forward
model. The main restraint being the L-band channel which has a
high influence on the retrieval, on ocean because of SSS and sea ice due to SIT. Higher up sampling of L-band than to the C-band resolution would require accurate knowledge about the antenna pattern of the instrument and a trade-off between accuracy and resolution.


## Level-2 end to end algorithm functional flow diagram

The diagram below shows the functional flow of the algorithm.

```{mermaid}
graph TD
	subgraph Input data
		subgraph CIMR L1b
			TBs[TOA TBs]
			TBe[TOA TB error]
		end
		subgraph External
			ECMWF[ECMWF analysis]
			ERA5["Historical ERA5"]
		end
	end
	subgraph Inversion
		F[Forward model]
		cost[Cost function]
		apriori[A priori] 
		OS[Optimal state]

	end

	subgraph Output data L1R
		geo[Geophysical variables]
		unc[Geophsysical uncertainties]
		outtb[Brightness tempreature residuals]
		flags[Quality flags]
	end

	resampling[Resampling processor L1R]

	TBs --> resampling
	TBe --> resampling
	ECMWF --> resampling
	resampling -- Measurement vector --> cost
	resampling -- Measurement uncertainty --> cost
	resampling -- State vector --> apriori 
	resampling -- State uncertainty --> apriori

	ERA5 --Covariance estimate--> apriori
	apriori --> cost
	F-->cost
	cost-->F

	cost --> OS
	OS --> geo
	OS --> unc 
	OS --> outtb
	OS --> flags
```


## Algorithm assumptions and simplifications

Compared to individual state-of-the-art algorithms for the retrieval of individual quantities, the
multi parameter retrieval uses a simplified approach. Particular trade-offs
include:
* Some ocean parameters contain many empirical values, developed and validated
  for different instruments which might have to be adjusted to match CIMR
	  instrument characteristics.
* The sea ice thickness is retrieved as a combination of first-year ice
  emissivity and the ocean emissivity. As this, it is prone to noise in the
  ocean emission, which can lead to erroneous occurrence of sea ice of low
  thickness, as the first-year ice emissivity is identical to the ocean
  emissivity at 0 cm thickness.
* The empirical parametrization of sea ice thickness stems from uncertain
  atmospheric conditions, so that the ice thickness dependence of higher
  frequency channels might be confounded with the dependence on the atmospheric
  conditions.
* The sensitivity of the CIMR channels to CLW, in particular over sea ice is
  low, so that the uncertainty of CLW is high. In addition, the effect of
  decreased CLW over first-year ice is similar to the effect of an increased
  MYI fraction, so that these two signals cannot be adequately separated in the
  current retrieval.
* In the present version of the algorithm, the forward model is assumed
  near-linear in its characteristic and the retrieval is unconstrained. This
  might lead to a bias in the retrieval of the parameters when other parameters
  might converge to nonphysical values. In test cases, this happened mostly
  with
	- IST > 273.15 K
	- SST < 273.15 K
	- TWV < 0 kg/m$^2$
* There are more modern approaches for the forward model of the ocean emission,
  in particular at L-band. This will influence the retrieval of the ocean
  parameters, in particular {term}`SSS` and {term}`SST`. The current algorithm is the approach
  of {cite}`Ruf2003` and {cite}`Scarlat2020`, but a switch to another
  parametrization like {cite}`Meissner2018` can be considered. However, the SSS is
  the parameter with the weakest signal in the forward model, so that other
  parametrization have to be improved, to justify this effort.
* The forward model is restricted to winter conditions over sea ice. The summer
  conditions are highly variable across frequencies and cannot be adequately
  described by the current forward model. In particular the definition of first-year
  ice and multiyear ice is ambiguous in melting conditions. In addition,
  the {term}`SIT` parametrization is only valid for winter conditions. 



## Functional description of each algorithm step

The retrieval is following a typical scheme with the objective to minimize

$$
\begin{align}
χ^2(\mathbf y,\mathbf x,\mathbf S_e,\mathbf S_a,\mathbf x_a) = \left(\mathbf y-F(\mathbf x)\right)^\mathbf{T}\mathbf S_e^{-1}(\mathbf
y-F(\mathbf x))+(\mathbf x_a-\mathbf x)^{\mathbf T}\mathbf S_a^{-1}(\mathbf x_a
- \mathbf x) 
\end{align}
$$ (eq:chi2)


where 

```{math}
:label: eqxy
\begin{aligned}
\mathbf y= \begin{bmatrix}
T_{b,h,1.4}\\
T_{b,v,1.4}\\
T_{b,h,6.9}\\
T_{b,v,6.9}\\
T_{b,h,10.7}\\
T_{b,v,10.7}\\
T_{b,h,18.7}\\
T_{b,v,18.7}\\
T_{b,h,36.5}\\
T_{b,v,36.5}\\
\end{bmatrix}, \quad
\mathbf x= \begin{bmatrix}
\text{WS}\\
\text{TWV}\\
\text{CLW}\\
\text{SST}\\
\text{IST}\\
\text{SIC}\\
\text{MYI}\\
\text{SIT}\\
\text{SSS}\\
\end{bmatrix}
\end{aligned}
```

i.e., $\mathbf{y}$ is the vector of input brightness temperatures (measurement vector) and
$\mathbf x$ is the vector with the corresponding physical quantities (state
vector), $\mathbf S_e$ is the error covariance matrix of the input brightness
temperatures, $\mathbf S_a$ is the covariance matrix of the a priori values,
$\textbf x_a$ is the a priori value, and $F$ is the forward operator, which is
the compositional forward model described in {ref}`fw-model`.  The optimal
solution $\hat {\mathbf y}$ which minimizes equation {eq}`eq:chi2` (and maximizes the posterior probability density function) is found by
iteration. The uncertainty of $\mathbf y$ can then be obtained as 

```{math}
:label: eq:error
\mathbf {\hat S}_a= (\mathbf S_a^{-1}+\mathbf{ \hat M}^\mathbf T \mathbf S_e^{-1}\mathbf {\hat M})^{-1}, 
```

with $\mathbf {\hat M}$ is the Jacobian of the optimal solution of
{eq}`eq:chi2`. The uncertainty of the individual parameters is then given by the
diagonal elements of $\mathbf{\hat S}_a$.


```{note}
$\mathbf S_e$ is often set higher than the pure radiometric uncertainty of the brightness temperatures, to account for the uncertainty of the forward model. This can be named an *effective error covariance of the measurements*. The effective contribution from the forward model error is not evaluated yet, but estimated to be below 2 K for all channels.
```
(fw-model)=
## Mathematical description of the Forward Model

The forward model is the compositional forward model. It consists of individual
components, namely the ocean surface, the sea ice, and the atmosphere. At the
low frequency channels of the {term}`CIMR` satellite, the sensitivity to
atmospheric parameters is relatively small, but nevertheless the atmosphere
needs to be considered.
The forward model for ocean and atmosphere for frequencies 6.9, 10.7, 18.7, and 36.5&nbsp;GHz is used from the {term}`AMSR2` {term}`ATBD` from {cite}`Wentz2000`. 
The surface contribution to the brightness temperature is given by
```{math}
:label: eq:surface
T_{b,s} = C_{\text{ow}}ε_{\text{ow}}T_{\text{ow}}+C_{\text{fyi}}ε_{\text{fyi}}T_{\text{fyi}}+C_{\text{myi}}ε_{\text{myi}}T_{\text{myi}},
```


where $C_{\text{ow}}$, $C_{\text{fyi}}$, and $C_{\text{myi}}$ are the area
fraction of ocean water, first-year ice, and multiyear ice, respectively.
$ε_{\text{ow}}$, $ε_{\text{fyi}}$, and $ε_{\text{myi}}$ are the emissivity of
ocean water, first-year ice, and multiyear ice, respectively. $T_{\text{ow}}$,
$T_{\text{fyi}}$, and $T_{\text{myi}}$ are the brightness temperature of ocean
water, first-year ice, and multiyear ice, respectively. $C_{\text{ow}}$,
$C_{\text{fyi}}$, and $C_{\text{myi}}$ are adding up to one. 


### Sea ice type
For sea ice, the frequency-dependent emissivities
for first-year ice and multiyear ice are derived from {cite}`Mathew2009` as 
```{math}
:label: eq:Nizy
U_{T, t, p}=(a_{t}T_{C}+b_{t}+273.15)ε_{t, p}
```
with $U_{t,h}$ being the temperature-corrected upwelling brightness temperature for the polarization $p$
at the air temperature at the surface $T_{C}$ (in °C), $a_{t}$ and $b_{t}$ are the frequency
dependent coefficients from {cite}`Mathew2009` (see {numref}`tab:emtemp`), and
$ε_{t,p}$ is the frequency-dependent emissivity for the polarization $p$ and ice type $t$ from {numref}`tab:c_ice`.
{cite}`Mathew2009` did not provide parametrization for L-band, but to match a
temperature dependence of emitted radiation at L-band we derive the effective
temperature as with {eq}`eq:Nizy` using $a_t=0.1$ and $b_t=0$ for both, first-year and multiyear ice. The emissivities used for L-band are taken from {cite}`Scarlat2020` with slight modifications and are given in {numref}`tab:c_ice`.

```{note}
$T_{C}$ in this parametrization is the IST in the retrieval as the only temperature dependence. 
The assumption is that in thermal equilibrium the IST is equal to the air temperature at the surface.
```

### Sea ice thickness
The ice thickness dependence at all frequencies is introduced via modification of the emission from the ice surface for the first-year ice fraction.
The ice emissivity is modified by the ice thickness $\text{SIT}$ according to
```{math}
:label: eq:ice_thickness
T_{b,p} = a_{p}-(a_{p}-b_{p})\exp\left(-\frac{\text{SIT}}{c_{p}}\right)
```
with the index $p$ indicate polarization ($h$ or $v$) the coefficients $a_{p}$,
$b_{p}$, and $c_{p}$ from {cite}`Scarlat2020` (see {numref}`tab:fy_thick`). To
combine ice temperature and ice thickness dependence, the ice emissivity is
modified with the ice thickness by substituting $a_{p}$ with the ice temperature sensitivity
$U_{T, t, p}$ from {eq}`eq:Nizy` for $t=\text{FYI}$. This was not performed in
{cite}`Scarlat2020` but is essential for the minimization of the cost function
to not introduce discontinuities in the forward model. The MYI emissivity is not
affected by the ice thickness in this forward model, as it was not part of the
thickness sensitivity study in {cite}`Scarlat2020`. As a consequence, the
retrieval of ice thickness accounts only for the FYI fraction and thus is very
noisy when the FYI fraction is small. However, the winter season during freeze-up conditions, this is not much of a restrictions as most of the ice is first-year ice thin multiyear ice is scarce.

### Ocean surface model

While for the sea ice an empirical model is used for the emissivity, for the
emission from the ocean the model uses the Fresnel reflection coefficient as a
basis, which relies on the dielectric constant of the sea water. The emission from calm sea water after {cite}`Meissner2012` is given by

```{math}
:label: eq:em_ocean
E_{0p} &= 1-|r_p|^2\\ 
r_v &= \frac{ε\cos(θ_i)-\sqrt{ε-\sin^2(θ_i)}}{ε\cos(θ_i)+\sqrt{ε-\sin^2(θ_i)}}\\
r_h &= \frac{\cos(θ_i)-\sqrt{ε-\sin^2(θ_i)}}{\cos(θ_i)+\sqrt{ε-\sin^2(θ_i)}}
```

with $ε$ being the dielectric constant of the sea water.
To account for the roughness and other disturbances on the ocean surface, the power reflectivity at each polarization $R_{0p}=|r_p|^2$ needs to be adjusted. The contribution of different ocean surface types are modeled by different parameters again following {cite}`Wentz2000`. The reflectivity is a composition of foam covered ocean and clear ocean. With a reduction of the ocean surface emission through the foam by a factor κ, the composite reflectivity is given by
```{math} 
R = (1-f_{\text{foam}})\cdot R_{\text{clear}}+f_{\text{foam}}\cdot κ\cdot R_{\text{clear}}
``` 
with $f_{\text{foam}}$ being the fraction of the ocean surface covered by foam and $R_{\text{clear}}$ being the reflectivity of the clear ocean. With a small loss from diffraction expressed in the term $β$ we can express the reflectivity as 
```{math}
:label: eq:genroughness
R_{\text{clear}} = (1-\beta)R_{\text{geo}}
```
with $R_{\text{geo}}$ being the reflectivity from a standard geometric optics
model{cite}`Wentz1975` defined later in {eq}`eq:roughness`. The combination can then be expressed by combining the foam and the
diffraction term into one quantity $F=f_{\text{foam}} + \beta -f_{foam}\cdot β - f_{foam} κβ$, which is a monotonic function of wind speed and is addressed by
{cite}`Wentz2000` as *catch-all* term. They determined $F$ empirically from
experiments with various radiometers. A fit for F is given by 
```{math}
:label: eq:catchall
\begin{aligned}
F & = m_1W  &(W<W_1)\\
F & = m_1W + \frac{(m_2-m_1)(W-W_1)^2}{2(W_2-W_1)}  &(W_1 \leq W \leq W_2)\\
F & = m_2W - (m_2-m_1)(W_2+W_1)  &(W>W_2),
\end{aligned}
```
a quadratic spline with knots at $W_1=3~\text{m/s}$ and $W_2=12~\text{m/s}$ for v- and $W_1=7~\text{m/s}$ and $W_2=12~\text{m/s}$ for h-polarization. The coefficients $m_1$ and $m_2$ are given in {numref}`tab:mc_m` in the appendix.

The geometrical optics surface roughness ($R_{geo}$ from {eq}`eq:genroughness`  is used as modeled by {cite}`Wentz2000` as
```{math}
:label: eq:roughness
R_{geo} = R_0 - (r_0 + r1(θ_i-53) + r_2(T_S-288) + r_3(θ_i-53)(T_S-288))W.
```
Here $R_0$ is the specular reflection, θ$_i$ is the incidence angle, T$_S$ is
the sea surface temperature, and W is the wind speed. The coefficients for each
polarization and frequency are in the appendix in {numref}`tab:mc_geo`.

### Dielectric constant of sea water
The dielectric constant of sea water depends on salinity and temperature {cite}`Meissner2004` as 
```{math}
:label: eq:eps_seawater
\begin{aligned}
\varepsilon(T, S)=\frac{\varepsilon_{\mathrm{S}}(T, S)-\varepsilon_1(T, S)}{1+i \nu / \nu_1(T, S)} & +\frac{\varepsilon_1(T, S)-\varepsilon_{\infty}(T, S)}{1+i \nu / \nu_2(T, S)} \\
& +\varepsilon_{\infty}(T, S)-i \frac{\sigma(T, S)}{\left(2 \pi \varepsilon_0\right) \nu}
\end{aligned}
```
with ε$_1$ as the intermediate frequency dielectric constant, ε$_S$ as the static dielectric constant, ε$_\infty$ as the high-frequency dielectric constant, $\nu_1$ and $\nu_2$ as first and second Debye relaxation frequencies, respectively and σ as the conductivity. Fits for S=0 are given by {cite}`Meissner2004` as

```{math}
:label: eq:eps_fresh
\begin{aligned}
\varepsilon_{\mathrm{S}}(T, S=0)&=\frac{3.70886 \cdot 10^4-8.2168 \cdot 10^1
T}{4.21854 \cdot 10^2+T}\\
 \varepsilon_1(T, S=0)&=a_0+a_1 T+a_2 T^2 \\
 \nu_1(T, S=0)&=\frac{45+T}{a_3+a_4 T+a_5 T^2} \\
 \varepsilon_{\infty}(T, S=0)&=a_6+a_7 T \\
 \nu_2(T, S=0)&=\frac{45+T}{a_8+a_9 T+a_{10} T^2} \\
\end{aligned}
```
where the salinity dependence is modeled {cite}`Meissner2004` as
```{math}
:label: eq:eps_sal
\begin{aligned}
\varepsilon_{\mathrm{S}}(T, S) & =\varepsilon_{\mathrm{S}}(T, S=0) \cdot \exp \left[b_0 S+b_1 S^2+b_2 T S\right] \\
\nu_1(T, S) & =\nu_1(T, S=0) \cdot\left[1+S \cdot\left(b_3+b_4 T+b_5 T^2\right)\right] \\
\varepsilon_1(T, S) & =\varepsilon_1(T, S=0) \cdot \exp \left[b_6 S+b_7 S^2+b_8 T S\right] \\
\nu_2(T, S) & =\nu_2(T, S=0) \cdot\left[1+S \cdot\left(b_9+b_{10} T\right)\right] \\
\varepsilon_{\infty}(T, S) & =\varepsilon_{\infty}(T, S=0) \cdot\left[1+S \cdot\left(b_{11}+b_{12} T\right)\right].
\end{aligned}
```

However, {cite}`Meissner2012` published an update to the salinity dependence of ε$_S$ and ν$_1$ as
```{math}
:label: eq:eps_water_sal_update
\begin{aligned}
ε_S(T_S,S) & = ε_S(T_S,S=0) \cdot \exp [b_0S +b_1S^2+b_2T_SS]\\
\nu_1\left(T_S, S\right) & =\nu_1\left(T_S, S=0\right) \cdot  {\left[1+S \cdot\left(d_0+d_1 T_S+d_2 T_S^2+d_3 T_S^3+d_4 t_s^4\right)\right] }.
\end{aligned}
```
We use this most recent update in the forward model.

For conductivity we follow {cite}`Meissner2004` and {cite}`Meissner2012` who used the earlier work by {cite}`Stogryn1995` as
```{math}
:label: eq:sigma
\sigma(T,S) = \sigma(T,S=35) \cdot R_{15}(S)\frac{R_T(S)}{R_{15}(S)},
```

with

```{math}
:label: eq:sigma_support
\begin{aligned}
 \sigma(T, S=35) & = 2.903602+8.607 \cdot 10^{-2} \cdot T+4.738817 \cdot 10^{-4} \cdot T^2\\
& -2.991 \cdot 10^{-6} \cdot T^3+4.3047 \cdot 10^{-9} \cdot T^4 \\
 R_{15}(S) & = S \cdot \frac{37.5109+5.45216 \cdot S+1.4409 \cdot 10^{-2} \cdot S^2}{1004.75+182.283 \cdot S+S^2} \\
 \frac{R_T(S)}{R_{15}(S)} & =1+\frac{\alpha_0(T-15)}{\alpha_1+T} \\
 \alpha_0 & = \frac{6.9431+3.2841 \cdot S-9.9486 \cdot 10^{-2} \cdot S^2}{84.850+69.024 \cdot S+S^2} \\
 \alpha_1 & = 49.843-0.2276 \cdot S+0.198 \cdot 10^{-2} \cdot S^2.
\end{aligned}
```

### Atmosphere
To get back to emissivity in order to account for the atmospheric contribution
to the brightness temperatures at surface level, the brightness temperature is just devided by the ice and ocean surface temperatures
```{math}
:label: eq:ist
ε_p = \frac{T_{b,p}}{T}
```


With $C_{\text{myi}}+C_{\text{fyi}}=\text{SIC}$ and $\text{SIC}\cdot C_\text{myi} =
\text{MYI}$ being part of the state vector {eq}`eqxy`, the state defines
the surface area fraction of all three considered surface types.

The reflectivity of
the surface, $R_{\text{surf}}$, is given by 
```{math}
:label: eq:surfref
R_{\text{surf}} = 1 -
ε_{\text{ow}}C_{\text{ow}} - ε_{\text{fyi}}C_{\text{fyi}} -
ε_{\text{myi}}C_{\text{myi}}.
```

The atmospheric contribution to the brightness temperature is calculated from
the parametrization from {cite}`Wentz2000`. They fitted the downwelling and
upwelling effective temperature using a least-squares fit to the {term}`TWV`.
The fit is given by 
```{math}
:label: eq:atm
\begin{aligned}
T_D &=b_0+b_1V+b_2V^2+b_3V^3+b_4V^4+b_5\zeta(T_s-T_v)\\
T_U &=T_D+b_6+b7V\\
\end{aligned}
```
where $T_v = 273.16+0.8337 V - 3.029\cdot 10^{-5}V^{3.33}$ for V<48 and
$T_v=301.16$ for V>48, $\zeta(x)=1.05x(1-x^2)/1200$ for $|x|<20\ \text{K}$ and
$\zeta(x)=\text{sign}(x)\cdot 14\ \text{K}$ for $|x|>20\ \text{K}$. 

The absorption by oxygen is given by 
```{math}
:label: eq:abs_oxy
A_o = a_{O1}+a_{O2}(T_{D}-270).
```

The vapor absorption is given by
```{math}
:label: eq:abs_vap
A_V = a_{V1}V+a_{V2}V^2.
```

The liquid water absorption is given by
```{math}
:label: eq:abs_liq
A_L = a_{L1}(1-a_{L2}(T_L - 283))L,
```

with $L$ being the liquid water path.

The total atmospheric attenuation is given by a combination of the individual terms from {eq}`eq:abs_oxy`, {eq}`eq:abs_vap`, and {eq}`eq:abs_liq`:

```{math}
:label: eq:abs_tot
\tau = \exp\left(-\frac{A_o+A_V+A_L}{\cos\theta}\right)
```
with $\theta$ being the incidence angle.

An appropriate approximation of the atmospheric integration for up- and downwelling brightness temperature used by {cite}`Wentz2000` is given by
```{math}
:label: eq:tbud
\begin{aligned}
T_{b,u} &= (1-\tau)(T_U)\\
T_{b,d} &= (1-\tau)(T_D)\\
\end{aligned}
```
This makes the upwelling and downwelling brightness temperature differ only slightly by some Kelvin.


### Addition of the L-band forward model
The addition for 1.4&nbsp;GHz was done by {cite}`Scarlat2020` and is based on {cite}`Ruf2003`. 

The atmospheric attenuation for L-band is
```{math}
:label: eq:tau_l
\tau = \exp\left(-\frac{0.009364 + 0.000024127V }{\cos(\theta)}\right),
```
with $V$ being the {term}`TWV` in mm and $\theta$ being the incidence angle.

The up- and down-welling brightness temperature for L-band is given by
```{math}
:label: eq:tbud_l
\begin{aligned}
T_{b,u} &= (1-\tau)(\text{ST}+258.15)\\
T_{b,d} &= (1-\tau)(\text{ST}+263.15)
\end{aligned}
``` 
with ST being the surface temperature in K. Over open ocean, this is
{term}`SST` and over ice it is {term}`IST`. The composite upwelling and
down-welling brightness temperature is then
calculated with the weight of the ice concentration.

The effect of wind roughening for L-band over open ocean is given by
```{math}
:label: eq:rough_L
\begin{aligned}
\epsilon_h &= E_{0,h} + u(0.0007 + 0.000015\theta) \\
\epsilon_v &= E_{0,v} + 0.0007u, 
\end{aligned}
```
with $u$ being the wind speed in m/s, $\theta$ the incidence angle in degrees, and $E_{0,h}$ and $E_{0,v}$ being the emissivity of the surface in the horizontal and vertical polarization from Eq. {eq}`eq:em_ocean`.

### Brightness temperature at instrument level 
Bringing together the different models for all frequencies, the brightness temperature at the instrument level is given by
```{math}
T_{b,p} = T_{b,u} + \left((T_{c} 𝜏 + T_{b,d})(1 − \epsilon_p) + T_{b,s}\right) \tau,
```
with $p$ being the polarization and $T_c$ being the cosmic background temperature. $T_{b,s}$ is the surface emitted brightness temperature for the combined surface and $\epsilon_p$ is the emissivity of the combined surface for the polarization $p$.

(OZA)=
### Correction for variation in incidence angle
While the algorithm is designed for a fixed incidence angle, the incidence
angle is varying slightly depending on the feed horn which is used for the
observation which, in turn, changes the brightness temperatures. We attempt a straightforward correction for the incidence angle dependence of the
brightness temperatures by estimating the surface reflectivity using the
polarization information for a given observation. The assumption here is a
simple Fresnel reflection model without an atmospheric contribution. First, we
define the Fresnel reflection coefficient of the surface with given effective
permittivity $\varepsilon$ and incidence angle $\theta_i$ as described in {eq}`eq:em_ocean`. The effective permittivity is used here which is a simplification calculated with the assumption that the atmosphere is transparent. We use an effective surface temperature $T_{\text{eff}}$ which we set to $T_{\text{eff}}=240$ if $T_{b,v}<220$ K otherwise it is set to $T_{\text{eff}}=T_{b,v}+20$ K as an approximation. Then we have
```{math}
:label: eq:eff_temp
T_{\text{eff}} = \frac{T_{b,h}(θ_i) \cdot (1-r_h(ε,θ_i)) + 
T_{b,v}(θ_i)\cdot(1-r_v(ε,θ_i))}{2}
```
which we can solve for a given $θ_i$ to get the effective permittivity as
```{math}
:label: eq:eff_eps
ε_{\text{eff}} = \operatorname{arg\,min}_{ε \in \mathbb{R}} \left| \left( \frac{T_{b,h}(θ_i)}{T_{\text{eff}}} - (1-r_h(ε,θ_i))\right)^2 +
\left( \frac{T_{b,v}(θ_i)}{T_{\text{eff}}} - (1-r_v(ε,θ_i))\right)^2 \right|
```
The corrected brightness temperature for the incidence angle $θ_{\text{ref}}$ at vertical polarization is then given by
```{math}
:label: eq:incidence
T_{b,p}(θ_{\text{ref}}) =  T_{b,p}(\theta) + \left(r_p(ε_{\text{eff}},θ_{\text{ref}}) - r_p(ε_{\text{eff}},θ)\right) \cdot ζ
```

with $θ$ being the incidence angle of the observation. In theory, we expect the correction factor $ζ$ to be equal to the effective temperature $T_{\text{eff}}$.
Our only source for the incidence angle dependence is a simulated scene from the SCEPS project, referred to as the SCEPS polar scene in this document. In this scene, the small differences in incidence angle between different scans are included and modeled in the brightness temperatures. 
For $T_{b,v}$ this correction works for all CIMR frequencies on the SCEPS polar
scene with a correction factor $ζ=T_{\text{eff}}$ as expected. For $T_{b,h}$ a slightly different correction factor was found in order
to remove the incidence angle dependence of the brightness temperatures in the SCEPS polar scene. The origin of this factor is not
yet understood and is subject to further investigation. In this version of the
algorithm, the correction factor $ζ$ for $T_{b,h}$ is shown in the table {numref}`tab:tbh_corr`
and is applied to the $T_{b,h}$ for each channel.

```{list-table} "Correction factor for $T_{b,h}$ for different channels"
:name: "tab:tbh_corr"
:header-rows: 1

* - Band
  - Correction factor $ζ$
* - L_BAND
  - $T_{\text{eff}}$
* - C_BAND
  - $T_{\text{eff}} - T_{b,h}(θ)$
* - X_BAND
  - $(T_{\text{eff}} - T_{b,h}(θ))/0.9$
* - KU_BAND
  - $(T_{\text{eff}} - T_{b,h}(θ))/1.4$
* - KA_BAND
  - $(T_{\text{eff}} - T_{b,h}(θ))/1.8$
```

This method requires a non-linear optimization for each observation, which can be computationally costly. It is expected to work without modification only on the SCEPS simulated scene as it is probably similar to the method used to create the incidence angle dependence in the scene. A real CIMR scene might require a different approach.

```{note}
The proposed method here is solely to correct the incidence angle dependence, the involved effective permittivity $ε_{\text{eff}}$ and effective temperature $T_{\text{eff}}$ should not be confused with the physical permittivity of the surface.
```


## Input data
The input to the retrieval is using the L1B data product from the CIMR instrument, which includes the brightness temperatures of the channels 1.4, 6.9, 10.7, 18.7, and 36.5&nbsp;GHz and their uncertainties.
Technically, missing values are allowed, in case of malfunctioning channels, but the retrieval
uncertainties will be larger. In addition, for a better constrain of the
solution space, {term}`ECMWF` analysis data are highly recommended as additional input (see {ref}`sec:auxiliary_data` below). 
$\mathbf{y}$ and $\mathbf{S}_e$ from {eq}`eq:chi2` are set with the brightness temperatures and their uncertainties from the L1B data product.


## Output data

The output data will include the retrieved geophysical parameters listed in equation
{eq}`eqxy` and their posterior uncertainties. In addition, quality flags
are derived for each quantity and the retrieval procedure in general.

(sec:auxiliary_data)=
## Auxiliary data
{term}`ECMWF` surface analysis data is used as background values for the retrieval. The
variables used are {term}`WS`, {term}`TWV`, {term}`CLW`, {term}`T2M`,
{term}`TSK`. They are used to fill the $\mathbf{S}_a$ matrix and the
$\mathbf{x}_a$ in equation {eq}`eq:chi2`. For near real-time retrieval, the
$\mathbf{S}_a$ and $\mathbf{x}_a$ can use monthly or seasonal values, as
the retrieval is not sensitive to the exact values of the background variables with the rather wide covariance matrix.
```{note}
In this document a retrieval with fixed background values and TB error is used. This may alter the results of the retrieval on the test cards compared to actual data where the background values are from ECMWF analysis data.
```

## Validation process

For validation of the multi parameter retrieval, the comparison to the
individual specialized products are planned, with the addition of their
validation datasets. For a complete retrieval with currently operational
instruments, however, the time period of the validation data must be aligned
with the operational period of the {term}`AMSR2` and {term}`SMOS` sensors.