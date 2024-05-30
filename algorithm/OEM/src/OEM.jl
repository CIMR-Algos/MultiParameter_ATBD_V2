module OEM
#export functions declared in this file
using PythonCall
 
#some workarounds for python imports 
const pr = Ref{Py}()
const np = Ref{Py}()

function __init__()
    np[] = pyimport("numpy")
    pr[] = pyimport("pyresample")
end

nsidcg() = pr[].geometry.AreaDefinition("nsidc_n", "NSIDC polar stereographic north", "nsidc_n", "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 608, 896, [-3850000, -5350000, 3750000, 5850000])


#all imports
using LinearAlgebra
using ForwardDiff
using StaticArrays
using Accessors
using NCDatasets
using NearestNeighbors
using Proj


LLtoEASE2n = Proj.Transformation("epsg:4326", "epsg:6931")
EASE2NtoLL = inv(LLtoEASE2n)
LLtoEASE2 = Proj.Transformation("epsg:4326", "epsg:6933")
EASE2toLL = inv(LLtoEASE2)
LLtoEASE2s = Proj.Transformation("epsg:4326", "epsg:6932")
EASE2StoLL = inv(LLtoEASE2s)


#functions defined in this file
export fw_fnct_amsre
export lm_retrieval
export run_oem

include("ncstore.jl")
export store_output

include("CIMR_IO_tools.jl")
export get_data_for_band 
export project_data
export read_testscene
export read_testscene_sceps
export project_data_to_area_nn

include("IOutils.jl")
export get_amsr2_data_on_NSIDC
export get_amsr2_tbs
export get_smos_data_on_NSIDC
export prepare_oem_tbs
export read_ecmwf
export get_ecmwf_on_NSIDC


function run_oem(band_data::T, band_nedt::Union{Nothing,T}=nothing) where {T<:AbstractMatrix}

    sqrtS_e_fix = SA[5.0, 5.0,  #1.4v,h
        2.356, 4.832,  #6.9v,h
        1.609, 5.460,  #10.6v,h
        0.977, 4.932,  #18.7v,h
        #200+1.042, 200+2.661,  #23.8v,h #no weight in CIMR
        2.540, 2.65, # 36.5v,h
        #200+4.903, 200+6.274] #89v,h #no weight in CIMR
        ]
    S_p = Diagonal(SA[10, #wsp
        10, #twv
        #0.0204, #clw
        0.5, #clw
        10.0, #sst
        10.0, #ist
        0.3, #ic #good guess when initial value from TBs
        0.3, #myif # good guess when initial value from TB
        100.0, #sit
        5 # sss
    ] .^ 2) #sit [cm]

    # startguess&background for geophysical values
    Sg = SA[
        5.0, #wsp
        2.0, #twv
        0.0204, #clw
        275.0, #sst
        270.0, #ist
        0.5, #ic #good guess when initial value from TBs
        0.5, #myif # good guess when initial value from TBs
        100.0, #sit
        30.0 # sss
    ]


    l = size(band_data, 1)
    dat = Array{Float64}(undef, l, 9)
    err = Array{Float64}(undef, l, 9)
    res = Array{Float64}(undef, l, 10)
    its = Array{Float64}(undef, l)

    cimridx = [1:8; 11:12]
    @sync Threads.@threads for i = 1:l #multi-threaded loop comment out if segfaults
    #for i = 1:l #single threaded loop, uncomment if multi-threading causes segfaults
        Tbs = SVector{10,Float64}(band_data[i, j] for j in 1:10)
        sqrtS_ee = isnothing(band_nedt) ? sqrtS_e_fix : SVector{10,Float64}(band_nedt[i, j] for j in 1:10)
        if Tbs[1] < 50
            dat[i, :] .= NaN
            err[i, :] .= NaN
            res[i, :] .= NaN
            continue
        end

        Tbs = SVector{14,Float64}(Tbs[1], Tbs[2], Tbs[3], Tbs[4], Tbs[5], Tbs[6], Tbs[7],
            Tbs[8], 200, 200, Tbs[9], Tbs[10], 200, 200)
        S_e = Diagonal(SVector{14,Float64}(sqrtS_ee[1], sqrtS_ee[2], sqrtS_ee[3], sqrtS_ee[4],
            sqrtS_ee[5], sqrtS_ee[6], sqrtS_ee[7], sqrtS_ee[8],
            200, 200, sqrtS_ee[9], sqrtS_ee[10], 200, 200) .^ 2) #weight of non-CIMR bands is very low

        a, c, b,it = lm_retrieval(Tbs, S_e, S_p, Sg, fw_fnct_amsre)
        c = diag(c)
        dat[i, :] = a
        err[i, :] = sqrt.(c)
        res[i, :] = SVector{10,Float64}(getindex(b', j) for j in cimridx)
        its[i] = it
    end
    return dat, err, res, its
end


function fresnel(eps1, θ)
    #calculates the reflection coefficient of a medium with dielectric constant eps1 to air (eps2=1) under incidence angle θ (in deg)
    ct = cosd(θ)
    ctt = sqrt(eps1 - sind(θ)^2)
    ρᵥ = (eps1 * ct - ctt) / (eps1 * ct + ctt)
    ρₕ = (ct - ctt) / (ct + ctt)
    return (abs2(ρₕ), abs2(ρᵥ))
end


function lm_retrieval(Ta, Sₑ, Sₐ, xₐ, F; kwargs...)
    return lm_retrieval(Ta, Sₑ, Sₐ, xₐ, F, xₐ; kwargs...)
end

"""
    dielectric_meissner_wentz(sst_in,s)

    Dielectric constant of sea water (pure water) as a function of temperature and salinity
"""
function dielectric_meissner_wentz(sst_in::T, s::T, freq) where {T}

    X_pure = SA[5.7230E+00, 2.2379E-02, -7.1237E-04, 5.0478E+00, -7.0315E-02, 6.0059E-04, 3.6143E+00, 2.8841E-02, 1.3652E-01, 1.4825E-03, 2.4166E-04]
    X_sea = SA[-3.56417E-03, 4.74868E-06, 1.15574E-05, 2.39357E-03, -3.13530E-05, 2.52477E-07, -6.28908E-03, 1.76032E-04, -9.22144E-05, -1.99723E-02, 1.81176E-04, -2.04265E-03, 1.57883E-04]
    # fit salinity dependence of the static dielectric constant
    a0coef = SA[-0.33330E-02, 4.74868E-06, 0.0E+00] #wentz 2012 Table VI
    #fit salinity  dependence of the first Debye relaxation frequency
    b1coef = SA[0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08]#Wentz 2012 Table VII

    sst = sst_in - 273.15 #to have sst in degree celsius
    if sst < -30.16
        sst = -30.16  #protects against n1 and n2 going to zero for very cold water
    end



    f0 = 17.97510

    ###     pure water, Debye relaxation fits, fit parameter for Eq. 6 (MW 2004)
    e0 = (3.70886e4 - 8.2168e1 * sst) / (4.21854e2 + sst) # stogryn et al.
    e1 = X_pure[1] + X_pure[2] * sst + X_pure[3] * sst^2
    n1 = (45.00 + sst) / (X_pure[4] + X_pure[5] * sst + X_pure[6] * sst^2)
    e2 = X_pure[7] + X_pure[8] * sst
    n2 = (45.00 + sst) / (X_pure[9] + X_pure[10] * sst + X_pure[11] * sst^2)

    ###     saline water
    #     conductivity [s/m] taken from stogryn et al.
    sig35 = 2.903602 + 8.60700E-2 * sst + 4.738817E-4 * sst^2 - 2.9910E-6 * sst^3 + 4.3047e-9 * sst^4
    r15 = s * (37.5109 + 5.45216 * s + 1.4409E-2 * s^2) / (1004.75 + 182.283 * s + s^2)
    alpha0 = (6.9431 + 3.2841 * s - 9.9486E-2 * s^2) / (84.850 + 69.024 * s + s^2)
    alpha1 = 49.843 - 0.2276 * s + 0.198E-2 * s^2
    rtr15 = 1.0 + (sst - 15.0) * alpha0 / (alpha1 + sst)

    sig = sig35 * r15 * rtr15

    #   permittivity
    a0 = exp(a0coef[1] * s + a0coef[2] * s^2 + a0coef[3] * s * sst) # static dielectric constant, Eq 29 MW 2012
    e0s = a0 * e0

    if sst < 30
        #b1 = 1.0 + s*(b1coef[1] + b1coef[2]*sst + b1coef[3]*sst^2 + b1coef[4]*sst^3 + b1coef[5]*sst^4)
        b1 = 1.0 + s * (b1coef[1] + sst * (b1coef[2] + sst * (b1coef[3] + sst * (b1coef[4] + sst * b1coef[5]))))
    else
        b1 = 1.0 + s * (9.1873715E-04 + 1.5012396E-04 * (sst - 30))
    end

    n1s = n1 * b1

    a1 = exp(X_sea[7] * s + X_sea[8] * s^2 + X_sea[9] * s * sst)
    e1s = e1 * a1

    b2 = 1.0 + s * (X_sea[10] + 0.5 * X_sea[11] * (sst + 30))
    n2s = n2 * b2

    a2 = 1.0 + s * (X_sea[12] + X_sea[13] * sst)
    e2s = e2 * a2

    permit = (e0s - e1s) / (1.0 - im * (freq / n1s)) + (e1s - e2s) / (1.0 - im * (freq / n2s)) + e2s + im * sig * f0 / freq #new: equation 6 from MW 2004
    permit = conj(permit) #to get the sign in the convention that the imaginary part is negative
    return real(permit), imag(permit)
end



function lm_retrieval(Ta, Sₑ, Sₐ, xₐ, F::T, x₀; verbose=false, debug=false) where {T<:Function}
    #lm method after Rodgers 2000
    #target: find x so that F(x)=Ta, given
    #Ta: measurement vector
    #Sₑ: error covariance of measurement
    #Sₐ: error covariance of physical state 
    #xₐ: expected physical state (also used as start, i.e. first guess)
    #F: the forward model translating measument space into state space
    # keywords:
    #   verbose, if true prints cost, d2 and iterations
    #   debug, is more verbose than verbose and prints these at each iterations, including γ
    # returns 
    #    x̂: retrieved parameters
    #    Ŝ: covariance after retrieval
    #    res: residuals, i.e. (Ta - F(x̂)) 
    Sₐ⁻¹ = inv(Sₐ)
    Sₑ⁻¹ = inv(Sₑ)
    #function to minimize with changing input x
    J(y, x, Sₑ⁻¹, Sₐ⁻¹, xₐ, Fx) = (y .- Fx)' * (Sₑ⁻¹ * (y .- Fx)) + (xₐ .- x)' * (Sₐ⁻¹ * (xₐ .- x)) #first two terms of eq 5.3
    xᵢ = x₀

    #inferene does not work well with DiffResults, so we are enforcing the type
    result::ForwardDiff.DiffResults.ImmutableDiffResult{1,SVector{14,Float64},Tuple{SMatrix{14,9,Float64,126}}} = ForwardDiff.DiffResults.JacobianResult(Ta, xₐ)

    γ = 1e-5 #set to 0 for gauss newton
    it = 0

    result = ForwardDiff.jacobian!(result, F, xᵢ)

    Kᵢ = result.derivs[1]
    Fxᵢ = result.value
    Jᵢ = J(Ta, xᵢ, Sₑ⁻¹, Sₐ⁻¹, xₐ, Fxᵢ)
    while true
        it += 1

        Ŝ⁻¹ = Sₐ⁻¹ + Kᵢ' * Sₑ⁻¹ * Kᵢ #eq 5.13
        xᵢ₊₁ = xᵢ + ((1 + γ) * Sₐ⁻¹ + Kᵢ' * Sₑ⁻¹ * Kᵢ) \ (Kᵢ' * Sₑ⁻¹ * (Ta - Fxᵢ) - Sₐ⁻¹ * (xᵢ - xₐ)) #eq 5.36

        result = ForwardDiff.jacobian!(result, F, xᵢ₊₁)
        Fxᵢ₊₁ = result.value
        Jᵢ₊₁ = J(Ta, xᵢ₊₁, Sₑ⁻¹, Sₐ⁻¹, xₐ, Fxᵢ₊₁)
        d² = (xᵢ - xᵢ₊₁)' * Ŝ⁻¹ * (xᵢ - xᵢ₊₁) #eq 5.29
        if (Jᵢ₊₁ > Jᵢ) #rejection criterion from our local OEM
            γ *= 10
        else #step accepted, reduce step size
            γ /= 10
            xᵢ = xᵢ₊₁
            Fxᵢ = Fxᵢ₊₁
            Jᵢ = Jᵢ₊₁
            Kᵢ = result.derivs[1]

        end

        debug && @show it, Jᵢ, d², γ
        if d² < length(xₐ) || it > 50
            verbose && @show it, Jᵢ, d²
            break
        end
    end
    Ŝ = inv(Sₐ⁻¹ + Kᵢ' * Sₑ⁻¹ * Kᵢ) # eq 5.38

    return xᵢ, Ŝ, Ta - Fxᵢ, it
end


function fw_fnct_amsre(x::SVector{9,T}) where {T} ##to use ForwardDiff the function must take an array
    Mc_atm = SA[239.50E+0 213.92E-2 -460.60E-4 457.11E-6 -16.84E-7 0.50E+0 -0.11E+0 -0.21E-2 8.34E-3 -0.48E-4 0.07E-3 0.00E-5
        239.51E+0 225.19E-2 -446.86E-4 391.82E-6 -12.20E-7 0.54E+0 -0.12E+0 -0.34E-2 9.08E-3 -0.47E-4 0.18E-3 0.00E-5
        240.24E+0 298.88E-2 -725.93E-4 814.50E-6 -36.07E-7 0.61E+0 -0.16E+0 -1.69E-2 12.15E-3 -0.61E-4 1.73E-3 -0.05E-5
        241.69E+0 310.32E-2 -814.29E-4 998.93E-6 -48.37E-7 0.20E+0 -0.20E+0 -5.21E-2 15.75E-3 -0.87E-4 5.14E-3 0.19E-5
        239.45E+0 254.41E-2 -512.84E-4 452.02E-6 -14.36E-7 0.58E-0 -0.57E+0 -2.38E-2 40.06E-3 -2.00E-4 1.88E-3 0.09E-5
        242.58E+0 302.33E-2 -749.76E-4 880.66E-6 -40.88E-7 0.62E+0 -0.57E+0 -8.07E-2 53.35E-3 -1.18E-4 8.78E-3 0.80E-5]

    Mc_abs = SA[0.0078 0.0303 0.0007 0.0 1.2216
        0.0183 0.0298 0.0027 0.0060 1.1795
        0.0556 0.0288 0.0113 0.0040 1.0636
        0.0891 0.0281 0.0188 0.0020 1.0220
        0.2027 0.0261 0.0425 -0.0020 0.9546
        0.9693 0.0146 0.1506 -0.0020 0.7961]

    Mc_geo = SA[-0.00027 0.00054 -0.000021 0.000032 -0.000021 -0.00002526 0.0 0.0
        -0.00032 0.00072 -0.000029 0.000044 -0.000021 -0.00002894 0.00000008 -0.00000002
        -0.00049 0.00113 -0.000053 0.000070 -0.000021 -0.00003690 0.00000031 -0.00000012
        -0.00063 0.00139 -0.000070 0.000085 -0.000021 -0.00004195 0.00000041 -0.00000020
        -0.00101 0.00191 -0.000105 0.000112 -0.000021 -0.00005451 0.00000045 -0.00000036
        -0.00153 0.00202 -0.000116 0.000130 -0.000021 -0.0000550 -0.0000009 -0.00000046]

    Mc_m = SA[0.00020 0.00200 0.00690 0.00600
        0.00020 0.00200 0.00690 0.00600
        0.00140 0.00293 0.00736 0.00656
        0.00178 0.00308 0.00730 0.00660
        0.00257 0.00329 0.00701 0.00660
        0.00260 0.00330 0.00700 0.00660]


    Emissivitet_FY_V = SA[0.958, 0.960, 0.965, 0.960, 0.946, 0.615]
    Emissivitet_FY_H = SA[0.868, 0.879, 0.887, 0.882, 0.864, 0.488]
    Emissivitet_MY_V = SA[0.972, 0.948, 0.885, 0.839, 0.731, 0.666]
    Emissivitet_MY_H = SA[0.866, 0.845, 0.799, 0.763, 0.675, 0.630]



    Freq = SA[6.93, 10.65, 18.70, 23.80, 36.50, 89.0]


    theta_d = 55.0
    theta_r = 55.0 * pi / 180.0
    T_C = 2.7

    # coefficients from Mathews et al 2008, regressions for emission layer temperature calculation (for months DJFM)

    a_fyi = SA[0.23 0.26 0.29 0.29 0.30 0.37]
    b_fyi = SA[-5.5 -5.2 -5.0 -4.9 -4.9 -4.2]
    a_myi = SA[0.27 0.34 0.42 0.43 0.45 0.49]
    b_myi = SA[-11.5 -10.5 -9.5 -9.2 -8.9 -8.4]


    #allocating output vector with correct type for Dual numbers to be filled later
    T_B = @SVector zeros(T, 14)

    #T_B=Array{Float64,1}(undef,12)
    varW, varV, varL, T_ow, T_is, C_is1, F_MY1, SI_T, s = x

    F_MY = F_MY1
    C_is = C_is1

    #     else
    W = varW
    V = varV
    L = varL
    #     end if

    if (F_MY < 0)
        F_MY = 0.0
    end

    if (F_MY > 1)
        F_MY = 1.0
    end

    IC_MY = C_is1 * F_MY
    IC_FY = C_is1 - IC_MY
    C_ow = 1 - C_is


    ############ Nizy's emission layer temperature calculation starts here

    T_2m_C = T_is - 273.15



    #eq_a_v = SA[250.170, 250.1, 255.06, 255.34, 254.338, 253.481] #unused due to thin ice 
    eq_b_v = SA[157.9397756705356, 163.4896743050617, 175.18191921420814, 185.91398347505353, 206.667948239106, 242.01900700540162]
    eq_c_v = SA[8.956615582385908,
        8.52431046569442,
        7.734378368070544,
        7.473947263410875,
        7.6681838621416425,
        2.8041750136460695]


    #eq_a_h = SA[220.703, 228.617, 238.249, 240.193, 242.620, 241.396] #unused due to thin ice
    eq_b_h = SA[74.22054192317813, 78.40477088765373, 90.60124656474672, 99.89685544430296, 128.36037504633393, 175.158277073163]
    eq_c_h = SA[11.894364589753874,
        11.645437092902535,
        10.165054195271782,
        9.128550147675183,
        8.98643673857489,
        5.686319297256588]

    #some constants needed lin the loop

    T_owC = T_ow - 273.15                 # Surface temperature T_ow [deg Celcius]



    #a = real((1im)^(1-ny))
    #b = imag((1im)^(1-ny))


    for i in eachindex(Freq) #loop over frequencies except L-band
        # ! -----------------------------------------------------------------
        # !!!!!!!!!! First year ice !!! Vertical polarisation
        # ! -----------------------------------------------------------------
        UL_FYI_H = (a_fyi[i] * T_2m_C + b_fyi[i] + 273.15) * Emissivitet_FY_H[i]
        UL_FYI_V = (a_fyi[i] * T_2m_C + b_fyi[i] + 273.15) * Emissivitet_FY_V[i]

        T_BV_FY = UL_FYI_V - (UL_FYI_V - eq_b_v[i]) * exp(-1 * SI_T / eq_c_v[i])
        T_BH_FY = UL_FYI_H - (UL_FYI_H - eq_b_h[i]) * exp(-1 * SI_T / eq_c_h[i])

        thinC_emi_v = T_BV_FY / T_is #effective emissivity including Nizys method (for ice temperauter calculation)
        thinC_emi_h = T_BH_FY / T_is


        #if (SI_T < 51)
        #T_BV_FY = eq_a_v[i]-(eq_a_v[i]-eq_b_v[i])*exp(-1*SI_T/eq_c_v[i])
        #T_BH_FY = eq_a_h[i]-(eq_a_h[i]-eq_b_h[i])*exp(-1*SI_T/eq_c_h[i])
        #else
        #T_BV_FY = (a_fyi[i]*T_2m_C +b_fyi[i] + 273.15)*Emissivitet_FY_V[i] 
        #T_BH_FY = (a_fyi[i]*T_2m_C +b_fyi[i] + 273.15)*Emissivitet_FY_H[i]
        #end

        # ! =================================================================
        # !!!!!!!!!! Multi year ice !!! Vertical polarisation

        T_BV_MY = (a_myi[i] * T_2m_C + b_myi[i] + 273.15) * Emissivitet_MY_V[i]

        # ! -----------------------------------------------------------------
        # !!!!!!!!!! Multi year ice !!! Horizontal polarisation
        # ! -----------------------------------------------------------------

        T_BH_MY = (a_myi[i] * T_2m_C + b_myi[i] + 273.15) * Emissivitet_MY_H[i]

        T_S_mix = C_is1 * T_is + C_ow * T_ow    #Temperature of a mixed surface (ice and open water)



        T_L = (T_S_mix + 273) / 2  #Temperature of the water droplets [Kelvin]; approximation: mean temp of surface and freezing level 

        ##########################################################3
        #Model for the Atmosphere

        if V <= 48 && V > 0
            # T_V = real(273.16+0.8337*V-3.029E-5*Complex(V)^3.33) #equation (27a) #protect against errors thrown by negative V
            T_V = 273.16 + 0.8337 * V - 3.029E-5 * V^3.33 #equation (27a)
        elseif V <= 0
            T_V = 273.16 + 0.8337 * V #deal with negative cases
        else
            T_V = 301.16  #equation (27b)
        end

        temp_test = abs(T_S_mix - T_V)

        if (temp_test <= 20)
            sig_TS_TV = 1.05 * (T_S_mix - T_V) * (1 - ((T_S_mix - T_V)^2) / 1200)  #equation (27c)
        else
            sig_TS_TV = sign(T_S_mix - T_V) * 14.0
        end

        if (V <= 58.0)
            T_D = Mc_atm[i, 1] + Mc_atm[i, 2] * V + Mc_atm[i, 3] * V^2 + Mc_atm[i, 4] * V^3 + Mc_atm[i, 5] * V^4 + Mc_atm[i, 6] * sig_TS_TV     #equation (26a)
        else
            T_D_54 = Mc_atm[i, 1] + Mc_atm[i, 2] * 54 + Mc_atm[i, 3] * 54^2 + Mc_atm[i, 4] * 54^3 + Mc_atm[i, 5] * 54^4 + Mc_atm[i, 6] * sig_TS_TV
            T_D_58 = Mc_atm[i, 1] + Mc_atm[i, 2] * 58 + Mc_atm[i, 3] * 58^2 + Mc_atm[i, 4] * 58^3 + Mc_atm[i, 5] * 58^4 + Mc_atm[i, 6] * sig_TS_TV
            T_D = T_D_58 + (V - 58) * (T_D_58 - T_D_54) / 4
        end

        T_U = T_D + Mc_atm[i, 7] + Mc_atm[i, 8] * V     #equation (26b)

        A_0 = Mc_atm[i, 9] + Mc_atm[i, 10] * (T_D - 270)        #equation (28)

        A_V = Mc_atm[i, 11] * V + Mc_atm[i, 12] * V^2           #equation (29)

        A_L = Mc_abs[i, 1] * (1 - Mc_abs[i, 2] * (T_L - 283)) * L   #equation (33) (no rain)

        tau = exp((-1 / cos(theta_r)) * (A_0 + A_V + A_L))    #Total transmittance; equation (22)
        T_BU = T_U * (1 - tau)   #The upwelling effective air temperature; equation (24a)
        T_BD = T_D * (1 - tau)   #The downwelling effective air temperature; equation (24b)

        #################################
        #Dielectric Constant of Sea-water

        #lambd = (light_speed / (Freq[i] * 1.0E9))    # Wavelength, [cm]

        ϵ1, ϵ2 = dielectric_meissner_wentz(T_ow, s, Freq[i])
        R_0H, R_0V = fresnel(Complex(ϵ1, ϵ2), theta_d)
        R_0V += (4.887E-8 - 6.108E-8 * (T_owC)^3)


        ###############################^
        # The wind-roughened Sea Surface

        R_geoH = R_0H - (Mc_geo[i, 2] + Mc_geo[i, 4] * (theta_d - 53) + Mc_geo[i, 6] * (T_ow - 288) + Mc_geo[i, 8] * (theta_d - 53) * (T_ow - 288)) * W #equation(57); horizontal pol.
        R_geoV = R_0V - (Mc_geo[i, 1] + Mc_geo[i, 3] * (theta_d - 53) + Mc_geo[i, 5] * (T_ow - 288) + Mc_geo[i, 7] * (theta_d - 53) * (T_ow - 288)) * W #equation(57); vertical pol.

        #	@show size(R_geoV)


        #horizontal polarisation
        W_1 = 7
        W_2 = 12

        if (W < W_1)
            F_H = Mc_m[i, 2] * W

        elseif ((W >= W_1) & (W <= W_2))
            F_H = Mc_m[i, 2] * W + 0.5 * (Mc_m[i, 4] - Mc_m[i, 2]) * (W - W_1)^2 / (W_2 - W_1)   #equation (60b)
        else
            F_H = Mc_m[i, 4] * W - 0.5 * (Mc_m[i, 4] - Mc_m[i, 2]) * (W_2 + W_1) #equation (60c)
        end

        W_1 = 3
        W_2 = 12

        if (W < W_1)
            F_V = Mc_m[i, 1] * W   #equation (60a)
        elseif ((W >= W_1) & (W <= W_2))
            F_V = Mc_m[i, 1] * W + 0.5 * (Mc_m[i, 3] - Mc_m[i, 1]) * (W - W_1)^2 / (W_2 - W_1)   #equation (60b)
        else
            F_V = Mc_m[i, 3] * W - 0.5 * (Mc_m[i, 3] - Mc_m[i, 1]) * (W_2 + W_1) #equation (60c)
        end


        R_H = (1 - F_H) * R_geoH   #equation (49) Composite reflektivity for open water; horisontal pol
        R_V = (1 - F_V) * R_geoV   #equation (49) Composite reflektivity for open water; vertical pol

        E_H = 1 - R_H       # surface emissivity for open water; Horizontal pol.; Equation (8); 
        E_V = 1 - R_V     # Surface emissivity for open water; Vertical pol; Equation (8); 

        #probably we have to reavaluate ice contribution after water, 
        #since lower tie point of thin ice should be open water tiepoint to avoid artifacts


        #     #Emissivity for mixed surface
        E_eff_H = C_ow * E_H + IC_FY * thinC_emi_h + IC_MY * Emissivitet_MY_H[i]
        E_eff_V = C_ow * E_V + IC_FY * thinC_emi_v + IC_MY * Emissivitet_MY_V[i]

        #     #Emissivity for mixed ice surface
        E_ice_H = (IC_FY * thinC_emi_h + IC_MY * Emissivitet_MY_H[i]) / (IC_FY + IC_MY)
        E_ice_V = (IC_FY * thinC_emi_v + IC_MY * Emissivitet_MY_V[i]) / (IC_FY + IC_MY)


        R_eff_H = 1 - E_eff_H   #Reflection coefficient for mixed surface
        R_eff_V = 1 - E_eff_V


        if i > 4 # 37 GHz and larger, equation (56a)
            Delta_S2 = 5.22E-3 * W
        else
            Delta_S2 = 5.22E-3 * (1 - 0.00748 * (37 - Freq[i])^1.3) * W # lower than 37 GHz, equ (56b)
        end
        #the following is from the text after equ (62b)
        if Delta_S2 > 0.069
            Delta_S2 = 0.069
        end
        term = Delta_S2 - 70 * Delta_S2^3 #Term for equation (62a+b)

        OmegaH = (6.2 - 0.001 * (37 - Freq[i])^2) * term * tau^2.0  #equation (62a); horizontal pol.
        OmegaV = (2.5 + 0.018 * (37 - Freq[i])) * term * tau^3.4     #equation (62b); vertical pol.
        T_BOmegaH = ((1 + OmegaH) * (1 - tau) * (T_D - T_C) + T_C) * R_eff_H    #equation (61) horizontal pol.
        T_BOmegaV = ((1 + OmegaV) * (1 - tau) * (T_D - T_C) + T_C) * R_eff_V    #equation (61) vertical pol.

        #     ##########################
        #     ## Output
        #     ##########################
        T_BH_ow = E_H * T_ow    #Horizontal brightnesstemperature from open water
        T_BV_ow = E_V * T_ow    #Vertikal brightnesstemperature from opent water

        T_BH_overflade = C_ow * T_BH_ow + IC_FY * T_BH_FY + IC_MY * T_BH_MY  #Horizontal brigthnesstemperature from mixed surface
        T_BV_overflade = C_ow * T_BV_ow + IC_FY * T_BV_FY + IC_MY * T_BV_MY  #Vertikal brigthnesstemperature from mixed surface

        T_BH1 = tau * T_BH_overflade
        T_BH2 = T_BU
        T_BH3 = C_is * tau * (T_BD + tau * T_C) * (1 - E_ice_H)
        T_BH4 = (1 - C_is) * T_BOmegaH * tau

        T_BV1 = tau * T_BV_overflade + T_BU
        T_BV2 = C_is * tau * (T_BD + tau * T_C) * (1 - E_ice_V)
        T_BV3 = (1 - C_is) * T_BOmegaV * tau

        T_BH = T_BH1 + T_BH2 + T_BH3 + T_BH4
        T_BV = T_BV1 + T_BV2 + T_BV3

        T_B = @set T_B[i*2+2] = T_BH
        T_B = @set T_B[i*2-1+2] = T_BV
    end


    #now deal with L-band
    #lambd_L = (light_speed / (1.413 * 1E9))
    Ltheta = 53.0

    #ugly repetition since it is calculated in the loop for the other frequency,
    #this is worth putting in a function to avoid  doubling
    #ocean properties first 

    ϵ1, ϵ2 = dielectric_meissner_wentz(T_ow, s, 1.413)
    R_0H_L, R_0V_L = fresnel(Complex(ϵ1, ϵ2), Ltheta)
    #equation (46)+correction vertical pol. has to be figured out if applicable at L-band
    R_0V_L += (4.887E-8 - 6.108E-8 * (T_owC)^3)


    E_H_L = 1 - R_0H_L#! surface emissivity for open water; L band horizontal pol.; Equation (8); 
    E_V_L = 1 - R_0V_L#! surface emissivity for open water; L band vertical pol.; Equation (8);


    Vcm = V * 0.1 #! Transfrom from mm (V) to cm for use in L band tau

    #now E_H_L and E_V_L are the emissivity for open water
    #correcting now for ocean roughening
    Lemi_v = W * 0.0007 + E_V_L   #!wind influenced ocean surface emissivity
    Lemi_h = W * (0.0007 + 0.000015 * Ltheta) + E_H_L

    LT_BV_OW = Lemi_v * T_ow #!ocean emitted 
    LT_BH_OW = Lemi_h * T_ow #!ocean emitted 

    # now deal with ice
    #temperature calculations like Mathews et al. 2008, with a=0.1 and b=0.0 
    a = 0.1
    T_sk = 273.15 - (273.15 - T_is) * a

    L_UL_FYI_H = 0.86 * T_sk #original
    L_UL_FYI_V = 0.92 * T_sk

    #LT_BH_FYI = L_UL_FYI_H + (72.1293 - L_UL_FYI_H) * exp(-1 * m_SIT / 0.182969)
    #LT_BV_FYI = L_UL_FYI_V + (147.361 - L_UL_FYI_V) * exp(-1 * m_SIT / 0.102263)
    #adapted from SIT ATBD for CIMR
    LT_BH_FYI = L_UL_FYI_H + (75.524 - L_UL_FYI_H) * exp(-1 * SI_T / 21.021)
    LT_BV_FYI = L_UL_FYI_V + (145.170 - L_UL_FYI_V) * exp(-1 * SI_T / 12.509)

    #replaced with higher, more relaistic values
    #(original values are in Scarlat 2020)
    LT_BV_MYI = 0.94 * T_sk#!v pol surface component over MYI
    LT_BH_MYI = 0.85 * T_sk#!h pol surface component over MYI

    #effective mixed emissivities for ice surface
    LE_ice_V = (IC_FY * LT_BV_FYI / T_sk + IC_MY * LT_BV_MYI / T_sk) / (IC_FY + IC_MY)#! effective emissivity of a mixed SI types scene, v pol
    LE_ice_H = (IC_FY * LT_BH_FYI / T_sk + IC_MY * LT_BH_MYI / T_sk) / (IC_FY + IC_MY)#! effective emissivity of a mixed SI types scene, h pol


    Vcm = V * 0.1 #! Transfrom from mm (V) to cm for use in L band tau

    Ltau = exp(-1 * (0.009364 + 0.000024127 * Vcm) * (1 / cosd(Ltheta))) #!atmospheric attenuation factor

    T_U = (T_ow - 273.16) + 258.15
    T_D = (T_ow - 273.16) + 263.15

    LT_U = (1 - Ltau) * T_U
    LT_D = (1 - Ltau) * T_D         #!downwelling component

    LT_BV_m_surf = C_ow * LT_BV_OW + IC_FY * LT_BV_FYI + IC_MY * LT_BV_MYI
    LT_BH_m_surf = C_ow * LT_BH_OW + IC_FY * LT_BH_FYI + IC_MY * LT_BH_MYI

    Tice_U = ((T_is - 273.16) + 258.15)           #!upwelling component over sea ice
    Tice_D = ((T_is - 273.16) + 263.15)           #!downwelling component over sea ice
    LTice_U = (1 - Ltau) * Tice_U
    LTice_D = (1 - Ltau) * Tice_D

    LT_U_tot = C_is * LTice_U + (1.0 - C_is) * LT_U#!combined upwelling based on surf temp of composite surf		
    #LT_D_tot = C_is * LTice_D + (1.0 - C_is) * LT_D#!combined downwelling based on surf temp of composite surf		

    #! LT_BV_SI = ((6*Ltau + LTice_D)*(1-LE_ice_V) + LE_ice_V*T_is)*Ltau #!sea ice emitted + reflected component v-pol
    #! LT_BH_SI = ((6*Ltau + LTice_D)*(1-LE_ice_H) + LE_ice_H*T_is)*Ltau #!sea ice emitted + reflected component h-pol

    #! LT_BT_BV_SI = ((6*Ltau + LTice_D)*(1-LE_ice_V) + LE_ice_V*T_is)*Ltau #!sea ice emitted + reflected component v-pol
    #! LT_BH_SI = ((6*Ltau + LTice_D)*(1-LE_ice_H) + LE_ice_H*T_is)*Ltau #!sea ice emitted + reflected component h-pol


    #! T_BV_overflade=C_ow*T_BV_ow+IC_FY*T_BV_FY+IC_MY*T_BV_MY


    #! some clarifications:
    # 'LE_ice_p' with p∈[H,V] is the emissivity of the ice surface
    # 'LT_Bp_m_surf' with p∈[H,V] is the surface brightness temperature of the mixed surface
    # 

    #L_TC = 6
    L_TC = 2.7

    #! T_BV1=tau*T_BV_overflade+T_BU;			#! surface component + atmo_up at TOA
    LT_BV1 = Ltau * LT_BV_m_surf + LT_U_tot
    LT_BV2 = Ltau * (LTice_D + Ltau * L_TC) * (1 - LE_ice_V) * C_is #ice component
    LT_BV3 = Ltau * (LT_D + Ltau * L_TC) * (1 - Lemi_v) * C_ow#! ocean component


    LT_BH1 = Ltau * LT_BH_m_surf + LT_U_tot
    LT_BH2 = Ltau * (LTice_D + Ltau * L_TC) * (1 - LE_ice_H) * C_is
    LT_BH3 = Ltau * (LT_D + Ltau * L_TC) * (1 - Lemi_h) * C_ow

    #! no ocean scattered component for L-band
    LT_BV = LT_BV1 + LT_BV2 + LT_BV3
    LT_BH = LT_BH1 + LT_BH2 + LT_BH3

    T_B = @set T_B[1] = LT_BV
    T_B = @set T_B[2] = LT_BH

    return T_B
end

end #module
