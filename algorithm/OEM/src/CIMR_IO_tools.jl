

bands = "L_BAND", "C_BAND", "X_BAND", "KU_BAND", "KA_BAND"
#reads h and v pol data for a given band
function get_data_for_band(fn, band, direction="forward";oza_corr=false,verbose=false)
    ncd = NCDataset(fn)
    d = ncd.group[band]
    flag = d["instrument_status"][:]
    
    ddir = direction == "forward" ? 0 : 1
    idx = flag .== ddir
    lon = vcat([d["lon"][i, :, :][idx] for i in 1:d.dim["n_horns"]]...)
    lat = vcat([d["lat"][i, :, :][idx] for i in 1:d.dim["n_horns"]]...)
    h = vcat([d["brightness_temperature_h"][i, :, :][idx] for i in 1:d.dim["n_horns"]]...)
    v = vcat([d["brightness_temperature_v"][i, :, :][idx] for i in 1:d.dim["n_horns"]]...)
    oza = vcat([d["OZA"][i, :, :][idx] for i in 1:d.dim["n_horns"]]...)
    

    #at the moment nedt is only one value for the whole swath per band, but to be consistent with the rest of the code we will project it as well
    if "nedt_150K_kelvin" in keys(d["brightness_temperature_h"].attrib)
        nedt_h=d["brightness_temperature_h"].attrib["nedt_150K_kelvin"]
        nedt_v=d["brightness_temperature_v"].attrib["nedt_150K_kelvin"]
        herr= fill(nedt_h, size(h))
        verr= fill(nedt_v, size(v))
        verbose && println("nedt found for $band using nedt_150K_kelvin attribute, h=$nedt_h, v=$nedt_v as uncertainty.")
    else
        verbose && println("nedt not found, using 1.0K as uncertainty for h and v.")
        herr = fill(1.0, size(h))
        verr = fill(1.0, size(v))
    end
    if oza_corr
        h,v = oza_correction(h,v,oza,band)
    end

    return (; h, v, herr, verr, lat, lon, oza)
end


"""
    oza_correction(h,v,oza)

Applies the OZA correction to tbh and tbv data
h and v are the brightness temperatures in K
oza is the observation zenith angle in degrees
"""
function oza_correction(h,v,oza,band)
    ttemp=240.0 #K 
    if band=="L_BAND"
        θref=52.0
    else
        θref=55.0
    end

    Tref=max.(ttemp,v.+20) #reference temperature (see ATBD)
    #estimate surface emissivity 

    rv= (v ./Tref)
    rh= (h ./Tref)

    function fcost(eps,rv1,rh1,oza1)
        refh,refv=OEM.fresnel(eps ,oza1)
        (((1.0-refv)-rv1)^2 + ((1.0-refh)-rh1)^2) |> sqrt 
    end

    #optimize the emissivity for each pixel
    eps_opt = [Optim.minimizer(optimize(x->fcost(x,rv[i],rh[i],oza[i]),1.0,1000.0)) for i in eachindex(rv)]

    rh_c,rv_c=OEM.fresnel.(eps_opt ,θref) |>x->(1 .- first.(x), 1 .- last.(x))
    rh_co,rv_co=OEM.fresnel.(eps_opt ,oza) |>x->(1 .- first.(x), 1 .- last.(x))

    #empirical correction factor (based on SCEPS polar scene)
    if band=="L_BAND"
        emp_corr_v= Tref
        emp_corr_h= Tref
    elseif band=="C_BAND"
        emp_corr_v= Tref
        emp_corr_h=(Tref .- h)
    elseif band=="X_BAND"
        emp_corr_v= Tref
        emp_corr_h= (Tref .- h)/0.9
    elseif band=="KU_BAND"
        emp_corr_v= Tref 
        emp_corr_h= (Tref .- h)/1.4
    elseif band=="KA_BAND"
        emp_corr_v= Tref
        emp_corr_h= (Tref .- h)/1.8
    end
    
    v_c = v + (rv_c - rv_co).*emp_corr_v
    h_c = h + (rh_c - rh_co).*emp_corr_h
    return h_c,v_c
end

#projects all data into a common grid of a given band
"""
    project_data(fn, target_band)

Projects all data into a common grid of a given band from a given l1b file
target_band is one of "L_BAND", "C_BAND", "X_BAND", "KU_BAND", "KA_BAND"
The output is a 2d array of size (npoints, nbands*2) where 
npoints is the number of data points in the target band. 
Also returns the lat and lon of the target band
"""
function project_data(fn, target_band, direction="forward";oza_corr=false)
    gm = pr[].geometry #pyresample.geometry
    kdt = pr[].kd_tree #pyresample.kd_tree 
    data = (; (Symbol(b) => get_data_for_band(fn, b, direction,oza_corr=oza_corr) for b in bands)...)
    target_lat = data[Symbol(target_band)].lat
    target_lon = data[Symbol(target_band)].lon
    target_lon[target_lon .> 180] .-= 360
    target_grid = gm.SwathDefinition(lons=target_lon, lats=target_lat)
    outarr = Array{Float64,2}(undef, length(target_lat), length(bands) * 2)
    #at the moment any error is only one value for the whole swath per band, #but to be consistent with the rest of the code we will project it as well
    outerr = Array{Float64,2}(undef, length(target_lat), length(bands) * 2)


    
    i = 0
    for band in bands
        dat = data[Symbol(band)]
        SD = gm.SwathDefinition(lons=dat.lon, lats=dat.lat)
        h = dat.h |> np[].array
        v = dat.v |> np[].array
        herr = dat.herr |> np[].array
        verr = dat.verr |> np[].array
        vproj = kdt.resample_nearest(SD, v, target_grid, 25000,reduce_data=false) |> PyArray |> Array
        hproj = kdt.resample_nearest(SD, h, target_grid, 25000,reduce_data=false) |> PyArray |> Array
        
        vprojerr = kdt.resample_nearest(SD, verr, target_grid, 25000,reduce_data=false) |> PyArray |> Array
        hprojerr = kdt.resample_nearest(SD, herr, target_grid, 25000,reduce_data=false) |> PyArray |> Array
        
        outarr[:, 1+2*i] = vproj
        outarr[:, 1+2*i+1] = hproj
        outerr[:, 1+2*i] = vprojerr
        outerr[:, 1+2*i+1] = hprojerr
        i += 1
    end
    return outarr, outerr, target_lat, target_lon
end


"""
    project_data_to_area_nn(data, lat,lon,area)

projects the data to a given area using nearest neighbour interpolation
data is a 1d array of data values
lat and lon are 1d arrays of latitudes and longitudes
area is a pyresample area definition which includes extent and number of x and y pixels
"""
function project_data_to_area_nn(data, lat, lon, area)
    gm = pr[].geometry #pyresample.geometry
    kdt = pr[].kd_tree #pyresample.kd_tree
    SD = gm.SwathDefinition(lons=lon, lats=lat)
    outarr = kdt.resample_nearest(SD, np[].array(data), area, 25000) |> PyArray |> Array
    return outarr
end


"""
    read_testscene(fn)

reads one of the test scenes and returns the data, the area definition and the extent of the image
"""
function read_testscene(fn)
    ods = NCDataset(fn)
    c = ods["crs"]
    epsg = c.attrib["epsg_code"]
    nrow = ods.dim["n_row"]
    ncol = ods.dim["n_col"]
    xstart, xstep, _, ystart, _, ystep = split(c.attrib["GeoTransform"]) .|> x -> parse(Float64, x)
    mea = pr[].geometry.AreaDefinition("x", "x", "x", epsg, nrow, ncol, [xstart, ystart, xstart + nrow * xstep, ystart + ncol * ystep])
    imextent = [mea.area_extent[0], mea.area_extent[2], mea.area_extent[1], mea.area_extent[3]]
    testscenebands = split("L_band_V L_band_H C_band_V C_band_H X_band_V X_band_H Ku_band_V Ku_band_H Ka_band_V Ka_band_H")
    datastack = cat(((ods[b][:, :]')[:] for b in testscenebands)..., dims=2)
    surfaces = ods["surfaces"][:,:]'
    return (datastack, mea, imextent, surfaces)
end

function read_testscene_sceps(fn)
    #fixing projection for simplicity
    #incidence angle is 55 at index 5
    incidx = 5
    lea = pr[].geometry.AreaDefinition("ease2_nh_testscene", "ease2_nh_testscene", "ease2_nh_testscene", "EPSG:6931", 1400, 1400, [0, -1.5e6, 1.4e6, -1e5])
    ods = NCDataset(fn)
    testbands = split("toa_tbs_L_Vpo toa_tbs_L_Hpo toa_tbs_C_Vpo toa_tbs_C_Hpo toa_tbs_X_Vpo toa_tbs_X_Hpo toa_tbs_Ku_Vpo toa_tbs_Ku_Hpo toa_tbs_Ka_Vpo toa_tbs_Ka_Hpo")
    leaextent = [lea.area_extent[0], lea.area_extent[2], lea.area_extent[1], lea.area_extent[3]]
    datastack = cat(((ods[b][:, :, incidx, 1]')[:] for b in testbands)..., dims=2)
    return (datastack, lea, leaextent)
end

function read_testscene_sceps_geo(fn)
  #  lea = pr.geometry.AreaDefinition("ease2_nh_testscene", "ease2_nh_testscene", "ease2_nh_testscene", "EPSG:6931", 1400, 1400, [0, -1.5e6, 1.4e6, -1e5])
    ods = NCDataset(fn)
    sic = ods["asi_sea_ice_concentration_nh"][:,:,1]
    return sic
end