using Dates
using NCDatasets
"""
    store_output(outfn, alldata, grid)

Store the output of the algorithm in a netcdf file
    with groups for forward, backward and joined data
outfn is the name of the output file
alldata is a namedtuple containing forward and backward data,
    each consisting of data, err, lat, lon
    data is a 2d array of size (npoints, 9) containing the data
    err is a 2d array of size (npoints, 9) containing the error
    lat is a 1d array containing the latitude
    lon is a 1d array containing the longitude
    its is a 1d array containing the number of iterations #not used yet
grid is a string containing the name of the grid used for the data
"""
function store_output(outfn, alldata, grid)
    
    @assert grid in ["native", "ease2_nh", "ease2_sh", "testgrid"]

    if grid == "ease2_nh"
        out_geo=pr[].geometry.AreaDefinition("ease2_nh", "ease2_nh", "ease2_nh", "EPSG:6931", 1440,1440, [-9.0e6, -9.0e6, 9.0e6, 9.0e6])
    elseif grid == "ease2_sh"
        out_geo=pr[].geometry.AreaDefinition("ease2_sh", "ease2_sh", "ease2_sh", "EPSG:6932", 1440, 1440, [-9.0e6, -9.0e6, 9.0e6, 9.0e6])
    elseif grid == "testgrid"
        out_geo=pr[].geometry.AreaDefinition("ease2_nh_testscene","ease2_nh_testscene","ease2_nh_testscene","EPSG:6931",1400,1400,[0,-1.5e6,1.4e6,-0.1e6]) #1400x1400 pixels, 1km resolution
    end
    


    standard_names = ["wind_speed", "total_water_vapor", "cloud_liquid_water", "sea_surface_temperature", "ice_surface_temperature", "sea_ice_concentration", "multi_year_ice_fraction", "first_year_ice_thickness", "sea_surface_salinity"]
    units = ["m s^-1", "kg m^-2", "kg m^-2", "K", "K", "1", "1", "m", "ppt"]
    longnames = ["Wind Speed", "Total Water Vapor", "Cloud Liquid Water", "Sea Surface Temperature", "Ice Surface Temperature", "Sea Ice Concentration", "Multi Year Ice Fraction", "First Year Ice Thickness", "Sea Surface Salinity"]

    # create the output file
    NCDataset(outfn, "c") do ncid
        ds = (;forward = defGroup(ncid,"forward"),
              backward = defGroup(ncid,"backward"),
              combined = defGroup(ncid,"combined"))
    
        # create the dimensions and variables , native grid is 1d, ease2 grids are 2d
        
        if grid == "native"
            defDim(ncid, "npoints", length(lat))
            latid = defVar(ncid,"lat", lat, ("npoints",),deflatelevel=5, shuffle=true, chunksizes=(length(lat),))

            lonid = defVar(ncid,"lon", lon, ("npoints",),deflatelevel=5, shuffle=true, chunksizes=(length(lon),))
        else
            xsize=pyconvert(Int,out_geo.x_size)
            ysize=pyconvert(Int,out_geo.y_size)
            defDim(ncid, "x", xsize)
            defDim(ncid, "y", ysize)
            lons,lats=out_geo.get_lonlats()
            x,y=out_geo.get_proj_vectors()
            @show xsize,ysize,length(x),length(y)
            defVar(ncid, "x" , pyconvert(Array{Int},x), ("x",),deflatelevel=5, shuffle=true, chunksizes=(xsize,), attrib=Dict("units"=>"m", "standard_name"=>"projection_x_coordinate"))
            defVar(ncid, "y", pyconvert(Array{Int},y), ("y",),deflatelevel=5, shuffle=true, chunksizes=(ysize,), attrib=Dict("units"=>"m", "standard_name"=>"projection_y_coordinate`"))
            latid = defVar(ncid,"lat", pyconvert(Array, lats)', ("x", "y"),deflatelevel=5, shuffle=true, chunksizes=(xsize,ysize), attrib=Dict("units"=>"degrees_north", "standard_name"=>"latitude","grid_mapping"=>"crs"))
            lonid = defVar(ncid,"lon", pyconvert(Array, lons)', ("x", "y"),deflatelevel=5, shuffle=true, chunksizes=(xsize,ysize), attrib=Dict("units"=>"degrees_east", "standard_name"=>"longitude","grid_mapping"=>"crs"))
        end
        
        # add data variables
        
        # add projection information
        cfdict=pyconvert(Dict,out_geo.to_cartopy_crs().to_cf())
        defVar(ncid, "crs", 1,() , attrib=cfdict)
        alldata=Dict("forward"=>alldata.forward, "backward"=>alldata.backward, "combined"=> let f=alldata.forward, b=alldata.backward
            (;data=[f.data;b.data],
            err=[f.err;b.err], 
            lat=[f.lat;b.lat], 
            lon=[f.lon;b.lon],
            it=[f.it;b.it],
            res=[f.res;b.res])
        end
        )
        #resampled and stored in group is the output file
        for dir in ("forward","backward","combined")
            datablock=alldata[dir].data
            errblock=alldata[dir].err
            lat=alldata[dir].lat
            lon=alldata[dir].lon
            itblock=alldata[dir].it
            resblock=alldata[dir].res

            npblock=np[].array(datablock, dtype=np[].float32)
            nperrblock=np[].array(errblock, dtype=np[].float32)
            in_geo=pr[].geometry.SwathDefinition(lons=lon, lats=lat)
            npitblock=np[].array(itblock, dtype=np[].int64)
            npresblock=np[].array(resblock, dtype=np[].float32)

            npblockresampled=pr[].kd_tree.resample_nearest(in_geo, npblock, out_geo, radius_of_influence=15000, fill_value=NaN) |> x -> pyconvert(Array, x)
            nperrblockresampled=pr[].kd_tree.resample_nearest(in_geo, nperrblock, out_geo, radius_of_influence=15000, fill_value=NaN) |> x -> pyconvert(Array, x)
            npitblockresampled=pr[].kd_tree.resample_nearest(in_geo, npitblock, out_geo, radius_of_influence=15000, fill_value=0) |> x -> pyconvert(Array, x)
            npresblockresampled=pr[].kd_tree.resample_nearest(in_geo, npresblock, out_geo, radius_of_influence=15000, fill_value=NaN) |> x -> pyconvert(Array, x)

            for i= 1:9
                if grid == "native"
                    defVar(getfield(ds, Symbol(dir)) , standard_names[i], 
                    Float64, ("npoints",), attrib=Dict("standard_name"=>standard_names[i], "units"=>units[i], "long_name"=>longnames[i]), deflatelevel=5, shuffle=true, chunksizes=(length(lat),))[:]=npblockresampled[:,i]
                    defVar(getfield(ds, Symbol(dir)) , standard_names[i]*"_error", Float64, ("npoints",) , attrib=Dict("standard_name"=>standard_names[i]*"_error", "units"=>units[i], "long_name"=>longnames[i]*" error"), deflatelevel=5, shuffle=true, chunksizes=(length(lat),))[:]=nperrblockresampled[:,i]
                    
                else
                    defVar(getfield(ds, Symbol(dir)) , standard_names[i],
                    Float64 , ("x", "y"), attrib=Dict("standard_name"=>standard_names[i], "units"=>units[i], "long_name"=>longnames[i],"grid_mapping"=>"crs"), deflatelevel=5, shuffle=true, chunksizes=(xsize,ysize))[:]=npblockresampled[:,:,i]'
                    defVar(getfield(ds, Symbol(dir)) , standard_names[i]*"_error", Float64, ("x","y"),attrib=Dict("standard_name"=>standard_names[i]*"_error", "units"=>units[i], "long_name"=>longnames[i]*" error","grid_mapping"=>"crs"), deflatelevel=5, shuffle=true, chunksizes=(xsize,ysize))[:] = nperrblockresampled[:,:,i]'
                end
            end
            #add flags according to the L2 product definition
            v= defVar(getfield(ds, Symbol(dir)) ,"quality_flag" , Int64, ("x", "y"),attrib=Dict("standard_name"=>"quality_flag", "units"=>"1", "long_name"=>"Quality flag","grid_mapping"=>"crs"), deflatelevel=5, shuffle=true, chunksizes=(xsize,ysize))
            #create the flags from the it and res blocks and the variables
            flag=zeros(Int64,size(npitblockresampled))
            #valid solution if all variables are finite
            flag .|= dropdims(any(isfinite.(npblockresampled[:,:,:]), dims=3),dims=3)' * 2^0  
            #convergence flag if it is less than 50 and larger than 0
            flag .|= (any(npitblockresampled[:,:].<50, dims=3) .| any(npitblockresampled[:,:].>0, dims=3))' * 2^1 
            #fallback solver used if it is less than 50 (not implemented, so set to false)
            flag .|= false * 2^2
            #fallback solver converged if it is less than 50 (not implemented, so set to false)
            flag .|= false * 2^3
            #no convergence if it is = 50 
            flag .|= npitblockresampled[:,:].==50 * 2^4
            #anomaly detected if it is = 50 #not implemented, so set to false
            flag .|= false * 2^5
            #invalid wind_speed if wind_speed is NaN
            flag .|= isnan.(npblockresampled[:,:,1]') * 2^6
            #invalid total_water_vapor if total_water_vapor is NaN
            flag .|= isnan.(npblockresampled[:,:,2]') * 2^7
            #invalid cloud_liq_water if cloud_liq_water is NaN
            flag .|= isnan.(npblockresampled[:,:,3]') * 2^8
            #invalid sea_surface_temperature if sea_surface_temperature is NaN
            flag .|= isnan.(npblockresampled[:,:,4]') * 2^9
            #invalid ice_surface_temperature if ice_surface_temperature is NaN
            flag .|= isnan.(npblockresampled[:,:,5]') * 2^10
            #invalid sea_ice_fraction if sea_ice_fraction is NaN
            flag .|= isnan.(npblockresampled[:,:,6]') * 2^11
            #invalid multi_year_ice_fraction if multi_year_ice_fraction is NaN
            flag .|= isnan.(npblockresampled[:,:,7]') * 2^12
            #invalid sea_ice_thickness if sea_ice_thickness is NaN
            flag .|= isnan.(npblockresampled[:,:,8]') * 2^13
            #anomaly in residual if residual is NaN or larger than 50 Kelvin, i.e. very far from forward model result
            flag .|= dropdims(any(npresblockresampled[:,:,:].>50, dims=3),dims=3)' * 2^14 

            #anomaly in L-band #not really defined, so set to false
            flag .|= false * 2^15
            #anomaly in C-band #not really defined, so set to false
            flag .|= false * 2^16
            #anomaly in X-band #not really defined, so set to false
            flag .|= false * 2^17
            #anomaly in Ku-band #not really defined, so set to false
            flag .|= false * 2^18
            #anomaly in K-band #not really defined, so set to false
            flag .|= false * 2^19
            #landmask if landmask is true
            lons,lats=out_geo.get_lonlats()
            landmask = basemap[].maskoceans(lons,lats,lats,resolution="i").mask |> PyArray |> Array |> x->.~x' 
            flag[landmask].=2^50
            #ice shelf mask if ice shelf mask is true #not implemented, so set to false
            flag .|= false * 2^51
            v[:]=flag
        end

        # add metadata
   
        ncid.attrib["project"]="ESA CIMR DEVALGO (contract 4000137493)"
        ncid.attrib["project_lead"]="Thomas Lavergne"
        ncid.attrib["project_lead_email"]="thomas.lavergne@met.no"
        ncid.attrib["date_created"]=Dates.now() |> string
        ncid.attrib["processing_level"]="Level-2"
        ncid.attrib["standard_name_vocabulary"]="CF Standard Name Table (Version 83, 17 October 2023)"
        ncid.attrib["spacecraft"]="CIMR"
        ncid.attrib["instrument"]="CIMR"
        ncid.attrib["product_level"]="2"
       
        ncid.attrib["product_name"]="multi parameter retrieval"
        ncid.attrib["variable_list"]="wind speed, total water vapor, cloud liquid water, sea surface temperature, ice surface temperature, sea ice concentration, multi year ice fraction, first year ice thickness, sea surface salinity"
        ncid.attrib["product_version"]="0.1"
        ncid.attrib["author"]="Marcus Huntemann"
        ncid.attrib["author_email"]="macrus.huntemann@uni-bremen.de"
    end
   nothing 
end

