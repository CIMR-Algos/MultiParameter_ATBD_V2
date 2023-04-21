export get_amsr2_data_on_NSIDC
export get_amsr2_tbs
export get_smos_data_on_NSIDC
export prepare_oem_tbs
export read_ecmwf
export get_ecmwf_on_NSIDC
export nsidcg


using HDF5
using NCDatasets
using PythonCall
using Statistics
using Dates
using Glob

pr=pyimport("pyresample")
np=pyimport("numpy")
gm=pr.geometry
kdt=pr.kd_tree
gm=pr.geometry
nsidcg=gm.AreaDefinition("nsidc_n","NSIDC polar stereographic north","nsidc_n","+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",608,896, [-3850000,-5350000,3750000,5850000])


#### ERA 5 data reading

function read_ecmwf(date)
    #date is of instance Date or DateTime from Dates from julia Base
    @assert date isa Union{Date,DateTime}
    
    y,m,d=Dates.format.(Ref(date),["yyyy","mm","dd"])

  
    #fname = "/mnt/spaces/Models/ERA5/incoming/$(y)/era5_$(y*m*d)_250_1hr_arc.nc"  #ERA 5 data
    @show fname = "/Users/huntemann/spaces/Models/ERA5/incoming/$(y)/era5_$(y*m*d)_250_1hr_arc.nc"  #ERA 5 data

    #    print(fname)
    ecmwf = NCDataset(fname);
    #Date = DateTime(parse(Int32,year),parse(Int32,month),parse(Int32,day),parse(Int32,hour),parse(Int32,minute))
    #     println("Date ",date)
    #target = round(Date, Dates.Hour(3)) #ERAInterim
    variables = ["tcwv","u10","v10","skt","t2m","tclw"]
    ec=Array{Float64}(undef,size(ecmwf["tcwv"])[1:2]...,6)


    if (date isa DateTime)
        target = floor(Date, Dates.Hour(1)) #ERA5  (hourly)   # floor to avoid next day
        timeind = findall(ecmwf["time"] .== target)

        for (i,var) in enumerate(variables)
            ec[:,:,i]=ecmwf[var][:,:,timeind]
        end
    else
        for (i,var) in enumerate(variables)
            ec[:,:,i]=mean(ecmwf[var][:,:,:],dims=3)
        end
    end
    lon=ecmwf["longitude"][:].+0.0
    lat=ecmwf["latitude"][:].+0.0
    olat=transpose(repeat(lat,outer=(1,length(lon))))
    olon=repeat(lon,outer=(1,length(lat)))
    #olon[:,1].+=360
    
    dt_out = Dict("wsp" => sqrt.(ec[:,:,2].^2 .+ ec[:,:,3].^2)[:],
    "twv" => ec[:,:,1][:],
    "clw" => ec[:,:,6][:],
    "tsk" => ec[:,:,4][:],
    "t2m" => ec[:,:,5][:],
    "lon" => olon[:],
    "lat" => olat[:])
    
    return dt_out
end

function get_ecmwf_on_NSIDC(date)
    ecmwf=read_ecmwf(date)
    
    lats=ecmwf["lat"]
    lons=ecmwf["lon"]
    SD=gm.SwathDefinition(lons,lats)
    data=np.array([ecmwf["wsp"] ecmwf["twv"] ecmwf["clw"] ecmwf["t2m"] ecmwf["tsk"]])
    out=kdt.resample_nearest(SD,data,nsidcg,25000) |> x->pyconvert(Array,x) |>x-> eachslice(x,dims=3) .|> collect
    return out
end


#### SMOS data reading

function get_smos_data_on_NSIDC(date)
    y,m,d=Dates.format.(Ref(date),["yyyy","mm","dd"])
    fn="../data/$(y*m*d)_TBs.nc"
    ncd=NCDataset(fn)
    outtbh=variable(ncd,"TB_val_h_inc_50_56")[:] |> x-> permutedims(x,[2,1]) |> x->reverse(x,dims=1)
    outtbv=variable(ncd,"TB_val_v_inc_50_56")[:] |> x-> permutedims(x,[2,1]) |> x->reverse(x,dims=1)

    return outtbv,outtbh
end

function get_amsr2_tbs(fn)
    dd=h5open(fn)
    V89=read(dd,"Brightness Temperature (89.0GHz-A,V)")[1:2:end,:][:]/100f0
    H89=read(dd,"Brightness Temperature (89.0GHz-A,H)")[1:2:end,:][:]/100f0
    V6=read(dd,"Brightness Temperature (6.9GHz,V)")[:]/100f0
    H6=read(dd,"Brightness Temperature (6.9GHz,H)")[:]/100f0
    V10=read(dd,"Brightness Temperature (10.7GHz,V)")[:]/100f0
    H10=read(dd,"Brightness Temperature (10.7GHz,H)")[:]/100f0
    V18=read(dd,"Brightness Temperature (18.7GHz,V)")[:]/100f0
    H18=read(dd,"Brightness Temperature (18.7GHz,H)")[:]/100f0
    V23=read(dd,"Brightness Temperature (23.8GHz,V)")[:]/100f0
    H23=read(dd,"Brightness Temperature (23.8GHz,H)")[:]/100f0
    V36=read(dd,"Brightness Temperature (36.5GHz,V)")[:]/100f0
    H36=read(dd,"Brightness Temperature (36.5GHz,H)")[:]/100f0
    lat=read(dd,"Latitude of Observation Point for 89A")[1:2:end,:][:]
    lon=read(dd,"Longitude of Observation Point for 89A")[1:2:end,:][:]
    Ql=read(dd,"Pixel Data Quality 6 to 36")[:]
    Qh=read(dd,"Pixel Data Quality 89")[:]
    #println(size(read(dd,"Brightness Temperature (89.0GHz-A,V)")/100f0))
    #println(size(read(dd,"Brightness Temperature (6.9GHz,V)")/100f0))
    slope = [1.008, #6V                       
    1.000, #6H                        
    1.005, #10V                       
    0.993, #10H                       
    1.031, #18V                     
    1.001, #18H                        
    0.999, #23V                        
    1.002, #23H                       
    0.997, #36V                               
    0.996, #36H                            
    0.989, #89V  
    0.977 #89H 
    ]

intercept = [-2.360,  #6V      #from Janna             
    -1.538,  #6H                          
    -4.551,  #10V                       
    -2.585,  #10H                       
    -9.710,  #18V                         
    -1.104,  #18H                        
    -1.706,  #23V                       
    -2.638,  #23H                       
    -2.610,  #36V                       
    -2.687,  #36H                       
    0.677,  #89V   
    3.184  #89H
    ]
    outvec=[V6,H6,V10,H10,V18,H18,V23,H23,V36,H36,V89,H89]
    for i in 1:12
        outvec[i]=outvec[i].*slope[i].+intercept[i]
    end

    return lat,lon, [V6,H6,V10,H10,V18,H18,V23,H23,V36,H36,V89,H89]
end



#### AMSR2 data reading

function get_amsr2_tbs(fns::Vector{String})
    alllat=Float32[]
    alllon=Float32[]
    #ttb=Float32[]
    alltbs=[Vector{Float32}() for i=1:12]
    for fn in fns
        lat,lon,tbs=get_amsr2_tbs(fn)
        append!(alllat,lat)
        append!(alllon,lon)
        #append!(ttb,tbs[4])
        append!.(alltbs,tbs)
    end
    return alllat,alllon,alltbs
end


function get_amsr2_data_on_NSIDC(date)
    y,m,d=Dates.format.(Ref(date),["yyyy","mm","dd"])
    
    #    G=SatResample.NSIDCN(12500)
    #    fnsamsr=glob("mnt/spaces/Radiometers/AMSR2/L1B/$(y)/$(y).$(m)/GW1AM2_$(y*m*d)*","/")
    fnsamsr=glob("Users/huntemann/spaces/Radiometers/AMSR2/L1B/$(y)/$(y).$(m)/GW1AM2_$(y*m*d)*","/")
    lat,lon,tbs=get_amsr2_tbs(fnsamsr)
    idx=lat.>50;
    rlat=lat[idx]
    rlon=lon[idx]
    #rtbs=tbs[idx]
    
    rtbs=np.array(hcat([tb[idx] for tb in tbs]...))
    sd=gm.SwathDefinition(rlon,rlat)
    
    outg=kdt.resample_nearest(sd,rtbs,nsidcg,10000,1.0)
    
    #    outg=SatResample.llresample(rlon,rlat,rtbs,G,12500,10)
    return pyconvert(Array,outg)
end

#### General routines
function prepare_oem_tbs(date)
    outg=[copy(get_amsr2_data_on_NSIDC(date)[:,:,i]) for i=1:12]
    @show typeof(outg)
    outtbv,outtbh=get_smos_data_on_NSIDC(date)
    insert!(outg,1,outtbh)
    insert!(outg,1,outtbv)
    return outg
end