#This is a command line script running the MultiParameter algorithm on a given CIMR L1b input file.

#The script reads the L1b file, projects the data to a common grid, runs the MultiParameter algorithm and stores the output in a netcdf file.

#The script is called with the following arguments:
#filename of the L1b file
#filename of the output netcdf file
#optional: outgrid (default is "ease2_nh")

#The script is called like this:
#julia run_mp.jl l1b_file.nc output_file.nc 

##
using Revise
includet(joinpath(@__DIR__(),"../book/common_functions.jl"))
#using OEM

inpath = get(ENV, "DEVALGO_INPUT_DATA_PATH", "/mnt/spaces/Projects/2022_CIMR-DEVALGO/DATA")
if isempty(ARGS)
    fn=joinpath(inpath,"SCEPS/SCEPS_l1b_sceps_geo_polar_scene_1_unfiltered_tot_minimal_nom_nedt_apc_tot_v2p1.nc")
else
    fn = ARGS[1]
    @assert isfile(fn)
end

if length(ARGS) < 2
   outfn = "out_polar.nc"
else
    outfn = ARGS[2]
end

if length(ARGS) < 3
    outgrid = "ease2_nh"
else
    outgrid = ARGS[3]
    @assert outgrid in ["ease2_nh", "ease2_sh"]
end


#forward
findata,finerr,flat,flon=project_data(fn,"C_BAND","forward")

fdat,ferr,fres,fit=run_oem(findata,band_error=finerr)
#backward
bindata,binerr,blat,blon=project_data(fn,"C_BAND","backward")
bdat,berr,bres,bit=run_oem(bindata,band_error=binerr)
##

thedata=(;forward=(;data=fdat,err=ferr,lat=flat,lon=flon,res=fres,it=fit),
         backward=(;data=bdat,err=berr,lat=blat,lon=blon,res=bres,it=bit))


println("storing output ", outfn, ", for grid ", outgrid)
store_output(outfn, thedata, outgrid)






