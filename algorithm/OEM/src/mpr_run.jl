#This is a command line script running the MultiParameter algorithm on a given CIMR L1b input file.

#The script reads the L1b file, projects the data to a common grid, runs the MultiParameter algorithm and stores the output in a netcdf file.

#The script is called with the following arguments:
#filename of the L1b file
#filename of the output netcdf file
#optional: outgrid (default is "ease2_nh")

#The script is called like this:
#julia run_mp.jl l1b_file.nc output_file.nc 

##
include("OEM_with_SIT.jl")
using .OEM


if isempty(ARGS)
   fn = "/mnt/spaces/Projects/2022_CIMR-DEVALGO/SCEPS_l1b_sceps_geo_polar_scene_1_unfiltered_tot_minimal_nom_nedt_apc_tot_v2p1.nc"
   outfn = "out_polar.nc"
else
    fn = ARGS[1]
    outfn = ARGS[2]
end

outgrid = length(ARGS) >= 3 ? ARGS[3] : "ease2_nh"
direction = "forward"

#forward
findata,finerr,flat,flon=project_data(fn,"C_BAND","forward")
println(sizeof(findata))
println(sizeof(finerr))

fdat,ferr,_=run_oem(findata,finerr)
#backward
bindata,binerr,blat,blon=project_data(fn,"C_BAND","backward")
bdat,berr,_=run_oem(bindata,binerr)
##

thedata=(;forward=(;data=fdat,err=ferr,lat=flat,lon=flon),
         backward=(;data=bdat,err=berr,lat=blat,lon=blon))


println("storing output ", outfn, ", for grid ", outgrid)
store_output(outfn, thedata, outgrid)







