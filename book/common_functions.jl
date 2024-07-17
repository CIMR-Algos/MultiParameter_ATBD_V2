##imports_start
using Pkg
Pkg.activate("../algorithm/algoenv")
using OEM
using IJulia
using Printf
using Markdown
using Proj
using LinearAlgebra
using NCDatasets
using PythonCall # to call python from julia
using PythonPlot # basically the matplotlib.pyplot module
using Statistics
using Markdown
using OrderedCollections
pygc = pyimport("gc")
pr = pyimport("pyresample")
#cfeature = pyimport("cartopy.feature")
basemap= pyimport("mpl_toolkits.basemap")
LLtoEASE2n = Proj.Transformation("epsg:4326", "epsg:6931")
EASE2NtoLL = inv(LLtoEASE2n)

lltoxy=LLtoEASE2n #conversion from lat/lon to EASE2 grid coordinates
xytoll=EASE2NtoLL #conversion from EASE2 grid coordinates to lat/lon
lea = pr.geometry.AreaDefinition("ease2_nh_testscene","ease2_nh_testscene","ease2_nh_testscene","EPSG:6931",1400,1400,[0,-1.5e6,1.4e6,-1e5]) # uses pyresample to define an area definition for the test region in the EASE2 grid
leaextent=[lea.area_extent[0],lea.area_extent[2],lea.area_extent[1],lea.area_extent[3]] #extent for the image, since order differs from the order in the area definition

#getting the landmask for the region
lealons,lealats=lea.get_lonlats()
landmask = basemap.maskoceans(lealons,lealats,lealats,resolution="i").mask |> PyArray |> Array |> x->.~x 

#defining a cartopy crs
ccrs=lea.to_cartopy_crs()
##imports_end



##plot_all_params_start
"""
    plot_all_params(dat, x, y, swath=true,mask=true,errors=false)

Plot all parameters from the OEM output in a 3x3 grid.
parameter:
* dat: the output of the OEM algorithm (or error)
* swath: if true, plot the data as a scatter plot, if false, project the data 
to the EASE2 grid and plot as an image
* mask: if true, mask the data according to the sea ice concentration
* errors: if :true (or :std or stderr), plot the error data instead of the data, changing the color scale, limits, etc. accordingly
"""
function plot_all_params(data,direction; swath=true,mask=(:sic,),errors=:false,save=false)
    if direction=="combined"
        ldata= let dat=vcat(data["forward"].dat,data["backward"].dat),
             x=vcat(data["forward"].x,data["backward"].x),y=vcat(data["forward"].y,data["backward"].y),lat=vcat(data["forward"].lat,data["backward"].lat),lon=vcat(data["forward"].lon,data["backward"].lon),errr=vcat(data["forward"].errr,data["backward"].errr),its=vcat(data["forward"].its,data["backward"].its)
            (;dat,x,y,lat,lon,errr,its) |> deepcopy
        end
    else
        ldata=data[direction] |> deepcopy
    end


    
    lealons,lealats=lea.get_lonlats()
    landmask = basemap.maskoceans(lealons,lealats,lealats,resolution="i").mask |> PyArray |> Array |> x->.~x 

    ccrs=lea.to_cartopy_crs()


    fig, ax = pyplot.subplots(nrows=3, ncols=3, figsize=[15, 15], subplot_kw=PyDict(Dict("projection" => ccrs)))
    ax = pyconvert(Array, ax)
    x=ldata.x
    y=ldata.y
    lat=ldata.lat
    lon=ldata.lon
    if errors in (:true,:std,:stderr)
        dat = ldata.errr
        mins = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        #maxs = [10, 10, 0.5, 10, 10, 0.3, 0.3, 100, 5]
        maxs = diag(sqrt.(OEM.S_p))'

        labs = ["wind speed error [m/s]", "total water vapor error [kg/m²]", "cloud liquid water error [kg/m²]",
            "sea surface temperature error [K]", "ice surface temperature error [K]", "sea ice concentration error [1]",
            "multi year ice fraction error [1]", "first year ice thickness error [cm]", "sea surface salinity error [ppt]"]
    elseif errors==:ratio
        dat =  1.0 .- (ldata.errr ./ diag(sqrt.(OEM.S_p))')
        mins = zeros(9)
        maxs = ones(9)
        labs = ["wind speed uncertainty ratio", "total water vapor uncertainty ratio", "cloud liquid water uncertainty ratio",
            "sea surface temperature uncertainty ratio", "ice surface temperature uncertainty ratio", "sea ice concentration uncertainty ratio",
            "multi year ice fraction uncertainty ratio", "first year ice thickness uncertainty ratio", "sea surface salinity uncertainty ratio"]
    elseif errors==:false    
        dat = ldata.dat
        mins = [0, 0, 0, 250, 250, 0, 0, 0, 25]
        maxs = [20, 20, 0.5, 280, 280, 1, 1, 150, 40]
        labs = ["wind speed [m/s]", "total water vapor [kg/m²]", "cloud liquid water [kg/m²]",
        "sea surface temperature [K]", "ice surface temperature [K]", "sea ice concentration [1]",
        "multi year ice fraction [1]", "first year ice thickness [cm]", "sea surface salinity [ppt]"]
    end
    
    cm=pyplot.cm.viridis
    cm.set_bad("grey")
    
    if :sic ∈ mask
        sic=ldata.dat[:,6] #sea ice concentration from original data input, independent if errors is set or not
        dat[sic.>0.9,1].=NaN
        dat[sic.>0.9,4].=NaN
        dat[sic.<0.1,5].=NaN
        dat[sic.<0.1,7].=NaN
        dat[sic.<0.1,8].=NaN
        dat[sic.>0.9,9].=NaN
    end


    for i = 1:9
        if swath
            pos = ax[i].scatter(x, y, c=dat[:, i], s=0.8, vmin=mins[i], vmax=maxs[i])
        else
            cimg = project_data_to_area_nn(dat[:, i], lat, lon, lea)
            cimg = reshape(cimg, 1400, 1400)
            leaextent = [lea.area_extent[0], lea.area_extent[2], lea.area_extent[1], lea.area_extent[3]]
            if :land ∈ mask
                cimg[landmask] .= NaN
            end
            pos = ax[i].imshow(cimg, zorder=1, extent=leaextent, vmin=mins[i], vmax=maxs[i], origin="upper",cmap=cm)
        end
        ax[i].set_title(labs[i])
        cb = fig.colorbar(pos, ax=ax[i], shrink=0.7, orientation="vertical", pad=0.01)
        ax[i].coastlines()
    end
    pyplot.subplots_adjust(wspace=0.01, hspace=0.01)

    if save!=false
        pyplot.savefig(save,dpi=200,bbox_inches="tight")
    end

    return fig
end
##plot_all_params_end

##calculate_overall_performance_start
"""
calculate_overall_performance(scenes,gdatas,scene_labels,borders)
calculates the overall performance of several scenes and returns an array of values. If direction is set to "combined" it combines the forward and backward retrieval during the reprojection.
scenes is just a array of scenes, gdatas is an array of the truths for corresponding to the scenes
"""
function calculate_overall_performances(scenes,gdatas,directions,borders)
    thevars= ["ws","twv","clw","sst","ist","sic","myi","sit","sss"]
    for direction in directions 
        @assert direction in ["forward","backward","combined"]
    end

    @assert length(scenes)==length(gdatas)==length(borders)
    results=OrderedDict()
    for j in axes(thevars,1)
        results[thevars[j]]=OrderedDict()
        for direction in directions
            diffs = []
            truthes = []
            for i in axes(scenes,1)
                if direction=="combined"
                    cimg = project_data_to_area_nn([scenes[i]["forward"].dat[:, j]; scenes[i]["backward"].dat[:,j]], [scenes[i]["forward"].lat; scenes[i]["backward"].lat], [scenes[i]["forward"].lon; scenes[i]["backward"].lon], lea)
                else
                    cimg = project_data_to_area_nn(scenes[i][direction].dat[:, j], scenes[i][direction].lat, scenes[i][direction].lon, lea)
                end
                cimg = reshape(cimg, 1400, 1400)
                cimg_ref=gdatas[i][Symbol(thevars[j])][:]
                cimg_ref=reshape(cimg_ref,1400,1400)
                diff=cimg.-cimg_ref
                
                diff=diff[borders[i]:end-borders[i],borders[i]:end-borders[i]][:] |> x->filter(!isnan,x)
                push!(diffs,diff)
                push!(truthes,cimg_ref[borders[i]:end-borders[i],borders[i]:end-borders[i]][:] |> x->filter(!isnan,x))
            end
            results[thevars[j]][direction]=(;meanvalue=mean(vcat(truthes...)),mean=mean(vcat(diffs...)),std=std(vcat(diffs...)))
        end
    end
     return results
end

function markdown_table(results,title,label;include_mean=false,show_units=true)
    thevars= ["ws","twv","clw","sst","ist","sic","myi","sit","sss"]
    units=Dict("ws"=>"m/s","twv"=>"kg/m²","clw"=>"kg/m²","sst"=>"K","ist"=>"K","sic"=>"1","myi"=>"1","sit"=>"cm","sss"=>"ppt")
    #create markdown table, MyST markdown if in jupyter book compilation, normal otherwise
    outstring="```{list-table} Performance metric for $title\n"
    outstring *=":name: \"performance-$label\"\n"
    outstring *=":header-rows: 1\n"
    
    directions= results[thevars[1]] |> keys

    sdirs=Dict("forward"=>"fw","backward"=>"bw","combined"=>"cb")
    if get(ENV,"JUPYTER_BOOK_BUILD","false") == "true"
        if include_mean
            outstring *="* - Variable\n "* mapreduce(x->"  - $(sdirs[x]) mean value\n  - $(sdirs[x]) mean(diff)\n  - $(sdirs[x]) std(diff)\n",*,directions)
        else
            outstring *="* - Variable\n "* mapreduce(x->"  - $(sdirs[x]) mean(diff)\n  - $(sdirs[x]) std(diff)\n",*,directions)
        end
        for var in thevars
            vals = results[var]
            if show_units
                outstring *= "* - $(uppercase(var)) [$(units[var])]\n"
            else
                outstring *= "* - $(uppercase(var))\n"
            end
            if include_mean
                outstring *= mapreduce(x->"  - $(@sprintf "%.3f" x.meanvalue)\n  - $(@sprintf "%.3f" x.mean)\n  - $(@sprintf "%.3f" x.std)\n",*,(vals[direction] for direction in directions))
            else
                outstring *= mapreduce(x->"  - $(@sprintf "%.3f" x.mean)\n  - $(@sprintf "%.3f" x.std)\n",*,(vals[direction] for direction in directions)) 
            end
            #outstring *=  (@sprintf "  - %.3f\n  - %.3f\n" vals.mean vals.std)
            #
        end
        outstring *= "```"
    else
        outstring= "~~~markdown\n"*outstring*"\n# table content (skipped in notebook) ```\n~~~\n"
        outstring *="| Variable |"
        if include_mean
            outstring *= mapreduce(x->" $(sdirs[x]) mean value | $(sdirs[x]) mean(diff) | $(sdirs[x]) std(diff) |",*,directions) *"\n"
        else
            outstring *= mapreduce(x->" $(sdirs[x]) mean(diff) | $(sdirs[x]) std(diff) |",*,directions) *"\n"
        end
        #outstring *="| Variable |"*" mean(diff) | std(diff) |"^(length(directions))*"\n"
        if include_mean
            outstring *="| --- "
        end
        outstring *="| --- |"*" --- | --- |"^(length(directions) )*"\n"
        for var in thevars
            vals = results[var]
            if show_units
                outstring *= "| $(uppercase(var)) [$(units[var])] |"
            else
                outstring *= "| $(uppercase(var)) |"
            end
            if include_mean
                outstring *= mapreduce(x->" $(@sprintf "%.3f" x.meanvalue) | $(@sprintf "%.3f" x.mean) | $(@sprintf "%.3f" x.std) |\n",*,(vals[direction] for direction in directions))
            else
                outstring *= mapreduce(x->" $(@sprintf "%.3f" x.mean) | $(@sprintf "%.3f" x.std) |\n",*,(vals[direction] for direction in directions))
            end
            #outstring *= "| $(uppercase(var)) | $(@sprintf "%.3f" vals.mean) | $(@sprintf "%.3f" vals.std) |\n"
        end
    end

    return outstring
end
##calculate_overall_performance_end

##calculate_forward_backward_performance_start
# calculate the mean and std of the difference between the CIMR input and the CIMR output for all variables, give output as makdown table (MYST format)

"""
    markdown_performance_table(scene, lat, lon, gdata;scene_label="1",border=100)
calculates the mean and rmsd difference between the CIMR input and the CIMR output for all variables, and returns the results as a markdown table
* scene: the output of the OEM algorithm in l1b geometry
* gdata: retrieval on the input TBs to the CIMR simulator
"""
function markdown_performance_table(scene, gdata;scene_label="1",border=100)
    outstring="```{list-table} performance metric for scene $scene_label retrieval parameters\n"
    outstring *=":name: \"performance-$scene_label\"\n"
    outstring *=":header-rows: 1\n"
    outstring *="* - variable\n  - fw mean difference\n  - fw spread\n  - bw mean difference\n  - bw spread\n"
    thevars= ["ws","twv","clw","sst","ist","sic","myi","sit","sss"]
    for (i,var) in enumerate(thevars)
        outstring *= "* - $(uppercase(var))\n"
        for dir in ["forward","backward"]
            cimg = project_data_to_area_nn(scene[dir].dat[:, i], scene[dir].lat, scene[dir].lon, lea)
            cimg = reshape(cimg, 1400, 1400)
            cimg_ref=gdata[Symbol(var)][:]
            cimg_ref=reshape(cimg_ref,1400,1400)
            diff=cimg.-cimg_ref
            diff=diff[border:end-border,border:end-border]
            diff=filter(!isnan,diff)
            outstring *=  (@sprintf "  - %.2f\n  - %.2f\n" mean(diff) std(diff))
        end
    end
    outstring *= "```"
    return outstring
end
##calculate_forward_backward_performance_end


##fns_start
#the input files can be found on the ESA DEVALGO FTP server,
#this is a named tuple with the paths to the input files
fns = (; scene1=
        (; l1b="/mnt/spaces/Projects/2022_CIMR-DEVALGO/CIMR_L1X/W_PT-DME-Lisbon-SAT-CIMR-1B_C_DME_20230420T103323_LD_20280110T114800_20280110T115700_TN.nc",
            input="/mnt/spaces/Projects/2022_CIMR-DEVALGO/Test_scenes_downscaled_projected_20230923/test_scene_1_compressed_lowres.nc"),
        scene2=
        (; l1b="/mnt/spaces/Projects/2022_CIMR-DEVALGO/CIMR_L1X/W_PT-DME-Lisbon-SAT-CIMR-1B_C_DME_20230417T105425_LD_20280110T114800_20280110T115700_TN.nc",
            input="/mnt/spaces/Projects/2022_CIMR-DEVALGO/Test_scenes_downscaled_projected_20230923/test_scene_2_compressed_lowres.nc"),
        scenep=
        (;
            l1b="/mnt/spaces/Projects/2022_CIMR-DEVALGO/SCEPS_l1b_sceps_geo_polar_scene_1_unfiltered_tot_minimal_nom_nedt_apc_tot_v2p1.nc",
            input="/mnt/spaces/Projects/2022_CIMR-DEVALGO/cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.nc",
            geo="/mnt/spaces/Projects/2022_CIMR-DEVALGO/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_surface.nc"
        )
)
##fns_end


function plot_with_caption(figure_instance, caption, label)
    # Generate the plot
    
    # Save the plot to a file
    mkpath("figures")
    filename = "figures/figure_$(label).png"
    figure_instance.savefig(filename, bbox_inches="tight")
    pyplot.close(figure_instance)
    
    # Generate MyST Markdown
    myst_markdown = """
    ```{figure} $(filename)
    :name: $(label)

    $(caption)
    ```
    """
    
    if get(ENV, "JUPYTER_BOOK_BUILD", "false") == "true"
         display("text/markdown",myst_markdown)
         return nothing #Markdown.parse(myst_markdown)
    else # Display the plot in the notebook if inside a Jupyter notebook
        display("image/png",figure_instance)
        display("text/markdown", "~~~markdown\n$(myst_markdown)\n~~~")
        return nothing
    end
end