module FBK_GIS_LCOH

    export GISwind, GISsolar, GIShydro, GIStemp, GISturbines, makedistances, annualwindindex, GISlcoh, read_save_region

    using MAT, HDF5, ProgressMeter, Random, Interpolations, BenchmarkTools, Images,
        Statistics, DelimitedFiles, Dates, NCDatasets, JLD, Parameters, ImageSegmentation,
        StatsBase, CSV, Distances, Printf, TimeZones, DataFrames, StatsBase, XLSX


    include("rasterize_shapefiles.jl")
    include("make_auxiliary_datasets.jl")
    include("getconfig.jl")
    include("era5download.jl")  #Occhio che qua ho commentato perché dentro ha PyCall (che no nva bene con sysimage)
    include("helperfunctions.jl")
    include("make_regions.jl")
    include("regiondefinitions.jl")
    include("solarposition.jl")
    include("makewindera5.jl")
    include("makesolarera5.jl")
    include("maketempera5.jl")
    include("makedistances.jl")
    include("GISwind.jl")
    include("GISsolar.jl")
    include("GIShydro.jl")
    include("GIStemp.jl")
    include("GISturbines.jl")
    include("downloaddatasets.jl")
    include("sspregiondefinitions.jl")
    include("mapping.jl")
    include("syntheticdemand_inputdata.jl")
    include("syntheticdemand_training.jl")
   #=  include("funz_test_binary.jl") =#
    include("GISlcoh_VSCode.jl")

    # function saveregions_x_exe()
    #     saveregions("Brunei", brunei)
    # end
    #
    # function julia_main()
    #     saveregions_x_exe()
    # end

    function GISlcoh_x_exe()
        path = homedir()
        excel_file = XLSX.readxlsx("$path\\save_regions.xlsx") ##Leggo Excel
        config_settings = excel_file["Settings"]
        datafolder_path = config_settings[:][4,2]
        copernicus_UID = config_settings[:][5,2]
        copernicus_API_Key = config_settings[:][6,2]
        download_convert_ERA5 = config_settings[:][1,2]
        download_convert_AUX = config_settings[:][3,2]
        ERA5_year = config_settings[:][2,2]
        saveconfig(datafolder_path, copernicus_UID, copernicus_API_Key; agree_terms=true)

        if download_convert_ERA5=="yes"
            println("Downloading and converting ERA5 data for the year $ERA5_year...")
        elseif download_convert_ERA5=="no"
            println("Using already existing data for the year $ERA5_year...")
        else
            println("Input for Download_convert_ERA5 not correct: possible options 'yes' or 'no'.")
        end

        if download_convert_AUX=="yes"
            println("Downloading and converting auxiliary datasets.")
        elseif download_convert_AUX=="no"
            println("Using already existing auxiliary datasets.")
        else
            println("Input for Download_convert_Auxiliary not correct: possible options 'yes' or 'no'.")
        end

        region_name = excel_file["Region Name"]
        data_path = "$path\\GISData"
        region_file = "$data_path\\regions_$(region_name[:][1]).jld"
        if isfile("$region_file") == false
            println("Creating region: $(region_name[:][1])")
            read_save_region()
        else
            println("Region: $(region_name[:][1]) already exists. \n")
        end
        println("Calculating profiles for region: $(region_name[:][1]).")
        GISlcoh(gisregion="$(region_name[:][1])")
    end
    
    function julia_main()
        GISlcoh_x_exe()
    end

# function GISsequence(regionname, regionmat)
#     saveregions(regionname, regionmat)
#     makedistances(regionname)
#     createmaps(regionname)
#     GISsolar(gisregion=regionname, plotmasks=true)
#     GISwind(gisregion=regionname, plotmasks=true)
#     GIShydro(gisregion=regionname)
#     predictdemand(gisregion=regionname, sspscenario="ssp2-26", sspyear=2050, era_year=2018)
# end

end