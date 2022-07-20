using BinDeps

export download_datasets

function download_datasets(startfile=1)
    datafolder = getconfig("datafolder")
    if !isfile(datafolder)
        mkpath(datafolder)
    end

    WDPA_filename = "WDPA_WDOECM_search_feadb1b4f30799a6dc3ad0b16116d3530ec4a477898f4e10e097e2030e167128_shp"

    # tuples of (dataset_name, filename, url)
    datasets = [
        # https://globalwindatlas.info/api/gis/global/wind-speed/10
        # https://globalwindatlas.info/api/gis/global/wind-speed/50
        # https://globalwindatlas.info/api/gis/global/wind-speed/100
        # https://globalwindatlas.info/api/gis/global/wind-speed/150
        # https://globalwindatlas.info/api/gis/global/wind-speed/200
        ("Global Wind Atlas", "Global Wind Atlas v3 - 100m wind speed.tif",
            "https://chalmersuniversity.box.com/shared/static/wfr6dm9bcmj0mcqtdn0uimhg0otd4ht1.tif"),
        ("WDPA (protected areas):", "WDPA.zip",
            "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_search_feadb1b4f30799a6dc3ad0b16116d3530ec4a477898f4e10e097e2030e167128_shp.zip"),
            # "https://chalmersuniversity.box.com/shared/static/wn1kznvy7qh1issqcxdlsq64kgtkaayi.zip"),
        ("GADM (global administrative areas)", "gadm36.zip",
            "https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_shp.zip"),
        ("NUTS (administrative areas in Europe)", "nuts-2016-01m.shp.zip",
            "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2016-01m.shp.zip"),
        ("Land Cover", "Landcover - USGS MODIS.tif",
            "https://chalmersuniversity.box.com/shared/static/zm8zdkv1mz0wns0u77afkna95wbo7ggo.tif"),
        ("ETOPO1 Topography", "ETOPO1_Ice_c_geotiff.zip",
            "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/georeferenced_tiff/ETOPO1_Ice_c_geotiff.zip"),
        ("1-km downscaled spatial population scenarios (Gao et al), file 1/3", "temp_ssp1.tar",
            "http://www.cgd.ucar.edu/iam/modeling/data/SSP1_1km_netcdf.tar"),
        ("1-km downscaled spatial population scenarios (Gao et al), file 2/3", "temp_ssp2.tar",
            "http://www.cgd.ucar.edu/iam/modeling/data/SSP2_1km_netcdf.tar"),
        ("1-km downscaled spatial population scenarios (Gao et al), file 3/3", "temp_ssp3.tar",
            "http://www.cgd.ucar.edu/iam/modeling/data/SSP3_1km_netcdf.tar"),
        ("Global population & GDP", "temp_popgdp.zip",
            "https://osf.io/7jv3n/download"),
        ("WRI Global Power Plant Database", "WRI - Global Power Plant Database v1.10.zip",
            "https://chalmersuniversity.box.com/shared/static/ss6gycw7hf10e1fiicbxl5rgk08q5xr9.zip"),
            # switch to the official v1.2 link later (some plant cleanup is hardcoded for v1.1) 
            # "http://datasets.wri.org/dataset/540dcf46-f287-47ac-985d-269b04bea4c6/resource/c240ed2e-1190-4d7e-b1da-c66b72e08858/download/globalpowerplantdatabasev120"),
        ("Time zone shape file", "timezones-with-oceans.shapefile.zip",
            "https://github.com/evansiroky/timezone-boundary-builder/releases/download/2019b/timezones-with-oceans.shapefile.zip"),
        ("Average monthly wind speeds 1979-2019", "era5monthlywind.h5",
            "https://chalmersuniversity.box.com/shared/static/otvb5nq0lz5ntqx0e65kocos2afo7gas.h5"),
        ("Average monthly solar insolation 1979-2019", "era5monthlysolar.h5",
            "https://chalmersuniversity.box.com/shared/static/jlpspmp9ou96hk7xno46rf79wg3d5hqc.h5"),
        ("Various smaller datasets", "Various_smaller_datasets.zip",
            "https://chalmersuniversity.box.com/shared/static/w3pmx4xhgorgd6jejv23gn4ycsnza8s6.zip")
    ]

    for (i, datasetinfo) in enumerate(datasets)
        i < startfile && continue
        name, filename, url = datasetinfo

        println("\nDownloading dataset $i: $name")
        download_progressbar(url, joinpath(datafolder, filename))
    end

    println("\nDownloads complete.")

    for datasetinfo in datasets
        name, filename, url = datasetinfo
        foldername, extension = splitext(filename)
        fullpath = joinpath(datafolder, filename)

        if extension in [".zip", ".tar"] && isfile(fullpath)
            println("\nUnpacking archive: $filename")
            unpack(fullpath, joinpath(datafolder, foldername), extension)
        end
    end

    function renameWDPAfiles(WDPAfolder)
        for filename in readdir(WDPAfolder)
            newname = replace(filename, WDPA_filename => "WDPA-shapefile")
            if newname != filename 
                mv(joinpath(WDPAfolder, filename), joinpath(WDPAfolder, newname))
            end
        end
    end

    if startfile <= 2
        renameWDPAfiles(joinpath(datafolder, "WDPA"))
        for i = 0:2
            foldername = "WDPA-shapefile$i"
            println("\nUnpacking archive: $foldername.zip")
            unpack(joinpath(datafolder, "WDPA", "$foldername.zip"),
                    joinpath(datafolder, "WDPA", foldername), ".zip")
            renameWDPAfiles(joinpath(datafolder, "WDPA", foldername))
        end
    end

    if startfile <= 4
        filename = "NUTS_RG_01M_2016_4326_LEVL_3.shp.zip"
        println("\nUnpacking archive: $filename")
        unpack(joinpath(datafolder, "nuts-2016-01m.shp", filename), joinpath(datafolder, "nuts2016-level3"), ".zip")
        rm(joinpath(datafolder, "nuts-2016-01m.shp"), force=true, recursive=true)
    end

    println("\n\n\nCleaning up...")

    if startfile <= 6
        mv(joinpath(datafolder, "ETOPO1_Ice_c_geotiff", "ETOPO1_Ice_c_geotiff.tif"), joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif"))
        rm(joinpath(datafolder, "ETOPO1_Ice_c_geotiff"))
    end
    if startfile <= 9
        sspfolder = joinpath(datafolder, "SSP_1km")
        !isdir(sspfolder) && mkdir(sspfolder)
        for ssp = 1:3
            for y = 2010:10:2100
                sspfile = joinpath(datafolder, "temp_ssp$ssp", "SSP$(ssp)_1km", "ssp$(ssp)_total_$y.nc4")
                isfile(sspfile) && mv(sspfile, joinpath(sspfolder, "ssp$(ssp)_total_$y.nc4"))
            end
            rm(joinpath(datafolder, "temp_ssp$ssp"), force=true, recursive=true)
        end
    end
    if startfile <= 10
        mv(joinpath(datafolder, "temp_popgdp", "data"), joinpath(datafolder, "global_population_and_gdp"))
        rm(joinpath(datafolder, "temp_popgdp"))
    end
    if startfile <= 11
        mv(joinpath(datafolder, "WRI - Global Power Plant Database v1.10"), joinpath(datafolder, "tempWRI"))
        mv(joinpath(datafolder, "tempWRI", "WRI - Global Power Plant Database v1.10"),
            joinpath(datafolder, "WRI - Global Power Plant Database v1.10"))
        rm(joinpath(datafolder, "tempWRI"))
    end
    if startfile <= 15
        mixeddir = joinpath(datafolder, "Various_smaller_datasets")
        for file in readdir(mixeddir)
            mv(joinpath(mixeddir, file), joinpath(datafolder, file))
        end
        rm(mixeddir)
    end

    for datasetinfo in datasets
        name, filename, url = datasetinfo
        foldername, extension = splitext(filename)

        if extension in [".zip", ".tar"]
            fullpath = joinpath(datafolder, filename)
            isfile(fullpath) && rm(fullpath)
        end
    end
    println("\n\nDone.")
end

function download_progressbar(url::AbstractString, filename::AbstractString)
    println("Downloading to $filename...")
    curl = Base.find_curl()
    if curl == nothing
        # if curl is not on the system, fall back to a download without progress bar
        download(url, filename)
    else
        download_curl(curl, url, filename)
    end
end

function download_curl(curl_exe::AbstractString, url::AbstractString, filename::AbstractString)
    run(`$curl_exe --progress-bar -g -L -f --retry 5 -o $filename $url`, wait=true)
    return filename
end

function unpack(inputfilename, outputpath, extension)
    mkdir(outputpath)
    run(BinDeps.unpack_cmd(inputfilename, outputpath, extension, ""))
end
