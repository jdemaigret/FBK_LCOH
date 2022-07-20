using DelimitedFiles, JuMP, GLPK, NPFinancial, Ipopt, NLopt, StatsBase, Plots

lcohoptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat"

    :pv_density => 45,                  # Solar PV land use 45 Wp/m2 = 45 MWp/km2 (includes PV efficiency & module spacing, add latitude dependency later)

    :plant_area => .05,                 # area available for PV or CSP plants after the masks have been applied

    :onshore_density => 5,              # about 30% of existing farms have at least 5 W/m2, will become more common
    :offshore_density => 8,             # varies a lot in existing parks (4-18 W/m2)
                                        # For reference: 10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2
    :area_onshore => .08,               # area available for onshore wind power after the masks have been applied
    :area_offshore => .33,              # area available for offshore wind power after the masks have been applied
    #:distance_elec_access => 150,       # max distance to grid [km] (for wind classes of category B and offshore)


    #:distance_elec_access => 150,       # max distance to grid [km] (for solar classes of category B)
    :pv_plant_persons_per_km2 => 150,      # not too crowded, max X persons/km2 (both PV and CSP plants)
    :wind_persons_per_km2 => 150,            # not too crowded, max X persons/km2


    :pv_exclude_landtypes => [0,1,2,3,4,5,8],       # exclude water and forests. See codes in table below.
    :pv_protected_codes => [1,2,3,4,5,8],          # IUCN codes to be excluded as protected areas. See codes in table below.

    :max_depth => 40,                   # max depth for offshore wind [m]
    :min_shore_distance => 5,           # minimum distance to shore for offshore wind [km]

    :wind_exclude_landtypes => [0,11,13],    # exclude water, wetlands and urban areas. See codes in table below.
    :wind_protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use
    :rescale_to_wind_atlas => false,     # rescale the ERA5 time series to fit annual wind speed averages from the Global Wind Atlas

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :downsample_masks => 1     # set to 2 or higher to scale down mask sizes to avoid GPU errors in Makie plots for large regions
    #:classB_threshold => 0.001  # minimum share of pixels within distance_elec_access km that must have grid access
                                # for a pixel to be considered for solar class B.
)
    # Land types
    #     0      'Water'
    #     1      'Evergreen Needleleaf Forests'
    #     2      'Evergreen Broadleaf Forests'
    #     3      'Deciduous Needleleaf Forests'
    #     4      'Deciduous Broadleaf Forests'
    #     5      'Mixed Forests'
    #     6      'Closed Shrublands'
    #     7      'Open Shrublands'
    #     8      'Woody Savannas'
    #     9      'Savannas'
    #    10      'Grasslands'
    #    11      'Permanent Wetlands'
    #    12      'Croplands'
    #    13      'Urban'
    #    14      'Cropland/Natural'
    #    15      'Snow/Ice'
    #    16      'Barren'

    # Protected areas (IUCN codes from the WDPA)
    #    1      'Ia'                'Strict Nature Reserve'
    #    2      'Ib'                'Wilderness Area'
    #    3      'II'                'National Park'
    #    4      'III'               'Natural Monument'
    #    5      'IV'                'Habitat/Species Management'
    #    6      'V'                 'Protected Landscape/Seascape'
    #    7      'VI'                'Managed Resource Protected Area'
    #    8      'Not Reported'      'Not Reported'
    #    9      'Not Applicable'    'Not Applicable'
    #    10     'Not Assigned'      'Not Assigned'

mutable struct LcohOptions
    gisregion                   ::String
    filenamesuffix              ::String
    pv_density                  ::Float64           # W/m2
    plant_area                  ::Float64           # share [0-1]
    onshore_density             ::Float64           # W/m2
    offshore_density            ::Float64           # W/m2
    area_onshore                ::Float64           # share [0-1]
    area_offshore               ::Float64           # share [0-1]
    #distance_elec_access       ::Float64           # km
    pv_plant_persons_per_km2    ::Float64           # persons/km2
    wind_persons_per_km2        ::Float64           # persons/km2
    pv_exclude_landtypes        ::Vector{Int}
    pv_protected_codes          ::Vector{Int}
    max_depth                   ::Float64           # m
    min_shore_distance          ::Float64           # km
    wind_exclude_landtypes      ::Vector{Int}
    wind_protected_codes        ::Vector{Int}
    scenarioyear                ::String
    era_year                    ::Int
    rescale_to_wind_atlas       ::Bool
    res                         ::Float64
    erares                      ::Float64
    downsample_masks            ::Int
end


LcohOptions() = LcohOptions("","",0,0,0,0,0,0,0,0,[],[],0,0,[],[],"",0,false,0,0,0)

function LcohOptions(d::Dict{Symbol,Any})
    loptions = LcohOptions()
    for (key,val) in d
        setproperty!(loptions, key, val)
    end
    return loptions
end


function GISlcoh(; savetodisk=false, plotmasks=false, optionlist...)

    loptions = LcohOptions(merge(lcohoptions(), optionlist))
    @unpack   gisregion, era_year, filenamesuffix, downsample_masks = loptions

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(loptions)


    mask_pv_lcoh = create_solar_mask_lcoh(loptions, regions, popdens, land, protected, lonrange, latrange; plotmasks=plotmasks, downsample=1)


    meanGTI, solarGTI = read_solar_datasets_lcoh(loptions, lonrange, latrange)
    windatlas, meanwind, windspeed = read_wind_datasets_lcoh(loptions, lonrange, latrange)

    mask_onshore_lcoh, mask_offshore_lcoh = create_wind_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, plotmasks=plotmasks, downsample=downsample_masks)

    pv_area, onshore_wind_area, offshore_wind_area = calc_area(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh, lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh)
    max_capacity_pv, CF_pvplant = calc_pv_cf(loptions, meanGTI, solarGTI, pv_area, regions, offshoreregions, regionlist, mask_pv_lcoh, lonrange, latrange)

    mask_lcoh = create_lcoh_mask(mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh)

    #plotmasks == :onlymasks && return nothing
    max_capacity_onshore, max_capacity_offshore, CF_onshore, CF_offshore = calc_wind_cf(loptions, windatlas, meanwind, windspeed, onshore_wind_area, offshore_wind_area, regions, offshoreregions, regionlist, mask_onshore_lcoh, mask_offshore_lcoh, lonrange, latrange)

    opt_capacity_pv, opt_capacity_wind, opt_capacity_el, LCOH = optimize_LCOH(loptions, pv_area, onshore_wind_area, offshore_wind_area,
                                                                                CF_pvplant, max_capacity_pv, max_capacity_onshore,
                                                                                max_capacity_offshore, CF_onshore, CF_offshore,
                                                                                lonrange, latrange, regionlist, solarGTI)
    #xs, ys = plot_lcoh(loptions, regions, mask_pv_lcoh,  mask_onshore_lcoh, mask_lcoh, lonrange, latrange; plotmasks=plotmasks, downsample=1)

    # if savetodisk
    #     mkpath(in_datafolder("output"))
    #     matopen(in_datafolder("output", "GISdata_lcoh$(era_year)_$gisregion$filenamesuffix.mat"), "w", compress=true) do file
    #         write(file, "CF_pvplant",  CF_pvplant)
    #         write(file, "max_apacity_pv", max_capacity_pv)
    #         write(file, "max_capacity_onshore", max_capacity_onshore)
    #         write(file, "max_capacity_offshore", max_capacity_offshore)
    #         write(file, "CF_onshore", CF_onshore)
    #         write(file, "CF_offshore", CF_offshore)
    #     end
    # end
    return regions, xs, ys
    nothing
    #  capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

function read_solar_datasets_lcoh(loptions, lonrange, latrange)
    @unpack res, erares, era_year = loptions

    println("Reading ERA5 solar datasets...")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)

    @time meanGTI, solarGTI = h5open(in_datafolder("era5solar$era_year.h5"), "r") do file
        if length(eralonranges) == 1
            file["meanGTI"][eralonranges[1], eralatrange],
                file["GTI"][:,eralonranges[1], eralatrange]
                #file["meanDNI"][eralonranges[1], eralatrange],
                #file["DNI"][:,eralonranges[1], eralatrange]
        else
            [file["meanGTI"][eralonranges[1], eralatrange]; file["meanGTI"][eralonranges[2], eralatrange]],
                [file["GTI"][:, eralonranges[1], eralatrange] file["GTI"][:, eralonranges[2], eralatrange]]
                #[file["meanDNI"][eralonranges[1], eralatrange]; file["meanDNI"][eralonranges[2], eralatrange]],
                #[file["DNI"][:, eralonranges[1], eralatrange] file["DNI"][:, eralonranges[2], eralatrange]]
        end
    end
    return meanGTI, solarGTI #meanDNI, solarDNI
end

function read_wind_datasets_lcoh(loptions, lonrange, latrange)
    @unpack res, erares, era_year = loptions
    windatlas = getwindatlas()[lonrange,latrange]

    println("Reading ERA5 wind speeds and calculating capacity factors...")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)

    @time meanwind, windspeed = h5open(in_datafolder("era5wind$era_year.h5"), "r") do file
        if length(eralonranges) == 1
            file["meanwind"][eralonranges[1], eralatrange],
                file["wind"][:, eralonranges[1], eralatrange]
        else
            [file["meanwind"][eralonranges[1], eralatrange]; file["meanwind"][eralonranges[2], eralatrange]],
                [file["wind"][:, eralonranges[1], eralatrange] file["wind"][:, eralonranges[2], eralatrange]]
        end
    end
    return windatlas, meanwind, windspeed
end

#Nuova funzione
function create_solar_mask_lcoh(loptions, regions, popdens, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, pv_exclude_landtypes, pv_protected_codes, pv_plant_persons_per_km2,
            filenamesuffix = loptions

    println("Creating PV masks...")

    goodland = (regions .> 0)
    for i in pv_exclude_landtypes
        goodland[land .== i] .= false
    end
    pv_protected_area = zeros(Bool, size(protected))
    for i in pv_protected_codes
        pv_protected_area[protected .== i] .= true
    end

    km_per_degree = π*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)

    mask_pv_lcoh = (popdens .< pv_plant_persons_per_km2) .& goodland .& .!pv_protected_area

    # if plotmasks != false   # can == :onlymasks as well
    #     # drawmap(land)
    #     isregion = (regions .> 0) .& (regions .!= NOREGION)
    #
    #     # mask values refer to colors in ColorBrewer Set2_7:
    #     # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
    #     masks = zeros(Int16, size(regions))
    #     # masks[(masks .== 0) .& (popdens .> pv_plant_persons_per_km2)] .= 2
    #     # masks[(masks .== 0) .& pv_protected_area] .= 3
    #     # #masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
    #     # masks[(masks .== 0) .& .!goodland] .= 1
    #     # #masks[(masks .== 0) .& .!gridA .& gridB] .= 6
    #     masks[(masks .== 0) .& isregion] .= 1
    #     masks[regions .== 0] .= 0
    #     masks[regions .== NOREGION] .= 2
    #     legendtext = ["", "", "","","","","",""]
    #     maskmap("$(gisregion)_masks_solAr$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    # end

    return mask_pv_lcoh
end

#Nuova funzione
function create_wind_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, wind_exclude_landtypes, wind_protected_codes, wind_persons_per_km2,
                min_shore_distance, max_depth, filenamesuffix = loptions

    println("Creating wind masks...")

    goodland = (regions .> 0)
    for i in wind_exclude_landtypes
        goodland[land .== i] .= false
    end
    wind_protected_area = zeros(Bool, size(protected))
    for i in wind_protected_codes
        wind_protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind A
    #gridA = (gridaccess .> 0)

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = π*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    #gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    #mask_onshoreA = gridA .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    #mask_onshoreB = (gridB .& .!gridA) .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    mask_onshore_lcoh = (popdens .< wind_persons_per_km2) .& goodland .& .!wind_protected_area
    # shoreline mask for offshore wind
    disk = diskfilterkernel(min_shore_distance/km_per_degree/res)
    shore = (imfilter(regions .> 0, disk) .> 1e-6)

    # all mask conditions
    mask_offshore_lcoh = .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!wind_protected_area

    # if plotmasks != false   # can == :onlymasks as well
    #     # drawmap(land)
    #     isregion = (regions .> 0) .& (regions .!= NOREGION)
    #
    #     # mask values refer to colors in ColorBrewer Set2_7:
    #     # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
    #     masks = zeros(Int16, size(regions))
    #     masks[(masks .== 0) .& (popdens .> wind_persons_per_km2)] .= 2
    #     masks[(masks .== 0) .& wind_protected_area] .= 3
    #     #masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
    #     masks[(masks .== 0) .& .!goodland] .= 1
    #     #masks[(masks .== 0) .& .!gridA .& gridB] .= 8
    #     masks[(masks .== 0) .& isregion] .= 4
    #     masks[(masks .== 0) .& .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!wind_protected_area] .= 5
    #     #masks[(masks .== 0) .& .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!protected_area] .= 5
    #     masks[regions .== 0] .= 0
    #     masks[regions .== NOREGION] .= NOREGION
    #     legendtext = ["bad land type", "high population", "protected area", "onshore", "offshore", "", "", ""]
    #     maskmap("$(gisregion)_masks_wind$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    #
    #     # drawmap(.!shore .& (offshoreregions .> 0))
    #     # drawmap((topo .> -max_depth) .& (offshoreregions .> 0))
    #     # drawmap(.!protected_area .& (offshoreregions .> 0))
    #     # drawmap(mask_offshore)
    # end

    return mask_onshore_lcoh, mask_offshore_lcoh
end

#Nuova funzione
function plot_lcoh(loptions, regions, mask_pv_lcoh,  mask_onshore_lcoh, mask_lcoh, lonrange, latrange; plotmasks=false, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = rand(Float16, size(regions))
        masks2 = ones(size(regions))*2

        masks2[mask_lcoh .== 1] .= 3
        # masks2[mask_pv_lcoh .== 0] .= 3        # Maschera dove non si puo' fare LCOH
        # masks2[mask_onshore_lcoh .== 0] .= 4   # Maschera dove non si puo' fare LCOH
        masks[regions .== 0] .= 0               # Ploting LCOH zona mare
        masks[regions .== NOREGION] .= 1        # Ploting LCOH zona NOREGION
        masks2[regions .== 0] .= 0              # Maschera zona mare
        masks2[regions .== NOREGION] .= 1       # Maschera zona NOREGION
        legendtext = ["", "", "","","","","",""]

        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]
        #nreg = length(regionlist)

        # colors = [
        #     RGB([174,140,114]/255...),  # bad land type
        #     RGB([214,64,64]/255...),    # high population
        #     RGB([255,100,255]/255...),  # protected area
        #     RGB([120,170,80]/255...),   # no grid
        #     RGB([255,217,47]/255...),   # solar plant A
        #     RGB([255,150,0]/255...),    # solar plant B
        #     RGB([140,156,209]/255...),  # wind plant A
        #     RGB([110,136,169]/255...),  # wind plant B
        # ]

        println("Mapping colors to regions (avoid same color in adjacent regions)...")
        #connected = connectedregions(regions, nreg)
        # colorindices = (nreg > 7) ? greedycolor(connected, 1:7, 1:nreg, randseed=randseed) : collect(1:nreg)
        # colors = colorschemes[Symbol("Set2_$nreg")].colors[colorindices]
        #onshorecolors = [RGB(0.3,0.3,0.45); colors; RGB(0.4,0.4,0.4)]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        println("\nOnshore map...")

        # surface(xs, ys; color=float.(masks), colormap= :blues, colorrange= (0,),
        #                         shading=false, show_axis=false, scale_plot=false, interpolate=false)
         scene = createmap("$(gisregion)_masks_lcoh$filenamesuffix", masks, masks2, [], lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
        #maskmap("$(gisregion)_masks_lcoh$filenamesuffix", masks, legendtext, lonrange, latrange; legend=false, downsample=downsample)
    end
    return xs, ys
end

function create_lcoh_mask(mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh)

    mask_lcoh =  mask_pv_lcoh .| mask_onshore_lcoh #&. mask_offshore_lcoh

    # println(mask_pv_lcoh[1:10,1:10])
    # println(mask_onshore_lcoh[1:10,1:10])
    # println(mask_lcoh[1:10,1:10])
    return mask_lcoh
end

#Funzione nuova
function calc_pv_cf(loptions, meanGTI, solarGTI, pv_area, regions, offshoreregions, regionlist,
                mask_pv_lcoh, lonrange, latrange)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    println("Calculating PV hourly capacity factors for each pixel (era5 resolution)")

    @unpack era_year, res, erares, pv_density, plant_area = loptions

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    max_capacity_pv = zeros(nlons,nlats)
    CF_pvplant = zeros(yearlength,nlons,nlats)
    #max_capacity_pv = zeros(size(rowrange,1),size(colrange,1)) #SIZE FORSE SBAGLIATE, prova rowrange*nlons,colrange*nlats
    #CF_pvplant = zeros(yearlength,size(rowrange,1),size(colrange,1)) #SIZE FORSE SBAGLIATE prova rowrange*nlons,colrange*nlats
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
        #colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            ##meanGTI[i,j] == 0 && continue #verificare se non necessario per escludere alcune zone
            GTI = solarGTI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            #rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]
            # for c in colrange, r in rowrange
            #     (c == 0 || r == 0) && continue
            #     reg = regions[r,c]
            #     (reg == 0 || reg == NOREGION) && continue
            #
            #     area = cellarea[c]
            #
            #     # @views is needed to make sure increment_windCF!() works with matrix slices
            #     # also faster since it avoids making copies
            #     @views
            #     if mask_pv[r,c] > 0
            #         max_capacity_pv[r,c] = 1/1000 * pv_density * plant_area * area
            #         CF_pvplant[:,r,c] = GTI
            #     end
            #
            # end
            max_capacity_pv[i,j] = 1/1000 * pv_density * pv_area[i,j]
            CF_pvplant[:,i,j] = GTI
        end
        next!(updateprogress)
    end

    return  max_capacity_pv, CF_pvplant
end

function assign_windCF!(cf::AbstractVector{<:AbstractFloat}, speed_or_cf::AbstractVector{<:AbstractFloat}, factor, rescale::Bool)
    if rescale
        @inbounds for i = 1:length(cf)
            cf[i] = speed2capacityfactor(speed_or_cf[i]*factor)
        end
    else
        @inbounds for i = 1:length(cf)
            cf[i] += speed_or_cf[i]
        end
    end
end

#Funzione Nuova
function calc_wind_cf(loptions, windatlas, meanwind, windspeed, onshore_wind_area, offshore_wind_area, regions, offshoreregions, regionlist,
            mask_onshore_lcoh, mask_offshore_lcoh, lonrange, latrange)

    println("Calculating wind hourly capacity factors for each pixel (era5 resolution)")

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    _, _, meanwind_allyears, _ = annualwindindex(loptions)
    @assert size(meanwind_allyears) == size(meanwind)

    @unpack rescale_to_wind_atlas, res, erares, onshore_density, offshore_density = loptions

    numreg = length(regionlist)
    #nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)

    max_capacity_onshore = zeros(nlons,nlats)
    max_capacity_offshore = zeros(nlons,nlats)
    CF_wind_onshore = zeros(yearlength,nlons,nlats)
    CF_wind_offshore = zeros(yearlength,nlons,nlats)

    if rescale_to_wind_atlas
        println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
        println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    end

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    @inbounds for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanwind[i,j] == 0 && continue
            wind = rescale_to_wind_atlas ? windspeed[:, i, j] : speed2capacityfactor.(windspeed[:, i, j])
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            #for c in colrange, r in rowrange
            #    (c == 0 || r == 0) && continue
                #@views if  mask_onshore_lcoh[r,c] > 0
            #        assign_windCF!(CF_wind_onshore[:,i,j], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                #elseif mask_offshore_lcoh[r,c] > 0
            #        assign_windCF!(CF_wind_offshore[:,i,j], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                #end
            #end
            CF_wind_onshore[:,i,j] = wind
            CF_wind_offshore[:,i,j] = wind
            #assign_windCF!(CF_wind_onshore[:,i,j], wind, 1, rescale_to_wind_atlas)
            #assign_windCF!(CF_wind_offshore[:,i,j], wind, 1, rescale_to_wind_atlas)
            max_capacity_onshore[i,j] = 1/1000 * onshore_density * onshore_wind_area[i,j]
            max_capacity_offshore[i,j] = 1/1000 * offshore_density * offshore_wind_area[i,j]

        end
        next!(updateprogress)
    end
    #println(speed2capacityfactor.(windspeed[1:10, 1, 1]))
    #println(CF_wind_onshore[1:10,1,1])
    #println(size(CF_wind_onshore[:,:,:]))
    #println( CF_wind_offshore[3000,1,1])
    return max_capacity_onshore, max_capacity_offshore, CF_wind_onshore, CF_wind_offshore
end

#Funzione nuova
function calc_area(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh, lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)
    #onshoreclass, offshoreclass = makewindclasses(loptions, windatlas)

    println("Calculating usable area for each pixel (using high res datasets to be used for low res computations)")
    # println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, res, erares, plant_area, rescale_to_wind_atlas, onshore_density, area_onshore, offshore_density, area_offshore = loptions

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)
    #nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)

    onshoremap = zeros(Int16, size(regions))
    offshoremap = zeros(Int16, size(regions))

    if rescale_to_wind_atlas
        println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
        println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    end

    #added for later -> #capacity_pv = zeros(nlons,nlats)
    #CF_pvplant = zeros(yearlength,nlons,nlats)
    pv_area = zeros(nlons,nlats)
    onshore_wind_area = zeros(nlons,nlats)
    offshore_wind_area = zeros(nlons,nlats)
    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        #println(size(colrange))
        for i = 1:nlons
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                (reg == 0 || reg == NOREGION) && continue

                area = cellarea[c]

                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if mask_pv_lcoh[r,c] > 0
                    pv_area[i,j] += plant_area * area
                end
                @views if mask_onshore_lcoh[r,c]>0
                    onshore_wind_area[i,j]  += area_onshore  * area
                end
                @views if  mask_offshore_lcoh[r,c]>0
                    offshore_wind_area[i,j] += area_offshore * area
                end

            end
        end
        next!(updateprogress)
    end

    return pv_area, onshore_wind_area, offshore_wind_area
end


#Funzione nuova
function optimize_LCOH(loptions, pv_area, onshore_wind_area, offshore_wind_area,
                        CF_pvplant, max_capacity_pv, max_capacity_onshore,
                        max_capacity_offshore, CF_onshore, CF_offshore,
                        lonrange, latrange, regionlist, solarGTI)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    println("Calculating optimal values (era5 resolution)")

    @unpack era_year, res, erares = loptions

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    CAPEX_el = 750;                      # USD/kW
    CAPEX_wind = 1500;                   # USD/kW
    CAPEX_pv = 1000;                     # USD/kW
    WACC = 0.08;                         # %

    opt_capacity_pv   = zeros(nlats, nlons)
    opt_capacity_wind = zeros(nlats, nlons)
    opt_capacity_el   = zeros(nlats, nlons)
    LCOH              = zeros(nlats, nlons)

    println(size(CF_pvplant))

    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
            for i = 1:nlons
                 opt_capacity_pv[i,j], opt_capacity_wind[i,j], opt_capacity_el[i,j], LCOH[i,j] = opt_model(loptions, pv_area[i,j], CF_pvplant[:,i,j],
                                                                                                            onshore_wind_area[i,j], CF_onshore[:,i,j],
                                                                                                            CAPEX_pv, CAPEX_wind, CAPEX_el, WACC)
            end
            next!(updateprogress)
    end

    return opt_capacity_pv, opt_capacity_wind, opt_capacity_el, LCOH
end

#Funzione nuova
function opt_model(loptions, pv_area, CF_pvplant,
                            onshore_wind_area, CF_onshore,
                            CAPEX_pv, CAPEX_wind, CAPEX_el, WACC)

    @unpack pv_density, plant_area, onshore_density, area_onshore  = loptions

    lifetime_stack = 80000;              # Hours
    lifetime_el = 15;                    # Years
    lifetime_wind = 25;                  # Years
    lifetime_pv = 25;                    # Years
    OPEX_el = 0.03;                      # %CAPEX
    OPEX_wind = 0.015;                   # %CAPEX
    OPEX_pv = 0.015;
    eff_el = 0.7;                       #% only electrolyzer
    min_CF_electrolyzer = 0.3;
    HHV_H2 = 39.4;
    H2_energy_consumption = HHV_H2/eff_el;

    #setting up model
    model = Model(Ipopt.Optimizer)

    #defining variables constraints
    @variable(model, 1 <= x <= 1000e+6);       #capacity pv kW
    @variable(model, 1 <= y <= 1000e+6);       #capacity wind kW
    @variable(model, 1 <= z <= 1000e+6);       #capacity electrolyzer kW

    max_cap_pv   = plant_area*(pv_area)*pv_density*1000;
    max_cap_wind = area_onshore*(onshore_wind_area)*onshore_density*1000;
    @constraint(model, x <= max_cap_pv);
    @constraint(model, y <= max_cap_wind);
    @constraint(model, x + y - z >= 0)                          #La somma delle cap. di PV e wind deve essere >= della cap. electrolyzer
    @NLexpression(model, CFe[n = 1:8760], (y*CF_onshore[n] + x*CF_pvplant[n])/z);
    @NLconstraint(model, [n = 1:8760], CFe[n] >= 0);
    @NLconstraint(model, [n = 1:8760], CFe[n] <= 1);

    function f_CFe(CFe::T...) where {T<:Real}
        a = zeros(T, length(CFe))
        for n = 1:length(CFe)
            if CFe[n] >= min_CF_electrolyzer
                a[n] = CFe[n]
            else
                a[n] = 0
            end
        end
        return a
    end

    function f_annual_hydrogen_prod_1(z,CFe::T...) where {T<:Real}
        a = zeros(T, length(CFe))
        a = f_CFe(CFe...)
        b = sum(a.*z/H2_energy_consumption)
        return b
    end

    f_annual_invest(x,y,z)=-pmt(WACC,lifetime_pv,CAPEX_pv*x)-pmt(WACC,lifetime_wind,CAPEX_wind*y)-pmt(WACC,lifetime_el,CAPEX_el*z)
                                +CAPEX_pv*x*OPEX_pv+CAPEX_wind*y*OPEX_wind+CAPEX_el*z*OPEX_el;

    f_LCOH(x,y,z,CFe...)=f_annual_invest(x,y,z)/f_annual_hydrogen_prod_1(z,CFe...);

    register(model, :f_CFe, 8760, f_CFe, autodiff=true);
    register(model, :f_annual_hydrogen_prod_1, 8763, f_annual_hydrogen_prod_1, autodiff=true);
    register(model, :f_annual_invest, 3, f_annual_invest, autodiff=true);
    register(model, :f_LCOH, 8763, f_LCOH, autodiff=true);

    @NLobjective(model, Min, f_LCOH(x,y,z,CFe...));#non-linear

    JuMP.optimize!(model) # JuMP v0.21

    #@show termination_status(model)
    opt_capacity_pv     = value(x)
    opt_capacity_wind   = value(y)
    opt_capacity_el     = value(z)
    LCOH                = objective_value(model)

    return opt_capacity_pv, opt_capacity_wind, opt_capacity_el, LCOH
end
