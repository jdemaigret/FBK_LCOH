solaroptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat"

    :pv_density => 45,                  # Solar PV land use 45 Wp/m2 = 45 MWp/km2 (includes PV efficiency & module spacing, add latitude dependency later)
    :csp_density => 35,                 # CSP land use 35 W/m2

    :pvroof_area => .05,                # area available for rooftop PV after the masks have been applied
    :plant_area => .05,                 # area available for PV or CSP plants after the masks have been applied

    :distance_elec_access => 150,       # max distance to grid [km] (for solar classes of category B)
    :plant_persons_per_km2 => 150,      # not too crowded, max X persons/km2 (both PV and CSP plants)
    :pvroof_persons_per_km2 => 200,     # only in populated areas, so AT LEAST x persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :exclude_landtypes => [0,1,2,3,4,5,8],       # exclude water and forests. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :pvclasses_min => [0.08,0.14,0.18,0.22,0.26],   # lower bound on annual PV capacity factor for class X    [0:0.01:0.49;]
    :pvclasses_max => [0.14,0.18,0.22,0.26,1.00],   # upper bound on annual PV capacity factor for class X    [0.01:0.01:0.50;]
    :cspclasses_min => [0.10,0.18,0.24,0.28,0.32],  # lower bound on annual CSP capacity factor for class X
    :cspclasses_max => [0.18,0.24,0.28,0.32,1.00],  # upper bound on annual CSP capacity factor for class X

    :downsample_masks => 1,     # set to 2 or higher to scale down mask sizes to avoid GPU errors in Makie plots for large regions
    :classB_threshold => 0.001  # minimum share of pixels within distance_elec_access km that must have grid access
                                # for a pixel to be considered for solar class B.
)

windoptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat"

    :onshore_density => 5,              # about 30% of existing farms have at least 5 W/m2, will become more common
    :offshore_density => 8,             # varies a lot in existing parks (4-18 W/m2)
                                        # For reference: 10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2
    :area_onshore => .08,               # area available for onshore wind power after the masks have been applied
    :area_offshore => .33,              # area available for offshore wind power after the masks have been applied

    :distance_elec_access => 150,       # max distance to grid [km] (for wind classes of category B and offshore)
    :persons_per_km2 => 150,            # not too crowded, max X persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :max_depth => 40,                   # max depth for offshore wind [m]
    :min_shore_distance => 5,           # minimum distance to shore for offshore wind [km]
    :exclude_landtypes => [0,11,13],    # exclude water, wetlands and urban areas. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use
    :rescale_to_wind_atlas => true,     # rescale the ERA5 time series to fit annual wind speed averages from the Global Wind Atlas

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :onshoreclasses_min => [2,5,6,7,8],     # lower bound on annual onshore wind speeds for class X    [0:0.25:12.25;]
    :onshoreclasses_max => [5,6,7,8,99],    # upper bound on annual onshore wind speeds for class X    [0.25:0.25:12.5;]
    :offshoreclasses_min => [3,6,7,8,9],    # lower bound on annual offshore wind speeds for class X
    :offshoreclasses_max => [6,7,8,9,99],   # upper bound on annual offshore wind speeds for class X

    :downsample_masks => 1,     # set to 2 or higher to scale down mask sizes to avoid GPU errors in Makie plots for large regions
    :classB_threshold => 0.001  # minimum share of pixels within distance_elec_access km that must have grid access
                                # for a pixel to be considered for wind class B.
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

mutable struct SolarOptions
    gisregion               ::String
    filenamesuffix          ::String
    pv_density              ::Float64           # W/m2
    csp_density             ::Float64           # W/m2
    pvroof_area             ::Float64           # share [0-1]
    plant_area              ::Float64           # share [0-1]
    distance_elec_access    ::Float64           # km
    plant_persons_per_km2   ::Float64           # persons/km2
    pvroof_persons_per_km2  ::Float64           # persons/km2
    exclude_landtypes       ::Vector{Int}
    protected_codes         ::Vector{Int}
    scenarioyear            ::String
    era_year                ::Int
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    pvclasses_min           ::Vector{Float64}
    pvclasses_max           ::Vector{Float64}
    cspclasses_min          ::Vector{Float64}
    cspclasses_max          ::Vector{Float64}
    downsample_masks        ::Int
    classB_threshold        ::Float64
end

mutable struct WindOptions
    gisregion               ::String
    filenamesuffix          ::String
    onshore_density         ::Float64           # W/m2
    offshore_density        ::Float64           # W/m2
    area_onshore            ::Float64           # share [0-1]
    area_offshore           ::Float64           # share [0-1]
    distance_elec_access    ::Float64           # km
    persons_per_km2         ::Float64           # persons/km2
    max_depth               ::Float64           # m
    min_shore_distance      ::Float64           # km
    exclude_landtypes       ::Vector{Int}
    protected_codes         ::Vector{Int}
    scenarioyear            ::String
    era_year                ::Int
    rescale_to_wind_atlas   ::Bool
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    onshoreclasses_min      ::Vector{Float64}
    onshoreclasses_max      ::Vector{Float64}
    offshoreclasses_min     ::Vector{Float64}
    offshoreclasses_max     ::Vector{Float64}
    downsample_masks        ::Int
    classB_threshold        ::Float64
end

SolarOptions() = SolarOptions("","",0,0,0,0,0,0,0,[],[],"",0,0,0,[],[],[],[],0,0.0)
WindOptions() = WindOptions("","",0,0,0,0,0,0,0,0,[],[],"",0,false,0,0,[],[],[],[],0,0.0)

function SolarOptions(Sd::Dict{Symbol,Any})
    Soptions = SolarOptions()
    for (key,val) in Sd
        setproperty!(Soptions, key, val)
    end
    return Soptions
end

function WindOptions(Wd::Dict{Symbol,Any})
    Woption = WindOptions()
    for (key,val) in Wd
        setproperty!(Woption, key, val)
    end
    return Woption
end

function GISsolar(; savetodisk=true, plotmasks=false, optionlist...)

    # IMPORTANT!! The function makesolarera5() uses ERA5 solar datasets to
    # calculate Global Tilted Irradiance (GTI) for solar PV and Direct Normal
    # Irradiance (DNI) for CSP solar towers. If we want CSP parabolic troughs
    # then we need to add that dataset in makesolarera5() (and capacity factors
    # will become somewhat lower).

    # NOTE ON SOLAR UNIT: the solar irradiance sets are in kW/m2. Since the
    # irradiance value used to represent "standard testing conditions" for PV
    # is 1000 W/m2, the solar datasets also directly give the capacity factor.
    # Actual insolation can occasionally go above 1000 W/m2.

    # Ideally, we should make direct assumptions of PV module efficiency as a
    # function of air temperature (see e.g. Bett & Thornton appendix A2), but
    # for now efficiency is included in our assumption of :pv_density. Wind
    # speed also affects PV module temperature and efficiency. However, the
    # uncertainties in :pv_density and :plant_area are so large that efficiency
    # variations as a function of temperature don't matter.

    Soptions = SolarOptions(merge(solaroptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix, pv_density, csp_density, downsample_masks = Soptions

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(Soptions)

    mask_rooftop, mask_plantA, mask_plantB =
        create_solar_masks(Soptions, regions, gridaccess, popdens, land, protected, lonrange, latrange,
                            plotmasks=plotmasks, downsample=downsample_masks)

    plotmasks == :onlymasks && return nothing

    meanGTI, solarGTI, meanDNI, solarDNI = read_solar_datasets(Soptions, lonrange, latrange)

    CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB, solar_overlap_areaA, solar_overlap_areaB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB =
        calc_solar_vars(Soption, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    if savetodisk
        mkpath(in_datafolder("output"))
        matopen(in_datafolder("output", "GISdata_solar$(era_year)_$gisregion$filenamesuffix.mat"), "w", compress=true) do file
            write(file, "CFtime_pvrooftop", CF_pvrooftop)
            write(file, "CFtime_pvplantA", CF_pvplantA)
            write(file, "CFtime_pvplantB", CF_pvplantB)
            write(file, "CFtime_cspplantA", CF_cspplantA)
            write(file, "CFtime_cspplantB", CF_cspplantB)
            write(file, "capacity_pvrooftop", capacity_pvrooftop)
            write(file, "capacity_pvplantA", capacity_pvplantA)
            write(file, "capacity_pvplantB", capacity_pvplantB)
            write(file, "capacity_cspplantA", capacity_cspplantA)
            write(file, "capacity_cspplantB", capacity_cspplantB)
            write(file, "solar_overlap_areaA", solar_overlap_areaA)
            write(file, "solar_overlap_areaB", solar_overlap_areaB)
            write(file, "pv_density", pv_density)
            write(file, "csp_density", csp_density)
        end
    end

    nothing
    # return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
    #         capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

function GISwind(; savetodisk=true, plotmasks=false, optionlist...)
    Woption = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix, downsample_masks = Woption

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(Woption)

    mask_onshoreA, mask_onshoreB, mask_offshore =
        create_wind_masks(Woption, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange,
                            plotmasks=plotmasks, downsample=downsample_masks)

    plotmasks == :onlymasks && return nothing

    windatlas, meanwind, windspeed = read_wind_datasets(Woption, lonrange, latrange)

    windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore =
        calc_wind_vars(Woption, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    if savetodisk
        mkpath(in_datafolder("output"))
        matopen(in_datafolder("output", "GISdata_wind$(era_year)_$gisregion$filenamesuffix.mat"), "w", compress=true) do file
            write(file, "CFtime_windonshoreA", windCF_onshoreA)
            write(file, "CFtime_windonshoreB", windCF_onshoreB)
            write(file, "CFtime_windoffshore", windCF_offshore)
            write(file, "capacity_onshoreA", capacity_onshoreA)
            write(file, "capacity_onshoreB", capacity_onshoreB)
            write(file, "capacity_offshore", capacity_offshore)
        end
    end

    nothing
    # return windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore
end

function read_solar_datasets(Soption, lonrange, latrange)
    @unpack res, erares, era_year = Soption

    println("Reading ERA5 solar datasets...")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)

    @time meanGTI, solarGTI, meanDNI, solarDNI = h5open(in_datafolder("era5solar$era_year.h5"), "r") do file
        if length(eralonranges) == 1
            file["meanGTI"][eralonranges[1], eralatrange],
                file["GTI"][:,eralonranges[1], eralatrange],
                file["meanDNI"][eralonranges[1], eralatrange],
                file["DNI"][:,eralonranges[1], eralatrange]
        else
            [file["meanGTI"][eralonranges[1], eralatrange]; file["meanGTI"][eralonranges[2], eralatrange]],
                [file["GTI"][:, eralonranges[1], eralatrange] file["GTI"][:, eralonranges[2], eralatrange]],
                [file["meanDNI"][eralonranges[1], eralatrange]; file["meanDNI"][eralonranges[2], eralatrange]],
                [file["DNI"][:, eralonranges[1], eralatrange] file["DNI"][:, eralonranges[2], eralatrange]]
        end
    end
    return meanGTI, solarGTI, meanDNI, solarDNI
end

function read_datasets(Woption)
    @unpack res, gisregion, scenarioyear = Woption

    println("\nReading auxiliary datasets...")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    gridaccess = JLD.load(in_datafolder("gridaccess_$scenarioyear.jld"), "gridaccess")[lonrange,latrange]
    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")[lonrange,latrange]
    topo = JLD.load(in_datafolder("topography.jld"), "topography")[lonrange,latrange]
    land = JLD.load(in_datafolder("landcover.jld"), "landcover")[lonrange,latrange]
    protected = JLD.load(in_datafolder("protected.jld"), "protected")[lonrange,latrange]

    popdens = pop ./ cellarea'

    # drawmap(log.(pop))

    return regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange
end

function read_wind_datasets(Woption, lonrange, latrange)
    @unpack res, erares, era_year = Woption
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
function create_pv_mask(Soption, regions, popdens, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, plant_persons_per_km2,
            filenamesuffix = Soption

    println("Creating masks...")

    goodland = (regions .> 0)
    for i in exclude_landtypes
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in protected_codes
        protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind A
    #gridA = (gridaccess .> 0)

    # Pixels with electricity access for onshore wind B and offshore wind
    ##km_per_degree = π*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    #gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    #mask_rooftop = gridA .& (popdens .> pvroof_persons_per_km2) .& .!protected_area
    #mask_plantA = gridA .& (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area
    #mask_plantB = (gridB .& .!gridA) .& (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area
    mask_pv = (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = zeros(Int16, size(regions))
        masks[(masks .== 0) .& (popdens .> plant_persons_per_km2)] .= 2
        masks[(masks .== 0) .& protected_area] .= 3
        #masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
        masks[(masks .== 0) .& .!goodland] .= 1
        #masks[(masks .== 0) .& .!gridA .& gridB] .= 6
        #masks[(masks .== 0) .& isregion] .= 5
        masks[regions .== 0] .= 0
        masks[regions .== NOREGION] .= NOREGION
        legendtext = ["bad land type", "high population", "protected area", "", ""]
        maskmap("$(gisregion)_masks_solar$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    end

    return mask_pv
end

#Nuova funzione
function create_wind_maskS(Woption, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, distance_elec_access, persons_per_km2,
                min_shore_distance, max_depth, classB_threshold, filenamesuffix = Woption

    println("Creating masks...")

    goodland = (regions .> 0)
    for i in exclude_landtypes
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in protected_codes
        protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind A
    #gridA = (gridaccess .> 0)

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    #gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    #mask_onshoreA = gridA .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    #mask_onshoreB = (gridB .& .!gridA) .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    mask_onshore = (popdens .< persons_per_km2) .& goodland .& .!protected_area
    # shoreline mask for offshore wind
    disk = diskfilterkernel(min_shore_distance/km_per_degree/res)
    shore = (imfilter(regions .> 0, disk) .> 1e-6)

    # all mask conditions
    mask_offshore = .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!protected_area

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = zeros(Int16, size(regions))
        masks[(masks .== 0) .& (popdens .> persons_per_km2)] .= 2
        masks[(masks .== 0) .& protected_area] .= 3
        #masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
        masks[(masks .== 0) .& .!goodland] .= 1
        #masks[(masks .== 0) .& .!gridA .& gridB] .= 8
        masks[(masks .== 0) .& isregion] .= 4
        masks[regions .== 0] .= 0
        masks[regions .== NOREGION] .= NOREGION
        legendtext = ["bad land type", "high population", "protected area", "onshore"]
        maskmap("$(gisregion)_masks_wind$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)

        # drawmap(.!shore .& (offshoreregions .> 0))
        # drawmap((topo .> -max_depth) .& (offshoreregions .> 0))
        # drawmap(.!protected_area .& (offshoreregions .> 0))
        # drawmap(mask_offshore)
    end

    return mask_onshore, mask_offshore
end

#=function create_solar_masks(Soption, regions, gridaccess, popdens, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, distance_elec_access, plant_persons_per_km2,
            pvroof_persons_per_km2, classB_threshold, filenamesuffix = Soption

    println("Creating masks...")

    goodland = (regions .> 0)
    for i in exclude_landtypes
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in protected_codes
        protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind A
    gridA = (gridaccess .> 0)

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    mask_rooftop = gridA .& (popdens .> pvroof_persons_per_km2) .& .!protected_area
    mask_plantA = gridA .& (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area
    mask_plantB = (gridB .& .!gridA) .& (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = zeros(Int16, size(regions))
        masks[(masks .== 0) .& (popdens .> plant_persons_per_km2)] .= 2
        masks[(masks .== 0) .& protected_area] .= 3
        masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
        masks[(masks .== 0) .& .!goodland] .= 1
        masks[(masks .== 0) .& .!gridA .& gridB] .= 6
        masks[(masks .== 0) .& isregion] .= 5
        masks[regions .== 0] .= 0
        masks[regions .== NOREGION] .= NOREGION
        legendtext = ["bad land type", "high population", "protected area", "no grid", "solar plant A", "solar plant B", "", ""]
        maskmap("$(gisregion)_masks_solar$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    end

    return mask_rooftop, mask_plantA, mask_plantB
end=#

# 0 - 29 m/s
const windparkcurve = [
    0.0, 0.0014, 0.0071, 0.0229, 0.0545, 0.1067, 0.1831, 0.2850, 0.4085, 0.5434,
    0.6744, 0.7847, 0.8614, 0.9048, 0.9266, 0.9353, 0.9373, 0.9375, 0.9375, 0.9375,
    0.9375, 0.9375, 0.9375, 0.9311, 0.8683, 0.6416, 0.2948, 0.0688, 0.0063, 0.0
]

function speed2capacityfactor(windspeed)
    if windspeed >= 29 || windspeed < 0
        return 0.0
    end
    fw = floor(Int, windspeed)
    frac = windspeed - fw
    return (1-frac).*windparkcurve[fw+1] + frac.*windparkcurve[ceil(Int, windspeed)+1]
end

function increment_solarCF!(cf::AbstractVector{<:AbstractFloat}, solardata::AbstractVector{<:AbstractFloat})
    @inbounds for i = 1:length(cf)
        cf[i] += solardata[i]
    end
end

function increment_windCF!(cf::AbstractVector{<:AbstractFloat}, speed_or_cf::AbstractVector{<:AbstractFloat}, factor, rescale::Bool)
    if rescale
        @inbounds for i = 1:length(cf)
            cf[i] += speed2capacityfactor(speed_or_cf[i]*factor)
        end
    else
        @inbounds for i = 1:length(cf)
            cf[i] += speed_or_cf[i]
        end
    end
end

function makesolarclasses(Soption, meanGTI, meanDNI)
    println("Allocating pixels to classes using ERA5 annual means...")

    @unpack pvclasses_min, pvclasses_max, cspclasses_min, cspclasses_max = Soption

    pvclass = getclasses(meanGTI, pvclasses_min, pvclasses_max)
    cspclass = getclasses(meanDNI, cspclasses_min, cspclasses_max)

    return pvclass, cspclass
end

function makewindclasses(Woption, windatlas)
    println("Allocating pixels to classes using the Global Wind Atlas...")
    # println("CHANGE TO CAPACITY FACTOR LATER!")

    @unpack onshoreclasses_min, onshoreclasses_max, offshoreclasses_min, offshoreclasses_max = Woption

    onshoreclass = getclasses(windatlas, onshoreclasses_min, onshoreclasses_max)
    offshoreclass = getclasses(windatlas, offshoreclasses_min, offshoreclasses_max)

    return onshoreclass, offshoreclass
end

function getclasses0(windatlas, classes_min, classes_max)
    class = zeros(Int16, size(windatlas))
    for c = 1:length(classes_min)
        class[(windatlas .>= classes_min[c]) .& (windatlas .< classes_max[c])] .= c
    end
    return class
end

function getclasses(windatlas, classes_min, classes_max)
    class = similar(windatlas, Int16)
    nclasses = length(classes_min)
    @inbounds for i = 1:length(windatlas)
        wa = windatlas[i]
        if wa < classes_min[1] || wa > classes_max[end]
            class[i] = 0
            continue
        end
        for c = 1:nclasses
            if wa < classes_max[c]
                class[i] = c
                break
            end
        end
    end
    return class
end

function calc_solar_vars(Soption, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    pvclass, cspclass = makesolarclasses(Soption, meanGTI, meanDNI)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Soption, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and class...")
    # println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, pvclasses_min, cspclasses_min, res, erares, pv_density, csp_density, pvroof_area, plant_area = Soption

    numreg = length(regionlist)
    npvclasses, ncspclasses = length(pvclasses_min), length(cspclasses_min)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    capacity_pvrooftop = zeros(numreg,npvclasses)
    capacity_pvplantA = zeros(numreg,npvclasses)
    capacity_pvplantB = zeros(numreg,npvclasses)
    capacity_cspplantA = zeros(numreg,ncspclasses)
    capacity_cspplantB = zeros(numreg,ncspclasses)
    CF_pvrooftop = zeros(yearlength,numreg,npvclasses)
    CF_pvplantA = zeros(yearlength,numreg,npvclasses)
    CF_pvplantB = zeros(yearlength,numreg,npvclasses)
    CF_cspplantA = zeros(yearlength,numreg,ncspclasses)
    CF_cspplantB = zeros(yearlength,numreg,ncspclasses)
    count_pvrooftop = zeros(Int,numreg,npvclasses)
    count_pvplantA = zeros(Int,numreg,npvclasses)
    count_pvplantB = zeros(Int,numreg,npvclasses)
    count_cspplantA = zeros(Int,numreg,ncspclasses)
    count_cspplantB = zeros(Int,numreg,ncspclasses)
    solar_overlap_areaA = zeros(numreg,npvclasses,ncspclasses)
    solar_overlap_areaB = zeros(numreg,npvclasses,ncspclasses)

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanGTI[i,j] == 0 && meanDNI[i,j] == 0 && continue
            GTI = solarGTI[:, i, j]
            DNI = solarDNI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                (reg == 0 || reg == NOREGION) && continue

                area = cellarea[c]
                class = pvclass[i,j]
                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if class > 0
                    if mask_rooftop[r,c] > 0
                        capacity_pvrooftop[reg,class] += 1/1000 * pv_density * pvroof_area * area
                        increment_solarCF!(CF_pvrooftop[:,reg,class], GTI)
                        count_pvrooftop[reg,class] += 1
                    elseif mask_plantA[r,c] > 0
                        capacity_pvplantA[reg,class] += 1/1000 * pv_density * plant_area * area
                        increment_solarCF!(CF_pvplantA[:,reg,class], GTI)
                        count_pvplantA[reg,class] += 1
                    elseif mask_plantB[r,c] > 0
                        capacity_pvplantB[reg,class] += 1/1000 * pv_density * 2 * plant_area * area
                        increment_solarCF!(CF_pvplantB[:,reg,class], GTI)
                        count_pvplantB[reg,class] += 1
                    end
                end

                class_pv = class
                class = cspclass[i,j]
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if class > 0
                    if mask_plantA[r,c] > 0
                        capacity_cspplantA[reg,class] += 1/1000 * csp_density * plant_area * area
                        increment_solarCF!(CF_cspplantA[:,reg,class], DNI)
                        count_cspplantA[reg,class] += 1
                    elseif mask_plantB[r,c] > 0
                        capacity_cspplantB[reg,class] += 1/1000 * csp_density * 2 * plant_area * area
                        increment_solarCF!(CF_cspplantB[:,reg,class], DNI)
                        count_cspplantB[reg,class] += 1
                    end
                end

                if class_pv > 0 && class > 0
                    if mask_plantA[r,c] > 0
                        solar_overlap_areaA[reg,class_pv,class] += 1/1000 * plant_area * area
                    elseif mask_plantB[r,c] > 0
                        solar_overlap_areaB[reg,class_pv,class] += 1/1000 * 2 * plant_area * area
                    end
                end
            end
        end
        next!(updateprogress)
    end

    for y = 1:yearlength
        CF_pvrooftop[y,:,:] ./= count_pvrooftop
        CF_pvplantA[y,:,:] ./= count_pvplantA
        CF_pvplantB[y,:,:] ./= count_pvplantB
        CF_cspplantA[y,:,:] ./= count_cspplantA
        CF_cspplantB[y,:,:] ./= count_cspplantB
    end

    return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB, solar_overlap_areaA, solar_overlap_areaB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

function calc_wind_vars(Woption, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and wind class...")
    println("Interpolate ERA5 wind speeds later (maybe 4x runtime).")

    onshoreclass, offshoreclass = makewindclasses(Woption, windatlas)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Woption, lonrange, latrange)

    _, _, meanwind_allyears, _ = annualwindindex(Woption)
    @assert size(meanwind_allyears) == size(meanwind)

    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                onshore_density, area_onshore, offshore_density, area_offshore = Woption

    numreg = length(regionlist)
    nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)

    capacity_onshoreA = zeros(numreg,nonshoreclasses)
    capacity_onshoreB = zeros(numreg,nonshoreclasses)
    capacity_offshore = zeros(numreg,noffshoreclasses)
    windCF_onshoreA = zeros(yearlength,numreg,nonshoreclasses)
    windCF_onshoreB = zeros(yearlength,numreg,nonshoreclasses)
    windCF_offshore = zeros(yearlength,numreg,noffshoreclasses)
    count_onshoreA = zeros(Int,numreg,nonshoreclasses)
    count_onshoreB = zeros(Int,numreg,nonshoreclasses)
    count_offshore = zeros(Int,numreg,noffshoreclasses)

    if rescale_to_wind_atlas
        println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
        # println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    end

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanwind[i,j] == 0 && continue
            wind = rescale_to_wind_atlas ? windspeed[:, i, j] : speed2capacityfactor.(windspeed[:, i, j])
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                offreg = offshoreregions[r,c]
                area = cellarea[c]

                # @views is needed to make sure increment_windCF!() works with matrix
                # slice. It's also faster since it avoids making copies.
                # We divide by meanwind_allyears instead of meanwind in order to
                # account for more and less windy years, while still rescaling
                # hourly ERA5 wind speeds to match annual averages from the Global
                # Wind Atlas.
                if reg > 0 && reg != NOREGION
                    class = onshoreclass[r,c]
                    @views if reg > 0 && class > 0 && mask_onshoreA[r,c] > 0
                        capacity_onshoreA[reg,class] += 1/1000 * onshore_density * area_onshore * area
                        increment_windCF!(windCF_onshoreA[:,reg,class], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                        count_onshoreA[reg,class] += 1
                    elseif reg > 0 && class > 0 && mask_onshoreB[r,c] > 0
                        capacity_onshoreB[reg,class] += 1/1000 * onshore_density * 2 * area_onshore * area
                        increment_windCF!(windCF_onshoreB[:,reg,class], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                        count_onshoreB[reg,class] += 1
                    end
                end
                if offreg > 0 && offreg != NOREGION
                    offclass = offshoreclass[r,c]
                    @views if offreg > 0 && offclass > 0 && mask_offshore[r,c] > 0
                        capacity_offshore[offreg,offclass] += 1/1000 * offshore_density * area_offshore * area
                        increment_windCF!(windCF_offshore[:,offreg,offclass], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                        count_offshore[offreg,offclass] += 1
                    end
                end
            end
        end
        next!(updateprogress)
    end

    for y = 1:yearlength
        windCF_onshoreA[y,:,:] ./= count_onshoreA
        windCF_onshoreB[y,:,:] ./= count_onshoreB
        windCF_offshore[y,:,:] ./= count_offshore
    end

    return windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore
end



#= Quick and ugly copy/paste hack to create resource maps for solar classes combined with masks.
function GISsolarmap(; optionlist...)
    Soption = SolarOptions(merge(solaroptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix = Soption

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(Soption)
    meanGTI, solarGTI, meanDNI, solarDNI = read_solar_datasets(Soption, lonrange, latrange)

    mask_rooftop, mask_plantA, mask_plantB =
        create_solar_masks(Soption, regions, gridaccess, popdens, land, protected, lonrange, latrange, plotmasks=true)

    pvmap, pvrooftopmap, cspmap =
        calc_solar_map(Soption, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    return pvmap, pvrooftopmap, cspmap
end=#

#= Quick and ugly copy/paste hack to create resource maps for wind classes combined with masks.
function GISwindmap(; optionlist...)
    Woption = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix = Woption

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(Woption)
    windatlas, meanwind, windspeed = read_wind_datasets(Woption, lonrange, latrange)

    mask_onshoreA, mask_onshoreB, mask_offshore =
        create_wind_masks(Woption, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, plotmasks=true)

    onshoremap, offshoremap = calc_wind_map(Woption, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    return onshoremap, offshoremap
end=#

function calc_solar_map(Soption, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    pvclass, cspclass = makesolarclasses(Soption, meanGTI, meanDNI)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Soption, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and class...")
    println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, pvclasses_min, cspclasses_min, res, erares, pv_density, csp_density, pvroof_area, plant_area = Soption

    numreg = length(regionlist)
    npvclasses, ncspclasses = length(pvclasses_min), length(cspclasses_min)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    # pvmap = zeros(size(regions))
    # pvrooftopmap = zeros(size(regions))
    # cspmap = zeros(size(regions))
    pvmap = zeros(Int16, size(regions))
    pvrooftopmap = zeros(Int16, size(regions))
    cspmap = zeros(Int16, size(regions))

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    @inbounds for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanGTI[i,j] == 0 && meanDNI[i,j] == 0 && continue
            GTI = solarGTI[:, i, j]
            DNI = solarDNI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                area = cellarea[c]

                class = pvclass[i,j]
                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0
                    if mask_rooftop[r,c] > 0
                        # pvrooftopmap[r,c] = meanGTI[i,j]
                        pvrooftopmap[r,c] = class
                    elseif mask_plantA[r,c] > 0
                        # pvmap[r,c] = meanGTI[i,j]
                        pvmap[r,c] = class
                    elseif mask_plantB[r,c] > 0
                        # pvmap[r,c] = class
                    end
                end

                class = cspclass[i,j]
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0
                    if mask_plantA[r,c] > 0
                        # cspmap[r,c] = meanGTI[i,j]
                        cspmap[r,c] = class
                    elseif mask_plantB[r,c] > 0
                        # cspmap[r,c] = class
                    end
                end
            end
        end
        next!(updateprogress)
    end

    return pvmap, pvrooftopmap, cspmap
end

function calc_wind_map(Woption, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and wind class...")
    # println("Interpolate ERA5 wind speeds later (maybe 4x runtime).")

    onshoreclass, offshoreclass = makewindclasses(Woption, windatlas)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Woption, lonrange, latrange)

    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                onshore_density, area_onshore, offshore_density, area_offshore = Woption

    numreg = length(regionlist)
    nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)

    onshoremap = zeros(Int16, size(regions))
    offshoremap = zeros(Int16, size(regions))

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

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                area = cellarea[c]
                class = onshoreclass[r,c]
                offreg = offshoreregions[r,c]
                offclass = offshoreclass[r,c]

                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0 && mask_onshoreA[r,c] > 0
                    onshoremap[r,c] = class
                elseif reg > 0 && class > 0 && mask_onshoreB[r,c] > 0
                    # onshoremap[r,c] = class
                elseif offreg > 0 && offclass > 0 && mask_offshore[r,c] > 0
                    offshoremap[r,c] = class
                end
            end
        end
        next!(updateprogress)
    end

    return onshoremap, offshoremap
end

#Funzione nuova
function calc_pv_cf(Soption, meanGTI, solarGTI, pv_area, regions, offshoreregions, regionlist,
                mask_pv, lonrange, latrange)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Soption, lonrange, latrange)

    println("Calculating hourly capacity factors for each pixel (era5 resolution)")
    # println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, res, erares, pv_density, plant_area = Soption

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
    rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

    #added for later -> #capacity_pv = zeros(nlons,nlats)
    max_capacity_pv = zeros(nlons,nlats)
    CF_pvplant = zeros(yearlength,nlons,nlats)
    #max_capacity_pv = zeros(size(rowrange,1),size(colrange,1)) #SIZE FORSE SBAGLIATE, prova rowrange*nlons,colrange*nlats
    #CF_pvplant = zeros(yearlength,size(rowrange,1),size(colrange,1)) #SIZE FORSE SBAGLIATE prova rowrange*nlons,colrange*nlats

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]

        for i = 1:nlons
            ##meanGTI[i,j] == 0 && continue #verificare se non necessario per escludere alcune zone
            GTI = solarGTI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell

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
            ###AGGIUNGERE QUI CALCOLO WIND e POI LCOH, prob necessario add funzioni wind eventualmente rinominando


        end
        next!(updateprogress)
    end

    return CF_pvplant, max_capacity_pv
end

#Funzione Nuova
function calc_wind_cf(Woption, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshore, mask_offshore, lonrange, latrange)

    println("Calculating hourly capacity factors for each pixel (era5 resolution)")
    # println("Interpolate ERA5 wind speeds later (maybe 4x runtime).")

    onshoreclass, offshoreclass = makewindclasses(Woption, windatlas)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Woption, lonrange, latrange)

    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                onshore_density, area_onshore, offshore_density, area_offshore = Woption

    numreg = length(regionlist)
    nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)

    onshoremap = zeros(Int16, size(regions))
    offshoremap = zeros(Int16, size(regions))

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

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                area = cellarea[c]
                class = onshoreclass[r,c]
                offreg = offshoreregions[r,c]
                offclass = offshoreclass[r,c]

                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0 && mask_onshoreA[r,c] > 0
                    onshoremap[r,c] = class
                elseif reg > 0 && class > 0 && mask_onshoreB[r,c] > 0
                    # onshoremap[r,c] = class
                elseif offreg > 0 && offclass > 0 && mask_offshore[r,c] > 0
                    offshoremap[r,c] = class
                end
            end
        end
        next!(updateprogress)
    end

    return onshoremap, offshoremap
end



#Funzione nuova
function calc_area(Soption, regions, offshoreregions, regionlist, mask_pv, lonrange, latrange
                    Woption, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                    mask_onshore, mask_offshore,)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Soption, Woption, lonrange, latrange)
    onshoreclass, offshoreclass = makewindclasses(Woption, windatlas)

    println("Calculating usable area for each pixel (using high res datasets to be used for low res computations)")
    # println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, res, erares, plant_area = Soption
    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                    onshore_density, area_onshore, offshore_density, area_offshore = Woption

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)
    nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)

    onshoremap = zeros(Int16, size(regions))
    offshoremap = zeros(Int16, size(regions))

    if rescale_to_wind_atlas
        println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
        println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    end

    colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
    rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

    #added for later -> #capacity_pv = zeros(nlons,nlats)
    #CF_pvplant = zeros(yearlength,nlons,nlats)
    pv_area = zeros(nlons,nlats)
    wind_area = zeros(nlons,nlats)
    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]

        for i = 1:nlons
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                (reg == 0 || reg == NOREGION) && continue

                area = cellarea[c]

                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views
                if mask_pv[r,c] > 0
                    pv_area[i,j] += plant_area * area
                end
                @views
                if mask_onshore[r,c]>0
                    onshore_area[i,j]  += area_onshore  * area
                end
                @views
                if  mask_offshore[r,c]>0
                    offshore_area[i,j] += area_offshore * area
                end

            end
        end
        next!(updateprogress)
    end

    return pv_area, onshore_area, offshore_area
end
