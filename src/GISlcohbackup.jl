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
    Woptions = WindOptions()
    for (key,val) in Wd
        setproperty!(Woptions, key, val)
    end
    return Woptions
end

function GISlcoh(; savetodisk=true, plotmasks=false, optionlist...)

    Soptions = SolarOptions(merge(solaroptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix, pv_density, csp_density, downsample_masks = Soptions

    Woptions = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix, downsample_masks = Woptions

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(Soptions)
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(Woptions)

    mask_pv = create_pv_mask(Soptions, regions, popdens, land, protected, lonrange,
                            latrange; plotmasks=false, downsample=1)

    mask_onshore, mask_offshore = create_wind_masks(Woptions,regions, offshoreregions, gridaccess, popdens,
                                                    topo, land, protected, lonrange, latrange, plotmasks=plotmasks,
                                                    downsample=downsample_masks)


    plotmasks == :onlymasks && return nothing

    meanGTI, solarGTI, meanDNI, solarDNI = read_solar_datasets(Soptions, lonrange, latrange)
    windatlas, meanwind, windspeed = read_wind_datasets(Woptions, lonrange, latrange)

    pv_area, onshore_area, offshore_area = calc_area(Soptions, regions, offshoreregions, regionlist, mask_pv, lonrange, latrange,
                        Woptions, windatlas, meanwind, windspeed, mask_onshore, mask_offshore)

    max_capacity_pv, CF_pvplant = calc_pv_cf(Soptions, meanGTI, solarGTI, pv_area, regions, offshoreregions, regionlist,
                                            mask_pv, lonrange, latrange)

    max_capacity_onshore, max_capacity_offshore, CF_onshore, CF_offshore = calc_wind_cf(Woptions, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshore, mask_offshore, lonrange, latrange)

    if savetodisk
        mkpath(in_datafolder("output"))
        matopen(in_datafolder("output", "GISdata_lcoh$(era_year)_$gisregion$filenamesuffix.mat"), "w", compress=true) do file
            write(file, "CF_pvplant",  CF_pvplant)
            write(file, "max_apacity_pv", max_capacity_pv)
            write(file, "pv_density", pv_density)
            write(file, "max_capacity_onshore", max_capacity_onshore)
            write(file, "max_capacity_offshore", max_capacity_offshore)
            write(file, "CF_onshore", CF_onshore)
            write(file, "CF_offshore", CF_offshore)
        end
    end

    nothing
    # return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
    #         capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

function read_solar_datasets(Soptions, lonrange, latrange)
    @unpack res, erares, era_year = Soptions

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

function read_datasets(Woptions)
    @unpack res, gisregion, scenarioyear = Woptions

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

function read_wind_datasets(Woptions, lonrange, latrange)
    @unpack res, erares, era_year = Woptions
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
function create_pv_mask(Soptions, regions, popdens, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, plant_persons_per_km2,
            filenamesuffix = Soptions

    println("Creating PV masks...")

    goodland = (regions .> 0)
    for i in exclude_landtypes
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in protected_codes
        protected_area[protected .== i] .= true
    end

    km_per_degree = π*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)

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
        legendtext = ["bad land type", "high population", "protected area"]
        maskmap("$(gisregion)_masks_solar$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    end

    return mask_pv
end

#Nuova funzione
function create_wind_masks(Woptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, distance_elec_access, persons_per_km2,
                min_shore_distance, max_depth, classB_threshold, filenamesuffix = Woptions

    println("Creating wind masks...")

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
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
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
        legendtext = ["bad land type", "high population", "protected area", "onshore", "offshore"]
        maskmap("$(gisregion)_masks_wind$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)

        # drawmap(.!shore .& (offshoreregions .> 0))
        # drawmap((topo .> -max_depth) .& (offshoreregions .> 0))
        # drawmap(.!protected_area .& (offshoreregions .> 0))
        # drawmap(mask_offshore)
    end

    return mask_onshore, mask_offshore
end

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


#Funzione nuova
function calc_pv_cf(Soptions, meanGTI, solarGTI, pv_area, regions, offshoreregions, regionlist,
                mask_pv, lonrange, latrange)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Soptions, lonrange, latrange)

    println("Calculating hourly capacity factors for each pixel (era5 resolution)")

    @unpack era_year, res, erares, pv_density, plant_area = Soptions

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
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            ##meanGTI[i,j] == 0 && continue #verificare se non necessario per escludere alcune zone
            GTI = solarGTI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]
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
function calc_wind_cf(Woptions, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshore, mask_offshore, lonrange, latrange)

    println("Calculating hourly capacity factors for each pixel (era5 resolution)")

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Woptions, lonrange, latrange)

    _, _, meanwind_allyears, _ = annualwindindex(options)
    @assert size(meanwind_allyears) == size(meanwind)

    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                onshore_density, area_onshore, offshore_density, area_offshore = Woptions

    numreg = length(regionlist)
    #nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)


    max_capacity_onshore = zeros(nlons,nlats)
    max_capacity_offshore = zeros(nlons,nlats)
    CF_onshore = zeros(yearlength,nlons,nlats)
    CF_offshore = zeros(yearlength,nlons,nlats)

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
            #    reg = regions[r,c]
            #    area = cellarea[c]
            #    offreg = offshoreregions[r,c]

            #    if reg > 0 && reg != NOREGION
            #        @views if reg > 0  && mask_onshore[r,c] > 0
            #        increment_windCF!(windCF_onshoreA[:,reg,class], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
            #        capacity_onshore[reg,class] += 1/1000 * onshore_density * area_onshore * area
            #    end
                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                max_capacity_onshore[i,j] = 1/1000 * onshore_density * onshore_wind_area[i,j]
                max_capacity_offshore[i,j] = 1/1000 * offshore_density * offshore_wind_area[i,j]
                #CF_wind_onshore[:,i,j]  =
                #CF_wind_offshore[:,i,j] =

                #end
            #end
        end
        next!(updateprogress)
    end

    return max_capacity_onshore, max_capacity_offshore, CF_onshore, CF_offshore
end

#Funzione nuova
function calc_area(Soptions, regions, offshoreregions, regionlist, mask_pv, lonrange, latrange,
                    Woptions, windatlas, meanwind, windspeed, mask_onshore, mask_offshore)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Soptions, Woptions, lonrange, latrange)
    onshoreclass, offshoreclass = makewindclasses(Woptions, windatlas)

    println("Calculating usable area for each pixel (using high res datasets to be used for low res computations)")
    # println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, res, erares, plant_area = Soptions
    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                    onshore_density, area_onshore, offshore_density, area_offshore = Woptions

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
    onshore_wind_area = zeros(nlons,nlats)
    offshore_wind_area = zeros(nlons,nlats)
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
                @views if mask_pv[r,c] > 0
                    pv_area[i,j] += plant_area * area
                end
                @views if mask_onshore[r,c]>0
                    onshore_wind_area[i,j]  += area_onshore  * area
                end
                @views if  mask_offshore[r,c]>0
                    offshore_wind_area[i,j] += area_offshore * area
                end

            end
        end
        next!(updateprogress)
    end

    return pv_area, onshore_wind_area, offshore_area
end
