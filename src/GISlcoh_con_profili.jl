using DelimitedFiles, StatsBase, Plots, FLoops

lcohoptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat"

    :pv_density => 45000,                  # Solar PV land use 45 Wp/m2 = 45 MWp/km2 (includes PV efficiency & module spacing, add latitude dependency later)

    :plant_area => 1,                 # area available for PV or CSP plants after the masks have been applied

    :onshore_density => 5000,              # about 30% of existing farms have at least 5 W/m2, will become more common
    :offshore_density => 7430,             # varies a lot in existing parks (4-18 W/m2)
                                        # For reference: 10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2
    :area_onshore => 1,               # area available for onshore wind power after the masks have been applied
    :area_offshore => 1,              # area available for offshore wind power after the masks have been applied
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

    #questa ?? dentro giswind
    regions,
    offshoreregions,
    regionlist,
    gridaccess,
    popdens,
    topo,
    land,
    protected,
    lonrange,
    latrange,
    aqueduct,
    slope =read_datasets(loptions);

    mask_pv_lcoh       = create_solar_mask_lcoh(loptions, regions, popdens, land, protected, lonrange, latrange; plotmasks=plotmasks, downsample=1)

    meanGTI,
    solarGTI           = read_solar_datasets_lcoh(loptions, lonrange, latrange)

    windatlas,
    meanwind,
    windspeed          = read_wind_datasets_lcoh(loptions, lonrange, latrange)

    mask_onshore_lcoh,
    mask_offshore_lcoh = create_wind_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, plotmasks=plotmasks, downsample=downsample_masks)

    mask_aqueduct_lcoh = create_aqueduct_mask_lcoh(regions, aqueduct)

    mask_desalination_lcoh,
    shore              = create_desalination_mask(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)

    mask_altitude_lcoh = create_altitude_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    #pv_area, onshore_wind_area, offshore_wind_area = calc_area(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh, lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh)

    #max_capacity_pv, CF_pvplant = calc_pv_cf(loptions, meanGTI, solarGTI, pv_area, regions, offshoreregions, regionlist, mask_pv_lcoh, lonrange, latrange)

    mask_lcoh = create_lcoh_mask(mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, mask_aqueduct_lcoh)

    #plotmasks == :onlymasks && return nothing
    #max_capacity_onshore, max_capacity_offshore, CF_onshore, CF_offshore = calc_wind_cf(loptions, windatlas, meanwind, windspeed, onshore_wind_area, offshore_wind_area, regions, offshoreregions, regionlist, mask_onshore_lcoh, mask_offshore_lcoh, lonrange, latrange)
    #WACC_value = 0.08
    WACC_value = "Variable_WACC"
    pv_region_matrix, onshore_region_matrix, offshore_region_matrix, el_region_matrix, lcoh_region_matrix, production_region_matrix,
    PV_pv_onshore_erares, ONSHORE_pv_onshore_erares, OFFSHORE_pv_onshore_erares, EL_pv_onshore_erares, LCOH_pv_onshore_erares, PROD_pv_onshore_erares,
    PV_pv_erares, ONSHORE_pv_erares, OFFSHORE_pv_erares, EL_pv_erares, LCOH_pv_erares, PROD_pv_erares,
    PV_onshore_erares, ONSHORE_onshore_erares, OFFSHORE_onshore_erares, EL_onshore_erares, LCOH_onshore_erares, PROD_onshore_erares,
    PV_offshore_erares, ONSHORE_offshore_erares, OFFSHORE_offshore_erares, EL_offshore_erares, LCOH_offshore_erares, PROD_offshore_erares,
    SUM_CF_PV, SUM_CF_ONSHORE, SUM_CF_OFFSHORE, GTI_export, CF_onshore_export, CF_offshore_export,
    sumcfpv_region_matrix, sumcfonshore_region_matrix, sumcfoffshore_region_matrix,
    pv_cf_distribution_matrix, onshore_cf_distribution_matrix, offshore_cf_distribution_matrix,
    area_idonea_pv, area_idonea_onshore, area_idonea_offshore, area_idonea_pv_onshore, area_tot  = calc_LCOH_CAPS_PROD(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh,
                                                                            lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh,
                                                                            meanGTI, WACC_value)

    println("area_idonea_pv")
    println(area_idonea_pv)

    println("area_idonea_onshore")
    println(area_idonea_onshore)

    println("area_idonea_offshore")
    println(area_idonea_offshore)

    println("area_idonea_pv_onshore")
    println(area_idonea_pv_onshore)

    println("area_tot")
    println(area_tot)

    vectorized_lcoh_matrix                 = vec(lcoh_region_matrix);
    vectorized_prod_matrix                 = vec(production_region_matrix);
    vectorized_pv_matrix                   = vec(pv_region_matrix)
    vectorized_onshore_matrix              = vec(onshore_region_matrix)
    vectorized_offshore_matrix             = vec(offshore_region_matrix)
    vectorized_el_matrix                   = vec(el_region_matrix)
    vectorized_sumcfpv_region_matrix       = vec(sumcfpv_region_matrix)
    vectorized_sumcfonshore_region_matrix  = vec(sumcfonshore_region_matrix)
    vectorized_sumcfoffshore_region_matrix = vec(sumcfoffshore_region_matrix)

    deleted_index = typemax(Int32)

    for s in 1:length(vectorized_lcoh_matrix)
        if vectorized_lcoh_matrix[s]==0 && vectorized_prod_matrix[s]==0 && vectorized_pv_matrix[s]==0 && vectorized_onshore_matrix[s]==0 && vectorized_offshore_matrix[s]==0 && vectorized_el_matrix[s]==0
            vectorized_lcoh_matrix[s]                 = deleted_index
            vectorized_prod_matrix[s]                 = deleted_index
            vectorized_pv_matrix[s]                   = deleted_index
            vectorized_onshore_matrix[s]              = deleted_index
            vectorized_offshore_matrix[s]             = deleted_index
            vectorized_el_matrix[s]                   = deleted_index
            vectorized_sumcfpv_region_matrix[s]       = deleted_index
            vectorized_sumcfonshore_region_matrix[s]  = deleted_index
            vectorized_sumcfoffshore_region_matrix[s] = deleted_index
        end
    end

    vectorized_lcoh_matrix_wo_0                 = filter(!isequal(deleted_index), vectorized_lcoh_matrix);
    vectorized_prod_matrix_wo_0                 = filter(!isequal(deleted_index), vectorized_prod_matrix);
    vectorized_pv_matrix_wo_0                   = filter(!isequal(deleted_index), vectorized_pv_matrix);
    vectorized_onshore_matrix_wo_0              = filter(!isequal(deleted_index), vectorized_onshore_matrix);
    vectorized_offshore_matrix_wo_0             = filter(!isequal(deleted_index), vectorized_offshore_matrix);
    vectorized_el_matrix_wo_0                   = filter(!isequal(deleted_index), vectorized_el_matrix);
    vectorized_sumcfpv_region_matrix_wo_0       = filter(!isequal(deleted_index), vectorized_sumcfpv_region_matrix);
    vectorized_sumcfonshore_region_matrix_wo_0  = filter(!isequal(deleted_index), vectorized_sumcfonshore_region_matrix);
    vectorized_sumcfoffshore_region_matrix_wo_0 = filter(!isequal(deleted_index), vectorized_sumcfoffshore_region_matrix);

    FLOH_vector              = zeros(length(vectorized_lcoh_matrix_wo_0));
    ENERGY_PROD_vector       = zeros(length(vectorized_lcoh_matrix_wo_0));
    ELECTROLYZER_CONS_vector = zeros(length(vectorized_lcoh_matrix_wo_0));
    CURTAILMENT_vector       = zeros(length(vectorized_lcoh_matrix_wo_0));
    eff_el                = 0.7;
    HHV_H2                = 39.4;
    H2_energy_consumption = HHV_H2/eff_el;

    for i in 1:length(FLOH_vector)
        FLOH_vector[i] = (vectorized_prod_matrix_wo_0[i]/(8760*vectorized_el_matrix_wo_0[i]/H2_energy_consumption))*8760;

        ENERGY_PROD_vector[i] = vectorized_sumcfpv_region_matrix_wo_0[i]*vectorized_pv_matrix_wo_0[i]+
                                vectorized_sumcfonshore_region_matrix_wo_0[i]*vectorized_onshore_matrix_wo_0[i]+
                                vectorized_sumcfoffshore_region_matrix_wo_0[i]*vectorized_offshore_matrix_wo_0[i];

        ELECTROLYZER_CONS_vector[i] = vectorized_prod_matrix_wo_0[i]*H2_energy_consumption;

        CURTAILMENT_vector[i] = ((ENERGY_PROD_vector[i]-ELECTROLYZER_CONS_vector[i])/ENERGY_PROD_vector[i])*8760;
    end

    results_matrix = [vectorized_lcoh_matrix_wo_0 vectorized_prod_matrix_wo_0 vectorized_pv_matrix_wo_0 vectorized_onshore_matrix_wo_0 vectorized_offshore_matrix_wo_0 vectorized_el_matrix_wo_0 FLOH_vector vectorized_sumcfpv_region_matrix_wo_0 vectorized_sumcfonshore_region_matrix_wo_0 vectorized_sumcfoffshore_region_matrix_wo_0 CURTAILMENT_vector]
    # for s in 1:length(size(vectorized_lcoh_matrix)[1])
    #     if vectorized_lcoh_matrix[s]==0 && vectorized_prod_matrix[s]==0 && vectorized_pv_matrix[s]==0 && vectorized_onshore_matrix[s]==0 && vectorized_offshore_matrix[s]==0 && vectorized_el_matrix[s]==0
    #         vectorized_lcoh_matrix_wo_0[s]     = 0
    #         vectorized_prod_matrix_wo_0[s]     = 0
    #         vectorized_pv_matrix_wo_0[s]       = 0
    #         vectorized_onshore_matrix_wo_0[s]  = 0
    #         vectorized_offshore_matrix_wo_0[s] = 0
    #         vectorized_el_matrix_wo_0[s]       = 0
    #     end
    # end

    #Levo tutti gli zeri, sono sicuro che adesso ho due vettori della stessa lunghezza, con i punti accoppiati bene
    # vectorized_lcoh_matrix_wo_0     = filter(!iszero, vectorized_lcoh_matrix)
    # vectorized_prod_matrix_wo_0     = filter(!iszero, vectorized_prod_matrix)
    # vectorized_pv_matrix_wo_0       = filter(!iszero, vectorized_pv_matrix)
    # vectorized_onshore_matrix_wo_0  = filter(!iszero, vectorized_onshore_matrix)
    # vectorized_offshore_matrix_wo_0 = filter(!iszero, vectorized_offshore_matrix)
    # vectorized_el_matrix_wo_0       = filter(!iszero, vectorized_el_matrix)

    # lcoh_prod_array = [vectorized_lcoh_matrix_wo_0 vectorized_prod_matrix_wo_0]
    # lcoh_prod_array_sorted = sort!(lcoh_prod_array, dims=1)
    # for i in 2:size(lcoh_prod_array_sorted)[1]
    #     lcoh_prod_array_sorted[i,2] = lcoh_prod_array_sorted[i-1,2] + lcoh_prod_array_sorted[i,2]
    # end

    # righe_matrice_prod = size(production_region_matrix)[1];
    # righe_matrice_prod = size(production_region_matrix)[2];
    # global sum_of_producions = 0;
    # for i in 1:righe_matrice_prod
    #     for j in 1:righe_matrice_prod
    #         global somma += production_region_matrix[i,j]
    #     end
    # end
    println("...saving data...");
    mkpath(in_datafolder("output"));
    header = ["LCOH (USD/kg)" "Prod (kg/year)" "PV (kW)" "Onshore (kW)" "Offshore (kW)" "Electrolyzer (kW)" "FLOH (h)" "Sum_CF_PV (h)" "Sum_CF_Onshore (h)" "Sum_CF_Offshore (h)" "Curtailment (h)"]
    filename0 = in_datafolder("output", "$gisregion results_matrix.txt");
    # filename1 = in_datafolder("output", "$gisregion vectorized_lcoh_matrix_wo_0.txt");
    # filename2 = in_datafolder("output", "$gisregion vectorized_prod_matrix_wo_0.txt");
    # filename3 = in_datafolder("output", "$gisregion vectorized_PV_matrix_wo_0.txt");
    # filename4 = in_datafolder("output", "$gisregion vectorized_onshore_matrix_wo_0.txt");
    # filename5 = in_datafolder("output", "$gisregion vectorized_offshore_matrix_wo_0.txt");
    # filename6 = in_datafolder("output", "$gisregion vectorized_el_matrix_wo_0.txt");
    filename7 = in_datafolder("output", "$gisregion GTI_CF_onshore_CF_offshore_export.txt")

    isfile(filename0) && rm(filename0);
    # isfile(filename1) && rm(filename1);
    # isfile(filename2) && rm(filename2);
    # isfile(filename3) && rm(filename3);
    # isfile(filename4) && rm(filename4);
    # isfile(filename5) && rm(filename5);
    # isfile(filename6) && rm(filename6);
    isfile(filename7) && rm(filename7);

    writedlm(filename0, [header ; results_matrix], '\t');
    # writedlm(filename1, vectorized_lcoh_matrix_wo_0);
    # writedlm(filename2, vectorized_prod_matrix_wo_0);
    # writedlm(filename3, vectorized_pv_matrix_wo_0);
    # writedlm(filename4, vectorized_onshore_matrix_wo_0);
    # writedlm(filename5, vectorized_offshore_matrix_wo_0);
    # writedlm(filename6, vectorized_el_matrix_wo_0);
    writedlm(filename7, [GTI_export CF_onshore_export CF_offshore_export], '\t');

    plot_lcoh(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, lcoh_region_matrix, regionlist, aqueduct, slope, mask_desalination_lcoh,shore, mask_altitude_lcoh; plotmasks=plotmasks, downsample=1)
    plot_pv(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, pv_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)
    plot_onshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, onshore_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)
    plot_offshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, offshore_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)
    plot_el(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, el_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)
    plot_prod(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, production_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)

    plot_PV_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_ONSHORE_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    #plot_OFFSHORE_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_pv_onshore_erares, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)
    plot_EL_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_LCOH_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_PROD_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    plot_PV_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    #plot_ONSHORE_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    #plot_OFFSHORE_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_EL_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_LCOH_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_PROD_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    #plot_PV_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_ONSHORE_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    #plot_OFFSHORE_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_EL_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_LCOH_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_PROD_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    #plot_PV_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    #plot_ONSHORE_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_OFFSHORE_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_EL_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_LCOH_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_PROD_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    #
    plot_GTI(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, SUM_CF_PV, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_wind_onshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, SUM_CF_ONSHORE, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
    plot_wind_offshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, SUM_CF_OFFSHORE, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    nothing
    #return pv_region_matrix, onshore_region_matrix, offshore_region_matrix, el_region_matrix,lcoh_region_matrix,production_region_matrix, aqueduct, slope, protected
end

function read_solar_datasets_lcoh(loptions, lonrange, latrange)
    @unpack res, erares, era_year = loptions

    println("Reading ERA5 solar datasets...\n")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
    #@time
     meanGTI, solarGTI = h5open(in_datafolder("era5solar$era_year.h5"), "r") do file
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

    println("Reading ERA5 wind speeds ...\n")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
    #@time
     meanwind, windspeed = h5open(in_datafolder("era5wind$era_year.h5"), "r") do file

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

    println("Creating PV masks...\n")

    goodland = (regions .> 0)
    for i in pv_exclude_landtypes
        goodland[land .== i] .= false
    end
    pv_protected_area = zeros(Bool, size(protected))
    for i in pv_protected_codes
        pv_protected_area[protected .== i] .= true
    end

    km_per_degree = ??*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)

    mask_pv_lcoh = (popdens .< pv_plant_persons_per_km2) .& goodland .& .!pv_protected_area

    return mask_pv_lcoh
end

#Nuova funzione
function create_wind_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, wind_exclude_landtypes, wind_protected_codes, wind_persons_per_km2,
                min_shore_distance, max_depth, filenamesuffix = loptions

    println("Creating wind masks...\n")

    goodland = (regions .> 0)
    for i in wind_exclude_landtypes
        goodland[land .== i] .= false
    end
    wind_protected_area = zeros(Bool, size(protected))
    for i in wind_protected_codes
        wind_protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = ??*2*6371/360
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

function create_aqueduct_mask_lcoh(regions, aqueduct)

    println("Creating aquedcut masks...\n")

    goodland = (regions .> 0)

    mask_aqueduct_lcoh = goodland .& (aqueduct.> 0.0) .& (aqueduct .<= 4.0)

    return mask_aqueduct_lcoh
end

function create_desalination_mask(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, wind_exclude_landtypes, wind_protected_codes, wind_persons_per_km2,
                min_shore_distance, max_depth, filenamesuffix = loptions

    println("Creating desalination masks...\n")
    goodland = (regions .> 0)
    for i in wind_exclude_landtypes
        goodland[land .== i] .= false
    end
    wind_protected_area = zeros(Bool, size(protected))
    for i in wind_protected_codes
        wind_protected_area[protected .== i] .= true
    end
    min_desalination_distance = 25
    km_per_degree = ??*2*6371/360
    disk = diskfilterkernel(min_desalination_distance/km_per_degree/res)
    shore = (imfilter(regions .== 0, disk) .> 1e-6)

    mask_desalination_lcoh = (shore) .& (regions .> 0) #.& .!wind_protected_area

    return mask_desalination_lcoh, shore
end

function create_altitude_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, wind_exclude_landtypes, wind_protected_codes, wind_persons_per_km2,
                min_shore_distance, max_depth, filenamesuffix = loptions

    println("Creating altitude masks...\n")

    max_height = 1500
    mask_altitude_lcoh = (topo .< max_height) .& (regions .> 0)

    return mask_altitude_lcoh
end

#Nuova funzione
function plot_lcoh(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, lcoh_region_matrix, regionlist, aqueduct, slope, mask_desalination_lcoh, shore, mask_altitude_lcoh; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare
        #to_plot = mask_altitude_lcoh
        to_plot = lcoh_region_matrix;       #Quello che voglio plottare

        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'

        # masks2[mask_lcoh .== 0] .= 3
        # # masks2[mask_pv_lcoh .== 0] .= 3        # Maschera dove non si puo' fare LCOH
        # # masks2[mask_onshore_lcoh .== 0] .= 4   # Maschera dove non si puo' fare LCOH
        # #masks[regions .== 0] .= 0               # Ploting LCOH zona mare
        # #masks[regions .== NOREGION] .= 1        # Ploting LCOH zona NOREGION
        # masks2[regions .== 0 ] .= 0              # Maschera zona mare
        # masks2[regions .== NOREGION] .= 1       # Maschera zona NOREGION

        # legendtext = ["", "", "","","","","",""]

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

        # println("Mapping colors to regions (avoid same color in adjacent regions)...")
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

        # println("\nOnshore map...")

        # surface(xs, ys; color=float.(masks), colormap= :blues, colorrange= (0,),
        #                         shading=false, show_axis=false, scale_plot=false, interpolate=false)
         # scene = createmap("$(gisregion)_masks_lcoh$filenamesuffix", masks, masks2, [], lons, lats, [], source, dest, xs, ys,
         #     [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)

        scene = createmap("$(gisregion)_masks_lcoh$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
        #maskmap("$(gisregion)_masks_lcoh$filenamesuffix", masks, legendtext, lonrange, latrange; legend=false, downsample=downsample)
    end
    nothing
end
function plot_pv(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, pv_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = pv_region_matrix;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_pv$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_onshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, onshore_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = onshore_region_matrix;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_onshore$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_offshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, offshore_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = offshore_region_matrix;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_offshore$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_el(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, el_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = el_region_matrix;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_el$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_prod(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, production_region_matrix, regionlist, aqueduct, slope; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = production_region_matrix;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_prod$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
#avevo fatto una funzione diversa per provare a plottare a risoluzione i j
function plot_PV_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    PV_pv_onshore_res = rescale_property(PV_pv_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1

        to_plot = PV_pv_onshore_res;       #Quello che voglio plottare

        regions  = regions[1:downsample:end, 1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)

        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)

        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_PV_pv_onshore_erares$filenamesuffix", to_plot, regions, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_ONSHORE_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    ONSHORE_pv_onshore_res = rescale_property(ONSHORE_pv_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)


    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = ONSHORE_pv_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_ONSHORE_pv_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_OFFSHORE_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    OFFSHORE_pv_onshore_res = rescale_property(OFFSHORE_pv_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = OFFSHORE_pv_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_OFFSHORE_pv_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_EL_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    EL_pv_onshore_res = rescale_property(EL_pv_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = EL_pv_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_EL_pv_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_LCOH_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    LCOH_pv_onshore_res = rescale_property(LCOH_pv_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = LCOH_pv_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_LCOH_pv_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_PROD_pv_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_pv_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    PROD_pv_onshore_res = rescale_property(PROD_pv_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = PROD_pv_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_PROD_pv_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end

function plot_PV_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    PV_pv_res = rescale_property(PV_pv_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = PV_pv_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_PV_pv_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
# function plot_ONSHORE_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
#
#     @unpack gisregion, filenamesuffix = loptions
#
#     if plotmasks != false   # can == :onlymasks as well
#         # drawmap(land)
#         isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa
#
#         # mask values refer to colors in ColorBrewer Set2_7:
#         # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
#         #
#         # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
#         # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
#         #to_plot = slope;       #Quello che voglio plottare
#
#         to_plot = ONSHORE_pv_erares;       #Quello che voglio plottare
#         to_mask = ones(size(regions))*4;     #Quello che voglio mascherare
#
#         to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
#                                                # Assegno a questa zona il colore '0'
#         #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#         to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#                                                # Assegno a questa zona il colore '1'
#         #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
#                                                # Assegno a questa zona il colore '2'
#
#
#         regions = regions[1:downsample:end, 1:downsample:end]
#         lonrange = lonrange[1:downsample:end]
#         latrange = latrange[1:downsample:end]
#
#         println("\nProjecting coordinates (Mollweide)...")
#         res = 0.01
#         res2 = res/2
#         lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
#         lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
#         source = LonLat()
#         dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
#         xs, ys = xygrid(lons, lats)
#         Proj4.transform!(source, dest, vec(xs), vec(ys))
#
#         scene = createmap("$(gisregion)_masks_ONSHORE_pv_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
#              [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
#     end
#     nothing
# end
# function plot_OFFSHORE_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
#
#     @unpack gisregion, filenamesuffix = loptions
#
#     if plotmasks != false   # can == :onlymasks as well
#         # drawmap(land)
#         isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa
#
#         # mask values refer to colors in ColorBrewer Set2_7:
#         # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
#         #
#         # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
#         # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
#         #to_plot = slope;       #Quello che voglio plottare
#
#         to_plot = OFFSHORE_pv_erares;       #Quello che voglio plottare
#         to_mask = ones(size(regions))*4;     #Quello che voglio mascherare
#
#         to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
#                                                # Assegno a questa zona il colore '0'
#         #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#         to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#                                                # Assegno a questa zona il colore '1'
#         #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
#                                                # Assegno a questa zona il colore '2'
#
#
#         regions = regions[1:downsample:end, 1:downsample:end]
#         lonrange = lonrange[1:downsample:end]
#         latrange = latrange[1:downsample:end]
#
#         println("\nProjecting coordinates (Mollweide)...")
#         res = 0.01
#         res2 = res/2
#         lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
#         lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
#         source = LonLat()
#         dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
#         xs, ys = xygrid(lons, lats)
#         Proj4.transform!(source, dest, vec(xs), vec(ys))
#
#         scene = createmap("$(gisregion)_masks_OFFSHORE_pv_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
#              [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
#     end
#     nothing
# end
function plot_EL_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    EL_pv_res = rescale_property(EL_pv_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = EL_pv_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_EL_pv_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_LCOH_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    LCOH_pv_res = rescale_property(LCOH_pv_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = LCOH_pv_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_LCOH_pv_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_PROD_pv_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_pv_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    PROD_pv_res = rescale_property(PROD_pv_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = PROD_pv_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_PROD_pv_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end

# function plot_PV_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
#
#     @unpack gisregion, filenamesuffix = loptions
#
#     if plotmasks != false   # can == :onlymasks as well
#         # drawmap(land)
#         isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa
#
#         # mask values refer to colors in ColorBrewer Set2_7:
#         # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
#         #
#         # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
#         # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
#         #to_plot = slope;       #Quello che voglio plottare
#
#         to_plot = PV_onshore_erares;       #Quello che voglio plottare
#         to_mask = ones(size(regions))*4;     #Quello che voglio mascherare
#
#         to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
#                                                # Assegno a questa zona il colore '0'
#         #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#         to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#                                                # Assegno a questa zona il colore '1'
#         #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
#                                                # Assegno a questa zona il colore '2'
#
#
#         regions = regions[1:downsample:end, 1:downsample:end]
#         lonrange = lonrange[1:downsample:end]
#         latrange = latrange[1:downsample:end]
#
#         println("\nProjecting coordinates (Mollweide)...")
#         res = 0.01
#         res2 = res/2
#         lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
#         lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
#         source = LonLat()
#         dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
#         xs, ys = xygrid(lons, lats)
#         Proj4.transform!(source, dest, vec(xs), vec(ys))
#
#         scene = createmap("$(gisregion)_masks_PV_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
#              [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
#     end
#     nothing
# end
function plot_ONSHORE_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    ONSHORE_onshore_res = rescale_property(ONSHORE_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = ONSHORE_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_ONSHORE_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
# function plot_OFFSHORE_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
#
#     @unpack gisregion, filenamesuffix = loptions
#
#     if plotmasks != false   # can == :onlymasks as well
#         # drawmap(land)
#         isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa
#
#         # mask values refer to colors in ColorBrewer Set2_7:
#         # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
#         #
#         # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
#         # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
#         #to_plot = slope;       #Quello che voglio plottare
#
#         to_plot = OFFSHORE_onshore_erares;       #Quello che voglio plottare
#         to_mask = ones(size(regions))*4;     #Quello che voglio mascherare
#
#         to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
#                                                # Assegno a questa zona il colore '0'
#         #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#         to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#                                                # Assegno a questa zona il colore '1'
#         #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
#                                                # Assegno a questa zona il colore '2'
#
#
#         regions = regions[1:downsample:end, 1:downsample:end]
#         lonrange = lonrange[1:downsample:end]
#         latrange = latrange[1:downsample:end]
#
#         println("\nProjecting coordinates (Mollweide)...")
#         res = 0.01
#         res2 = res/2
#         lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
#         lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
#         source = LonLat()
#         dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
#         xs, ys = xygrid(lons, lats)
#         Proj4.transform!(source, dest, vec(xs), vec(ys))
#
#         scene = createmap("$(gisregion)_masks_OFFSHORE_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
#              [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
#     end
#     nothing
# end
function plot_EL_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    EL_onshore_res = rescale_property(EL_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = EL_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_EL_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_LCOH_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    LCOH_onshore_res = rescale_property(LCOH_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = LCOH_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_LCOH_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_PROD_onshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_onshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    PROD_onshore_res = rescale_property(PROD_onshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = PROD_onshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_PROD_onshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end

# function plot_PV_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PV_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
#
#     @unpack gisregion, filenamesuffix = loptions
#
#     if plotmasks != false   # can == :onlymasks as well
#         # drawmap(land)
#         isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa
#
#         # mask values refer to colors in ColorBrewer Set2_7:
#         # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
#         #
#         # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
#         # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
#         #to_plot = slope;       #Quello che voglio plottare
#
#         to_plot = PV_offshore_erares;       #Quello che voglio plottare
#         to_mask = ones(size(regions))*4;     #Quello che voglio mascherare
#
#         to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
#                                                # Assegno a questa zona il colore '0'
#         #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#         to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#                                                # Assegno a questa zona il colore '1'
#         #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
#                                                # Assegno a questa zona il colore '2'
#
#
#         regions = regions[1:downsample:end, 1:downsample:end]
#         lonrange = lonrange[1:downsample:end]
#         latrange = latrange[1:downsample:end]
#
#         println("\nProjecting coordinates (Mollweide)...")
#         res = 0.01
#         res2 = res/2
#         lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
#         lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
#         source = LonLat()
#         dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
#         xs, ys = xygrid(lons, lats)
#         Proj4.transform!(source, dest, vec(xs), vec(ys))
#
#         scene = createmap("$(gisregion)_masks_PV_offshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
#              [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
#     end
#     nothing
# end
# function plot_ONSHORE_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, ONSHORE_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)
#
#     @unpack gisregion, filenamesuffix = loptions
#
#     if plotmasks != false   # can == :onlymasks as well
#         # drawmap(land)
#         isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa
#
#         # mask values refer to colors in ColorBrewer Set2_7:
#         # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
#         #
#         # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
#         # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
#         #to_plot = slope;       #Quello che voglio plottare
#
#         to_plot = ONSHORE_offshore_erares;       #Quello che voglio plottare
#         to_mask = ones(size(regions))*4;     #Quello che voglio mascherare
#
#         to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
#                                                # Assegno a questa zona il colore '0'
#         #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#         to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
#                                                # Assegno a questa zona il colore '1'
#         #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
#                                                # Assegno a questa zona il colore '2'
#
#
#         regions = regions[1:downsample:end, 1:downsample:end]
#         lonrange = lonrange[1:downsample:end]
#         latrange = latrange[1:downsample:end]
#
#         println("\nProjecting coordinates (Mollweide)...")
#         res = 0.01
#         res2 = res/2
#         lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
#         lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
#         source = LonLat()
#         dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
#         xs, ys = xygrid(lons, lats)
#         Proj4.transform!(source, dest, vec(xs), vec(ys))
#
#         scene = createmap("$(gisregion)_masks_ONSHORE_offshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
#              [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
#     end
#     nothing
# end
function plot_OFFSHORE_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, OFFSHORE_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    OFFSHORE_offshore_res = rescale_property(OFFSHORE_offshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = OFFSHORE_offshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_OFFSHORE_offshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_EL_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, EL_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    EL_offshore_res = rescale_property(EL_offshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = EL_offshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_EL_offshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_LCOH_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, LCOH_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    LCOH_offshore_res = rescale_property(LCOH_offshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = LCOH_offshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_LCOH_offshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_PROD_offshore_erares(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, PROD_offshore_erares, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    PROD_offshore_res = rescale_property(PROD_offshore_erares, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = PROD_offshore_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_masks_PROD_offshore_erares$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end

function plot_GTI(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, SUM_CF_PV, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    SUM_CF_PV_res = rescale_property(SUM_CF_PV,  loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = SUM_CF_PV_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_SUM_CF_PV_res$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_wind_onshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, SUM_CF_ONSHORE, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    SUM_CF_ONSHORE_res = rescale_property(SUM_CF_ONSHORE,  loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = SUM_CF_ONSHORE_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'


        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_SUM_CF_ONSHORE_res$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end
function plot_wind_offshore(loptions, regions, mask_lcoh, mask_offshore_lcoh, lonrange, latrange, SUM_CF_OFFSHORE, regionlist, aqueduct, slope, solarGTI; plotmasks=plotmasks, downsample=1)

    @unpack gisregion, filenamesuffix = loptions

    SUM_CF_OFFSHORE_res = rescale_property(SUM_CF_OFFSHORE,  loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION) #questo mi trova la terraferma della regione che mi interessa

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        #
        # masks = lcoh_region_matrix  #questa maschera mi individua la quantit?? che voglio plottare
        # masks2 = ones(size(regions))*2 #mi inizializzo una matrice di grandezza regions di tutti 2
        #to_plot = slope;       #Quello che voglio plottare

        to_plot = SUM_CF_OFFSHORE_res;       #Quello che voglio plottare
        to_mask = ones(size(regions))*4;     #Quello che voglio mascherare

        to_mask[regions .== NOREGION] .= 0;    # Dentro to_mask individuo la zona NOREGION,
                                               # Assegno a questa zona il colore '0'
        #to_mask[(regions .== 0)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
        to_mask[(regions .== 0) .& (.!mask_offshore_lcoh)]      .= 1;    # Dentro to_mask individuo la zona mare (TUTTA, inclusa la off_shore)
                                               # Assegno a questa zona il colore '1'
        #to_mask[(mask_lcoh .== false)]    .= 2;    # Dentro to_mask individuo la zona dove NON posso fare idrogeno
                                               # Assegno a questa zona il colore '2'

        regions = regions[1:downsample:end, 1:downsample:end]
        lonrange = lonrange[1:downsample:end]
        latrange = latrange[1:downsample:end]

        println("\nProjecting coordinates (Mollweide)...")
        res = 0.01
        res2 = res/2
        lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
        lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
        source = LonLat()
        dest = Projection("+proj=moll +lon_0=$(mean(lons)) +R=3185500")
        xs, ys = xygrid(lons, lats)
        Proj4.transform!(source, dest, vec(xs), vec(ys))

        scene = createmap("$(gisregion)_SUM_CF_OFFSHORE_res$filenamesuffix", to_plot, to_mask, regionlist, lons, lats, [], source, dest, xs, ys,
             [], [], [], [], lines=false, labels=false, resolutionscale=1, textscale=1, legend=false)
    end
    nothing
end


function create_lcoh_mask(mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, mask_aqueduct_lcoh)

    #mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh sono di tipo BitArray.
    #Space-efficient N-dimensional boolean array, using just one bit for each boolean value.
    # Hanno dentro tutti 0 e 1 con valore booleano
    # secondo me ha senso fargli fare le operazioni logiche di 'or'

    mask_lcoh =  mask_pv_lcoh .| mask_onshore_lcoh .| mask_offshore_lcoh

    return mask_lcoh
end

#Funzione nuova
function calc_LCOH_CAPS_PROD(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh,
                            lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh,
                            meanGTI, WACC_value)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    global yearlength, nlons, nlats = size(solarGTI)

    # for i in 1:size(solarGTI)[2]
    #     for j in 1:(solarGTI)[3]
    #         for t in 1:size(solarGTI)[1]
    #             if solarGTI[i]<1
    #                 solarGTI[t,i,j] = 0
    #             end
    #         end
    #     end
    # end
    # Matrici bassa risoluzione
        # PV and ONSHORE
        PV_pv_onshore_erares       = zeros(nlons, nlats);
        ONSHORE_pv_onshore_erares  = zeros(nlons, nlats);
        OFFSHORE_pv_onshore_erares = zeros(nlons, nlats);
        EL_pv_onshore_erares       = zeros(nlons, nlats);
        LCOH_pv_onshore_erares     = zeros(nlons, nlats);
        PROD_pv_onshore_erares     = zeros(nlons, nlats);
        # PV
        PV_pv_erares               = zeros(nlons, nlats);
        ONSHORE_pv_erares          = zeros(nlons, nlats);
        OFFSHORE_pv_erares         = zeros(nlons, nlats);
        EL_pv_erares               = zeros(nlons, nlats);
        LCOH_pv_erares             = zeros(nlons, nlats);
        PROD_pv_erares             = zeros(nlons, nlats);
        # ONSHORE
        PV_onshore_erares          = zeros(nlons, nlats);
        ONSHORE_onshore_erares     = zeros(nlons, nlats);
        OFFSHORE_onshore_erares    = zeros(nlons, nlats);
        EL_onshore_erares          = zeros(nlons, nlats);
        LCOH_onshore_erares        = zeros(nlons, nlats);
        PROD_onshore_erares        = zeros(nlons, nlats);
        # OFFSHORE
        PV_offshore_erares         = zeros(nlons, nlats);
        ONSHORE_offshore_erares    = zeros(nlons, nlats);
        OFFSHORE_offshore_erares   = zeros(nlons, nlats);
        EL_offshore_erares         = zeros(nlons, nlats);
        LCOH_offshore_erares       = zeros(nlons, nlats);
        PROD_offshore_erares       = zeros(nlons, nlats);

    # Matrice alta risoluzione
        pv_region_matrix            = zeros(size(regions));
        onshore_region_matrix       = zeros(size(regions));
        offshore_region_matrix      = zeros(size(regions));
        el_region_matrix            = zeros(size(regions));
        lcoh_region_matrix          = zeros(size(regions));
        production_region_matrix    = zeros(size(regions));
        sumcfpv_region_matrix       = zeros(size(regions));
        sumcfonshore_region_matrix  = zeros(size(regions));
        sumcfoffshore_region_matrix = zeros(size(regions));

        pv_cf_distribution_matrix       = zeros(yearlength, size(regions)[1],size(regions)[2]);
        onshore_cf_distribution_matrix  = zeros(yearlength, size(regions)[1],size(regions)[2]);
        offshore_cf_distribution_matrix = zeros(yearlength, size(regions)[1],size(regions)[2]);

        SUM_CF_PV       = reshape(sum(solarGTI, dims=1), (nlons, nlats));
        SUM_CF_ONSHORE  = reshape(sum(speed2capacityfactor_onshore.(windspeed), dims=1), (nlons, nlats));
        SUM_CF_OFFSHORE = reshape(sum(speed2capacityfactor_offshore.(windspeed), dims=1), (nlons, nlats));

    #MA SERVONO?
    _, _, meanwind_allyears, _ = annualwindindex(loptions)
    @assert size(meanwind_allyears) == size(meanwind)

    @unpack era_year, res, erares, rescale_to_wind_atlas,
            pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore = loptions

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    #yearlength, nlons, nlats = size(windspeed)
    firsttime = DateTime(era_year, 1, 1)

    #MA SERVONO?
    onshoremap = zeros(Int16, size(regions))
    offshoremap = zeros(Int16, size(regions))

    # Matrici r,c resolution
    LCOH_res = zeros(size(regions))

    global area_idonea_pv = 0;
    global area_idonea_onshore = 0;
    global area_idonea_offshore = 0;
    global area_idonea_pv_onshore = 0;
    global area_tot = 0;
    # if rescale_to_wind_atlas
    #     println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
    #     println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    # end
    # pv_area = zeros(nlons,nlats)
    # onshore_wind_area = zeros(nlons,nlats)
    # offshore_wind_area = zeros(nlons,nlats)

    global dim_year_lon_lats = zeros(3,1)
    global dim_year_lon_lats = size(solarGTI) #yearlength, nlons, nlats

    global contatore = 0;

    num_tot_pixel = nlats*nlons
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    #nlats
    #for j in randperm(nlats)
    #for j in 2:4
    #@floop ThreadedEx()
    #Threads.@threads
    @time   Threads.@threads for j in 1:nlats
                eralat = eralats[j]
                colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
                # println("colrange:")
                # println(colrange)
                calculation_pixel_column = latmap[lat2col(eralat, res)] # Si introduce un errore prendendo il pixel in centro alla sottomatrice

                if calculation_pixel_column == 0 && j == 1
                    calculation_pixel_column = 1
                elseif calculation_pixel_column == 0 && j == nlats
                    calculation_pixel_column = latmap[findlast(!isequal(0), latmap)]
                end

                area = cellarea[calculation_pixel_column]

                for i in 1:nlons
                    contatore = contatore + 1;
                    println("Pixel: $contatore/$num_tot_pixel")
                    eralon = eralons[i]

                    for s in 1:size(solarGTI)[1]
                        if solarGTI[s,i,j]<0
                            solarGTI[s,i,j] = 0
                        end
                    end

                    GTI = solarGTI[:, i, j]
                    wind_onshore  = speed2capacityfactor_onshore.(windspeed[:, i, j])
                    wind_offshore = speed2capacityfactor_offshore.(windspeed[:, i, j])
                    rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

                    PV_pv_onshore_erares[i,j],
                    ONSHORE_pv_onshore_erares[i,j],
                    OFFSHORE_pv_onshore_erares[i,j],
                    EL_pv_onshore_erares[i,j],
                    LCOH_pv_onshore_erares[i,j],
                    PROD_pv_onshore_erares[i,j] = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                    GTI, area, wind_onshore, wind_offshore, 1, WACC_value); # Il numero sarebbe il type_min = 1
                    PV_pv_erares[i,j],
                    ONSHORE_pv_erares[i,j],
                    OFFSHORE_pv_erares[i,j],
                    EL_pv_erares[i,j],
                    LCOH_pv_erares[i,j],
                    PROD_pv_erares[i,j]         = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                    GTI, area, wind_onshore, wind_offshore, 2, WACC_value); # Il numero sarebbe il type_min  = 2
                    PV_onshore_erares[i,j],
                    ONSHORE_onshore_erares[i,j],
                    OFFSHORE_onshore_erares[i,j],
                    EL_onshore_erares[i,j],
                    LCOH_onshore_erares[i,j],
                    PROD_onshore_erares[i,j]    = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                    GTI, area, wind_onshore, wind_offshore, 3, WACC_value); # Il numero sarebbe il type_min  = 3
                    PV_offshore_erares[i,j],
                    ONSHORE_offshore_erares[i,j],
                    OFFSHORE_offshore_erares[i,j],
                    EL_offshore_erares[i,j],
                    LCOH_offshore_erares[i,j],
                    PROD_offshore_erares[i,j]   = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                    GTI, area, wind_onshore, wind_offshore, 4, WACC_value); # Il numero sarebbe il type_min  = 4
                    for c in colrange
                        for r in rowrange
                            #if contatore%1000 == 0
                            # println("Pixel $contatore/$num_tot_pixel")
                            # contatore = contatore + 1;
                            #end
                            (c == 0 || r == 0) && continue
                            reg = regions[r,c] #region all'indice r, c
                            #(reg == 0 || reg == NOREGION) && continue
                            (reg == NOREGION) && continue

                            pv_cf_distribution_matrix[:,r,c]       = GTI;
                            onshore_cf_distribution_matrix[:,r,c]  = wind_onshore;
                            offshore_cf_distribution_matrix[:,r,c] = wind_offshore;

                            area = cellarea[c]

                            if mask_pv_lcoh[r,c]>0
                                area_idonea_pv += area
                            end
                            if mask_onshore_lcoh[r,c]>0
                                area_idonea_onshore += area
                            end
                            if mask_offshore_lcoh[r,c]>0
                                area_idonea_offshore += area
                            end
                            if mask_pv_lcoh[r,c]>0 && mask_onshore_lcoh[r,c]>0
                                area_idonea_pv_onshore += area
                            end

                                area_tot += area

                            #davanti all if c'era un @views
                            if mask_pv_lcoh[r,c]>0 && mask_onshore_lcoh[r,c]>0  #sia pv che onshore
                                pv_region_matrix[r,c]            = PV_pv_onshore_erares[i,j];
                                onshore_region_matrix[r,c]       = ONSHORE_pv_onshore_erares[i,j];
                                offshore_region_matrix[r,c]      = OFFSHORE_pv_onshore_erares[i,j];
                                el_region_matrix[r,c]            = EL_pv_onshore_erares[i,j];
                                lcoh_region_matrix[r,c]          = LCOH_pv_onshore_erares[i,j]#min(LCOH_pv_onshore_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_pv_onshore_erares[i,j];
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];

                                area_idonea_pv_onshore += area
                            end
                            if mask_pv_lcoh[r,c]>0 && mask_onshore_lcoh[r,c]==0 #solo pv
                                pv_region_matrix[r,c]            = PV_pv_erares[i,j];
                                onshore_region_matrix[r,c]       = ONSHORE_pv_erares[i,j];
                                offshore_region_matrix[r,c]      = OFFSHORE_pv_erares[i,j];
                                el_region_matrix[r,c]            = EL_pv_erares[i,j];
                                lcoh_region_matrix[r,c]          = LCOH_pv_erares[i,j]#min(LCOH_pv_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_pv_erares[i,j];
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];

                                area_idonea_pv += area

                            end
                            if mask_pv_lcoh[r,c]==0 && mask_onshore_lcoh[r,c]>0 #solo onshore
                                pv_region_matrix[r,c]            = PV_onshore_erares[i,j];
                                onshore_region_matrix[r,c]       = ONSHORE_onshore_erares[i,j];
                                offshore_region_matrix[r,c]      = OFFSHORE_onshore_erares[i,j];
                                el_region_matrix[r,c]            = EL_onshore_erares[i,j];
                                lcoh_region_matrix[r,c]          = LCOH_onshore_erares[i,j]#min(LCOH_onshore_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_onshore_erares[i,j];
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];

                                area_idonea_onshore += area

                            end
                            if mask_offshore_lcoh[r,c]>0 #solo offshore
                                pv_region_matrix[r,c]            = PV_offshore_erares[i,j];
                                onshore_region_matrix[r,c]       = ONSHORE_offshore_erares[i,j];
                                offshore_region_matrix[r,c]      = OFFSHORE_offshore_erares[i,j];
                                el_region_matrix[r,c]            = EL_offshore_erares[i,j];
                                lcoh_region_matrix[r,c]          = LCOH_offshore_erares[i,j]#min(LCOH_offshore_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_offshore_erares[i,j]
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];

                                area_idonea_offshore += area

                            end
                        end
                    end
                end
             next!(updateprogress)
            end

            GTI_export = solarGTI[:, 4, 3]
            CF_onshore_export = speed2capacityfactor_onshore.(windspeed[:, 4, 3])
            CF_offshore_export  = speed2capacityfactor_onshore.(windspeed[:, 4, 3])
    return pv_region_matrix, onshore_region_matrix, offshore_region_matrix, el_region_matrix, lcoh_region_matrix, production_region_matrix,
            PV_pv_onshore_erares, ONSHORE_pv_onshore_erares, OFFSHORE_pv_onshore_erares, EL_pv_onshore_erares, LCOH_pv_onshore_erares, PROD_pv_onshore_erares,
            PV_pv_erares, ONSHORE_pv_erares, OFFSHORE_pv_erares, EL_pv_erares, LCOH_pv_erares, PROD_pv_erares,
            PV_onshore_erares, ONSHORE_onshore_erares, OFFSHORE_onshore_erares, EL_onshore_erares, LCOH_onshore_erares, PROD_onshore_erares,
            PV_offshore_erares, ONSHORE_offshore_erares, OFFSHORE_offshore_erares, EL_offshore_erares, LCOH_offshore_erares, PROD_offshore_erares,
            SUM_CF_PV, SUM_CF_ONSHORE, SUM_CF_OFFSHORE, GTI_export, CF_onshore_export, CF_offshore_export,sumcfpv_region_matrix, sumcfonshore_region_matrix, sumcfoffshore_region_matrix,
            pv_cf_distribution_matrix, onshore_cf_distribution_matrix, offshore_cf_distribution_matrix,
            area_idonea_pv, area_idonea_onshore, area_idonea_offshore, area_idonea_pv_onshore, area_tot
end

function find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                           CF_pvplant, area, CF_onshore, CF_offshore, type_min, WACC_value)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    @unpack era_year, res, erares, pv_density, onshore_density, offshore_density = loptions
    numreg = length(regionlist)
    # yearlength, nlons, nlats = size(solarGTI)
    yearlength, nlons, nlats = dim_year_lon_lats
    firsttime = DateTime(era_year, 1, 1)

    CAPEX_el               = 250;
    CAPEX_pv               = 185;
    CAPEX_wind_onshore     = 896;
    CAPEX_wind_offshore    = 1746;
    WACC_el                = 0.08
    # WACC_pv                = 0.0543
    # WACC_onshore           = 0.0647
    # WACC_offshore          = 0.0785

    if WACC_value == "Variable_WACC"
        WACC_pv       = 0.0364
        WACC_onshore  = 0.0364
        WACC_offshore = 0.0462
    else
        WACC_pv       = WACC_value
        WACC_onshore  = WACC_value
        WACC_offshore = WACC_value
    end

    x_min_LCOH_rescaled                 = 0
    y_onshore_min_LCOH_rescaled         = 0
    y_offshore_min_LCOH_rescaled        = 0
    z_min_LCOH_rescaled                 = 0
    minLCOH                             = 0
    rescaled_annual_hydrogen_production = 0
    max_capacity_pv                     = 0
    max_capacity_onshore                = 0
    max_capacity_offshore               = 0


    if type_min == 1
        max_capacity_pv      = area * pv_density
        max_capacity_onshore = area * onshore_density

        x_min_LCOH_rescaled,
        y_onshore_min_LCOH_rescaled,
        z_min_LCOH_rescaled,
        minLCOH,
        rescaled_annual_hydrogen_production = iterative_find_min_3var(loptions, CF_pvplant, max_capacity_pv,
                                                                        CF_onshore, max_capacity_onshore,
                                                                        CAPEX_pv, CAPEX_wind_onshore, CAPEX_el,
                                                                        WACC_el, WACC_pv, WACC_onshore)
    else
        max_capacity_pv       = area * pv_density
        max_capacity_onshore  = area * onshore_density
        max_capacity_offshore = area * offshore_density

        x_min_LCOH_rescaled,
        y_onshore_min_LCOH_rescaled,
        y_offshore_min_LCOH_rescaled,
        z_min_LCOH_rescaled,
        minLCOH,
        rescaled_annual_hydrogen_production = iterative_find_min_2var(loptions, CF_pvplant, max_capacity_pv,
                                                                        CF_onshore, max_capacity_onshore,
                                                                        CF_offshore, max_capacity_offshore,
                                                                        CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, CAPEX_el,
                                                                        WACC_el, WACC_pv, WACC_onshore, WACC_offshore,
                                                                        type_min)
    end

    return x_min_LCOH_rescaled, y_onshore_min_LCOH_rescaled, y_offshore_min_LCOH_rescaled, z_min_LCOH_rescaled, minLCOH, rescaled_annual_hydrogen_production
end

#Funzione nuova
function iterative_find_min_3var(loptions, CF_pvplant, max_capacity_pv,
                                CF_onshore, max_capacity_onshore,
                                CAPEX_pv, CAPEX_wind, CAPEX_el,
                                WACC_el, WACC_pv, WACC_onshore)

    # println("Minimum di CF_pvplant")
    # println(minimum(CF_pvplant))
    # println("Minimum di CF_onshore")
    # println(minimum(CF_onshore))
    # println("Maximum di CF_pvplant")
    # println(maximum(CF_pvplant))
    # println("Maximum di CF_onshore")
    # println(maximum(CF_onshore))

    lifetime_stack        = 80000;              # Hours
    lifetime_el           = 15;                 # Years
    lifetime_pv           = 25;                 # Years
    lifetime_wind         = 25;                 # Years
    OPEX_el               = 0.03;               # %CAPEX
    OPEX_wind             = 0.015;              # %CAPEX

    OPEX_pv               = 0.015;
    eff_el                = 0.7;
    min_CF_electrolyzer   = 0.3;
    HHV_H2                = 39.4;
    H2_energy_consumption = HHV_H2/eff_el;

    cost_pv            =   (CAPEX_pv*WACC_pv)/(1-(1+WACC_pv)^(-lifetime_pv))+(CAPEX_pv*OPEX_pv);
    cost_wind          =   (CAPEX_wind*WACC_onshore)/(1-(1+WACC_onshore)^(-lifetime_wind))+(CAPEX_wind*OPEX_wind);
    cost_el            =   (CAPEX_el*WACC_el)/(1-(1+WACC_el)^(-lifetime_el))+(CAPEX_el*OPEX_el);

    global LCOH_min     = 0;
    global xmin_LCOH    = 0;
    global ymin_LCOH    = 0;
    global zmin_LCOH    = 0;
    global max_cap_pv   = max_capacity_pv;
    global max_cap_wind = max_capacity_onshore;

    # println("\nMax cap pv:")
    # println(max_cap_pv)
    # println("\nMax cap wind:")
    # println(max_cap_wind)

    global startvalue_pv   = 0;
    global startvalue_wind = 0;
    global endvalue_pv     = max_cap_pv;
    global endvalue_wind   = max_cap_wind;
    global startvalue_el   = 1;
    global endvalue_el     = max_cap_pv + max_cap_wind;

    if max_cap_pv<=10 && max_cap_wind<=10
            return 0,0,0,0,0;
        elseif max_cap_pv<=0 && max_cap_wind>0
            endvalue_pv=startvalue_pv;
        elseif max_cap_pv>0 && max_cap_wind<=0
            endvalue_wind=startvalue_wind;
    end

    global dd = (7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4);
    for l in 1:length(dd)

        d=dd[l];
        # println("\nIter: $l");
        # println("\nstartvalue_pv");
        # println(startvalue_pv);
        # println("endvalue_pv");
        # println(endvalue_pv);
        # println("\nstartvalue_wind");
        # println(startvalue_wind);
        # println("endvalue_wind");
        # println(endvalue_wind);
        # println("\nstartvalue_el");
        # println(startvalue_el);
        # println("endvalue_el");
        # println(endvalue_el);

        x=collect(range(startvalue_pv,   length=d, stop=endvalue_pv));
        y=collect(range(startvalue_wind, length=d, stop=endvalue_wind));
        z=collect(range(startvalue_el,   length=d, stop=endvalue_el));

        annual_invest_LCOE=zeros(length(x),length(y));
        annual_energy_production=zeros(length(x),length(y));
        annual_invest_LCOH=zeros(length(x),length(y),length(z));
        LCOE=zeros(length(x),length(y));
        LCOH=zeros(length(x),length(y),length(z));
        annual_hydrogen_production=zeros(length(x),length(y),length(z));

        sumCF_pv=sum(CF_pvplant);
        sumCF_wind=sum(CF_onshore);

        for i in 1:length(x)
    	    for j in 1:length(y)
    		    annual_invest_LCOE[i,j]       = cost_pv*x[i]  + cost_wind*y[j];
                annual_energy_production[i,j] = sumCF_pv*x[i] + sumCF_wind*y[j];
    		    for k in 1:length(z)
    	            annual_invest_LCOH[i,j,k] = annual_invest_LCOE[i,j] + cost_el*z[k];
                    for t in 1:8760
                        annual_hydrogen_production[i,j,k]+=min((CF_pvplant[t]*x[i]+CF_onshore[t]*y[j]),z[k])/H2_energy_consumption;
                    end
                end
            end
        end

        LCOE = 1000*annual_invest_LCOE./annual_energy_production;
        LCOH = annual_invest_LCOH./annual_hydrogen_production;
        LCOE_min, indicesLCOEmin=findmin(LCOE);

        xmin_LCOE = x[indicesLCOEmin[1]];
        ymin_LCOE = y[indicesLCOEmin[2]];

        minLCOH, indicesLCOHmin=findmin(LCOH);
        global LCOH_min  = minLCOH;
        global xmin_LCOH = x[indicesLCOHmin[1]];
        # println("xmin_LCOH: $xmin_LCOH")
        # println("lowerbound: $startvalue_pv")
        # println("upperbound: $endvalue_pv")
        global ymin_LCOH=y[indicesLCOHmin[2]];
        global zmin_LCOH=z[indicesLCOHmin[3]];
        delta_pv=(endvalue_pv-startvalue_pv)/(d-1);
        delta_wind=(endvalue_wind-startvalue_wind)/(d-1);
        delta_el=(endvalue_el-startvalue_el)/(d-1);

        global startvalue_pv=max(xmin_LCOH-delta_pv,1);
        global endvalue_pv=min(xmin_LCOH+delta_pv,max_cap_pv);
        global startvalue_wind=max(ymin_LCOH-delta_wind,1);
        global endvalue_wind=min(ymin_LCOH+delta_wind,max_cap_wind);
        global startvalue_el=zmin_LCOH-delta_el;
        global endvalue_el=zmin_LCOH+delta_el;

        if startvalue_pv < 0
            startvalue_pv = 0
        end

        if endvalue_pv < 0
            endvalue_pv = 0
        end

        if startvalue_wind < 0
            startvalue_wind = 0
        end

        if endvalue_wind < 0
            endvalue_wind = 0
        end

        if startvalue_el < 0
            startvalue_el = 0
        end

        if endvalue_el < 0
            endvalue_el = 0
        end

        # println("\nMinimum LCOH");
        # println(LCOH_min);
        #
        # println("\nMinimum X");
        # println(xmin_LCOH);
        #
        # println("\nMinimum Y");
        # println(ymin_LCOH);
        #
        # println("\nMinimum Z");
        # println(zmin_LCOH);
    end

    # global annual_hydrogen_production = Float64(0.0);
    #
    # for t in 1:8760
    #     global annual_hydrogen_production +=min((CF_pvplant[t]*xmin_LCOH+CF_onshore[t]*ymin_LCOH),zmin_LCOH)/H2_energy_consumption;
    # end

    if xmin_LCOH/max_cap_pv > ymin_LCOH/max_cap_wind
    	rescaling_factor = xmin_LCOH/max_cap_pv
    else
    	rescaling_factor = ymin_LCOH/max_cap_wind
    end

    x_min_LCOH_rescaled	        = xmin_LCOH/rescaling_factor;
    y_onshore_min_LCOH_rescaled	= ymin_LCOH/rescaling_factor;
    z_min_LCOH_rescaled		    = zmin_LCOH/rescaling_factor;

    global rescaled_annual_hydrogen_production = Float64(0.0);

    for t in 1:8760
        global rescaled_annual_hydrogen_production += min((CF_pvplant[t]*x_min_LCOH_rescaled + CF_onshore[t]*y_onshore_min_LCOH_rescaled),z_min_LCOH_rescaled)/H2_energy_consumption;
    end


    FLOH_el_pv_onshore = rescaled_annual_hydrogen_production/(z_min_LCOH_rescaled*8760/H2_energy_consumption)

    return x_min_LCOH_rescaled, y_onshore_min_LCOH_rescaled, z_min_LCOH_rescaled, LCOH_min, rescaled_annual_hydrogen_production
end

#Funzione nuova
function iterative_find_min_2var(loptions, CF_pvplant, max_capacity_pv,
                                CF_onshore, max_capacity_onshore,
                                CF_offshore, max_capacity_offshore,
                                CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, CAPEX_el,
                                WACC_el, WACC_pv, WACC_onshore, WACC_offshore,
                                type_min)


    @unpack pv_density, onshore_density  = loptions

    lifetime_stack = 80000;              # Hours
    lifetime_el = 15;                    # Years
    lifetime_wind_onshore = 25;          # Years
    lifetime_wind_offshore = 25;         # Years
    lifetime_pv = 25;                    # Years
    OPEX_el = 0.03;                      # %CAPEX
    OPEX_wind_onshore = 0.015;           # %CAPEX
    OPEX_wind_offshore = 0.015;          # %CAPEX
    OPEX_pv = 0.015;
    eff_el = 0.7;
    min_CF_electrolyzer = 0.3;
    HHV_H2 = 39.4;
    H2_energy_consumption = HHV_H2/eff_el;

    cost_pv            =   (CAPEX_pv*WACC_pv)/(1-(1+WACC_pv)^(-lifetime_pv))+(CAPEX_pv*OPEX_pv);
    cost_wind_onshore  =   (CAPEX_wind_onshore*WACC_onshore)/(1-(1+WACC_onshore)^(-lifetime_wind_onshore))+(CAPEX_wind_onshore*OPEX_wind_onshore);
    cost_wind_offshore =   (CAPEX_wind_offshore*WACC_offshore)/(1-(1+WACC_offshore)^(-lifetime_wind_offshore))+(CAPEX_wind_offshore*OPEX_wind_offshore);
    cost_el            =   (CAPEX_el*WACC_el)/(1-(1+WACC_el)^(-lifetime_el))+(CAPEX_el*OPEX_el);

    if type_min == 2 #solo pv
        global LCOH_min      = 0;
        global xmin_LCOH     = 0;
        global max_cap_pv    = max_capacity_pv;
        global zmin_LCOH=0;

        global startvalue_pv = 0;
        global endvalue_pv   = max_cap_pv;
        global startvalue_el = 1;
        global endvalue_el   = max_cap_pv;

        if max_cap_pv<=0
            return 0,0,0,0,0,0
        elseif max_cap_pv<=10
            return 0,0,0,0,0,0
        end

        dd=(7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4);
        for l in 1:length(dd)
            d=dd[l];

            x=collect(range(startvalue_pv, length=d, stop=endvalue_pv));
            z=collect(range(startvalue_el, length=d, stop=endvalue_el));

            annual_invest_LCOE=zeros(length(x));
            annual_energy_production=zeros(length(x));
            annual_invest_LCOH=zeros(length(x),length(z));
            annual_hydrogen_production=zeros(length(x),length(z));

            sumCF_pv=sum(CF_pvplant);

            for i in 1:length(x)
    	        annual_invest_LCOE[i]=cost_pv*x[i];
    	        annual_energy_production[i]=sumCF_pv*x[i];
    		        for k in 1:length(z)
    	                annual_invest_LCOH[i,k]=annual_invest_LCOE[i]+cost_el*z[k];
                            for t in 1:8760
                                annual_hydrogen_production[i,k]+=min((CF_pvplant[t]*x[i]),z[k])/H2_energy_consumption;
                            end
                    end
            end

            LCOE=1000*annual_invest_LCOE./annual_energy_production;
            LCOH=annual_invest_LCOH./annual_hydrogen_production;
            LCOE_min, indicesLCOEmin=findmin(LCOE);
            xmin_LCOE=x[indicesLCOEmin[1]];

            global minLCOH, indicesLCOHmin=findmin(LCOH);
            global LCOH_min= minLCOH;
            global xmin_LCOH=x[indicesLCOHmin[1]];
            global zmin_LCOH=z[indicesLCOHmin[2]];

            delta_pv=(endvalue_pv-startvalue_pv)/(d-1);
            delta_el=(endvalue_el-startvalue_el)/(d-1);

            global startvalue_pv=max(xmin_LCOH-delta_pv,1);
            global endvalue_pv=min(xmin_LCOH+delta_pv,max_cap_pv);
            global startvalue_el=zmin_LCOH-delta_el;
            global endvalue_el=zmin_LCOH+delta_el;

            if startvalue_pv < 0
                startvalue_pv = 0
            end

            if endvalue_pv < 0
                endvalue_pv = 0
            end

            if startvalue_el < 0
                startvalue_el = 0
            end

            if endvalue_el < 0
                endvalue_el = 0
            end

            # println("\nMinimum LCOH");
            # println(LCOH_min);
            #
            # println("\nMinimum X");
            # println(xmin_LCOH);
            #
            # println("\nMinimum Y");
            # println(ymin_LCOH);
            #
            # println("\nMinimum Z");
            # println(zmin_LCOH);
            if l==length(dd)
                global annual_hydrogen_production_final = annual_hydrogen_production
            end
        end

        # global annual_hydrogen_production = Float64(0.0);
        #
        # for t in 1:8760
        #     global annual_hydrogen_production +=min((CF_pvplant[t]*xmin_LCOH+CF_onshore[t]*ymin_LCOH),zmin_LCOH)/H2_energy_consumption;
        # end

    	rescaling_factor = xmin_LCOH/max_cap_pv

        x_min_LCOH_rescaled	         = xmin_LCOH / rescaling_factor;
        y_onshore_min_LCOH_rescaled  = 0;
        y_offshore_min_LCOH_rescaled = 0;
        z_min_LCOH_rescaled		     = zmin_LCOH / rescaling_factor;

        rescaled_annual_hydrogen_production = annual_hydrogen_production_final[indicesLCOHmin[1],indicesLCOHmin[2]]/rescaling_factor

        # global rescaled_annual_hydrogen_production = Float64(0.0);
        #
        # for t in 1:8760
        #     global rescaled_annual_hydrogen_production += min((CF_pvplant[t]*x_min_LCOH_rescaled),z_min_LCOH_rescaled)/H2_energy_consumption;
        # end

    end

    if type_min == 3 #solo wind onshore
        global LCOH_min      = 0;
        global yonshoremin_LCOH        = 0;
        global max_cap_wind_onshore    = max_capacity_onshore;
        global zmin_LCOH=0;

        global startvalue_wind_onshore = 0;
        global endvalue_wind_onshore   = max_cap_wind_onshore;
        global startvalue_el           = 1;
        global endvalue_el             = max_cap_wind_onshore;

        if max_cap_wind_onshore<=0
            return 0,0,0,0,0,0
        elseif max_cap_wind_onshore<=10
            return 0,0,0,0,0,0
        end

        dd=(7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4);
        for l in 1:length(dd)
            d=dd[l];

            y_onshore = collect(range(startvalue_wind_onshore, length=d, stop=endvalue_wind_onshore));
            z         = collect(range(startvalue_el, length=d, stop=endvalue_el));

            annual_invest_LCOE=zeros(length(y_onshore));
            annual_energy_production=zeros(length(y_onshore));
            annual_invest_LCOH=zeros(length(y_onshore),length(z));
            annual_hydrogen_production=zeros(length(y_onshore),length(z));

            sumCF_wind_onshore=sum(CF_onshore);

            for i in 1:length(y_onshore)
                annual_invest_LCOE[i]=cost_wind_onshore*y_onshore[i];
                annual_energy_production[i]=sumCF_wind_onshore*y_onshore[i];
                    for k in 1:length(z)
                        annual_invest_LCOH[i,k]=annual_invest_LCOE[i]+cost_el*z[k];
                            for t in 1:8760
                                annual_hydrogen_production[i,k]+=min((CF_onshore[t]*y_onshore[i]),z[k])/H2_energy_consumption;
                            end
                    end
            end

            LCOE =1000*annual_invest_LCOE./annual_energy_production;
            LCOH =annual_invest_LCOH./annual_hydrogen_production;
            LCOE_min, indicesLCOEmin=findmin(LCOE);
            yonshoremin_LCOE=y_onshore[indicesLCOEmin[1]];

            global minLCOH, indicesLCOHmin=findmin(LCOH);
            global LCOH_min = minLCOH;
            global yonshoremin_LCOH = y_onshore[indicesLCOHmin[1]];
            global zmin_LCOH = z[indicesLCOHmin[2]];
            delta_wind_onshore=(endvalue_wind_onshore-startvalue_wind_onshore)/(d-1);
            delta_el=(endvalue_el-startvalue_el)/(d-1);

            global startvalue_wind_onshore=max(yonshoremin_LCOH-delta_wind_onshore,1);
            global endvalue_wind_onshore=min(yonshoremin_LCOH+delta_wind_onshore,max_cap_wind_onshore);
            global startvalue_el=zmin_LCOH-delta_el;
            global endvalue_el=zmin_LCOH+delta_el;

            if startvalue_wind_onshore < 0
                startvalue_wind_onshore = 0
            end

            if endvalue_wind_onshore < 0
                endvalue_wind_onshore = 0
            end

            if startvalue_el < 0
                startvalue_el = 0
            end

            if endvalue_el < 0
                endvalue_el = 0
            end

            # println("\nMinimum LCOH");
            # println(LCOH_min);
            #
            # println("\nMinimum X");
            # println(xmin_LCOH);
            #
            # println("\nMinimum Y");
            # println(ymin_LCOH);
            #
            # println("\nMinimum Z");
            # println(zmin_LCOH);
            if l==length(dd)
                global annual_hydrogen_production_final = annual_hydrogen_production
            end
        end

        # global annual_hydrogen_production = Float64(0.0);
        #
        # for t in 1:8760
        #     global annual_hydrogen_production +=min((CF_pvplant[t]*xmin_LCOH+CF_onshore[t]*ymin_LCOH),zmin_LCOH)/H2_energy_consumption;
        # end

        rescaling_factor = yonshoremin_LCOH/max_cap_wind_onshore

        x_min_LCOH_rescaled          = 0;
        y_onshore_min_LCOH_rescaled	 = yonshoremin_LCOH/rescaling_factor;
        y_offshore_min_LCOH_rescaled = 0;
        z_min_LCOH_rescaled		     = zmin_LCOH/rescaling_factor;

        rescaled_annual_hydrogen_production = annual_hydrogen_production_final[indicesLCOHmin[1],indicesLCOHmin[2]]/rescaling_factor
        #global rescaled_annual_hydrogen_production = Float64(0.0);

        # for t in 1:8760
        #     global rescaled_annual_hydrogen_production += min((CF_onshore[t]*y_onshore_min_LCOH_rescaled),z_min_LCOH_rescaled)/H2_energy_consumption;
        # end

    end

    if type_min == 4 #solo wind offshore
        global LCOH_min              = 0;
        global yoffshoremin_LCOH     = 0;
        global max_cap_wind_offshore = max_capacity_offshore;
        global zmin_LCOH             = 0;

        global startvalue_wind_offshore = 0;
        global endvalue_wind_offshore   = max_cap_wind_offshore;
        global startvalue_el            = 1;
        global endvalue_el              = max_cap_wind_offshore;

        if max_cap_wind_offshore<=0
            return 0,0,0,0,0,0
        elseif max_cap_wind_offshore<=10
            return 0,0,0,0,0,0
        end

        dd=(7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4);
        for l in 1:length(dd)
            d=dd[l];

            y_offshore = collect(range(startvalue_wind_offshore, length=d, stop=endvalue_wind_offshore));
            z          = collect(range(startvalue_el, length=d, stop=endvalue_el));

            annual_invest_LCOE       = zeros(length(y_offshore));
            annual_energy_production = zeros(length(y_offshore));
            annual_invest_LCOH       = zeros(length(y_offshore),length(z));
            LCOE = zeros(length(y_offshore));
            LCOH = zeros(length(y_offshore),length(z));
            annual_hydrogen_production = zeros(length(y_offshore),length(z));

            sumCF_wind_offshore = sum(CF_offshore);

            for i in 1:length(y_offshore)
            annual_invest_LCOE[i]=cost_wind_offshore*y_offshore[i];
            annual_energy_production[i]=sumCF_wind_offshore*y_offshore[i];
                for k in 1:length(z)
                annual_invest_LCOH[i,k]=annual_invest_LCOE[i]+cost_el*z[k];
                    for t in 1:8760
                        annual_hydrogen_production[i,k]+=min((CF_offshore[t]*y_offshore[i]),z[k])/H2_energy_consumption;
                    end
                end
            end

            LCOE=1000*annual_invest_LCOE./annual_energy_production;
            LCOH=annual_invest_LCOH./annual_hydrogen_production;
            LCOE_min, indicesLCOEmin=findmin(LCOE);
            yonshoremin_LCOE=y_offshore[indicesLCOEmin[1]];

            global minLCOH, indicesLCOHmin=findmin(LCOH);
            global LCOH_min = minLCOH;
            global yoffshoremin_LCOH=y_offshore[indicesLCOHmin[1]];
            global zmin_LCOH=z[indicesLCOHmin[2]];
            delta_wind_offshore=(endvalue_wind_offshore-startvalue_wind_offshore)/(d-1);
            delta_el=(endvalue_el-startvalue_el)/(d-1);

            global startvalue_wind_offshore=max(yoffshoremin_LCOH-delta_wind_offshore,1);
            global endvalue_wind_offshore=min(yoffshoremin_LCOH+delta_wind_offshore,max_cap_wind_offshore);
            global startvalue_el=zmin_LCOH-delta_el;
            global endvalue_el=zmin_LCOH+delta_el;

            if startvalue_wind_offshore < 0
                startvalue_wind_offshore = 0
            end

            if endvalue_wind_offshore < 0
                endvalue_wind_offshore = 0
            end

            if startvalue_el < 0
                startvalue_el = 0
            end

            if endvalue_el < 0
                endvalue_el = 0
            end

            # println("\nMinimum LCOH");
            # println(LCOH_min);
            #
            # println("\nMinimum X");
            # println(xmin_LCOH);
            #
            # println("\nMinimum Y");
            # println(ymin_LCOH);
            #
            # println("\nMinimum Z");
            # println(zmin_LCOH);
            if l==length(dd)
                global annual_hydrogen_production_final = annual_hydrogen_production
            end
        end

        # global annual_hydrogen_production = Float64(0.0);
        #
        # for t in 1:8760
        #     global annual_hydrogen_production +=min((CF_pvplant[t]*xmin_LCOH+CF_onshore[t]*ymin_LCOH),zmin_LCOH)/H2_energy_consumption;
        # end

        rescaling_factor = yoffshoremin_LCOH/max_cap_wind_offshore

        x_min_LCOH_rescaled          = 0;
        y_onshore_min_LCOH_rescaled  = 0;
        y_offshore_min_LCOH_rescaled = yoffshoremin_LCOH/rescaling_factor;
        z_min_LCOH_rescaled		     = zmin_LCOH/rescaling_factor;

        rescaled_annual_hydrogen_production = annual_hydrogen_production_final[indicesLCOHmin[1],indicesLCOHmin[2]]/rescaling_factor

        # global rescaled_annual_hydrogen_production = Float64(0.0);
        #
        # for t in 1:8760
        #     global rescaled_annual_hydrogen_production += min((CF_offshore[t]*y_offshore_min_LCOH_rescaled),z_min_LCOH_rescaled)/H2_energy_consumption;
        # end

    end

    return x_min_LCOH_rescaled, y_onshore_min_LCOH_rescaled, y_offshore_min_LCOH_rescaled, z_min_LCOH_rescaled, LCOH_min, rescaled_annual_hydrogen_production
end

#Funzione Nuova
function rescale_property(low_res_property, loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    @unpack era_year, res, erares, pv_density, plant_area = loptions

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    high_res_property = zeros(size(regions))


    for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        #println(length(colrange))
        #println(size(colrange))
        for i = 1:nlons
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]
            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                (reg == 0 || reg == NOREGION) && continue
                #area = cellarea[c]
                #if exclusion_mask[r,c] > 0
                    high_res_property[r,c] = low_res_property[i,j]
                #end
             end
        end
    end
    return  high_res_property
end

# function calc_area()
#
#
#
# end
