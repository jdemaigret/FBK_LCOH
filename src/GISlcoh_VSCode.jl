
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

function create_loptions()
    path = homedir()
    excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
    region_name = excel_file["Region Name"]
    inputs = excel_file["Inputs"]
    PV_power_dens = inputs[:][1,2]
    PV_plant_area = inputs[:][2,2]
    Onshore_power_dens = inputs[:][3,2]
    Onshore_plant_area = inputs[:][4,2]
    Offshore_power_dens = inputs[:][5,2]
    Offshore_plant_area = inputs[:][6,2]
    exclusion_zones = excel_file["Exclusion Zones"]
    max_popdens_PV = exclusion_zones[:][2,2]
    max_popdens_Onshore = exclusion_zones[:][3,2]
    landtype_exclusion_PV_input = exclusion_zones[:][6:22,2]
    landtype_exclusion_Onshore_input = exclusion_zones[:][6:22,3]
    landtype_exclusion_PV = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    landtype_exclusion_Onshore = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    for i in 1:length(landtype_exclusion_PV)
        if  isequal(landtype_exclusion_PV_input[i],missing)
            landtype_exclusion_PV[i]=typemax(Int8)
        end
        if  isequal(landtype_exclusion_Onshore_input[i],missing)
            landtype_exclusion_Onshore[i]=typemax(Int8)
        end
    end
    protected_exclusion_PV_input = exclusion_zones[:][25:34,2]
    protected_exclusion_Onshore_input = exclusion_zones[:][25:34,3]
    protected_exclusion_PV = [1,2,3,4,5,6,7,8,9,10]
    protected_exclusion_Onshore = [1,2,3,4,5,6,7,8,9,10]
    for i in 1:length(protected_exclusion_PV)
        if  isequal(protected_exclusion_PV_input[i],missing)
            protected_exclusion_PV[i]=typemax(Int8)
        end
        if  isequal(protected_exclusion_Onshore_input[i],missing)
            protected_exclusion_Onshore[i]=typemax(Int8)
        end
    end
    filter(!isequal(typemax(Int8)),landtype_exclusion_PV)
    filter(!isequal(typemax(Int8)),landtype_exclusion_Onshore)
    filter(!isequal(typemax(Int8)),protected_exclusion_PV)
    filter(!isequal(typemax(Int8)),protected_exclusion_Onshore)
    max_depth_offshore = exclusion_zones[:][44,2]
    min_shore_distance_offshore = exclusion_zones[:][45,2]

    lcohoptions() = Dict(
        :gisregion => "$(region_name[:][1,1])",
        :filenamesuffix => "",

        :pv_density => PV_power_dens,
        :plant_area => PV_plant_area,

        :onshore_density => Onshore_power_dens,
        :offshore_density => Offshore_power_dens,

        :area_onshore => Onshore_plant_area,
        :area_offshore => Offshore_plant_area,

        :pv_plant_persons_per_km2 => max_popdens_PV,
        :wind_persons_per_km2 => max_popdens_Onshore,

        :pv_exclude_landtypes => landtype_exclusion_PV,
        :pv_protected_codes => protected_exclusion_PV,

        :max_depth => max_depth_offshore,
        :min_shore_distance => min_shore_distance_offshore,

        :wind_exclude_landtypes => landtype_exclusion_Onshore,
        :wind_protected_codes => protected_exclusion_Onshore,

        :scenarioyear => "ssp2_2050",
        :era_year => 2018,
        :rescale_to_wind_atlas => false,

        :res => 0.01,
        :erares => 0.28125,
        :downsample_masks => 1
    )

    return lcohoptions()
end



function GISlcoh(; savetodisk=false, plotmasks=false, optionlist...)

    imported_options = create_loptions() 

    loptions = LcohOptions(merge(imported_options, optionlist))

    @unpack   gisregion, era_year, filenamesuffix, downsample_masks = loptions

    #questa è dentro giswind
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange, aqueduct, slope = read_datasets(loptions);

    # 2050 Optimistic
        # scenario = "2050 Opt"
        scenario = "2050 Opt ws"
        eff_el   = 0.875;
        CAPEX_el = 133.77;

    # 2050 Pessimistic
        # scenario = "2050 Pess"
        # scenario = "2050 Pess ws"
        # eff_el   = 0.8204;
        # CAPEX_el = 326.27;

    # 2030 Optimistic
        # scenario = "2030 Opt"
        # scenario = "2030 Opt ws"
        # CAPEX_el = 383.7660952;
        # eff_el   = 0.812075691897299

    # 2030 Pessimistic
        # scenario = "2030 Pess"
        # scenario = "2030 Pess ws"
        # CAPEX_el = 688.4291624;
        # eff_el     = 0.754154718239413
        # scenario = "V2"


    path = homedir()
    excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
    exclusion_zones = excel_file["Exclusion Zones"]
    max_slope_pv = exclusion_zones[:][37,2]
    max_slope_onshore = exclusion_zones[:][37,3]
    #max_slope_pv = 5;
    #max_slope_onshore = 20;
    pv_system_efficiency = 0.774;
    min_cf_pv       = 0.10
    min_cf_onshore  = 0.20
    min_cf_offshore = 0.23
    crop_nat_excl_coeff = 0.6;
    WACC_value = "Variable_WACC"
    #WACC_value = 0.08

    mask_altitude_lcoh = create_altitude_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)

    mask_slope_pv_lcoh = create_slope_masks_pv_lcoh(loptions, regions, offshoreregions, slope, max_slope_pv)
    mask_slope_onshore_lcoh = create_slope_masks_onshore_lcoh(loptions, regions, offshoreregions, slope, max_slope_onshore)

    mask_aqueduct_lcoh = create_aqueduct_mask_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, aqueduct; plotmasks=false, downsample=1)

    meanGTI, solarGTI = read_solar_datasets_lcoh(loptions, lonrange, latrange)
    solarGTI = solarGTI*pv_system_efficiency;
    windatlas, meanwind, windspeed = read_wind_datasets_lcoh(loptions, lonrange, latrange)

    resource_quality_mask_pv       = create_resource_quality_mask_pv(loptions, lonrange, latrange, regions, regionlist, solarGTI, min_cf_pv)
    resource_quality_mask_onshore  = create_resource_quality_mask_onshore(loptions, lonrange, latrange, regions, regionlist, solarGTI, windspeed, min_cf_onshore)
    resource_quality_mask_offshore = create_resource_quality_mask_offshore(loptions, lonrange, latrange, regions, regionlist, solarGTI, windspeed, min_cf_offshore, topo, protected, offshoreregions)

    mask_pv_lcoh = create_solar_mask_lcoh(loptions, regions, popdens, land, protected, lonrange, latrange, mask_altitude_lcoh, mask_aqueduct_lcoh, mask_slope_pv_lcoh, resource_quality_mask_pv; plotmasks=plotmasks, downsample=1)
    mask_onshore_lcoh,
    mask_offshore_lcoh = create_wind_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange,mask_altitude_lcoh, mask_aqueduct_lcoh, mask_slope_onshore_lcoh, resource_quality_mask_onshore, resource_quality_mask_offshore; plotmasks=plotmasks, downsample=downsample_masks)
   
    # filename_mask_pv       = in_datafolder("output", "$gisregion mask_pv_matrix (prova plot).txt");
    # filename_mask_onshore  = in_datafolder("output", "$gisregion mask_onshore_matrix (prova plot).txt");
    # filename_mask_offshore = in_datafolder("output", "$gisregion mask_offshore_matrix (prova plot).txt");
    #
    # isfile(filename_mask_pv) && rm(filename_mask_pv);
    # isfile(filename_mask_onshore) && rm(filename_mask_onshore);
    # isfile(filename_mask_offshore) && rm(filename_mask_offshore);
    #
    # writedlm(filename_mask_pv, mask_pv_lcoh, ',');
    # writedlm(filename_mask_onshore, mask_onshore_lcoh, ',');
    # writedlm(filename_mask_offshore, mask_offshore_lcoh, ',');

    # mask_aqueduct_lcoh = create_aqueduct_mask_lcoh(regions, aqueduct)
    #
    # mask_desalination_lcoh,
    # shore              = create_desalination_mask(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)

    #pv_area, onshore_wind_area, offshore_wind_area = calc_area(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh, lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh)

    #mask_lcoh = create_lcoh_mask(mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, mask_aqueduct_lcoh)

    #plotmasks == :onlymasks && return nothing

    area_idonea_pv, area_idonea_onshore, area_idonea_offshore, area_idonea_pv_onshore_PV, area_idonea_pv_onshore_ONSHORE,
    installable_pv, installable_onshore, installable_offshore, installable_hybrid_pv, installable_hybrid_onshore,
    area_tot, area_esclusa = calculate_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, land, crop_nat_excl_coeff)

    mean_pv_dist_1, mean_pv_dist_2, mean_pv_dist_3, mean_pv_dist_4, mean_pv_dist_5,
    mean_onshore_dist_1, mean_onshore_dist_2, mean_onshore_dist_3, mean_onshore_dist_4, mean_onshore_dist_5,
    mean_offshore_dist_1, mean_offshore_dist_2, mean_offshore_dist_3, mean_offshore_dist_4, mean_offshore_dist_5,
    pv_cap_1, pv_cap_2, pv_cap_3, pv_cap_4, pv_cap_5,
    onshore_cap_1, onshore_cap_2, onshore_cap_3, onshore_cap_4, onshore_cap_5,
    offshore_cap_1, offshore_cap_2, offshore_cap_3, offshore_cap_4, offshore_cap_5,
    sum_pv_dist_1, sum_pv_dist_2, sum_pv_dist_3, sum_pv_dist_4, sum_pv_dist_5,
    sum_onshore_dist_1, sum_onshore_dist_2, sum_onshore_dist_3, sum_onshore_dist_4, sum_onshore_dist_5,
    sum_offshore_dist_1, sum_offshore_dist_2, sum_offshore_dist_3, sum_offshore_dist_4, sum_offshore_dist_5,
    avg_pv_dist_1, avg_pv_dist_2, avg_pv_dist_3, avg_pv_dist_4, avg_pv_dist_5,
    avg_onshore_dist_1, avg_onshore_dist_2, avg_onshore_dist_3, avg_onshore_dist_4, avg_onshore_dist_5,
    avg_offshore_dist_1, avg_offshore_dist_2, avg_offshore_dist_3, avg_offshore_dist_4, avg_offshore_dist_5 = calculate_mean_CF(loptions, regions, lonrange, latrange, solarGTI, windspeed, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh,pv_system_efficiency,land, crop_nat_excl_coeff)

    max_sum_cf_pv, min_sum_cf_pv,
    max_sum_cf_onshore, min_sum_cf_onshore,
    max_sum_cf_offshore, min_sum_cf_offshore = find_min_max_cf(loptions, lonrange, latrange, solarGTI, windspeed)

    # LCOE_pv_1, LCOE_pv_2, LCOE_pv_3, LCOE_pv_4, LCOE_pv_5,
    # LCOE_onshore_1, LCOE_onshore_2, LCOE_onshore_3, LCOE_onshore_4, LCOE_onshore_5,
    # LCOE_offshore_1, LCOE_offshore_2, LCOE_offshore_3, LCOE_offshore_4, LCOE_offshore_5 = calculate_LCOE(gisregion, WACC_value)

    ###COMMENTO PER BINARY

    #= pv_region_matrix, onshore_region_matrix, offshore_region_matrix, el_region_matrix, lcoh_region_matrix, production_region_matrix,
    PV_pv_onshore_erares, ONSHORE_pv_onshore_erares, OFFSHORE_pv_onshore_erares, EL_pv_onshore_erares, LCOH_pv_onshore_erares, PROD_pv_onshore_erares,
    PV_pv_erares, ONSHORE_pv_erares, OFFSHORE_pv_erares, EL_pv_erares, LCOH_pv_erares, PROD_pv_erares,
    PV_onshore_erares, ONSHORE_onshore_erares, OFFSHORE_onshore_erares, EL_onshore_erares, LCOH_onshore_erares, PROD_onshore_erares,
    PV_offshore_erares, ONSHORE_offshore_erares, OFFSHORE_offshore_erares, EL_offshore_erares, LCOH_offshore_erares, PROD_offshore_erares,
    SUM_CF_PV, SUM_CF_ONSHORE, SUM_CF_OFFSHORE, sumcfpv_region_matrix, sumcfonshore_region_matrix, sumcfoffshore_region_matrix,
    lcoe_nocurt_region_matrix = calc_LCOH_CAPS_PROD(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh,
                                                                            lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh,
                                                                            meanGTI, WACC_value, land, crop_nat_excl_coeff, eff_el, CAPEX_el, scenario)
    =#
    # println("Saving the big matrices...")
    filename_pv       = in_datafolder("output", "$gisregion pv_matrix ($scenario).txt");
    filename_onshore  = in_datafolder("output", "$gisregion onshore_matrix ($scenario).txt");
    filename_offshore = in_datafolder("output", "$gisregion offshore_matrix ($scenario).txt");
    filename_el       = in_datafolder("output", "$gisregion el_matrix ($scenario).txt");
    filename_lcoh     = in_datafolder("output", "$gisregion lcoh_matrix ($scenario).txt");
    filename_prod     = in_datafolder("output", "$gisregion prod_matrix ($scenario).txt");
    filename_lcoe     = in_datafolder("output", "$gisregion lcoe_matrix ($scenario).txt");
   
    isfile(filename_pv) && rm(filename_pv);
    isfile(filename_onshore) && rm(filename_onshore);
    isfile(filename_offshore) && rm(filename_offshore);
    isfile(filename_el) && rm(filename_el);
    isfile(filename_lcoh) && rm(filename_lcoh);
    isfile(filename_prod) && rm(filename_prod);
    isfile(filename_lcoe) && rm(filename_lcoe);
   
    # writedlm(filename_pv, pv_region_matrix, ',');
    # writedlm(filename_onshore, onshore_region_matrix, ',');
    # writedlm(filename_offshore, offshore_region_matrix, ',');
    # writedlm(filename_el, el_region_matrix, ',');
    # writedlm(filename_lcoh, lcoh_region_matrix, ',');
    # writedlm(filename_prod, production_region_matrix, ',');
    # writedlm(filename_lcoe, lcoe_nocurt_region_matrix, ',');

    # println("Reading big matrices....")
    #
    # pv_region_matrix          = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion pv_matrix ($scenario).txt" ,',')
    # onshore_region_matrix     = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion onshore_matrix ($scenario).txt" ,',')
    # offshore_region_matrix    = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion offshore_matrix ($scenario).txt" ,',')
    # el_region_matrix          = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion el_matrix ($scenario).txt" ,',')
    # lcoh_region_matrix        = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion lcoh_matrix ($scenario).txt" ,',')
    # production_region_matrix  = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion prod_matrix ($scenario).txt" ,',')
    # lcoe_nocurt_region_matrix = readdlm("D:\\Julia_Jacopo\\GISData\\Result Matrices\\$scenario\\$gisregion\\$gisregion prod_matrix ($scenario).txt" ,',')

    ###COMMENTO PER BINARY

    #= vectorized_lcoh_matrix                 = vec(lcoh_region_matrix);
    vectorized_prod_matrix                 = vec(production_region_matrix);
    vectorized_pv_matrix                   = vec(pv_region_matrix)
    vectorized_onshore_matrix              = vec(onshore_region_matrix)
    vectorized_offshore_matrix             = vec(offshore_region_matrix)
    vectorized_el_matrix                   = vec(el_region_matrix)
    vectorized_sumcfpv_region_matrix       = vec(sumcfpv_region_matrix)
    vectorized_sumcfonshore_region_matrix  = vec(sumcfonshore_region_matrix)
    vectorized_sumcfoffshore_region_matrix = vec(sumcfoffshore_region_matrix)
    vectorized_lcoe_nocurt_region_matrix   = vec(lcoe_nocurt_region_matrix)

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
            vectorized_lcoe_nocurt_region_matrix[s]   = deleted_index
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
    vectorized_lcoe_nocurt_region_matrix_wo_0   = filter(!isequal(deleted_index), vectorized_lcoe_nocurt_region_matrix);
    
    FLOH_vector              = zeros(length(vectorized_lcoh_matrix_wo_0));
    ENERGY_PROD_vector       = zeros(length(vectorized_lcoh_matrix_wo_0));
    ELECTROLYZER_CONS_vector = zeros(length(vectorized_lcoh_matrix_wo_0));
    CURTAILMENT_vector       = zeros(length(vectorized_lcoh_matrix_wo_0));

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

    results_matrix = [vectorized_lcoh_matrix_wo_0 vectorized_prod_matrix_wo_0 vectorized_pv_matrix_wo_0 vectorized_onshore_matrix_wo_0 vectorized_offshore_matrix_wo_0 vectorized_el_matrix_wo_0 FLOH_vector vectorized_sumcfpv_region_matrix_wo_0 vectorized_sumcfonshore_region_matrix_wo_0 vectorized_sumcfoffshore_region_matrix_wo_0 CURTAILMENT_vector] =#
    charachteristic_cf_curves = [mean_pv_dist_1 mean_pv_dist_2 mean_pv_dist_3 mean_pv_dist_4 mean_pv_dist_5 mean_onshore_dist_1 mean_onshore_dist_2 mean_onshore_dist_3 mean_onshore_dist_4 mean_onshore_dist_5 mean_offshore_dist_1 mean_offshore_dist_2 mean_offshore_dist_3 mean_offshore_dist_4 mean_offshore_dist_5]
    max_and_min_CF = [min_sum_cf_pv max_sum_cf_pv min_sum_cf_onshore max_sum_cf_onshore min_sum_cf_offshore max_sum_cf_offshore]
    capacities = [pv_cap_1 pv_cap_2 pv_cap_3 pv_cap_4 pv_cap_5 onshore_cap_1 onshore_cap_2 onshore_cap_3 onshore_cap_4 onshore_cap_5 offshore_cap_1 offshore_cap_2 offshore_cap_3 offshore_cap_4 offshore_cap_5]
    cf_sums = [sum_pv_dist_1 sum_pv_dist_2 sum_pv_dist_3 sum_pv_dist_4 sum_pv_dist_5 sum_onshore_dist_1 sum_onshore_dist_2 sum_onshore_dist_3 sum_onshore_dist_4 sum_onshore_dist_5 sum_offshore_dist_1 sum_offshore_dist_2 sum_offshore_dist_3 sum_offshore_dist_4 sum_offshore_dist_5]
    cf_avg =[avg_pv_dist_1 avg_pv_dist_2 avg_pv_dist_3 avg_pv_dist_4 avg_pv_dist_5 avg_onshore_dist_1 avg_onshore_dist_2 avg_onshore_dist_3 avg_onshore_dist_4 avg_onshore_dist_5 avg_offshore_dist_1 avg_offshore_dist_2 avg_offshore_dist_3 avg_offshore_dist_4 avg_offshore_dist_5]

    area_idonea_pv_totale      = area_idonea_pv + area_idonea_pv_onshore_PV;
    area_idonea_onshore_totale = area_idonea_onshore + area_idonea_pv_onshore_ONSHORE;
    installable_pv_totale      = installable_pv + installable_hybrid_pv;
    installable_onshore_totale =  installable_onshore + installable_hybrid_onshore;
    areas = [area_idonea_pv, area_idonea_onshore, area_idonea_offshore, area_idonea_pv_onshore_PV, area_idonea_pv_onshore_ONSHORE, area_idonea_pv_totale, area_idonea_onshore_totale, area_tot, 0, area_esclusa]
    potentials = [installable_pv, installable_onshore, installable_offshore, installable_hybrid_pv, installable_hybrid_onshore, installable_pv_totale, installable_onshore_totale,  0, 0, 0]
    # LCOE = [LCOE_pv_1 LCOE_pv_2 LCOE_pv_3 LCOE_pv_4 LCOE_pv_5 LCOE_onshore_1 LCOE_onshore_2 LCOE_onshore_3 LCOE_onshore_4 LCOE_onshore_5 LCOE_offshore_1 LCOE_offshore_2 LCOE_offshore_3 LCOE_offshore_4 LCOE_offshore_5 ]
    println("...saving data...");
    mkpath(in_datafolder("output"));
    #header0 = ["LCOH (USD/kg)" "Prod (kg/year)" "PV (kW)" "Onshore (kW)" "Offshore (kW)" "Electrolyzer (kW)" "FLOH (h)" "Sum_CF_PV (h)" "Sum_CF_Onshore (h)" "Sum_CF_Offshore (h)" "Curtailment (h)"]
    header1 = ["$gisregion PV 1" "$gisregion PV 2" "$gisregion PV 3" "$gisregion PV 4" "$gisregion PV 5" "$gisregion Onshore 1" "$gisregion Onshore 2" "$gisregion Onshore 3" "$gisregion Onshore 4" "$gisregion Onshore 5" "$gisregion Offshore 1" "$gisregion Offshore 2" "$gisregion Offshore 3" "$gisregion Offshore 4" "$gisregion Offshore 5" ]
    header2 = ["$gisregion Min CF PV" "$gisregion Max CF PV" "$gisregion Min CF Onshore" "$gisregion Max CF Onshore" "$gisregion Min CF Offshore" "$gisregion Max CF Offshore"]
    header4 = ["" "Eligible Area [km2]" "Installable Potential [GW]"]
    header4a = ["PV Only", "Onshore Only", "Offshore", "Hybrid PV", "Hybrid Onshore", "PV Total", "Onshore Total", "$gisregion Total Area", "", "Excluded Land Area"]
    # header5 = ["$gisregion LCOE PV 1" "$gisregion LCOE PV 2" "$gisregion LCOE PV 3" "$gisregion LCOE PV 4" "$gisregion LCOE PV 5" "$gisregion LCOE Onshore 1" "$gisregion LCOE Onshore 2" "$gisregion LCOE Onshore 3" "$gisregion LCOE Onshore 4" "$gisregion LCOE Onshore 5" "$gisregion LCOE Offshore 1" "$gisregion LCOE Offshore 2" "$gisregion LCOE Offshore 3" "$gisregion LCOE Offshore 4" "$gisregion LCOE Offshore 5" ]
    
    areas_potential_mat = [header4a areas potentials]
    #filename0 = in_datafolder("output", "$gisregion results_matrix ($scenario).txt");
    #= filename1 = in_datafolder("output", "$gisregion charachteristic_cf_curves.txt");
    filename2 = in_datafolder("output", "$gisregion Max and Min CF.txt");
    filename3 = in_datafolder("output", "$gisregion Caps CF sums and average.txt");
    filename4 = in_datafolder("output", "$gisregion Areas and Potential.txt"); =#

    

    filename1 = "$path\\GISData\\output\\$gisregion charachteristic_cf_curves.txt"
    filename2 = "$path\\GISData\\output\\$gisregion Max and Min CF.txt"
    filename3 = "$path\\GISData\\output\\$gisregion Caps CF sums and average.txt"
    filename4 = "$path\\GISData\\output\\$gisregion Areas and Potential.txt"
    # filename5 = in_datafolder("output", "$gisregion LCOEs (prima mask).txt");
    # filename7 = in_datafolder("output", "$gisregion GTI  (prova plot).txt")
   
    #isfile(filename0) && rm(filename0);
    isfile(filename1) && rm(filename1);
    isfile(filename2) && rm(filename2);
    isfile(filename3) && rm(filename3);
    isfile(filename4) && rm(filename4);
 
    # isfile(filename5) && rm(filename5);
    # isfile(filename7) && rm(filename7);
    # writedlm(filename0, [header0 ; results_matrix], '\t');
    writedlm(filename1, [header1 ; charachteristic_cf_curves], ',');
    writedlm(filename2, [header2 ; max_and_min_CF], '\t');
    writedlm(filename3, [header1 ; capacities ; cf_sums ; cf_avg], ',');
    writedlm(filename4, [header4 ; areas_potential_mat], '\t');
    # # writedlm(filename5, [header5 ; LCOE], ',')
    # writedlm(filename7, solarGTI[:,1:1,1:1]);
 
    
    Evergreen_Needleleaf_Forests, Evergreen_Broadleaf_Forests, Deciduous_Needleleaf_Forests, Deciduous_Broadleaf_Forests,
    Mixed_Forests, Closed_Shrublands, Open_Shrublands, Woody_Savannas, Savannas, Grasslands, Permanent_Wetlands, Croplands,
    Urban, Cropland_Natural,  Snow_Ice, Barren = calculate_landtype_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, land)
    land_types_names = ["Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests", "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests", "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands", "Permanent_Wetlands", "Croplands", "Urban", "Cropland_Natural", "Snow_Ice", "Barren"]
    land_types_values =  [Evergreen_Needleleaf_Forests, Evergreen_Broadleaf_Forests, Deciduous_Needleleaf_Forests, Deciduous_Broadleaf_Forests, Mixed_Forests, Closed_Shrublands, Open_Shrublands, Woody_Savannas, Savannas, Grasslands, Permanent_Wetlands, Croplands, Urban, Cropland_Natural, Snow_Ice, Barren]
    land_types_matrix = [land_types_names land_types_values]

    filename_landtype = in_datafolder("output", "$gisregion land type areas.txt");
    isfile(filename_landtype) && rm(filename_landtype);


    Strict_Nature_Reserve, Wilderness_Area, National_Park, Natural_Monument, Habitat_Species_Management,
    Protected_Landscape_Seascape, Managed_Resource_Protected_Area, Not_Reported, Not_Applicable, Not_Assigned = calculate_protected_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, protected)
    protected_names = ["Strict_Nature_Reserve", "Wilderness_Area", "National_Park", "Natural_Monument", "Habitat_Species_Management", "Protected_Landscape_Seascape", "Managed_Resource_Protected_Area", "Not_Reported", "Not_Applicable", "Not_Assigned"]
    protected_values =  [Strict_Nature_Reserve, Wilderness_Area, National_Park, Natural_Monument, Habitat_Species_Management, Protected_Landscape_Seascape, Managed_Resource_Protected_Area, Not_Reported, Not_Applicable, Not_Assigned]
    protected_matrix = [protected_names protected_values]

    filename_protected = in_datafolder("output", "$gisregion protected areas ($scenario).txt");
    isfile(filename_protected) && rm(filename_protected);


    @unpack pv_plant_persons_per_km2, wind_persons_per_km2 = loptions
    popdens_excl_area_PV, popdens_excl_area_Onshore = calculate_popdens_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, popdens)
    popdens_names = ["Excluded PV ($pv_plant_persons_per_km2)" "Excluded Onshore ($wind_persons_per_km2)"]
    popdens_values =  [popdens_excl_area_PV popdens_excl_area_Onshore]
    popdens_matrix = [popdens_names ; popdens_values]

    filename_popdens = in_datafolder("output", "$gisregion popdens areas ($scenario).txt");
    isfile(filename_popdens) && rm(filename_popdens);


    slope_excl_area_PV, slope_excl_area_Onshore =  calculate_slope_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, slope, max_slope_pv, max_slope_onshore)
    slope_names = ["Excluded PV ($(max_slope_pv)%)" "Excluded Onshore ($(max_slope_onshore)%)"]
    slope_values =  [slope_excl_area_PV slope_excl_area_Onshore]
    slope_matrix = [slope_names ; slope_values]

    filename_slope = in_datafolder("output", "$gisregion slope areas ($scenario).txt");
    isfile(filename_slope) && rm(filename_slope);


    aqueduct_excl_area = calculate_aqueduct_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, mask_aqueduct_lcoh)
    aqueduct_names = ["Excluded area"]
    aqueduct_values =  [aqueduct_excl_area]
    aqueduct_matrix = [aqueduct_names ; aqueduct_values]

    filename_aqueduct = in_datafolder("output", "$gisregion aqueduct areas ($scenario).txt");
    isfile(filename_aqueduct) && rm(filename_aqueduct);


    resource_excl_area_PV, resource_excl_area_Onshore, resource_excl_area_Offshore =  calculate_pv_resource_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, resource_quality_mask_pv, resource_quality_mask_onshore, resource_quality_mask_offshore)
    resource_names = ["Excluded PV ($(min_cf_pv)%)" "Excluded Onshore ($(min_cf_onshore)%)" "Excluded Onshore ($(min_cf_offshore)%)"]
    resource_values =  [resource_excl_area_PV resource_excl_area_Onshore resource_excl_area_Offshore]
    resource_matrix = [resource_names ; resource_values]

    filename_resource = in_datafolder("output", "$gisregion resource areas ($scenario).txt");
    isfile(filename_resource) && rm(filename_resource);

    # writedlm(filename_landtype, land_types_matrix, '\t')
    # writedlm(filename_protected, protected_matrix, '\t')
    # writedlm(filename_popdens, popdens_matrix, '\t')
    # writedlm(filename_slope, slope_matrix, '\t')
    # writedlm(filename_aqueduct, aqueduct_matrix, '\t')
    # writedlm(filename_resource, resource_matrix, '\t')

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
#Nuova funzione
function create_solar_mask_lcoh(loptions, regions, popdens, land, protected, lonrange, latrange, mask_altitude_lcoh, mask_aqueduct_lcoh, mask_slope_pv_lcoh, resource_quality_mask_pv; plotmasks=false, downsample=1)

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

    km_per_degree = π*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    path = homedir()
    excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
    exclusion_zones = excel_file["Exclusion Zones"]
    include_desal = exclusion_zones[:][40,2]

    if include_desal=="yes"
        mask_pv_lcoh = (popdens .< pv_plant_persons_per_km2) .& goodland .& .!pv_protected_area .& mask_slope_pv_lcoh .& mask_aqueduct_lcoh #.& mask_altitude_lcoh
    elseif include_desal=="no"
        mask_pv_lcoh = (popdens .< pv_plant_persons_per_km2) .& goodland .& .!pv_protected_area .& mask_slope_pv_lcoh
    else
        println("Statement for water stress inclusion unclear: use 'yes' or 'no'.")
    end
    #mask_pv_lcoh = (popdens .< pv_plant_persons_per_km2) .& goodland .& .!pv_protected_area .& mask_slope_pv_lcoh .& mask_aqueduct_lcoh #.& mask_altitude_lcoh
    
    return mask_pv_lcoh
end
#Nuova funzione
function create_wind_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, mask_altitude_lcoh, mask_aqueduct_lcoh, mask_slope_onshore_lcoh, resource_quality_mask_onshore, resource_quality_mask_offshore; plotmasks=false, downsample=1)
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
    km_per_degree = π*2*6371/360
    #disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    #gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    #mask_onshoreA = gridA .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    #mask_onshoreB = (gridB .& .!gridA) .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    path = homedir()
    excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
    exclusion_zones = excel_file["Exclusion Zones"]
    include_desal = exclusion_zones[:][40,2]

    if include_desal=="yes"
        mask_onshore_lcoh = (popdens .< wind_persons_per_km2) .& goodland .& .!wind_protected_area .& mask_slope_onshore_lcoh .& mask_aqueduct_lcoh
    elseif include_desal=="no"
        mask_onshore_lcoh = (popdens .< wind_persons_per_km2) .& goodland .& .!wind_protected_area .& mask_slope_onshore_lcoh
    end
    #mask_onshore_lcoh = (popdens .< wind_persons_per_km2) .& goodland .& .!wind_protected_area .& mask_slope_onshore_lcoh .& mask_aqueduct_lcoh
    
    # shoreline mask for offshore wind
    disk = diskfilterkernel(min_shore_distance/km_per_degree/res)
    shore = (imfilter(regions .> 0, disk) .> 1e-6)

    # all mask conditions
    # mask_offshore_lcoh = .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!wind_protected_area
    # println("typeof(resource_quality_mask_offshore)")
    # println(typeof(resource_quality_mask_offshore))
    #
    # println("unique(resource_quality_mask_offshore)")
    # println(unique(resource_quality_mask_offshore))

    mask_offshore_lcoh = .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!wind_protected_area #.& resource_quality_mask_offshore
    
    #= mask_offshore_lcoh =  (topo .> -max_depth) .& (offshoreregions .> 0) .& .!wind_protected_area #.& resource_quality_mask_offshore =#

    return mask_onshore_lcoh, mask_offshore_lcoh
end

function create_aqueduct_mask_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, aqueduct; plotmasks=false, downsample=1)
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
    path = homedir()
    excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
    exclusion_zones = excel_file["Exclusion Zones"]
    min_desalination_distance = exclusion_zones[:][48,2]
    #min_desalination_distance = 50
    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(min_desalination_distance/km_per_degree/res)
    shore = (imfilter(regions .== 0, disk) .> 1e-6)

    mask_desalination_lcoh = (shore) .& (regions .> 0) #.& .!wind_protected_area
    println("Creating aquedcut masks...\n")

    goodland = (regions .> 0)

    BWS_values = sort(unique(aqueduct))
    #println(BWS_values)

    #mask_BWS_lcoh = zeros(size(aqueduct))
    BWS_threshold = exclusion_zones[:][41,2]

    #BWS_threshold = 1
    mask_BWS_lcoh = goodland .& (aqueduct .<= BWS_threshold)
    mask_desalination_lcoh = (shore) .& (regions .> 0) #.& (aqueduct .> BWS_threshold)
    mask_aqueduct_lcoh = mask_BWS_lcoh .| mask_desalination_lcoh

    return mask_aqueduct_lcoh
end
function create_altitude_masks_lcoh(loptions, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, wind_exclude_landtypes, wind_protected_codes, wind_persons_per_km2,
                min_shore_distance, max_depth, filenamesuffix = loptions

    println("Creating altitude masks...\n")

    topo_land = (regions .> 0) .& topo
    # topo_land = zeros(size(topo))
    #
    # println("unique(regions)")
    # println(unique(regions))
    # for i in 1:size(topo)[1]
    #     for j in 1:size(topo)[2]
    #         if topo[i,j] > 0
    #             topo_land[i,j] = topo[i,j]
    #         end
    #     end
    # end


    max_height = 2000;
    mask_altitude_lcoh = (topo .< max_height) .& (regions .> 0) #Qua io posso costruire!

    return mask_altitude_lcoh
end
function create_slope_masks_pv_lcoh(loptions, regions, offshoreregions, slope, max_slope_pv)

    println("Creating PV slope masks...\n")

    mask_slope_pv_lcoh = (slope .< max_slope_pv) .& (regions .> 0)

    return mask_slope_pv_lcoh
end
function create_slope_masks_onshore_lcoh(loptions, regions, offshoreregions, slope, max_slope_onshore)

    println("Creating Onshore slope masks...\n")

    mask_slope_onshore_lcoh = (slope .< max_slope_onshore) .& (regions .> 0)

    return mask_slope_onshore_lcoh
end
function create_resource_quality_mask_pv(loptions, lonrange, latrange, regions, regionlist, solarGTI, min_cf_pv)
    # println("Creating PV resource quality masks...\n")
    global yearlength, nlons, nlats = size(solarGTI)

    SUM_CF_PV       = reshape(sum(solarGTI, dims=1)/8760, (nlons, nlats));
    SUM_CF_PV_res = rescale_property(SUM_CF_PV,  loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    resource_quality_mask_pv = (SUM_CF_PV_res .> min_cf_pv) .& (regions .> 0)

    return resource_quality_mask_pv
end
function create_resource_quality_mask_onshore(loptions, lonrange, latrange, regions, regionlist, solarGTI, windspeed, min_cf_onshore)
    # println("Creating Onshore resource quality masks...\n")
    global yearlength, nlons, nlats = size(solarGTI)

    SUM_CF_ONSHORE  = reshape(sum(speed2capacityfactor_onshore.(windspeed), dims=1)/8760, (nlons, nlats));
    SUM_CF_ONSHORE_res = rescale_property(SUM_CF_ONSHORE,  loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    resource_quality_mask_onshore = (SUM_CF_ONSHORE_res .> min_cf_onshore) .& (regions .> 0)

    return resource_quality_mask_onshore
end
function create_resource_quality_mask_offshore(loptions, lonrange, latrange, regions, regionlist, solarGTI, windspeed, min_cf_offshore, topo, protected, offshoreregions)

    @unpack res,min_shore_distance, max_depth, wind_protected_codes = loptions
    # println("Creating Offshore resource quality masks...\n")
    global yearlength, nlons, nlats = size(solarGTI)

    SUM_CF_OFFSHORE = reshape(sum(speed2capacityfactor_offshore.(windspeed), dims=1)/8760, (nlons, nlats));
    SUM_CF_OFFSHORE_res = rescale_property(SUM_CF_OFFSHORE,  loptions,  lonrange , latrange, regionlist, solarGTI, regions)

    wind_protected_area = zeros(Bool, size(protected))
    for i in wind_protected_codes
        wind_protected_area[protected .== i] .= true
    end

    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(min_shore_distance/km_per_degree/res)
    shore = (imfilter(regions .> 0, disk) .> 1e-6)

    resource_quality_mask_offshore = .!shore .& (SUM_CF_OFFSHORE_res .> min_cf_offshore) .& (offshoreregions .> 0)
    resource_quality_mask_offshore = (SUM_CF_OFFSHORE_res .> min_cf_offshore)
    # resource_quality_mask_offshore = .!shore .& (SUM_CF_OFFSHORE_res .> min_cf_offshore)
    # resource_quality_mask_offshore =  (SUM_CF_OFFSHORE_res .> min_cf_offshore) .& (offshoreregions .> 0)

    return resource_quality_mask_offshore
end

#Funzione nuova
function calc_LCOH_CAPS_PROD(loptions, regions, solarGTI, offshoreregions, regionlist, mask_pv_lcoh,
                            lonrange, latrange, windatlas, meanwind, windspeed, mask_onshore_lcoh, mask_offshore_lcoh,
                            meanGTI, WACC_value, land, crop_nat_excl_coeff, eff_el, CAPEX_el, scenario)
  
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    @unpack   gisregion = loptions

    global yearlength, nlons, nlats = size(solarGTI)

    # Matrici bassa risoluzione
    # PV and ONSHORE
    PV_pv_onshore_erares       = zeros(nlons, nlats);
    ONSHORE_pv_onshore_erares  = zeros(nlons, nlats);
    OFFSHORE_pv_onshore_erares = zeros(nlons, nlats);
    EL_pv_onshore_erares       = zeros(nlons, nlats);
    LCOH_pv_onshore_erares     = zeros(nlons, nlats);
    PROD_pv_onshore_erares     = zeros(nlons, nlats);
    LCOE_pv_onshore_erares     = zeros(nlons, nlats);

    # PV
    PV_pv_erares               = zeros(nlons, nlats);
    ONSHORE_pv_erares          = zeros(nlons, nlats);
    OFFSHORE_pv_erares         = zeros(nlons, nlats);
    EL_pv_erares               = zeros(nlons, nlats);
    LCOH_pv_erares             = zeros(nlons, nlats);
    PROD_pv_erares             = zeros(nlons, nlats);
    LCOE_pv_erares             = zeros(nlons, nlats);

    # ONSHORE
    PV_onshore_erares          = zeros(nlons, nlats);
    ONSHORE_onshore_erares     = zeros(nlons, nlats);
    OFFSHORE_onshore_erares    = zeros(nlons, nlats);
    EL_onshore_erares          = zeros(nlons, nlats);
    LCOH_onshore_erares        = zeros(nlons, nlats);
    PROD_onshore_erares        = zeros(nlons, nlats);
    LCOE_onshore_erares        = zeros(nlons, nlats);

    # OFFSHORE
    PV_offshore_erares         = zeros(nlons, nlats);
    ONSHORE_offshore_erares    = zeros(nlons, nlats);
    OFFSHORE_offshore_erares   = zeros(nlons, nlats);
    EL_offshore_erares         = zeros(nlons, nlats);
    LCOH_offshore_erares       = zeros(nlons, nlats);
    PROD_offshore_erares       = zeros(nlons, nlats);
    LCOE_offshore_erares       = zeros(nlons, nlats);

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
    lcoe_nocurt_region_matrix   = zeros(size(regions));
  

    SUM_CF_PV       = reshape(sum(solarGTI, dims=1), (nlons, nlats));
   
    SUM_CF_ONSHORE  = reshape(sum(speed2capacityfactor_onshore.(windspeed), dims=1), (nlons, nlats));
    SUM_CF_OFFSHORE = reshape(sum(speed2capacityfactor_offshore.(windspeed), dims=1), (nlons, nlats));
   

    CAPEX_and_WACC = readdlm("WACC_CAPEX Versione 9.csv" ,',')

    

    # nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
    # CAPEX_and_WACC_country_index = findfirst(isequal("$gisregion"), nomi_stringa)[2];

    # ReadMe for WACC_CAPEX Versione 9.csv: The numbers indicate the column.
    # OPTIMISTIC
    # 2020 ([1]PV CAPEX	   [2]Onshore CAPEX	[3]Offshore CAPEX)
    # 2030 ([4]PV CAPEX	   [5]Onshore CAPEX	[6]Offshore CAPEX)
    # 2050 ([7]PV CAPEX    [8]Onshore CAPEX	[9]Offshore CAPEX)
    # WACC ([10]PV WACC    [11]Onshore WACC [12]Offshore WACC)
    # PESSIMISTIC
    # 2020 ([13]PV CAPEX   [14]Onshore CAPEX [15]Offshore CAPEX)
    # 2030 ([16]PV CAPEX   [17]Onshore CAPEX [18]Offshore CAPEX)
    # 2050 ([19]PV CAPEX   [20]Onshore CAPEX [21]Offshore CAPEX)
    # WACC ([22]PV WACC    [23]Onshore WACC  [24]Offshore WACC)
    #
    # CAPEX_pv               = CAPEX_and_WACC[30,16];
    # CAPEX_wind_onshore     = CAPEX_and_WACC[30,17];
    # CAPEX_wind_offshore    = CAPEX_and_WACC[30,18];
    #
    # if WACC_value == "Variable_WACC"
    #     WACC_pv       = CAPEX_and_WACC[30,16];
    #     WACC_onshore  = CAPEX_and_WACC[30,17];
    #     WACC_offshore = CAPEX_and_WACC[30,18];
    # else
    #     WACC_pv       = WACC_value
    #     WACC_onshore  = WACC_value
    #     WACC_offshore = WACC_value
    # end

    CAPEX_pv               = 343.36
    CAPEX_wind_onshore     = 905.48
    CAPEX_wind_offshore    = 1494.22

    WACC_pv       = 0.0713
    WACC_onshore  = 0.0713
    WACC_offshore = 0.0713


    # println("CAPEX_pv: $CAPEX_pv")
    # println("CAPEX_wind_onshore: $CAPEX_wind_onshore")
    # println("CAPEX_wind_offshore: $CAPEX_wind_offshore")
    # println("WACC_pv: $WACC_pv")
    # println("WACC_onshore: $WACC_onshore")
    # println("WACC_offshore: $WACC_offshore")
    # println("scenario: $scenario")
    # println("CAPEX_el: $CAPEX_el")
    # println("eff_el: $eff_el")

    #MA SERVONO?
    #_, _, meanwind_allyears, _ = annualwindindex(loptions)
    #@assert size(meanwind_allyears) == size(meanwind)

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

    # if rescale_to_wind_atlas
    #     println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
    #     println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    # end
    # pv_area = zeros(nlons,nlats)
    # onshore_wind_area = zeros(nlons,nlats)
    # offshore_wind_area = zeros(nlons,nlats)
    hybrid_system_mask = mask_pv_lcoh .& mask_onshore_lcoh;
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
    #Basta scrivere nlats
    # l = Threads.ReentrantLock()
    # l = Threads.SpinLock()
   
    for j in 1:1
        # Threads.lock(l)
        # println((thread=Threads.threadid(), iteration=j))
        #
        # Threads.unlock(l)
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
        #Basta scrivere nlons
        for i in 1:1
            contatore = contatore + 1;

            println("Pixel opt: $contatore/$num_tot_pixel")
            # println("Pixel=$contatore/$num_tot_pixel thread=$(Threads.threadid())")
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
            ((length(rowrange)==0) || (length(colrange)==0)) && continue
            ((last(rowrange) == 0) || (last(colrange) == 0)) && continue
            rowrange_wo_zero = filter(!isequal(0), rowrange);
            colrange_wo_zero = filter(!isequal(0), colrange);

          
            if findfirst(!isequal(0), hybrid_system_mask[rowrange_wo_zero,colrange_wo_zero]) != nothing

                PV_pv_onshore_erares[i,j],
                ONSHORE_pv_onshore_erares[i,j],
                OFFSHORE_pv_onshore_erares[i,j],
                EL_pv_onshore_erares[i,j],
                LCOH_pv_onshore_erares[i,j],
                PROD_pv_onshore_erares[i,j],
                LCOE_pv_onshore_erares[i,j] = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                GTI, area, wind_onshore, wind_offshore, 1,
                                                CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, WACC_pv, WACC_onshore, WACC_offshore, eff_el,CAPEX_el); # Il numero sarebbe il type_min = 1
            end

            if  findfirst(!isequal(0), mask_pv_lcoh[rowrange_wo_zero,colrange_wo_zero]) != nothing

                PV_pv_erares[i,j],
                ONSHORE_pv_erares[i,j],
                OFFSHORE_pv_erares[i,j],
                EL_pv_erares[i,j],
                LCOH_pv_erares[i,j],
                PROD_pv_erares[i,j],
                LCOE_pv_erares[i,j]         = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                GTI, area, wind_onshore, wind_offshore, 2,
                                                CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, WACC_pv, WACC_onshore, WACC_offshore, eff_el, CAPEX_el); # Il numero sarebbe il type_min  = 2
            end

            if  findfirst(!isequal(0), mask_onshore_lcoh[rowrange_wo_zero,colrange_wo_zero]) != nothing

                PV_onshore_erares[i,j],
                ONSHORE_onshore_erares[i,j],
                OFFSHORE_onshore_erares[i,j],
                EL_onshore_erares[i,j],
                LCOH_onshore_erares[i,j],
                PROD_onshore_erares[i,j],
                LCOE_onshore_erares[i,j]   = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                GTI, area, wind_onshore, wind_offshore, 3,
                                                CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, WACC_pv, WACC_onshore, WACC_offshore, eff_el, CAPEX_el); # Il numero sarebbe il type_min  = 3
            end

            if  findfirst(!isequal(0), mask_offshore_lcoh[rowrange_wo_zero,colrange_wo_zero]) != nothing

                PV_offshore_erares[i,j],
                ONSHORE_offshore_erares[i,j],
                OFFSHORE_offshore_erares[i,j],
                EL_offshore_erares[i,j],
                LCOH_offshore_erares[i,j],
                PROD_offshore_erares[i,j],
                LCOE_offshore_erares[i,j]  = find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                                                GTI, area, wind_onshore, wind_offshore, 4,
                                                CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, WACC_pv, WACC_onshore, WACC_offshore, eff_el, CAPEX_el); # Il numero sarebbe il type_min  = 4
            end

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c] #region all'indice r, c
                    #(reg == 0 || reg == NOREGION) && continue
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    #davanti all if c'era un @views
                    if hybrid_system_mask[r,c]>0#mask_pv_lcoh[r,c]>0 && mask_onshore_lcoh[r,c]>0  #sia pv che onshore
                        if land[r,c]!=14 # NO RESCALE
                            pv_region_matrix[r,c]            = PV_pv_onshore_erares[i,j];
                            onshore_region_matrix[r,c]       = ONSHORE_pv_onshore_erares[i,j];
                            offshore_region_matrix[r,c]      = OFFSHORE_pv_onshore_erares[i,j];
                            el_region_matrix[r,c]            = EL_pv_onshore_erares[i,j];
                            lcoh_region_matrix[r,c]          = LCOH_pv_onshore_erares[i,j]#min(LCOH_pv_onshore_erares[i,j], 10);
                            production_region_matrix[r,c]    = PROD_pv_onshore_erares[i,j];
                            sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                            sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                            sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                            lcoe_nocurt_region_matrix[r,c]   = LCOE_pv_onshore_erares[i,j];
                        else
                            if (PV_pv_onshore_erares[i,j]==area*pv_density) && (ONSHORE_pv_onshore_erares[i,j]<0.01*area*onshore_density)

                                pv_region_matrix[r,c]            = PV_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                onshore_region_matrix[r,c]       = ONSHORE_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                offshore_region_matrix[r,c]      = OFFSHORE_pv_onshore_erares[i,j];
                                el_region_matrix[r,c]            = EL_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                lcoh_region_matrix[r,c]          = LCOH_pv_onshore_erares[i,j]#min(LCOH_pv_onshore_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                                lcoe_nocurt_region_matrix[r,c]   = LCOE_pv_onshore_erares[i,j];

                            elseif (ONSHORE_pv_onshore_erares[i,j]==area*onshore_density) && (PV_pv_onshore_erares[i,j]<0.01*area*pv_density)

                                pv_region_matrix[r,c]            = PV_pv_onshore_erares[i,j];
                                onshore_region_matrix[r,c]       = ONSHORE_pv_onshore_erares[i,j];
                                offshore_region_matrix[r,c]      = OFFSHORE_pv_onshore_erares[i,j];
                                el_region_matrix[r,c]            = EL_pv_onshore_erares[i,j];
                                lcoh_region_matrix[r,c]          = LCOH_pv_onshore_erares[i,j]#min(LCOH_pv_onshore_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_pv_onshore_erares[i,j];
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                                lcoe_nocurt_region_matrix[r,c]   = LCOE_pv_onshore_erares[i,j];

                            else

                                pv_region_matrix[r,c]            = PV_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                onshore_region_matrix[r,c]       = ONSHORE_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                offshore_region_matrix[r,c]      = OFFSHORE_pv_onshore_erares[i,j];
                                el_region_matrix[r,c]            = EL_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                lcoh_region_matrix[r,c]          = LCOH_pv_onshore_erares[i,j]#min(LCOH_pv_onshore_erares[i,j], 10);
                                production_region_matrix[r,c]    = PROD_pv_onshore_erares[i,j] * (1 - crop_nat_excl_coeff);;
                                sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                                sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                                sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                                lcoe_nocurt_region_matrix[r,c]   = LCOE_pv_onshore_erares[i,j];
                            end
                        end
                    end
                    if mask_pv_lcoh[r,c]>0 && !hybrid_system_mask[r,c]#mask_onshore_lcoh[r,c]==0 #solo pv
                        if land[r,c]!=14
                            pv_region_matrix[r,c]            = PV_pv_erares[i,j];
                            onshore_region_matrix[r,c]       = ONSHORE_pv_erares[i,j];
                            offshore_region_matrix[r,c]      = OFFSHORE_pv_erares[i,j];
                            el_region_matrix[r,c]            = EL_pv_erares[i,j];
                            lcoh_region_matrix[r,c]          = LCOH_pv_erares[i,j]#min(LCOH_pv_erares[i,j], 10);
                            production_region_matrix[r,c]    = PROD_pv_erares[i,j];
                            sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                            sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                            sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                            lcoe_nocurt_region_matrix[r,c]   = LCOE_pv_erares[i,j];
                        end
                        if land[r,c]==14
                            pv_region_matrix[r,c]            = PV_pv_erares[i,j]* (1 - crop_nat_excl_coeff);
                            onshore_region_matrix[r,c]       = ONSHORE_pv_erares[i,j];
                            offshore_region_matrix[r,c]      = OFFSHORE_pv_erares[i,j];
                            el_region_matrix[r,c]            = EL_pv_erares[i,j]* (1 - crop_nat_excl_coeff);
                            lcoh_region_matrix[r,c]          = LCOH_pv_erares[i,j]#min(LCOH_pv_erares[i,j], 10);
                            production_region_matrix[r,c]    = PROD_pv_erares[i,j]* (1 - crop_nat_excl_coeff);
                            sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                            sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                            sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                            lcoe_nocurt_region_matrix[r,c]   = LCOE_pv_erares[i,j];
                        end
                    end
                    if  mask_onshore_lcoh[r,c]>0 && !hybrid_system_mask[r,c] #solo onshore mask_pv_lcoh[r,c]==0 &&
                        pv_region_matrix[r,c]            = PV_onshore_erares[i,j];
                        onshore_region_matrix[r,c]       = ONSHORE_onshore_erares[i,j];
                        offshore_region_matrix[r,c]      = OFFSHORE_onshore_erares[i,j];
                        el_region_matrix[r,c]            = EL_onshore_erares[i,j];
                        lcoh_region_matrix[r,c]          = LCOH_onshore_erares[i,j]#min(LCOH_onshore_erares[i,j], 10);
                        production_region_matrix[r,c]    = PROD_onshore_erares[i,j];
                        sumcfpv_region_matrix[r,c]       = SUM_CF_PV[i,j];
                        sumcfonshore_region_matrix[r,c]  = SUM_CF_ONSHORE[i,j];
                        sumcfoffshore_region_matrix[r,c] = SUM_CF_OFFSHORE[i,j];
                        lcoe_nocurt_region_matrix[r,c]   = LCOE_onshore_erares[i,j];

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
                        lcoe_nocurt_region_matrix[r,c]   = LCOE_offshore_erares[i,j]
                    end
                end
            end
        end
        next!(updateprogress)
        # Threads.unlock(l)
    end

    return pv_region_matrix, onshore_region_matrix, offshore_region_matrix, el_region_matrix, lcoh_region_matrix, production_region_matrix,
            PV_pv_onshore_erares, ONSHORE_pv_onshore_erares, OFFSHORE_pv_onshore_erares, EL_pv_onshore_erares, LCOH_pv_onshore_erares, PROD_pv_onshore_erares,
            PV_pv_erares, ONSHORE_pv_erares, OFFSHORE_pv_erares, EL_pv_erares, LCOH_pv_erares, PROD_pv_erares,
            PV_onshore_erares, ONSHORE_onshore_erares, OFFSHORE_onshore_erares, EL_onshore_erares, LCOH_onshore_erares, PROD_onshore_erares,
            PV_offshore_erares, ONSHORE_offshore_erares, OFFSHORE_offshore_erares, EL_offshore_erares, LCOH_offshore_erares, PROD_offshore_erares,
            SUM_CF_PV, SUM_CF_ONSHORE, SUM_CF_OFFSHORE,sumcfpv_region_matrix, sumcfonshore_region_matrix, sumcfoffshore_region_matrix,
            lcoe_nocurt_region_matrix

end

function find_minimum_LCOH(loptions, lonrange, latrange, regionlist, dim_year_lon_lats,
                           CF_pvplant, area, CF_onshore, CF_offshore, type_min,
                           CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, WACC_pv, WACC_onshore, WACC_offshore, eff_el, CAPEX_el)
   
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    @unpack era_year, res, erares, pv_density, onshore_density, offshore_density = loptions
    numreg = length(regionlist)
    # yearlength, nlons, nlats = size(solarGTI)
    yearlength, nlons, nlats = dim_year_lon_lats
    firsttime = DateTime(era_year, 1, 1)

    WACC_el = WACC_pv # Tanto sono tutti uguali, ho gia' controllato

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
        rescaled_annual_hydrogen_production,
        LCOE_min = iterative_find_min_3var(loptions, CF_pvplant, max_capacity_pv,
                                                                        CF_onshore, max_capacity_onshore,
                                                                        CAPEX_pv, CAPEX_wind_onshore, CAPEX_el,
                                                                        WACC_el, WACC_pv, WACC_onshore, eff_el)
    else
        max_capacity_pv       = area * pv_density
        max_capacity_onshore  = area * onshore_density
        max_capacity_offshore = area * offshore_density

        x_min_LCOH_rescaled,
        y_onshore_min_LCOH_rescaled,
        y_offshore_min_LCOH_rescaled,
        z_min_LCOH_rescaled,
        minLCOH,
        rescaled_annual_hydrogen_production,
        LCOE_min = iterative_find_min_2var(loptions, CF_pvplant, max_capacity_pv,
                                                                        CF_onshore, max_capacity_onshore,
                                                                        CF_offshore, max_capacity_offshore,
                                                                        CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, CAPEX_el,
                                                                        WACC_el, WACC_pv, WACC_onshore, WACC_offshore,
                                                                        type_min, eff_el)
    end

    return x_min_LCOH_rescaled, y_onshore_min_LCOH_rescaled, y_offshore_min_LCOH_rescaled, z_min_LCOH_rescaled, minLCOH, rescaled_annual_hydrogen_production, LCOE_min
end

#Funzione nuova
function iterative_find_min_3var(loptions, CF_pvplant, max_capacity_pv,
                                CF_onshore, max_capacity_onshore,
                                CAPEX_pv, CAPEX_wind, CAPEX_el,
                                WACC_el, WACC_pv, WACC_onshore, eff_el)
    
    # println("Minimum di CF_pvplant")
    # println(minimum(CF_pvplant))
    # println("Minimum di CF_onshore")
    # println(minimum(CF_onshore))
    # println("Maximum di CF_pvplant")
    # println(maximum(CF_pvplant))
    # println("Maximum di CF_onshore")
    # println(maximum(CF_onshore))

    lifetime_stack         = 80000;              # Hours
    lifetime_el            = 25;                 # Years
    lifetime_pv            = 25;                 # Years
    lifetime_wind          = 25;                 # Years
    OPEX_el                = 0.03;               # %CAPEX
    OPEX_wind              = 0.03;              # %CAPEX
    global system_lifetime = min(lifetime_el, lifetime_pv, lifetime_wind);
    global generation_system_lifetime = min(lifetime_pv, lifetime_wind);

    OPEX_pv               = 0.01;
    # eff_el                = 0.875;
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
    # global el_over_pv_residual      = 0;
    # global el_over_onshore_residual = 0;

    #global dd = (7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4);
    #global dd = (7,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5);
    global dd = (7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6);
    #global dd = (10,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7);

    for l in 1:length(dd)

        d=dd[l];

        LCOH_min0 = LCOH_min;

        if l!=0
            el_over_pv0      = 0;
            el_over_onshore0 = 0;
        else
            el_over_pv0      = zmin_LCOH/xmin_LCOH;
            el_over_onshore0 = zmin_LCOH/ymin_LCOH;
        end

        # println("\nLCOH / el_over_pv0 / el_over_onshore0 (type_min=1) in=$l : $LCOH_min / $el_over_pv0 / $el_over_onshore0")

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

        annual_invest_LCOE                   = zeros(length(x),length(y));
        total_cost_LCOE                      = zeros(length(x),length(y));
        annual_energy_production             = zeros(length(x),length(y));
        total_energy_production_discounted   = zeros(length(x),length(y));

        annual_invest_LCOH                   = zeros(length(x),length(y),length(z));
        total_cost_LCOH                      = zeros(length(x),length(y),length(z));
        annual_hydrogen_production           = zeros(length(x),length(y),length(z));
        total_hydrogen_production_discounted = zeros(length(x),length(y),length(z));

        LCOE = zeros(length(x),length(y));
        LCOH = zeros(length(x),length(y),length(z));

        sumCF_pv=sum(CF_pvplant);
        sumCF_wind=sum(CF_onshore);

        for i in 1:length(x)
    	    for j in 1:length(y)

    		    annual_invest_LCOE[i,j]       = cost_pv*x[i]  + cost_wind*y[j];
                annual_energy_production[i,j] = sumCF_pv*x[i] + sumCF_wind*y[j];

                WACC_wave_gen  = (WACC_pv*cost_pv*x[i] + WACC_onshore*cost_wind*y[j]) / (cost_pv*x[i] + cost_wind*y[j]);
                D_gen          = sum([1 / (1+WACC_wave_gen)^s for s in 1:generation_system_lifetime]);

                total_cost_LCOE[i,j]                     = D_gen * annual_invest_LCOE[i,j]       * generation_system_lifetime;
                total_energy_production_discounted[i,j]  = D_gen * annual_energy_production[i,j] * generation_system_lifetime;

    		    for k in 1:length(z)

                    # for t in 1:8760
                    #     annual_hydrogen_production[i,j,k]+=min((CF_pvplant[t]*x[i]+CF_onshore[t]*y[j]),z[k])/H2_energy_consumption;
                    # end
                    annual_hydrogen_production[i,j,k] = sum([min((CF_pvplant[t]*x[i]+CF_onshore[t]*y[j]),z[k])/H2_energy_consumption for t in 1:8760]);

                    WACC_wave_H2 = (WACC_pv*cost_pv*x[i] + WACC_onshore*cost_wind*y[j] + WACC_el*cost_el*z[k]) / (cost_pv*x[i] + cost_wind*y[j] + cost_el*z[k]);
                    D_H2 = sum([1 / (1+WACC_wave_H2)^s for s in 1:system_lifetime]);

                    total_cost_LCOH[i,j,k]                      = D_H2 * (annual_invest_LCOE[i,j] + cost_el*z[k]) * system_lifetime;
                    total_hydrogen_production_discounted[i,j,k] = D_H2 * annual_hydrogen_production[i,j,k]        * system_lifetime;
                end
            end
        end

        LCOE = 1000*total_cost_LCOE./total_energy_production_discounted;
        LCOH = total_cost_LCOH./total_hydrogen_production_discounted;

        # min_LCOE, indicesLCOEmin=findmin(LCOE);
        # global LCOE_min = min_LCOE
        # xmin_LCOE = x[indicesLCOEmin[1]];
        # ymin_LCOE = y[indicesLCOEmin[2]];

        minLCOH, indicesLCOHmin=findmin(LCOH);
        global LCOH_min  = minLCOH;
        global xmin_LCOH = x[indicesLCOHmin[1]];
        global ymin_LCOH = y[indicesLCOHmin[2]];
        global zmin_LCOH = z[indicesLCOHmin[3]];

        global LCOE_in_minLCOH = LCOE[indicesLCOHmin[1],indicesLCOHmin[2]]

        # println("xmin_LCOH in l=$l : $xmin_LCOH")
        # println("ymin_LCOH in l=$l : $ymin_LCOH")
        # println("zmin_LCOH in l=$l : $zmin_LCOH")

        delta_pv   = (endvalue_pv   - startvalue_pv)   / (d-1);
        delta_wind = (endvalue_wind - startvalue_wind) / (d-1);
        delta_el   = (endvalue_el   - startvalue_el)   / (d-1);

        global startvalue_pv   = max(xmin_LCOH-delta_pv,1);
        global endvalue_pv     = min(xmin_LCOH+delta_pv,max_cap_pv);
        global startvalue_wind = max(ymin_LCOH-delta_wind,1);
        global endvalue_wind   = min(ymin_LCOH+delta_wind,max_cap_wind);
        global startvalue_el   = zmin_LCOH-delta_el;
        global endvalue_el     = zmin_LCOH+delta_el;

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

        # println("zmin_LCOH/xmin_LCOH = $(zmin_LCOH/xmin_LCOH)")
        # println("zmin_LCOH/ymin_LCOH = $(zmin_LCOH/ymin_LCOH)")
        el_over_pv_residual      = abs(el_over_pv0-zmin_LCOH/xmin_LCOH);
        # println("el_over_pv_residual: $el_over_pv_residual")
        el_over_onshore_residual = abs(el_over_onshore0-zmin_LCOH/ymin_LCOH);
        # println("el_over_onshore_residual: $el_over_onshore_residual")
        LCOH_residual = abs(LCOH_min0-LCOH_min)
        # (Precedente-Attuale)/Precedente
        # el_over_pv_residual_percent      = el_over_pv_residual      / el_over_pv0;
        # el_over_onshore_residual_percent = el_over_onshore_residual / el_over_onshore0;
        el_over_pv_residual_percent      = el_over_pv_residual      / (zmin_LCOH/xmin_LCOH);
        el_over_onshore_residual_percent = el_over_onshore_residual / (zmin_LCOH/ymin_LCOH);
        # println("el_over_pv_residual_percent = $el_over_pv_residual_percent")
        # println("el_over_onshore_residual_percent = $el_over_onshore_residual_percent")
        # println("LCOH_residual / el_over_pv_residual / el_over_onshore_residual (type_min=1) in=$l : $LCOH_residual / $el_over_pv_residual / $el_over_onshore_residual")
        ((el_over_pv_residual_percent < 0.01) && (el_over_onshore_residual_percent < 0.01)) && break
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

    return x_min_LCOH_rescaled, y_onshore_min_LCOH_rescaled, z_min_LCOH_rescaled, LCOH_min, rescaled_annual_hydrogen_production, LCOE_in_minLCOH
end

#Funzione nuova
function iterative_find_min_2var(loptions, CF_pvplant, max_capacity_pv,
                                CF_onshore, max_capacity_onshore,
                                CF_offshore, max_capacity_offshore,
                                CAPEX_pv, CAPEX_wind_onshore, CAPEX_wind_offshore, CAPEX_el,
                                WACC_el, WACC_pv, WACC_onshore, WACC_offshore,
                                type_min, eff_el)

 
    @unpack pv_density, onshore_density  = loptions

    lifetime_stack = 80000;              # Hours
    lifetime_el = 25;                    # Years
    global lifetime_wind_onshore = 25;          # Years
    global lifetime_wind_offshore = 25;         # Years
    global lifetime_pv = 25;                    # Years
    global system_lifetime_2 = min(lifetime_el, lifetime_pv);
    global system_lifetime_3 = min(lifetime_el, lifetime_wind_onshore);
    global system_lifetime_4 = min(lifetime_el, lifetime_wind_offshore);

    OPEX_el = 0.03;                      # %CAPEX
    OPEX_wind_onshore = 0.03;           # %CAPEX
    OPEX_wind_offshore = 0.025;          # %CAPEX
    OPEX_pv = 0.01;
    # eff_el = 0.875;
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
        global zmin_LCOH     = 0;

        global startvalue_pv = 0;
        global endvalue_pv   = max_cap_pv;
        global startvalue_el = 1;
        global endvalue_el   = max_cap_pv;

        if max_cap_pv<=0
            return 0,0,0,0,0,0
        elseif max_cap_pv<=10
            return 0,0,0,0,0,0
        end

        global dd = (7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6);
        for l in 1:length(dd)
            d=dd[l];

            if l==1
                el_over_pv0 = 0;
            else
                el_over_pv0 = zmin_LCOH/xmin_LCOH;
            end

            LCOH_min0   = LCOH_min

            # println("\nLCOH / el_over_pv0 (type_min=2) iter n.$l : $LCOH_min / $el_over_pv0")

            x=collect(range(startvalue_pv, length=d, stop=endvalue_pv));
            z=collect(range(startvalue_el, length=d, stop=endvalue_el));

            annual_invest_LCOE                   = zeros(length(x));
            total_cost_LCOE                      = zeros(length(x));
            annual_energy_production             = zeros(length(x));
            total_energy_production_discounted   = zeros(length(x));

            annual_invest_LCOH                   = zeros(length(x),length(z));
            total_cost_LCOH                      = zeros(length(x),length(z));
            annual_hydrogen_production           = zeros(length(x),length(z));
            total_hydrogen_production_discounted = zeros(length(x),length(z));

            LCOE = zeros(length(x));
            LCOH = zeros(length(x),length(z));

            sumCF_pv=sum(CF_pvplant);

            for i in 1:length(x)
    	        annual_invest_LCOE[i]       = cost_pv  * x[i];
    	        annual_energy_production[i] = sumCF_pv * x[i];

                D_gen = sum([1 / (1+WACC_pv)^s for s in 1:lifetime_pv]);

                total_cost_LCOE[i]                     = D_gen * annual_invest_LCOE[i]       * lifetime_pv;
                total_energy_production_discounted[i]  = D_gen * annual_energy_production[i] * lifetime_pv;

    		        for k in 1:length(z)

    	                annual_invest_LCOH[i,k] = annual_invest_LCOE[i] + cost_el * z[k];

                        # for t in 1:8760
                        #     annual_hydrogen_production[i,k] += min((CF_pvplant[t] * x[i]), z[k]) / H2_energy_consumption;
                        # end
                        annual_hydrogen_production[i,k] = sum([min(CF_pvplant[t]*x[i],z[k])/H2_energy_consumption for t in 1:8760]);

                        WACC_wave = (WACC_pv * cost_pv * x[i] + WACC_el * cost_el * z[k]) / (cost_pv * x[i] + cost_el * z[k]);
                        D = sum([1 / (1+WACC_wave)^s for s in 1:system_lifetime_2]);

                        total_cost_LCOH[i,k] = D * (annual_invest_LCOE[i] + cost_el*z[k]) * system_lifetime_2;
                        total_hydrogen_production_discounted[i,k] =  D * annual_hydrogen_production[i,k] * system_lifetime_2;

                    end
            end

            LCOE=1000*annual_invest_LCOE./annual_energy_production;
            LCOH=total_cost_LCOH./total_hydrogen_production_discounted;

            # global min_LCOE, indicesLCOEmin = findmin(LCOE);
            # global LCOE_min  = min_LCOE
            # global xmin_LCOE = x[indicesLCOEmin[1]];

            global minLCOH, indicesLCOHmin=findmin(LCOH);
            global LCOH_min= minLCOH;
            global xmin_LCOH=x[indicesLCOHmin[1]];
            global zmin_LCOH=z[indicesLCOHmin[2]];

            global LCOE_in_minLCOH = LCOE[indicesLCOHmin[1]];

            # println("xmin_LCOH in l=$l : $xmin_LCOH")
            # println("zmin_LCOH in l=$l : $zmin_LCOH")

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
            # if l==length(dd)
                global annual_hydrogen_production_final = annual_hydrogen_production
            # end

            # println("zmin_LCOH/xmin_LCOH = $(zmin_LCOH/xmin_LCOH)")

            el_over_pv_residual = abs(el_over_pv0-zmin_LCOH/xmin_LCOH);
            # println("el_over_pv_residual: $el_over_pv_residual")

            LCOH_residual = abs(LCOH_min0-LCOH_min)
            # (Precedente-Attuale)/Precedente
            el_over_pv_residual_percent   = el_over_pv_residual / el_over_pv0;
            # el_over_pv_residual_percent = el_over_pv_residual / (zmin_LCOH/xmin_LCOH);

            # println("el_over_pv_residual_percent = $el_over_pv_residual_percent")
            # println("LCOH_residual / el_over_pv_residual  (type_min=2) in=$l : $LCOH_residual / $el_over_pv_residual ")

            (el_over_pv_residual_percent < 0.01) && break

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
        global LCOH_min             = 0;
        global yonshoremin_LCOH     = 0;
        global max_cap_wind_onshore = max_capacity_onshore;
        global zmin_LCOH            = 0;

        global startvalue_wind_onshore = 0;
        global endvalue_wind_onshore   = max_cap_wind_onshore;
        global startvalue_el           = 1;
        global endvalue_el             = max_cap_wind_onshore;

        if max_cap_wind_onshore<=0
            return 0,0,0,0,0,0
        elseif max_cap_wind_onshore<=10
            return 0,0,0,0,0,0
        end

        global dd = (7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6);
        for l in 1:length(dd)
            d=dd[l];

            if l==1
                el_over_onshore0 = 0;
            else
                el_over_onshore0 = zmin_LCOH/yonshoremin_LCOH;
            end
            LCOH_min0=LCOH_min

            # println("\nLCOH / el_over_onshore0 (type_min=3) iter n.$l : $LCOH_min / $el_over_onshore0")

            y_onshore = collect(range(startvalue_wind_onshore, length=d, stop=endvalue_wind_onshore));
            z         = collect(range(startvalue_el, length=d, stop=endvalue_el));

            annual_invest_LCOE                   =zeros(length(y_onshore));
            total_cost_LCOE                      =zeros(length(y_onshore));
            annual_energy_production             =zeros(length(y_onshore));
            total_energy_production_discounted   =zeros(length(y_onshore));

            annual_invest_LCOH                   =zeros(length(y_onshore),length(z));
            total_cost_LCOH                      =zeros(length(y_onshore),length(z));
            annual_hydrogen_production           =zeros(length(y_onshore),length(z));
            total_hydrogen_production_discounted =zeros(length(y_onshore),length(z));

            LCOE = zeros(length(y_onshore));
            LCOH = zeros(length(y_onshore),length(z));

            sumCF_wind_onshore = sum(CF_onshore);

            for i in 1:length(y_onshore)
                annual_invest_LCOE[i]       = cost_wind_onshore  * y_onshore[i];
                annual_energy_production[i] = sumCF_wind_onshore * y_onshore[i];

                D_gen = sum([1 / (1+WACC_onshore)^s for s in 1:lifetime_wind_onshore]);

                total_cost_LCOE[i]                     = D_gen * annual_invest_LCOE[i]       * lifetime_wind_onshore;
                total_energy_production_discounted[i]  = D_gen * annual_energy_production[i] * lifetime_wind_onshore;

                    for k in 1:length(z)
                        annual_invest_LCOH[i,k] = annual_invest_LCOE[i] + cost_el * z[k];

                        # for t in 1:8760
                        #     annual_hydrogen_production[i,k] += min((CF_onshore[t] * y_onshore[i]), z[k]) / H2_energy_consumption;
                        # end
                        annual_hydrogen_production[i,k] = sum([min(CF_onshore[t]*y_onshore[i],z[k])/H2_energy_consumption for t in 1:8760]);

                        WACC_wave = (WACC_onshore * cost_wind_onshore * y_onshore[i] + WACC_el * cost_el * z[k]) / (cost_wind_onshore * y_onshore[i] + cost_el * z[k]);
                        D = sum([1 / (1+WACC_wave)^s for s in 1:system_lifetime_3]);

                        total_cost_LCOH[i,k] = D * (annual_invest_LCOE[i] + cost_el*z[k]) * system_lifetime_3;
                        total_hydrogen_production_discounted[i,k] =  D * annual_hydrogen_production[i,k] * system_lifetime_3;

                    end
            end

            LCOE =1000*total_cost_LCOE./total_energy_production_discounted;
            LCOH =total_cost_LCOH./total_hydrogen_production_discounted;

            # global min_LCOE, indicesLCOEmin = findmin(LCOE);
            # global LCOE_min = min_LCOE;
            # global yonshoremin_LCOE = y_onshore[indicesLCOEmin[1]];

            global minLCOH, indicesLCOHmin = findmin(LCOH);
            global LCOH_min = minLCOH;
            global yonshoremin_LCOH = y_onshore[indicesLCOHmin[1]];
            global zmin_LCOH = z[indicesLCOHmin[2]];

            global LCOE_in_minLCOH = LCOE[indicesLCOHmin[1]];


            # println("yonshoremin_LCOH in l=$l : $yonshoremin_LCOH")
            # println("zmin_LCOH in l=$l : $zmin_LCOH")

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
            # if l==length(dd)
                global annual_hydrogen_production_final = annual_hydrogen_production
            # end
            # println("zmin_LCOH/yonshoremin_LCOH = $(zmin_LCOH/yonshoremin_LCOH)")

            el_over_onshore_residual = abs(el_over_onshore0-zmin_LCOH/yonshoremin_LCOH);
            # println("el_over_onshore_residual: $el_over_onshore_residual")

            LCOH_residual = abs(LCOH_min0-LCOH_min)
            # (Precedente-Attuale)/Precedente
            el_over_onshore_residual_percent   = el_over_onshore_residual / el_over_onshore0;
            # el_over_onshore_residual_percent = el_over_onshore_residual / (zmin_LCOH/yonshoremin_LCOH);

            # println("el_over_onshore_residual_percent = $el_over_onshore_residual_percent")
            # println("LCOH_residual / el_over_onshore_residual  (type_min=3) in=$l : $LCOH_residual / $el_over_onshore_residual ")

            (el_over_onshore_residual_percent < 0.01) && break
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

        global dd = (7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6);
        for l in 1:length(dd)
            d=dd[l];

            if l==1
                el_over_offshore0 = 0;
            else
                el_over_offshore0 = zmin_LCOH/yoffshoremin_LCOH;
            end

            LCOH_min0=LCOH_min

            # println("\nLCOH / el_over_offshore0 (type_min=4) iter n.$l : $LCOH_min / $el_over_offshore0")

            y_offshore = collect(range(startvalue_wind_offshore, length=d, stop=endvalue_wind_offshore));
            z          = collect(range(startvalue_el, length=d, stop=endvalue_el));

            annual_invest_LCOE                   = zeros(length(y_offshore));
            total_cost_LCOE                      = zeros(length(y_offshore));
            annual_energy_production             = zeros(length(y_offshore));
            total_energy_production_discounted   = zeros(length(y_offshore));

            annual_invest_LCOH                   = zeros(length(y_offshore),length(z));
            total_cost_LCOH                      = zeros(length(y_offshore),length(z));
            annual_hydrogen_production           = zeros(length(y_offshore),length(z));
            total_hydrogen_production_discounted = zeros(length(y_offshore),length(z));

            LCOE = zeros(length(y_offshore));
            LCOH = zeros(length(y_offshore),length(z));

            sumCF_wind_offshore = sum(CF_offshore);

            for i in 1:length(y_offshore)
                annual_invest_LCOE[i]       = cost_wind_offshore  * y_offshore[i];
                annual_energy_production[i] = sumCF_wind_offshore * y_offshore[i];

                D_gen = sum([1 / (1+WACC_offshore)^s for s in 1:lifetime_wind_offshore]);

                total_cost_LCOE[i]                     = D_gen * annual_invest_LCOE[i]       * lifetime_wind_offshore;
                total_energy_production_discounted[i]  = D_gen * annual_energy_production[i] * lifetime_wind_offshore;

                for k in 1:length(z)
                    annual_invest_LCOH[i,k] = annual_invest_LCOE[i] + cost_el * z[k];
                    # for t in 1:8760
                    #     annual_hydrogen_production[i,k] += min((CF_offshore[t] * y_offshore[i]), z[k]) / H2_energy_consumption;
                    # end
                    annual_hydrogen_production[i,k] = sum([min(CF_offshore[t]*y_offshore[i],z[k])/H2_energy_consumption for t in 1:8760]);

                    WACC_wave = (WACC_offshore * cost_wind_offshore * y_offshore[i] + WACC_el * cost_el * z[k]) / (cost_wind_offshore * y_offshore[i] + cost_el * z[k]);
                    D = sum([1 / (1+WACC_wave)^s for s in 1:system_lifetime_4]);

                    total_cost_LCOH[i,k] = D * (annual_invest_LCOE[i] + cost_el * z[k]) * system_lifetime_4;
                    total_hydrogen_production_discounted[i,k] =  D * annual_hydrogen_production[i,k] * system_lifetime_4;

                end
            end

            LCOE=1000*total_cost_LCOE./total_energy_production_discounted;
            LCOH=total_cost_LCOH./total_hydrogen_production_discounted;

            # global min_LCOE, indicesLCOEmin = findmin(LCOE);
            # global LCOE_min = min_LCOE
            # global yoffshoremin_LCOE = y_offshore[indicesLCOEmin[1]];

            global minLCOH, indicesLCOHmin = findmin(LCOH);
            global LCOH_min = minLCOH;
            global yoffshoremin_LCOH = y_offshore[indicesLCOHmin[1]];
            global zmin_LCOH =z[indicesLCOHmin[2]];

            global LCOE_in_minLCOH = LCOE[indicesLCOHmin[1]];

            # println("yoffshoremin_LCOH in l=$l : $yoffshoremin_LCOH")
            # println("zmin_LCOH in l=$l : $zmin_LCOH")

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
            # if l==length(dd)
                global annual_hydrogen_production_final = annual_hydrogen_production
            # end
            # println("zmin_LCOH/yoffshoremin_LCOH = $(zmin_LCOH/yoffshoremin_LCOH)")

            el_over_offshore_residual = abs(el_over_offshore0-zmin_LCOH/yoffshoremin_LCOH);
            # println("el_over_offshore_residual: $el_over_offshore_residual")

            LCOH_residual = abs(LCOH_min0-LCOH_min)
            # (Precedente-Attuale)/Precedente
            el_over_offshore_residual_percent   = el_over_offshore_residual / el_over_offshore0;
            # el_over_offshore_residual_percent = el_over_offshore_residual / (zmin_LCOH/yoffshoremin_LCOH);

            # println("el_over_offshore_residual_percent = $el_over_offshore_residual_percent")
            # println("LCOH_residual / el_over_offshore_residual  (type_min=4) in=$l : $LCOH_residual / $el_over_offshore_residual ")

            (el_over_offshore_residual_percent < 0.01) && break
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

    return x_min_LCOH_rescaled, y_onshore_min_LCOH_rescaled, y_offshore_min_LCOH_rescaled, z_min_LCOH_rescaled, LCOH_min, rescaled_annual_hydrogen_production, LCOE_in_minLCOH
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


#Funzione Nuova
function calculate_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, land, crop_nat_excl_coeff)
   
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore = loptions

    println("PV density: $pv_density")
    println("Modified Onshore density: $onshore_density")
    println("Modified Offshore density: $offshore_density")

    hybrid_system_mask = mask_pv_lcoh .& mask_onshore_lcoh;

    global area_idonea_pv = 0;
    global area_idonea_onshore = 0;
    global area_idonea_offshore = 0;
    global area_idonea_pv_onshore = 0;
    global area_idonea_pv_onshore_PV = 0;
    global area_idonea_pv_onshore_ONSHORE = 0;
    global area_tot = 0;
    global area_esclusa = 0;

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue

                    area = cellarea[c]

                    if mask_pv_lcoh[r,c]>0 && !hybrid_system_mask[r,c] && land[r,c]!=14
                        area_idonea_pv += area * plant_area
                    end
                    if mask_pv_lcoh[r,c]>0 && !hybrid_system_mask[r,c] && land[r,c]==14
                        area_idonea_pv += area * plant_area * (1-crop_nat_excl_coeff)
                    end

                    if mask_onshore_lcoh[r,c]>0 && !hybrid_system_mask[r,c]
                        area_idonea_onshore += area * area_onshore
                    end
                    if mask_offshore_lcoh[r,c]>0
                        area_idonea_offshore += area * area_offshore
                    end

                    if hybrid_system_mask[r,c]>0 && land[r,c]!=14
                        area_idonea_pv_onshore_PV += area
                        area_idonea_pv_onshore_ONSHORE += area
                    end

                    if hybrid_system_mask[r,c]>0 && land[r,c]==14
                        area_idonea_pv_onshore_PV += area * (1-crop_nat_excl_coeff)
                        area_idonea_pv_onshore_ONSHORE += area
                    end

                    if (mask_pv_lcoh[r,c]==0 && mask_onshore_lcoh[r,c]==0) && regions[r,c]!=0 #&& regions[r,c]!=NOREGION
                        area_esclusa += area
                    end

                    if regions[r,c]!=0 && regions[r,c]!=NOREGION
                        area_tot += area
                    end

                end
            end
        end
    end



    installable_pv             = (area_idonea_pv                 * pv_density       * plant_area)    / 1000000;
    installable_onshore        = (area_idonea_onshore            * onshore_density  * area_onshore)  / 1000000;
    installable_offshore       = (area_idonea_offshore           * offshore_density * area_offshore) / 1000000;
    installable_hybrid_pv      = (area_idonea_pv_onshore_PV      * pv_density       * plant_area)    / 1000000;
    installable_hybrid_onshore = (area_idonea_pv_onshore_ONSHORE * onshore_density  * area_onshore)  / 1000000;


    return area_idonea_pv, area_idonea_onshore, area_idonea_offshore, area_idonea_pv_onshore_PV, area_idonea_pv_onshore_ONSHORE,
            installable_pv, installable_onshore, installable_offshore, installable_hybrid_pv, installable_hybrid_onshore,
            area_tot, area_esclusa
end

function calculate_mean_CF(loptions, regions, lonrange, latrange, solarGTI, windspeed, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, pv_system_efficiency, land, crop_nat_excl_coeff)
    
    @unpack era_year, res, erares,pv_density, plant_area,
    onshore_density, area_onshore, offshore_density, area_offshore  = loptions
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)
    global yearlength, nlons, nlats = size(solarGTI)

    SUM_CF_PV       = reshape(sum(solarGTI, dims=1)/8760, (nlons, nlats));
    SUM_CF_ONSHORE  = reshape(sum(speed2capacityfactor_onshore.(windspeed), dims=1)/8760, (nlons, nlats));
    SUM_CF_OFFSHORE = reshape(sum(speed2capacityfactor_offshore.(windspeed), dims=1)/8760, (nlons, nlats));

    global pv_dist_1 = zeros(1:8760)
    global pv_cap_1  = 0
    global pv_dist_2 = zeros(1:8760)
    global pv_cap_2  = 0
    global pv_dist_3 = zeros(1:8760)
    global pv_cap_3  = 0
    global pv_dist_4 = zeros(1:8760)
    global pv_cap_4  = 0
    global pv_dist_5 = zeros(1:8760)
    global pv_cap_5  = 0

    global onshore_dist_1 = zeros(1:8760)
    global onshore_cap_1  = 0
    global onshore_dist_2 = zeros(1:8760)
    global onshore_cap_2  = 0
    global onshore_dist_3 = zeros(1:8760)
    global onshore_cap_3  = 0
    global onshore_dist_4 = zeros(1:8760)
    global onshore_cap_4  = 0
    global onshore_dist_5 = zeros(1:8760)
    global onshore_cap_5  = 0

    global offshore_dist_1 = zeros(1:8760)
    global offshore_cap_1  = 0
    global offshore_dist_2 = zeros(1:8760)
    global offshore_cap_2  = 0
    global offshore_dist_3 = zeros(1:8760)
    global offshore_cap_3  = 0
    global offshore_dist_4 = zeros(1:8760)
    global offshore_cap_4  = 0
    global offshore_dist_5 = zeros(1:8760)
    global offshore_cap_5  = 0

    global cf_pv_matrix       = zeros(size(regions));
    global cf_onshore_matrix  = zeros(size(regions));
    global cf_offshore_matrix = zeros(size(regions));

    global counter_pv_dist_1=0;
    global counter_pv_dist_2=0;
    global counter_pv_dist_3=0;
    global counter_pv_dist_4=0;
    global counter_pv_dist_5=0;
    global counter_onshore_dist_1=0;
    global counter_onshore_dist_2=0;
    global counter_onshore_dist_3=0;
    global counter_onshore_dist_4=0;
    global counter_onshore_dist_5=0;
    global counter_offshore_dist_1=0;
    global counter_offshore_dist_2=0;
    global counter_offshore_dist_3=0;
    global counter_offshore_dist_4=0;
    global counter_offshore_dist_5=0;

    global contatore = 0;
    num_tot_pixel = nlats*nlons;

    #@floop ThreadedEx()
    #Threads.@threads
    # l = Threads.ReentrantLock()
    # l = Threads.SpinLock()

    for j in 1:nlats
        # Threads.lock(l)
        # #println((thread=Threads.threadid(), iteration=j))
        # Threads.unlock(l)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            contatore = contatore + 1;
            # if contatore%100 == 0
                println("Pixel mean: $contatore/$num_tot_pixel")
            # end
            eralon = eralons[i]

            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for s in 1:size(solarGTI)[1]
                if solarGTI[s,i,j]<0
                    solarGTI[s,i,j] = 0
                end
            end

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c] #region all'indice r, c
                    #(reg == 0 || reg == NOREGION) && continue
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if mask_pv_lcoh[r,c]>0

                        cf_pv_matrix[r,c] = SUM_CF_PV[i,j];

                        if cf_pv_matrix[r,c]>0.26*pv_system_efficiency
                            pv_dist_1 .+= solarGTI[:,i,j];
                            counter_pv_dist_1 +=  1;
                            if land[r,c]==14
                                pv_cap_1 += cellarea[c] * pv_density * plant_area * (1-crop_nat_excl_coeff);
                            else
                                pv_cap_1 += cellarea[c] * pv_density * plant_area;
                            end
                        end

                        if cf_pv_matrix[r,c]>0.22*pv_system_efficiency && cf_pv_matrix[r,c]<=0.26*pv_system_efficiency
                            pv_dist_2 .+= solarGTI[:,i,j];
                            counter_pv_dist_2 +=  1;
                            if land[r,c]==14
                                pv_cap_2 += cellarea[c] * pv_density * plant_area * (1-crop_nat_excl_coeff);
                            else
                                pv_cap_2 += cellarea[c] * pv_density * plant_area
                            end
                        end

                        if cf_pv_matrix[r,c]>0.18*pv_system_efficiency && cf_pv_matrix[r,c]<=0.22*pv_system_efficiency
                            pv_dist_3 .+= solarGTI[:,i,j];
                            counter_pv_dist_3 +=  1;
                            if land[r,c]==14
                                pv_cap_3 += cellarea[c] * pv_density * plant_area * (1-crop_nat_excl_coeff);
                            else
                                pv_cap_3 += cellarea[c] * pv_density * plant_area
                            end
                        end

                        if cf_pv_matrix[r,c]>0.14*pv_system_efficiency && cf_pv_matrix[r,c]<=0.18*pv_system_efficiency
                            pv_dist_4 .+= solarGTI[:,i,j];
                            counter_pv_dist_4 +=  1;
                            if land[r,c]==14
                                pv_cap_4 += cellarea[c] * pv_density * plant_area * (1-crop_nat_excl_coeff);
                            else
                                pv_cap_4 += cellarea[c] * pv_density * plant_area
                            end
                        end

                        if cf_pv_matrix[r,c]>0*pv_system_efficiency && cf_pv_matrix[r,c]<=0.14*pv_system_efficiency
                            pv_dist_5 .+= solarGTI[:,i,j];
                            counter_pv_dist_5 +=  1;
                            if land[r,c]==14
                                pv_cap_5 += cellarea[c] * pv_density * plant_area * (1-crop_nat_excl_coeff);
                            else
                                pv_cap_5 += cellarea[c] * pv_density * plant_area
                            end
                        end
                    end

                    if  mask_onshore_lcoh[r,c]>0

                        cf_onshore_matrix[r,c]  = SUM_CF_ONSHORE[i,j];

                        if cf_onshore_matrix[r,c]>0.60

                            onshore_dist_1 .+= speed2capacityfactor_onshore.(windspeed[:,i,j]);
                            counter_onshore_dist_1 +=  1;
                            onshore_cap_1 += cellarea[c] *  onshore_density * area_onshore;

                        end

                        if cf_onshore_matrix[r,c]>0.45 && cf_onshore_matrix[r,c]<=0.60

                            onshore_dist_2 .+= speed2capacityfactor_onshore.(windspeed[:,i,j]);
                            counter_onshore_dist_2 +=  1;
                            onshore_cap_2 += cellarea[c] *  onshore_density * area_onshore;

                        end

                        if cf_onshore_matrix[r,c]>0.30 && cf_onshore_matrix[r,c]<=0.45

                            onshore_dist_3 .+=  speed2capacityfactor_onshore.(windspeed[:,i,j]);
                            counter_onshore_dist_3 +=  1;
                            onshore_cap_3 += cellarea[c] *  onshore_density * area_onshore;

                        end

                        if cf_onshore_matrix[r,c]>0.15 && cf_onshore_matrix[r,c]<0.30

                            onshore_dist_4 .+=  speed2capacityfactor_onshore.(windspeed[:,i,j]);
                            counter_onshore_dist_4 +=  1;
                            onshore_cap_4 += cellarea[c] *  onshore_density * area_onshore;

                        end

                        if cf_onshore_matrix[r,c]>0 && cf_onshore_matrix[r,c]<=0.15

                            onshore_dist_5 .+=  speed2capacityfactor_onshore.(windspeed[:,i,j]);
                            counter_onshore_dist_5 +=  1;
                            onshore_cap_5 += cellarea[c] *  onshore_density * area_onshore;

                        end
                    end

                    if mask_offshore_lcoh[r,c]>0

                        cf_offshore_matrix[r,c] = SUM_CF_OFFSHORE[i,j];

                        if cf_offshore_matrix[r,c]>0.6

                            offshore_dist_1 .+= speed2capacityfactor_offshore.(windspeed[:,i,j]);
                            counter_offshore_dist_1 +=  1;
                            offshore_cap_1 += cellarea[c] *  offshore_density * area_offshore;

                        end

                        if cf_offshore_matrix[r,c]>0.45 && cf_offshore_matrix[r,c]<=0.60

                            offshore_dist_2 .+= speed2capacityfactor_offshore.(windspeed[:,i,j]);
                            counter_offshore_dist_2 +=  1;
                            offshore_cap_2 += cellarea[c] *  offshore_density * area_offshore;

                        end

                        if cf_offshore_matrix[r,c]>0.30 && cf_offshore_matrix[r,c]<=0.45

                            offshore_dist_3 .+= speed2capacityfactor_offshore.(windspeed[:,i,j]);
                            counter_offshore_dist_3 +=  1;
                            offshore_cap_3 += cellarea[c] *  offshore_density * area_offshore;

                        end

                        if cf_offshore_matrix[r,c]>0.15 && cf_offshore_matrix[r,c]<0.30

                            offshore_dist_4 .+= speed2capacityfactor_offshore.(windspeed[:,i,j]);
                            counter_offshore_dist_4 +=  1;
                            offshore_cap_4 += cellarea[c] *  offshore_density * area_offshore;

                        end

                        if cf_offshore_matrix[r,c]>0 && cf_offshore_matrix[r,c]<=0.15

                            offshore_dist_5 .+= speed2capacityfactor_offshore.(windspeed[:,i,j]);
                            counter_offshore_dist_5 +=  1;
                            offshore_cap_5 += cellarea[c] *  offshore_density * area_offshore;

                        end
                    end
                end
            end
        end
    end

    global mean_pv_dist_1       = zeros(8760);
    global mean_pv_dist_2       = zeros(8760);
    global mean_pv_dist_3       = zeros(8760);
    global mean_pv_dist_4       = zeros(8760);
    global mean_pv_dist_5       = zeros(8760);
    global mean_onshore_dist_1  = zeros(8760);
    global mean_onshore_dist_2  = zeros(8760);
    global mean_onshore_dist_3  = zeros(8760);
    global mean_onshore_dist_4  = zeros(8760);
    global mean_onshore_dist_5  = zeros(8760);
    global mean_offshore_dist_1 = zeros(8760);
    global mean_offshore_dist_2 = zeros(8760);
    global mean_offshore_dist_3 = zeros(8760);
    global mean_offshore_dist_4 = zeros(8760);
    global mean_offshore_dist_5 = zeros(8760);

    if counter_pv_dist_1>0
        mean_pv_dist_1 = pv_dist_1 ./ counter_pv_dist_1;
    end
    if counter_pv_dist_2>0
        mean_pv_dist_2 = pv_dist_2 ./ counter_pv_dist_2;
    end
    if counter_pv_dist_3>0
        mean_pv_dist_3 = pv_dist_3 ./ counter_pv_dist_3;
    end
    if counter_pv_dist_4>0
        mean_pv_dist_4 = pv_dist_4 ./ counter_pv_dist_4;
    end
    if counter_pv_dist_5>0
        mean_pv_dist_5 = pv_dist_5 ./ counter_pv_dist_5;
    end

    if counter_onshore_dist_1>0
        mean_onshore_dist_1 = onshore_dist_1 ./ counter_onshore_dist_1;
    end
    if counter_onshore_dist_2>0
        mean_onshore_dist_2 = onshore_dist_2 ./ counter_onshore_dist_2;
    end
    if counter_onshore_dist_3>0
        mean_onshore_dist_3 = onshore_dist_3 ./ counter_onshore_dist_3;
    end
    if counter_onshore_dist_4>0
        mean_onshore_dist_4 = onshore_dist_4 ./ counter_onshore_dist_4;
    end
    if counter_onshore_dist_5>0
        mean_onshore_dist_5 = onshore_dist_5 ./ counter_onshore_dist_5;
    end

    if counter_offshore_dist_1>0
        mean_offshore_dist_1 = offshore_dist_1 ./ counter_offshore_dist_1;
    end
    if counter_offshore_dist_2>0
        mean_offshore_dist_2 = offshore_dist_2 ./ counter_offshore_dist_2;
    end
    if counter_offshore_dist_3>0
        mean_offshore_dist_3 = offshore_dist_3 ./ counter_offshore_dist_3;
    end
    if counter_offshore_dist_4>0
        mean_offshore_dist_4 = offshore_dist_4 ./ counter_offshore_dist_4;
    end
    if counter_offshore_dist_5>0
        mean_offshore_dist_5 = offshore_dist_5 ./ counter_offshore_dist_5;
    end

    sum_pv_dist_1 = sum(mean_pv_dist_1);
    sum_pv_dist_2 = sum(mean_pv_dist_2);
    sum_pv_dist_3 = sum(mean_pv_dist_3);
    sum_pv_dist_4 = sum(mean_pv_dist_4);
    sum_pv_dist_5 = sum(mean_pv_dist_5);

    sum_onshore_dist_1 = sum(mean_onshore_dist_1);
    sum_onshore_dist_2 = sum(mean_onshore_dist_2);
    sum_onshore_dist_3 = sum(mean_onshore_dist_3);
    sum_onshore_dist_4 = sum(mean_onshore_dist_4);
    sum_onshore_dist_5 = sum(mean_onshore_dist_5);

    sum_offshore_dist_1 = sum(mean_offshore_dist_1);
    sum_offshore_dist_2 = sum(mean_offshore_dist_2);
    sum_offshore_dist_3 = sum(mean_offshore_dist_3);
    sum_offshore_dist_4 = sum(mean_offshore_dist_4);
    sum_offshore_dist_5 = sum(mean_offshore_dist_5);

    avg_pv_dist_1 = sum(mean_pv_dist_1)/8760;
    avg_pv_dist_2 = sum(mean_pv_dist_2)/8760;
    avg_pv_dist_3 = sum(mean_pv_dist_3)/8760;
    avg_pv_dist_4 = sum(mean_pv_dist_4)/8760;
    avg_pv_dist_5 = sum(mean_pv_dist_5)/8760;

    avg_onshore_dist_1 = sum(mean_onshore_dist_1)/8760;
    avg_onshore_dist_2 = sum(mean_onshore_dist_2)/8760;
    avg_onshore_dist_3 = sum(mean_onshore_dist_3)/8760;
    avg_onshore_dist_4 = sum(mean_onshore_dist_4)/8760;
    avg_onshore_dist_5 = sum(mean_onshore_dist_5)/8760;

    avg_offshore_dist_1 = sum(mean_offshore_dist_1)/8760;
    avg_offshore_dist_2 = sum(mean_offshore_dist_2)/8760;
    avg_offshore_dist_3 = sum(mean_offshore_dist_3)/8760;
    avg_offshore_dist_4 = sum(mean_offshore_dist_4)/8760;
    avg_offshore_dist_5 = sum(mean_offshore_dist_5)/8760;



    return  mean_pv_dist_1, mean_pv_dist_2, mean_pv_dist_3, mean_pv_dist_4, mean_pv_dist_5,
            mean_onshore_dist_1, mean_onshore_dist_2, mean_onshore_dist_3, mean_onshore_dist_4, mean_onshore_dist_5,
            mean_offshore_dist_1, mean_offshore_dist_2, mean_offshore_dist_3, mean_offshore_dist_4, mean_offshore_dist_5,
            pv_cap_1, pv_cap_2, pv_cap_3, pv_cap_4, pv_cap_5,
            onshore_cap_1, onshore_cap_2, onshore_cap_3, onshore_cap_4, onshore_cap_5,
            offshore_cap_1, offshore_cap_2, offshore_cap_3, offshore_cap_4, offshore_cap_5,
            sum_pv_dist_1, sum_pv_dist_2, sum_pv_dist_3, sum_pv_dist_4, sum_pv_dist_5,
            sum_onshore_dist_1, sum_onshore_dist_2, sum_onshore_dist_3, sum_onshore_dist_4, sum_onshore_dist_5,
            sum_offshore_dist_1, sum_offshore_dist_2, sum_offshore_dist_3, sum_offshore_dist_4, sum_offshore_dist_5,
            avg_pv_dist_1, avg_pv_dist_2, avg_pv_dist_3, avg_pv_dist_4, avg_pv_dist_5,
            avg_onshore_dist_1, avg_onshore_dist_2, avg_onshore_dist_3, avg_onshore_dist_4, avg_onshore_dist_5,
            avg_offshore_dist_1, avg_offshore_dist_2, avg_offshore_dist_3, avg_offshore_dist_4, avg_offshore_dist_5

end


function find_min_max_cf(loptions, lonrange, latrange, solarGTI, windspeed)
    
    # @unpack era_year, res, erares = loptions
    # eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)
    global yearlength, nlons, nlats = size(solarGTI)

    SUM_CF_PV       = reshape(sum(solarGTI, dims=1)/8760, (nlons, nlats));
    SUM_CF_ONSHORE  = reshape(sum(speed2capacityfactor_onshore.(windspeed), dims=1)/8760, (nlons, nlats));
    SUM_CF_OFFSHORE = reshape(sum(speed2capacityfactor_offshore.(windspeed), dims=1)/8760, (nlons, nlats));

    max_sum_cf_pv       = maximum(SUM_CF_PV);
    min_sum_cf_pv       = minimum(SUM_CF_PV);
    max_sum_cf_onshore  = maximum(SUM_CF_ONSHORE);
    min_sum_cf_onshore  = minimum(SUM_CF_ONSHORE);
    max_sum_cf_offshore = maximum(SUM_CF_OFFSHORE);
    min_sum_cf_offshore = minimum(SUM_CF_OFFSHORE);

    return max_sum_cf_pv, min_sum_cf_pv,
            max_sum_cf_onshore, min_sum_cf_onshore,
             max_sum_cf_offshore, min_sum_cf_offshore
end

function calculate_landtype_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, land)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore = loptions

    global Evergreen_Needleleaf_Forests = 0;
    global Evergreen_Broadleaf_Forests = 0;
    global Deciduous_Needleleaf_Forests = 0;
    global Deciduous_Broadleaf_Forests = 0;
    global Mixed_Forests = 0;
    global Closed_Shrublands = 0;
    global Open_Shrublands = 0;
    global Woody_Savannas = 0;
    global Savannas = 0;
    global Grasslands = 0;
    global Permanent_Wetlands = 0;
    global Croplands = 0;
    global Urban = 0;
    global Cropland_Natural = 0;
    global Snow_Ice = 0;
    global Barren = 0;

    goodland = (regions .> 0)

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if land[r,c]==1 && regions[r,c]!=0
                        Evergreen_Needleleaf_Forests += area
                    elseif land[r,c]==2 && regions[r,c]!=0
                        Evergreen_Broadleaf_Forests += area
                    elseif land[r,c]==3 && regions[r,c]!=0
                        Deciduous_Needleleaf_Forests += area
                    elseif land[r,c]==4 && regions[r,c]!=0
                        Deciduous_Broadleaf_Forests += area
                    elseif land[r,c]==5 && regions[r,c]!=0
                        Mixed_Forests += area
                    elseif land[r,c]==6 && regions[r,c]!=0
                        Closed_Shrublands += area
                    elseif land[r,c]==7 && regions[r,c]!=0
                        Open_Shrublands += area
                    elseif land[r,c]==8 && regions[r,c]!=0
                        Woody_Savannas += area
                    elseif land[r,c]==9 && regions[r,c]!=0
                        Savannas += area
                    elseif land[r,c]==10 && regions[r,c]!=0
                        Grasslands += area
                    elseif land[r,c]==11 && regions[r,c]!=0
                        Permanent_Wetlands += area
                    elseif land[r,c]==12 && regions[r,c]!=0
                        Croplands += area
                    elseif land[r,c]==13 && regions[r,c]!=0
                        Urban += area
                    elseif land[r,c]==14 && regions[r,c]!=0
                        Cropland_Natural += area
                    elseif land[r,c]==15 && regions[r,c]!=0
                        Snow_Ice += area
                    elseif land[r,c]==16 && regions[r,c]!=0
                        Barren += area
                    end

                end
            end
        end
    end

    return  Evergreen_Needleleaf_Forests,
     Evergreen_Broadleaf_Forests,
     Deciduous_Needleleaf_Forests,
     Deciduous_Broadleaf_Forests,
     Mixed_Forests,
     Closed_Shrublands,
     Open_Shrublands,
     Woody_Savannas,
     Savannas,
     Grasslands,
     Permanent_Wetlands,
     Croplands,
     Urban,
     Cropland_Natural,
     Snow_Ice,
     Barren
end
function calculate_protected_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, protected)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore = loptions

    global Strict_Nature_Reserve = 0;
    global Wilderness_Area = 0;
    global National_Park = 0;
    global Natural_Monument = 0;
    global Habitat_Species_Management = 0;
    global Protected_Landscape_Seascape = 0;
    global Managed_Resource_Protected_Area = 0;
    global Not_Reported = 0;
    global Not_Applicable = 0;
    global Not_Assigned = 0;

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if protected[r,c]==1 && regions[r,c]!=0
                        Strict_Nature_Reserve += area
                    elseif protected[r,c]==2 && regions[r,c]!=0
                        Wilderness_Area += area
                    elseif protected[r,c]==3 && regions[r,c]!=0
                        National_Park += area
                    elseif protected[r,c]==4 && regions[r,c]!=0
                        Natural_Monument += area
                    elseif protected[r,c]==5 && regions[r,c]!=0
                        Habitat_Species_Management += area
                    elseif protected[r,c]==6 && regions[r,c]!=0
                        Protected_Landscape_Seascape += area
                    elseif protected[r,c]==7 && regions[r,c]!=0
                        Managed_Resource_Protected_Area += area
                    elseif protected[r,c]==8 && regions[r,c]!=0
                        Not_Reported += area
                    elseif protected[r,c]==9 && regions[r,c]!=0
                        Not_Applicable += area
                    elseif protected[r,c]==10 && regions[r,c]!=0
                        Not_Assigned += area
                    end

                end
            end
        end
    end

    return  Strict_Nature_Reserve, Wilderness_Area, National_Park, Natural_Monument, Habitat_Species_Management,
        Protected_Landscape_Seascape, Managed_Resource_Protected_Area, Not_Reported, Not_Applicable, Not_Assigned
end
function calculate_popdens_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, popdens)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore, pv_plant_persons_per_km2, wind_persons_per_km2 = loptions

    global popdens_excl_area_PV = 0;
    global popdens_excl_area_Onshore = 0;

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if popdens[r,c]>=pv_plant_persons_per_km2 && regions[r,c]!=0
                        popdens_excl_area_PV += area
                    end

                    if popdens[r,c]>=wind_persons_per_km2 && regions[r,c]!=0
                        popdens_excl_area_Onshore += area
                    end
                end
            end
        end
    end
    return popdens_excl_area_PV, popdens_excl_area_Onshore
end
function calculate_slope_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, slope, max_slope_pv, max_slope_onshore)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore, pv_plant_persons_per_km2, wind_persons_per_km2 = loptions

    global slope_excl_area_PV = 0;
    global slope_excl_area_Onshore = 0;

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if slope[r,c]>=max_slope_pv && regions[r,c]!=0
                        slope_excl_area_PV += area
                    end

                    if slope[r,c]>=max_slope_onshore && regions[r,c]!=0
                        slope_excl_area_Onshore += area
                    end
                end
            end
        end
    end
    return slope_excl_area_PV, slope_excl_area_Onshore
end
function calculate_aqueduct_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, mask_aqueduct_lcoh)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore, pv_plant_persons_per_km2, wind_persons_per_km2 = loptions

    global aqueduct_excl_area = 0;

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if mask_aqueduct_lcoh[r,c]==0 && regions[r,c]>0 && regions[r,c]!=NOREGION
                        aqueduct_excl_area += area
                    end
                end
            end
        end
    end
    return aqueduct_excl_area
end
function calculate_pv_resource_areas(loptions, regions, solarGTI, lonrange, latrange, mask_pv_lcoh, mask_onshore_lcoh, mask_offshore_lcoh, resource_quality_mask_pv, resource_quality_mask_onshore, resource_quality_mask_offshore)

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(loptions, lonrange, latrange)

    nlons = size(solarGTI)[2]
    nlats = size(solarGTI)[3]

    @unpack era_year, res, erares, pv_density, plant_area,
            onshore_density, area_onshore, offshore_density, area_offshore, pv_plant_persons_per_km2, wind_persons_per_km2 = loptions

    global resource_excl_area_PV = 0;
    global resource_excl_area_Onshore = 0;
    global resource_excl_area_Offshore = 0;

    for j in 1:nlats

        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]

        for i in 1:nlons

            eralon = eralons[i]
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange
                for r in rowrange

                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    (reg == NOREGION) && continue
                    area = cellarea[c]

                    if resource_quality_mask_pv[r,c]==0 && regions[r,c]>0 && regions[r,c]!=NOREGION
                        resource_excl_area_PV += area
                    end

                    if resource_quality_mask_onshore[r,c]==0 && regions[r,c]>0 && regions[r,c]!=NOREGION
                        resource_excl_area_Onshore += area
                    end

                    if resource_quality_mask_offshore[r,c]==0 && regions[r,c]==0
                        resource_excl_area_Offshore += area
                    end
                end
            end
        end
    end
    return resource_excl_area_PV, resource_excl_area_Onshore, resource_excl_area_Offshore
end
