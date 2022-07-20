using DelimitedFiles, Dierckx#, LsqFit


nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
# nomi_stringa = ["AUS" "BRA" "CAN" "CHI" "CHL" "JAP" "RUS" "SAU" "SSA" "USA"]
nomi_stringa = ["IDN"]
# scenario = "Come PLEXOS con BWS"

nomi_stringa = ["COD" "EGY" "ETH" "KEN" "MRT" "NAM"]
# scenario_string = ["2030 Opt" "2030 Opt ws" "2030 Pess" "2030 Pess ws" "2050 Opt" "2050 Opt ws"  "2050 Pess" "2050 Pess ws"]
scenario_string = ["2050 Opt ws" "2050 Pess ws"]
for scenario in scenario_string
    for country in nomi_stringa

        println("Region: $country")
        println("Scenario: $scenario")
        # local scenario = "2050 V Pess"
        # local scenario = "2050 Pess ws"

        # local supply_cost_curve_data = readdlm("C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\File di testo\\$country results_matrix ($scenario).txt", skipstart=1)
        # local supply_cost_curve_data = readdlm("E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\File di testo\\$country results_matrix ($scenario).txt", skipstart=1)
        # local supply_cost_curve_data = readdlm("E:\\Users\\jdemaigret\\Documents\\EF\\Output\\Water stress considered\\$scenario\\File di testo\\$country results_matrix ($scenario).txt", skipstart=1)
        local supply_cost_curve_data = readdlm("C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\Report\\Africa BOX\\$scenario\\$country\\$country results_matrix ($scenario).txt", skipstart=1)
        local LCOH            = supply_cost_curve_data[:,1];
        local prod_kg         = supply_cost_curve_data[:,2];
        local pv              = supply_cost_curve_data[:,3];
        local onshore         = supply_cost_curve_data[:,4];
        local offshore        = supply_cost_curve_data[:,5];
        local electrolyzer    = supply_cost_curve_data[:,6];
        local FLOH            = supply_cost_curve_data[:,7];
        local sum_cf_pv       = supply_cost_curve_data[:,8];
        local sum_cf_onshore  = supply_cost_curve_data[:,9];
        local sum_cf_offshore = supply_cost_curve_data[:,10];
        local curtailment     = supply_cost_curve_data[:,11];

        local N = length(LCOH)

        local deleted_index = typemax(Int32)

        for i in 1:N
            if LCOH[i] == Inf
                local LCOH[i]            = deleted_index;
                local prod_kg[i]         = deleted_index;
                local pv[i]              = deleted_index;
                local onshore[i]         = deleted_index;
                local offshore[i]        = deleted_index;
                local electrolyzer[i]    = deleted_index;
                local FLOH[i]            = deleted_index;
                local sum_cf_pv[i]       = deleted_index;
                local sum_cf_onshore[i]  = deleted_index;
                local sum_cf_offshore[i] = deleted_index;
                local curtailment[i]     = deleted_index;
            end
        end

        local LCOH_clean            = zeros(length(LCOH));
        local prod_kg_clean         = zeros(length(prod_kg));
        local pv_clean              = zeros(length(pv));
        local onshore_clean         = zeros(length(onshore));
        local offshore_clean        = zeros(length(offshore));
        local electrolyzer_clean    = zeros(length(electrolyzer));
        local FLOH_clean            = zeros(length(FLOH));
        local sum_cf_pv_clean       = zeros(length(sum_cf_pv));
        local sum_cf_onshore_clean  = zeros(length(sum_cf_onshore));
        local sum_cf_offshore_clean = zeros(length(sum_cf_offshore));
        local curtailment_clean     = zeros(length(curtailment));

        local LCOH_clean            = filter(!isequal(deleted_index), LCOH);
        local prod_kg_clean         = filter(!isequal(deleted_index), prod_kg);
        local pv_clean              = filter(!isequal(deleted_index), pv);
        local onshore_clean         = filter(!isequal(deleted_index), onshore);
        local offshore_clean        = filter(!isequal(deleted_index), offshore);
        local electrolyzer_clean    = filter(!isequal(deleted_index), electrolyzer);
        local FLOH_clean            = filter(!isequal(deleted_index), FLOH);
        local sum_cf_pv_clean       = filter(!isequal(deleted_index), sum_cf_pv);
        local sum_cf_onshore_clean  = filter(!isequal(deleted_index), sum_cf_onshore);
        local sum_cf_offshore_clean = filter(!isequal(deleted_index), sum_cf_offshore);
        local curtailment_clean     = filter(!isequal(deleted_index), curtailment);

        local supply_cost_curve_data_clean = [ LCOH_clean prod_kg_clean pv_clean onshore_clean offshore_clean electrolyzer_clean FLOH_clean sum_cf_pv_clean sum_cf_onshore_clean sum_cf_offshore_clean curtailment_clean];

        local supply_cost_curve_data_clean_sorted = sortslices(supply_cost_curve_data_clean, dims=1);

        local N_clean_sorted = size(supply_cost_curve_data_clean_sorted)[1];

        local prod_PJ         = supply_cost_curve_data_clean_sorted[:,2] .* 120 ./ 10^9;
        local prog_prod_PJ    = zeros(length(prod_PJ));
        local prog_prod_PJ[1] = prod_PJ[1];

        for i in 2:N_clean_sorted
            local prog_prod_PJ[i] = prog_prod_PJ[i-1] + prod_PJ[i];
        end

        local supply_cost_curve_data_clean_sorted_w_prog = hcat(supply_cost_curve_data_clean_sorted,prog_prod_PJ);

        local pv_onshore_capacities = zeros(length(1:N_clean_sorted),2); #1 PV 2 ONSHORE
        local pv_capacities         = zeros(length(1:N_clean_sorted));
        local onshore_capacities    = zeros(length(1:N_clean_sorted));
        local offshore_capacities   = zeros(length(1:N_clean_sorted));
        local pv_onshore_cf_sum     = zeros(length(1:N_clean_sorted),2);
        local pv_cf_sum             = zeros(length(1:N_clean_sorted));
        local onshore_cf_sum        = zeros(length(1:N_clean_sorted));
        local offshore_cf_sum       = zeros(length(1:N_clean_sorted));


        for i in 1:N_clean_sorted
            if supply_cost_curve_data_clean_sorted_w_prog[i,3]!=0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]!=0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]==0
                local pv_onshore_cf_sum[i,1] = supply_cost_curve_data_clean_sorted_w_prog[i,8];
                local pv_onshore_cf_sum[i,2] = supply_cost_curve_data_clean_sorted_w_prog[i,9];
                local pv_onshore_capacities[i,1] = supply_cost_curve_data_clean_sorted_w_prog[i,3];
                local pv_onshore_capacities[i,2] = supply_cost_curve_data_clean_sorted_w_prog[i,4];
            elseif supply_cost_curve_data_clean_sorted_w_prog[i,3]!=0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]==0
                local pv_cf_sum[i] =  supply_cost_curve_data_clean_sorted_w_prog[i,8];
                #supply_cost_curve_data_clean_sorted_w_prog[i,9]=0;
                #supply_cost_curve_data_clean_sorted_w_prog[i,10]=0;
                local pv_capacities[i] = supply_cost_curve_data_clean_sorted_w_prog[i,3];
            elseif supply_cost_curve_data_clean_sorted_w_prog[i,3]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]!=0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]==0
                local onshore_cf_sum[i] =  supply_cost_curve_data_clean_sorted_w_prog[i,9];
                #supply_cost_curve_data_clean_sorted_w_prog[i,8]=0;
                #supply_cost_curve_data_clean_sorted_w_prog[i,10]=0;
                local onshore_capacities[i] = supply_cost_curve_data_clean_sorted_w_prog[i,4];
            elseif supply_cost_curve_data_clean_sorted_w_prog[i,3]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]!=0
                local offshore_cf_sum[i] =  supply_cost_curve_data_clean_sorted_w_prog[i,10];
                #supply_cost_curve_data_clean_sorted_w_prog[i,8]=0;
                #supply_cost_curve_data_clean_sorted_w_prog[i,9]=0;
                local offshore_capacities[i] = supply_cost_curve_data_clean_sorted_w_prog[i,5];
            end
        end

        # Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
        # x e y sarebbero gli array che voglio iterpolare
        # k sarebbe il grado della Spline
        # bc sarebbe il comportamento della spline agli estremi del dominio di supporto
        # i quali sarebbero minimo e massimo di x. I valori che può assumere bc sono:
        # nearest, zero,extrapolate, error. Secondo me è da lasciare nearest.
        # s ivnece indica lo smoothing factor.
        # posso fare la derivata di grado nu-esimo nel punto x con
        # derivative(spl, x; nu=1). Posso fare l'integrale da x=a a x=b con
        # integrate(spl, a, b).Posso trovare gli zeri con roots(spl; maxn=8). maxn mi
        # da il numero massimo di radici da ritornare.
        #  penso che w siano i weight che vengono asegnati.

        local x = supply_cost_curve_data_clean_sorted_w_prog[:,12];
        local y = supply_cost_curve_data_clean_sorted_w_prog[:,1];

        local fit1 = Spline1D(x,y,w=ones(length(x)), k=1, bc="nearest", s=0.0); #LCOH in funzione di PJ
        local fit2 = Spline1D(y,x,w=ones(length(y)), k=1, bc="nearest", s=0.0); #PJ in funzione di LCOH

        local xs = range(0; length= 1000, stop=last(x))
        local ys = fit1(xs)
        # local path = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\Interpolate"
        # local path = "E:\\Users\\jdemaigret\\Documents\\EF\\Output\\Water stress considered\\$scenario\\Interpolate"
        # local path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\Interpolate"
        local path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\Report\\Africa BOX\\$scenario\\Interpolate"
        local filename = ("$path\\$country Supply-Cost Curve Interpolate $scenario.txt")
        isfile(filename) && rm(filename)
        writedlm(filename, [xs  ys], '\t')
    end
end
    # LCOH_upper_lim = 2
    # Δ_LCOH = 0.05
    # prod_0 = 0
    #
    # yt = range(0.89; length=round(Int, (2.45-0.89)/0.05), stop=2.45)
    # xt = fit2(yt)
    # #prod_upper_lim = fit2(LCOH_upper_lim)
    #
    # supply_cost_step_PROD = zeros()
    # supply_cost_step_LCOH = zeros()
    #
    # supply_cost_step_PROD[1] = prod_0
    # supply_cost_step_LCOH[1] = (fit1(xs[1]) + fit1(xs[1]) + Δ_LCOH) / 2;
    #
    # for i in 2:size(supply_cost_step)[1]
    #
    #     if  i%2==0
    #         supply_cost_step[i,1] = prod0 + fit2(fit1(prod_0))
    #         supply_cost_step[i,2] = fit1(xs[i-1])
    #     end
    #
    #     supply_cost_step[i,1] = fit2(fit1(xs[i-1]))
    #     supply_cost_step[i,2] = fit1(xs[i])
    #
    # end
    #
    #
    # fit1(xs[1]) + Δ_LCOH #LCOH al punto dopo
    # fit2(fit1(xs[1]) + Δ_LCOH) #PJ corrispondente al punto dopo
