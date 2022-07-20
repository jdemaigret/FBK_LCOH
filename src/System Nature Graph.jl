using DelimitedFiles, Plots
region = "GER"
time_horizon = 2050
# time_horizon = 2030
scenario = "Optimistic"
# scenario = "Pessimistic"

# Country    Best  PV_Best  Onshore_Best
# AUS    ->   1-2    1-4     3-2
# CHL    ->   1-1    1-3     3-1
# CHI    ->   1-2    1-4     3-2
# GER    ->   3-2    3-4     4-2
# JAP    ->   2-2    2-4     4-2
# SAU    ->   1-3    1-5     2-3
# USA    ->   1-1    1-3     3-1

PV_quality = 4
Onshore_quality = 4
min_PV_CAPEX = 100
max_PV_CAPEX = 700
length_pv = 50
min_Onshore_CAPEX = 600
max_Onshore_CAPEX = 1100
length_onshore = 50
global max_cap_pv   = 45000;
global max_cap_wind = 5000;
global tot_points = length_pv * length_onshore;
global contatore = 0;
println("PV_quality = $PV_quality")
println("Onshore_quality = $Onshore_quality")
total_cf_data = readdlm("Characteristic CF Profiles\\$region charachteristic_cf_curves V9.txt",',', skipstart=1)
global CF_pvplant = total_cf_data[:,PV_quality]
global CF_onshore = total_cf_data[:,Onshore_quality+5]
if (scenario == "Optimistic") && (time_horizon == 2050)
    global CAPEX_el = 133.77;
    global eff_el = 0.875;
elseif (scenario == "Pessimistic") && (time_horizon == 2050)
    global CAPEX_el = 326.27;
    global eff_el = 0.8204;
elseif (scenario == "Optimistic") && (time_horizon == 2030)
    global CAPEX_el = 383.77;
    global eff_el = 0.812;
elseif (scenario == "Pessimistic") && (time_horizon == 2030)
    global CAPEX_el = 688.43;
    global eff_el = 0.754;
end
global HHV_H2 = 39.4;
global H2_energy_consumption = HHV_H2/eff_el;
global PV_CAPEXs = collect(range(min_PV_CAPEX, length = length_pv, stop=max_PV_CAPEX))
global Onshore_CAPEXs = collect(range(min_Onshore_CAPEX, length = length_onshore, stop=max_Onshore_CAPEX))
WACC_data = readdlm("WACC_CAPEX Versione 9.csv" ,',')
nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
WACC_country_index = findfirst(isequal("$region"), nomi_stringa)[2];
if scenario == "Optimistic"
    scenario_col = 11
else
    scenario_col = 22
end
global WACC_pv      = WACC_data[WACC_country_index,scenario_col];
global WACC_onshore = WACC_data[WACC_country_index,scenario_col+1];
global WACC_el = (WACC_pv+WACC_onshore)/2;
global gen_tech_ratio_mat = zeros(length(PV_CAPEXs), length(Onshore_CAPEXs))
global lcoh_mat = zeros(length(PV_CAPEXs), length(Onshore_CAPEXs))
for r in 1:length(PV_CAPEXs)
    for s in 1:length(Onshore_CAPEXs)
        global contatore = contatore + 1;
        println("Point: $contatore/$tot_points")
        lifetime_el   = 25;
        lifetime_pv   = 25;
        lifetime_wind = 25;
        OPEX_el       = 0.03;
        OPEX_wind     = 0.03;
        OPEX_pv       = 0.01;
        global system_lifetime = min(lifetime_el, lifetime_pv, lifetime_wind);
        global generation_system_lifetime = min(lifetime_pv, lifetime_wind);
        cost_pv   = (PV_CAPEXs[r]*WACC_pv)/(1-(1+WACC_pv)^(-lifetime_pv))+(PV_CAPEXs[r]*OPEX_pv);
        cost_wind = (Onshore_CAPEXs[s]*WACC_onshore)/(1-(1+WACC_onshore)^(-lifetime_wind))+(Onshore_CAPEXs[s]*OPEX_wind);
        cost_el   = (CAPEX_el*WACC_el)/(1-(1+WACC_el)^(-lifetime_el))+(CAPEX_el*OPEX_el);
        global LCOH_min     = 0;
        global xmin_LCOH    = 0;
        global ymin_LCOH    = 0;
        global zmin_LCOH    = 0;
        global startvalue_pv   = 0;
        global startvalue_wind = 0;
        global endvalue_pv     = max_cap_pv;
        global endvalue_wind   = max_cap_wind;
        global startvalue_el   = 1;
        global endvalue_el     = max_cap_pv + max_cap_wind;
        global dd = (7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6);
        # global dd = (5,4,4,4,4,4);
        @time for l in 1:length(dd)
            d=dd[l];
            LCOH_min0 = LCOH_min;
            if l==1
                el_over_pv0      = 0;
                el_over_onshore0 = 0;
            else
                el_over_pv0      = zmin_LCOH/xmin_LCOH;
                el_over_onshore0 = zmin_LCOH/ymin_LCOH;
            end
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
            # println("el_over_pv_residual: $el_over_pv_residual")
            # println("el_over_onshore_residual: $el_over_onshore_residual")
            # (Precedente-Attuale)/Precedente
            # el_over_pv_residual_percent      = el_over_pv_residual      / el_over_pv0;
            # el_over_onshore_residual_percent = el_over_onshore_residual / el_over_onshore0;
            # println("el_over_pv_residual_percent = $el_over_pv_residual_percent")
            # println("el_over_onshore_residual_percent = $el_over_onshore_residual_percent")
            # println("LCOH_residual / el_over_pv_residual / el_over_onshore_residual (type_min=1) in=$l : $LCOH_residual / $el_over_pv_residual / $el_over_onshore_residual")
            el_over_pv_residual      = abs(el_over_pv0-zmin_LCOH/xmin_LCOH);
            el_over_onshore_residual = abs(el_over_onshore0-zmin_LCOH/ymin_LCOH);
            LCOH_residual = abs(LCOH_min0-LCOH_min)
            el_over_pv_residual_percent      = el_over_pv_residual      / (zmin_LCOH/xmin_LCOH);
            el_over_onshore_residual_percent = el_over_onshore_residual / (zmin_LCOH/ymin_LCOH);
            ((el_over_pv_residual_percent < 0.01) && (el_over_onshore_residual_percent < 0.01)) && break
        end
        println("xmin_LCOH = $xmin_LCOH")
        println("ymin_LCOH = $ymin_LCOH")
        println("zmin_LCOH = $zmin_LCOH")
        println("LCOH_min = $LCOH_min")

        global gen_tech_ratio_mat[r,s] = xmin_LCOH/(ymin_LCOH+xmin_LCOH)
        global lcoh_mat[r,s] = LCOH_min
        # Uscito da il ciclo per la lunghezza di d ho un LCOH ottimizzato per:
        # Una combinazione di CAPEX di wind e pv
    end
end


gen_tech_ratio_mat_inv = zeros(size(gen_tech_ratio_mat))
lcoh_mat_inv = zeros(size(lcoh_mat))
global count=0
for colonna in 1:size(gen_tech_ratio_mat)[2]
    gen_tech_ratio_mat_inv[:,end-count]=gen_tech_ratio_mat[:,colonna]
    lcoh_mat_inv[:,end-count]=lcoh_mat[:,colonna]
    global count+=1
end
gen_tech_ratio_mat_inv_trans=transpose(gen_tech_ratio_mat_inv)
lcoh_mat_inv_trans=transpose(lcoh_mat_inv)

my_color = cgrad([:blue, :green ,:orange])
heatmap(PV_CAPEXs, Onshore_CAPEXs, gen_tech_ratio_mat_inv_trans; color=my_color, xlabel="PV CAPEX [USD/kW]", ylabel="Onshore CAPEX [USD/kW]",  title="PV Class $PV_quality - Onshore Class $Onshore_quality")
heatmap(PV_CAPEXs, Onshore_CAPEXs, lcoh_mat_inv_trans, xlabel="PV CAPEX [USD/kW]", ylabel="Onshore CAPEX [USD/kW]",  title="PV Class $PV_quality - Onshore Class $Onshore_quality")
filename_gen_tech_ratio_mat = "$region Gen Tech Ratio Matrix $time_horizon $scenario (pv$PV_quality, w$Onshore_quality).txt"
filename_lcoh_mat = "$region LCOH Matrix $time_horizon $scenario (pv$PV_quality, w$Onshore_quality).txt"
isfile(filename_gen_tech_ratio_mat) && rm(filename_gen_tech_ratio_mat);
isfile(filename_lcoh_mat) && rm(filename_lcoh_mat);
writedlm(filename_gen_tech_ratio_mat, gen_tech_ratio_mat_inv_trans, '\t');
writedlm(filename_lcoh_mat, lcoh_mat_inv_trans, '\t');

# pv1_w1 = readdlm("HEATMAPS\\CHL Gen Tech Ratio Matrix 2050 Optimistic (pv1, w1).txt")
# pv1_w3 = readdlm("HEATMAPS\\CHL Gen Tech Ratio Matrix 2050 Optimistic (pv1, w3).txt")
# pv3_w1 = readdlm("HEATMAPS\\CHL Gen Tech Ratio Matrix 2050 Optimistic (pv3, w1).txt")
#
# pv1_w1_lcoh = readdlm("HEATMAPS\\CHL LCOH Matrix 2050 Optimistic (pv1, w1).txt")
# pv1_w3_lcoh = readdlm("HEATMAPS\\CHL LCOH Matrix 2050 Optimistic (pv1, w3).txt")
# pv3_w1_lcoh = readdlm("HEATMAPS\\CHL LCOH Matrix 2050 Optimistic (pv3, w1).txt")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv1_w1; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="pv1_w1")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv1_w3; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="pv1_w3")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv3_w1; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="pv3_w1")
#
#
# pv1_w1_inverted = zeros(size(pv1_w1))
# pv1_w1_lcoh_inverted = zeros(size(pv1_w1_lcoh))
# count=0
# for colonna in 1:size(pv1_w1)[2]
#     pv1_w1_inverted[:,end-count]=pv1_w1[:,colonna]
#     pv1_w1_lcoh_inverted[:,end-count]=pv1_w1_lcoh[:,colonna]
#     count+=1
# end
# pv1_w1_inverted_transposed=transpose(pv1_w1_inverted)
# pv1_w1_lcoh_inverted_transposed=transpose(pv1_w1_lcoh_inverted)
#
# pv1_w3_inverted = zeros(size(pv1_w3))
# pv1_w3_lcoh_inverted = zeros(size(pv1_w3_lcoh))
# count=0
# for colonna in 1:size(pv1_w3)[2]
#     pv1_w3_inverted[:,end-count]=pv1_w3[:,colonna]
#     pv1_w3_lcoh_inverted[:,end-count]=pv1_w3_lcoh[:,colonna]
#     count+=1
# end
# pv1_w3_inverted_transposed=transpose(pv1_w3_inverted)
# pv1_w3_lcoh_inverted_transposed=transpose(pv1_w3_lcoh_inverted)
#
# pv3_w1_inverted = zeros(size(pv3_w1))
# pv3_w1_lcoh_inverted = zeros(size(pv3_w1_lcoh))
# count=0
# for colonna in 1:size(pv3_w1)[2]
#     pv3_w1_inverted[:,end-count]=pv3_w1[:,colonna]
#     pv3_w1_lcoh_inverted[:,end-count]=pv3_w1_lcoh[:,colonna]
#     count+=1
# end
# pv3_w1_inverted_transposed=transpose(pv3_w1_inverted)
# pv3_w1_lcoh_inverted_transposed=transpose(pv3_w1_lcoh_inverted)
#
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv1_w1_inverted_transposed; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="CHL ratio pv1_w1")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv1_w3_inverted_transposed; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="CHL ratio pv1_w3")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv3_w1_inverted_transposed; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="CHL ratio pv3_w1")
#
#
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv1_w1_lcoh_inverted; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="CHL LCOH pv1_w1")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv1_w3_lcoh_inverted; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="CHL LCOH pv1_w3")
# heatmap(Onshore_CAPEXs, PV_CAPEXs, pv3_w1_lcoh_inverted; color=my_color, xlabel="Onshore CAPEX [USD/kW]", ylabel="PV CAPEX [USD/kW]", title="CHL LCOH pv3_w1")
