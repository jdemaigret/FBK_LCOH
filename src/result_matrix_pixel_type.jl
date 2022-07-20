using DelimitedFiles

country = "ITA"
scenario = "Come Plexos"
# scenario = "Come PLEXOS con BWS"

supply_cost_curve_data = readdlm("C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Da mostrare\\V13\\File di testo\\$country results_matrix (come PLEXOS).txt", skipstart=1)

LCOH            = supply_cost_curve_data[:,1];
prod_kg         = supply_cost_curve_data[:,2];
pv              = supply_cost_curve_data[:,3];
onshore         = supply_cost_curve_data[:,4];
offshore        = supply_cost_curve_data[:,5];
electrolyzer    = supply_cost_curve_data[:,6];
FLOH            = supply_cost_curve_data[:,7];
sum_cf_pv       = supply_cost_curve_data[:,8];
sum_cf_onshore  = supply_cost_curve_data[:,9];
sum_cf_offshore = supply_cost_curve_data[:,10];
curtailment     = supply_cost_curve_data[:,11];
LCOE_no_curt    = supply_cost_curve_data[:,12];

N = length(LCOH)

deleted_index = typemax(Int32)

for i in 1:N
    if LCOH[i] == Inf
        LCOH[i]            = deleted_index;
        prod_kg[i]         = deleted_index;
        pv[i]              = deleted_index;
        onshore[i]         = deleted_index;
        offshore[i]        = deleted_index;
        electrolyzer[i]    = deleted_index;
        FLOH[i]            = deleted_index;
        sum_cf_pv[i]       = deleted_index;
        sum_cf_onshore[i]  = deleted_index;
        sum_cf_offshore[i] = deleted_index;
        curtailment[i]     = deleted_index;
        LCOE_no_curt[i]    = deleted_index
    end
end

LCOH_clean            = zeros(length(LCOH));
prod_kg_clean         = zeros(length(prod_kg));
pv_clean              = zeros(length(pv));
onshore_clean         = zeros(length(onshore));
offshore_clean        = zeros(length(offshore));
electrolyzer_clean    = zeros(length(electrolyzer));
FLOH_clean            = zeros(length(FLOH));
sum_cf_pv_clean       = zeros(length(sum_cf_pv));
sum_cf_onshore_clean  = zeros(length(sum_cf_onshore));
sum_cf_offshore_clean = zeros(length(sum_cf_offshore));
curtailment_clean     = zeros(length(curtailment));
LCOE_no_curt_clean    = zeros(length(LCOE_no_curt));

LCOH_clean            = filter(!isequal(deleted_index), LCOH);
prod_kg_clean         = filter(!isequal(deleted_index), prod_kg);
pv_clean              = filter(!isequal(deleted_index), pv);
onshore_clean         = filter(!isequal(deleted_index), onshore);
offshore_clean        = filter(!isequal(deleted_index), offshore);
electrolyzer_clean    = filter(!isequal(deleted_index), electrolyzer);
FLOH_clean            = filter(!isequal(deleted_index), FLOH);
sum_cf_pv_clean       = filter(!isequal(deleted_index), sum_cf_pv);
sum_cf_onshore_clean  = filter(!isequal(deleted_index), sum_cf_onshore);
sum_cf_offshore_clean = filter(!isequal(deleted_index), sum_cf_offshore);
curtailment_clean     = filter(!isequal(deleted_index), curtailment);
LCOE_no_curt_clean    = filter(!isequal(deleted_index), LCOE_no_curt);

supply_cost_curve_data_clean = [ LCOH_clean prod_kg_clean pv_clean onshore_clean offshore_clean electrolyzer_clean FLOH_clean sum_cf_pv_clean sum_cf_onshore_clean sum_cf_offshore_clean curtailment_clean LCOE_no_curt_clean];

supply_cost_curve_data_clean_sorted = sortslices(supply_cost_curve_data_clean, dims=1);

N_clean_sorted = size(supply_cost_curve_data_clean_sorted)[1];

prod_PJ         = supply_cost_curve_data_clean_sorted[:,2] .* 120 ./ 10^9;
prog_prod_PJ    = zeros(length(prod_PJ));
prog_prod_PJ[1] = prod_PJ[1];

for i in 2:N_clean_sorted
    prog_prod_PJ[i] = prog_prod_PJ[i-1] + prod_PJ[i];
end

supply_cost_curve_data_clean_sorted_w_prog = hcat(supply_cost_curve_data_clean_sorted,prog_prod_PJ);
res_mat = supply_cost_curve_data_clean_sorted_w_prog

pixel_type_vec = Vector{String}(undef, size(res_mat)[1])

for i in 1:length(pixel_type_vec)
    if (res_mat[i,3]!=0) && (res_mat[i,4]==0) && (res_mat[i,5]==0)
        pixel_type_vec[i] = "PV"
    elseif (res_mat[i,3]==0) && (res_mat[i,4]!=0) && (res_mat[i,5]==0)
        pixel_type_vec[i] = "Onshore"
    elseif (res_mat[i,3]==0) && (res_mat[i,4]==0) && (res_mat[i,5]!=0)
        pixel_type_vec[i] = "Offshore"
    elseif (res_mat[i,3]!=0) && (res_mat[i,4]!=0) && ((res_mat[i,4]/res_mat[i,3])<0.05)
        pixel_type_vec[i] = "Hybrid -> PV"
    elseif (res_mat[i,3]!=0) && (res_mat[i,4]!=0) && ((res_mat[i,3]/res_mat[i,4])<0.05)
        pixel_type_vec[i] = "Hybrid -> Onshore"
    else
        pixel_type_vec[i] = "Hybrid"
    end
end

new_PV_pot_vec      = zeros(length(pixel_type_vec))
new_Onshore_pot_vec = zeros(length(pixel_type_vec))

for i in 1:length(pixel_type_vec)
    if pixel_type_vec[i] == "PV"
        new_PV_pot_vec[i] = res_mat[i,3];
    elseif pixel_type_vec[i] == "Onshore"
        new_Onshore_pot_vec[i] = res_mat[i,4];
    elseif pixel_type_vec[i] == "Hybrid -> PV"
        new_PV_pot_vec[i]  = res_mat[i,3];
    elseif pixel_type_vec[i] == "Hybrid -> Onshore"
        new_Onshore_pot_vec[i] = res_mat[i,3];
    elseif pixel_type_vec[i] == "Hybrid"
        new_PV_pot_vec[i]  = res_mat[i,3];
        new_Onshore_pot_vec[i] = res_mat[i,3];
    end
end
type_mat = [pixel_type_vec new_PV_pot_vec new_Onshore_pot_vec];
fin_res_mat = hcat(res_mat, type_mat)
