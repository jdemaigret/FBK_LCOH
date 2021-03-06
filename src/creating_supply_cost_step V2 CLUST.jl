using DelimitedFiles, Clustering,  Plots

# nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
nomi_stringa = ["OCE"]

for reg in nomi_stringa

    println("reg")
    println(reg)

    local country = reg
    local scenario = "Come Plexos"
    # scenario = "Come PLEXOS con BWS"

    local supply_cost_curve_data = readdlm("C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Da mostrare\\V13\\File di testo\\$country results_matrix (come PLEXOS).txt", skipstart=1)


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

    local supply_cost_curve_data_clean = [ LCOH_clean prod_kg_clean pv_clean onshore_clean offshore_clean electrolyzer_clean FLOH_clean sum_cf_pv_clean sum_cf_onshore_clean sum_cf_offshore_clean curtailment_clean];

    local supply_cost_curve_data_clean_sorted = sortslices(supply_cost_curve_data_clean, dims=1);

    local N_clean_sorted = size(supply_cost_curve_data_clean_sorted)[1];

    local prod_PJ         = supply_cost_curve_data_clean_sorted[:,2] .* 120 ./ 10^9;
    local prog_prod_PJ    = zeros(length(prod_PJ));
    local prog_prod_PJ[1] = prod_PJ[1];

    for i in 2:N_clean_sorted
        prog_prod_PJ[i] = prog_prod_PJ[i-1] + prod_PJ[i];
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
            pv_onshore_cf_sum[i,1] = supply_cost_curve_data_clean_sorted_w_prog[i,8];
            pv_onshore_cf_sum[i,2] = supply_cost_curve_data_clean_sorted_w_prog[i,9];
            pv_onshore_capacities[i,1] = supply_cost_curve_data_clean_sorted_w_prog[i,3];
            pv_onshore_capacities[i,2] = supply_cost_curve_data_clean_sorted_w_prog[i,4];
        elseif supply_cost_curve_data_clean_sorted_w_prog[i,3]!=0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]==0
            pv_cf_sum[i] =  supply_cost_curve_data_clean_sorted_w_prog[i,8];
            #supply_cost_curve_data_clean_sorted_w_prog[i,9]=0;
            #supply_cost_curve_data_clean_sorted_w_prog[i,10]=0;
            pv_capacities[i] = supply_cost_curve_data_clean_sorted_w_prog[i,3];
        elseif supply_cost_curve_data_clean_sorted_w_prog[i,3]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]!=0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]==0
            onshore_cf_sum[i] =  supply_cost_curve_data_clean_sorted_w_prog[i,9];
            #supply_cost_curve_data_clean_sorted_w_prog[i,8]=0;
            #supply_cost_curve_data_clean_sorted_w_prog[i,10]=0;
            onshore_capacities[i] = supply_cost_curve_data_clean_sorted_w_prog[i,4];
        elseif supply_cost_curve_data_clean_sorted_w_prog[i,3]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,4]==0 && supply_cost_curve_data_clean_sorted_w_prog[i,5]!=0
            offshore_cf_sum[i] =  supply_cost_curve_data_clean_sorted_w_prog[i,10];
            #supply_cost_curve_data_clean_sorted_w_prog[i,8]=0;
            #supply_cost_curve_data_clean_sorted_w_prog[i,9]=0;
            offshore_capacities[i] = supply_cost_curve_data_clean_sorted_w_prog[i,5];
        end
    end

    # Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
    # x e y sarebbero gli array che voglio iterpolare
    # k sarebbe il grado della Spline
    # bc sarebbe il comportamento della spline agli estremi del dominio di supporto
    # i quali sarebbero minimo e massimo di x. I valori che pu?? assumere bc sono:
    # nearest, zero,extrapolate, error. Secondo me ?? da lasciare nearest.
    # s ivnece indica lo smoothing factor.
    # posso fare la derivata di grado nu-esimo nel punto x con
    # derivative(spl, x; nu=1). Posso fare l'integrale da x=a a x=b con
    # integrate(spl, a, b).Posso trovare gli zeri con roots(spl; maxn=8). maxn mi
    # da il numero massimo di radici da ritornare.
    #  penso che w siano i weight che vengono asegnati.

    local x = supply_cost_curve_data_clean_sorted_w_prog[:,12];
    local y = supply_cost_curve_data_clean_sorted_w_prog[:,1];

    local max_LCOH = 2

    for i in 1:length(x)
        if y[i]>max_LCOH
            y[i]=typemax(Int32)
            x[i]=typemax(Int32)
        end
    end

    filter!(!isequal(typemax(Int32)), y);
    filter!(!isequal(typemax(Int32)), x);

    local points = [transpose(x); transpose(y)]

    local n_clusters = 5

    local R = kmeans(points, n_clusters; maxiter=500, display=:iter);

    @assert nclusters(R) == n_clusters # verify the number of clusters

    local a = assignments(R) # get the assignments of points to clusters (Array{Int64,1})
    local c = counts(R) # get the cluster sizes (Array{Int64,1})
    local M = R.centers # get the cluster centers (Array{Float64,2})

    local centers = [M[2,:] M[1,:] zeros(size(M)[2]) zeros(size(M)[2]) c];

    local centers_sorted = sortslices(centers, dims=1)

    for i in 1:n_clusters
        if i==1
            centers_sorted[i,3] = 1
            centers_sorted[i,4] = centers_sorted[i,5];
            continue
        elseif i==n_clusters
            centers_sorted[i,3] = centers_sorted[i-1,4] + 1;
            centers_sorted[i,4] = centers_sorted[i,3] + centers_sorted[i,5] - 1
            break
        end
        centers_sorted[i,3] = centers_sorted[i-1,4];
        centers_sorted[i,4] = centers_sorted[i,3] + centers_sorted[i,5];
    end

    local step_matrix = zeros(size(centers_sorted))
    local step_matrix[:,1:2] = centers_sorted[:,1:2]

    for i in 1:n_clusters
        step_matrix[i,3] = x[trunc(Int, centers_sorted[i,3])]
        step_matrix[i,4] = x[trunc(Int, centers_sorted[i,4])]
        step_matrix[i,5] = x[trunc(Int, centers_sorted[i,4])] - x[trunc(Int, centers_sorted[i,3])]
    end


    local region_vector = Vector{String}(undef, n_clusters)

    for i in 1:n_clusters
        region_vector[i] = "$country"
    end

    local LCOH_class_vector = Vector{String}(undef, n_clusters)

    for i in 1:n_clusters
        LCOH_class_vector[i] = "$i"
    end
    local header = ["Region" "Class" "LCOH" "Center Potential" "Lower Lim Potential" "Upper Lim Potential" "Potential Step"]
    local step_matrix = [region_vector LCOH_class_vector step_matrix]
    local step_matrix_w_header = [header ; step_matrix]

    local path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Final\\$scenario\\Step"
    local filename = ("$path\\$country Step Supply-Cost Curve $scenario.txt")
    isfile(filename) && rm(filename)
    writedlm(filename,step_matrix_w_header , '\t')

    local step_function_prod = zeros(n_clusters * 2);

    global j = 0

    for i in 1:n_clusters
        global j += 1
        if i==1
            step_function_prod[i] = step_matrix[i,5]
            continue
        elseif i==n_clusters
            step_function_prod[i*2] = step_matrix[i,6]
        end
        step_function_prod[j] = step_matrix[i-1,6]
        j += 1
        step_function_prod[j] = step_matrix[i-1,6]
    end

    # println(step_function_prod)

    local step_function_LCOH = zeros(n_clusters * 2);

    global j = 0

    for i in 1:n_clusters
        global j += 1
        step_function_LCOH[j] = step_matrix[i,3]
        j += 1
        step_function_LCOH[j] = step_matrix[i,3]
    end

    # println(step_function_LCOH)

    local step_supply_cost_coordinates = [step_function_prod step_function_LCOH]
    local header = ["Potential" "LCOH"]
    local step_supply_cost_coordinates_w_header = [header ; step_supply_cost_coordinates]

    local path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Final\\$scenario\\Step"
    local filename = ("$path\\$country Step Supply-Cost Curve Coordinate $scenario.txt")
    isfile(filename) && rm(filename)
    writedlm(filename,step_supply_cost_coordinates_w_header , '\t')

    local filename = ("$path\\$country Step Supply-Cost Curve $scenario")

    local plt = plot(x,y, title = "Supply Cost Curve $country", label="Raw data", xlabel="PJ")
    plot!(step_function_prod, step_function_LCOH, label="Clustered", ylabel="LCOH")

    savefig(filename)
end
