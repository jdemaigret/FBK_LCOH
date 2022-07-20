using DelimitedFiles, Clustering,  Plots

# nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "POR" "ROA" "ROE" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
nomi_stringa = ["CHI" "CHL" "GER" "IDN" "ITA" "JAP"]
country = "IDN"
scenario = "2050 V"
# scenario = "Come PLEXOS con BWS"

supply_cost_curve_data = readdlm("C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\File di testo\\$country results_matrix ($scenario).txt", skipstart=1)

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

supply_cost_curve_data_clean = [ LCOH_clean prod_kg_clean pv_clean onshore_clean offshore_clean electrolyzer_clean FLOH_clean sum_cf_pv_clean sum_cf_onshore_clean sum_cf_offshore_clean curtailment_clean];

supply_cost_curve_data_clean_sorted = sortslices(supply_cost_curve_data_clean, dims=1);

N_clean_sorted = size(supply_cost_curve_data_clean_sorted)[1];

prod_PJ         = supply_cost_curve_data_clean_sorted[:,2] .* 120 ./ 10^9;
prog_prod_PJ    = zeros(length(prod_PJ));
prog_prod_PJ[1] = prod_PJ[1];

for i in 2:N_clean_sorted
    prog_prod_PJ[i] = prog_prod_PJ[i-1] + prod_PJ[i];
end

supply_cost_curve_data_clean_sorted_w_prog = hcat(supply_cost_curve_data_clean_sorted,prog_prod_PJ);

pv_onshore_capacities = zeros(length(1:N_clean_sorted),2); #1 PV 2 ONSHORE
pv_capacities         = zeros(length(1:N_clean_sorted));
onshore_capacities    = zeros(length(1:N_clean_sorted));
offshore_capacities   = zeros(length(1:N_clean_sorted));
pv_onshore_cf_sum     = zeros(length(1:N_clean_sorted),2);
pv_cf_sum             = zeros(length(1:N_clean_sorted));
onshore_cf_sum        = zeros(length(1:N_clean_sorted));
offshore_cf_sum       = zeros(length(1:N_clean_sorted));


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
# i quali sarebbero minimo e massimo di x. I valori che può assumere bc sono:
# nearest, zero,extrapolate, error. Secondo me è da lasciare nearest.
# s ivnece indica lo smoothing factor.
# posso fare la derivata di grado nu-esimo nel punto x con
# derivative(spl, x; nu=1). Posso fare l'integrale da x=a a x=b con
# integrate(spl, a, b).Posso trovare gli zeri con roots(spl; maxn=8). maxn mi
# da il numero massimo di radici da ritornare.
#  penso che w siano i weight che vengono asegnati.

x = supply_cost_curve_data_clean_sorted_w_prog[:,12];
y = supply_cost_curve_data_clean_sorted_w_prog[:,1];


max_LCOH = 4

for i in 1:length(x)
    if y[i]>max_LCOH
        y[i]=typemax(Int32)
        x[i]=typemax(Int32)
    end
end

filter!(!isequal(typemax(Int32)), y);
filter!(!isequal(typemax(Int32)), x);

points = [transpose(x); transpose(y)]

n_clusters = 5

R = kmeans(points, n_clusters; maxiter=200, display=:iter);

@assert nclusters(R) == n_clusters # verify the number of clusters

a = assignments(R) # get the assignments of points to clusters (Array{Int64,1})
c = counts(R) # get the cluster sizes (Array{Int64,1})
M = R.centers # get the cluster centers (Array{Float64,2})

centers = [M[2,:] M[1,:] zeros(size(M)[2]) zeros(size(M)[2]) c];

centers_sorted = sortslices(centers, dims=1)

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

step_matrix = zeros(size(centers_sorted))
step_matrix[:,1:2] = centers_sorted[:,1:2]

for i in 1:n_clusters
    step_matrix[i,3] = x[trunc(Int, centers_sorted[i,3])]
    step_matrix[i,4] = x[trunc(Int, centers_sorted[i,4])]
    step_matrix[i,5] = x[trunc(Int, centers_sorted[i,4])] - x[trunc(Int, centers_sorted[i,3])]
end


region_vector = Vector{String}(undef, n_clusters)

for i in 1:n_clusters
    region_vector[i] = "$country"
end

LCOH_class_vector = Vector{String}(undef, n_clusters)

for i in 1:n_clusters
    LCOH_class_vector[i] = "$i"
end
header = ["Region" "Class" "LCOH (<$max_LCOH)" "Center Potential" "Lower Lim Potential" "Upper Lim Potential" "Potential Step"]
step_matrix = [region_vector LCOH_class_vector step_matrix]
step_matrix_w_header = [header ; step_matrix]

path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\K-Means Cluster"
filename = ("$path\\$country Step Supply-Cost Curve $scenario.txt")
isfile(filename) && rm(filename)
writedlm(filename,step_matrix_w_header , '\t')

step_function_prod = zeros(n_clusters * 2);

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

println(step_function_prod)

step_function_LCOH = zeros(n_clusters * 2);

global j = 0

for i in 1:n_clusters
    global j += 1
    step_function_LCOH[j] = step_matrix[i,3]
    j += 1
    step_function_LCOH[j] = step_matrix[i,3]
end

println(step_function_LCOH)

step_supply_cost_coordinates = [step_function_prod step_function_LCOH]
header = ["Potential" "LCOH"]
step_supply_cost_coordinates_w_header = [header ; step_supply_cost_coordinates]

path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\K-Means Cluster"
filename = ("$path\\$country Step Supply-Cost Curve Coordinate $scenario.txt")
isfile(filename) && rm(filename)
writedlm(filename,step_supply_cost_coordinates_w_header , '\t')

filename = ("$path\\$country Step Supply-Cost Curve $scenario")

plt = plot(x,y, title = "Supply Cost Curve $country", label="Raw data", xlabel="PJ")
plot!(step_function_prod, step_function_LCOH, label="Clustered", ylabel="LCOH")


savefig(filename)
