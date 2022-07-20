using DelimitedFiles, Dierckx, Plots#, LsqFit

country = "KOR"
scenario = "Come PLEXOS"
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

plot(x,y)

check_vector=zeros(length(y))

for i in 2:length(y)
    (y[i]>y[i-1]) && continue
    if  y[i]<y[i-1]
        check_vector[i]=1
    end
end
println("")
if sum(check_vector)==0
    println("\n...Raw Data is ordered...")
elseif sum(check_vector)!=0
    println("\n...Raw Data isn't ordered...")
end

fit1 = Spline1D(x,y,w=ones(length(x)), k=1, bc="nearest", s=0.0); #LCOH in funzione di PJ

max_num_step = 1000;
find_injective_vector = zeros(length(1:max_num_step))

for i in 2:max_num_step
    local xfit1 = collect(range(0; length= i, stop=last(x)))
    local yfit1 = fit1(xfit1)
    find_injective_vector[i] = (length(yfit1)==length(unique(yfit1)))
end

global step_fit1 = findlast(isequal(1),find_injective_vector)

println("\nMax number of steps to have an injective fit1 : $step_fit1")

xfit1 = collect(range(0; length= step_fit1, stop=last(x)));
yfit1 = fit1(xfit1);

plot(x,y)
scatter!(xfit1,yfit1)

check_injcetive = length(yfit1)==length(unique(yfit1))
println("\nIs fit1 injective? -> $check_injcetive")

# max_LCOH = 3
# for i in 1:length(yfit1)
#     if yfit1[i] >= max_LCOH
#         xfit1[i] = typemax(Int32);
#         yfit1[i] = typemax(Int32);
#     end
# end
# for i in 1:length(y)
#     if y[i] >= max_LCOH
#         x[i] = typemax(Int32);
#         y[i] = typemax(Int32);
#     end
# end
# filter!(!isequal(typemax(Int32)), xfit1)
# filter!(!isequal(typemax(Int32)), yfit1)
# filter!(!isequal(typemax(Int32)), x)
# filter!(!isequal(typemax(Int32)), y)

prod_funct_of_LCOH = [yfit1 xfit1]
prod_funct_of_LCOH_sorted = sortslices(prod_funct_of_LCOH, dims=1);
yfit1 = prod_funct_of_LCOH_sorted[:,1];
xfit1 = prod_funct_of_LCOH_sorted[:,2];

check_vector=zeros(length(yfit1))
for i in 2:length(yfit1)
    (yfit1[i]>yfit1[i-1]) && continue
    if  yfit1[i]<yfit1[i-1]
        check_vector[i]=1
    end
end
println("")
if sum(check_vector)==0
    println("\n...fit1 Data is ordered...")
elseif sum(check_vector)!=0
    println("\n...fit1 Data isn't ordered...")
end

fit2 = Spline1D(yfit1,xfit1;w=ones(length(yfit1)), k=1, bc="nearest", s=0.0); #PJ in funzione di LCOH

step = 0.01
lower_lim = round(minimum(yfit1), digits=2);
upper_lim = round(maximum(yfit1), digits=2);
# lower_lim = minimum(yfit1);
# upper_lim = maximum(yfit1);

println("\nLower LCOH limit : $lower_lim")
println("\nUpper LCOH limit : $upper_lim")

yfit2 = collect(range(lower_lim; step=step, stop=upper_lim));
xfit2 = fit2(yfit2);

check_vector=zeros(length(xfit2))
for i in 2:length(xfit2)
    (xfit2[i]>xfit2[i-1]) && continue
    if  xfit2[i]<xfit2[i-1]
        check_vector[i]=1
    end
end
println("")
if sum(check_vector)==0
    println("\n...fit2 Data is ordered...")
elseif sum(check_vector)!=0
    println("\n...fit2 Data isn't ordered...")
end

println("\n...calculating LCOH averages...")
global step_num = 20;
Δ_LCOH = step*step_num;
LCOH_avg_vector_length = trunc(Int, ceil((upper_lim-lower_lim)/Δ_LCOH) - 1);
LCOH_avg_vector = zeros(length(1:LCOH_avg_vector_length))
global j=1
for i in 1:LCOH_avg_vector_length
    global j
    LCOH_avg_vector[i] = (yfit2[j] + yfit2[j+step_num]) / 2;
    j += step_num
end

LCOH_step_vector = collect(range(lower_lim; step=Δ_LCOH, stop=upper_lim));
PROD_step_vector = fit2(LCOH_step_vector)

println("\n...creating step matrix...")
lower_lim_pot_vec = zeros(length(LCOH_avg_vector))
for i in 1:length(lower_lim_pot_vec)
    lower_lim_pot_vec[i] = PROD_step_vector[i]
end
upper_lim_pot_vec = zeros(length(LCOH_avg_vector))
for i in 1:length(lower_lim_pot_vec)
    upper_lim_pot_vec[i] = PROD_step_vector[i+1]
end
pot_step_vec = zeros(length(LCOH_avg_vector))
for i in 1:length(lower_lim_pot_vec)
    pot_step_vec[i] = upper_lim_pot_vec[i] - lower_lim_pot_vec[i]
end

step_matrix = [LCOH_avg_vector lower_lim_pot_vec upper_lim_pot_vec pot_step_vec]

region_vector = Vector{String}(undef, length(LCOH_avg_vector))
for i in 1:length(LCOH_avg_vector)
    region_vector[i] = "$country"
end
LCOH_class_vector = Vector{String}(undef, length(LCOH_avg_vector))
for i in 1:length(LCOH_avg_vector)
    LCOH_class_vector[i] = "$i"
end

header = ["Region" "Class" "LCOH" "Lower Lim Potential" "Upper Lim Potential" "Potential Step"]
step_matrix_w_region_class = [region_vector LCOH_class_vector step_matrix]
step_matrix_w_region_class_header = [header ; step_matrix_w_region_class]
println("\n...saving step matrix...")
path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Final\\$scenario\\Step"
filename = ("$path\\$country Step SC SPLINE $scenario.txt")
isfile(filename) && rm(filename)
writedlm(filename,step_matrix_w_region_class_header , '\t')

println("\n...creating step function coordinates matrix...")
global prod_coord = zeros(length(LCOH_avg_vector) * 2);
global j = 0
for i in 1:length(LCOH_avg_vector)
    global j += 1
    if i==1
        prod_coord[i] = step_matrix[i,2]
        continue
    elseif i==length(LCOH_avg_vector)
        prod_coord[i*2] = step_matrix[i,3]
    end
    prod_coord[j] = step_matrix[i-1,3]
    j += 1
    prod_coord[j] = step_matrix[i-1,3]
end
# println(prod_coord)
global lcoh_coord = zeros(length(LCOH_avg_vector) * 2);
global j = 0
for i in 1:length(LCOH_avg_vector)
    global j += 1
    lcoh_coord[j] = step_matrix[i,1]
    j += 1
    lcoh_coord[j] = step_matrix[i,1]
end
println("\n...saving step function coordinates matrix...")
step_supply_cost_coordinates = [prod_coord lcoh_coord]
header = ["Potential" "LCOH"]
step_supply_cost_coordinates_w_header = [header ; step_supply_cost_coordinates]
path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Final\\$scenario\\Step"
filename = ("$path\\$country Step SC Coord SPLINE Coordinate $scenario.txt")
isfile(filename) && rm(filename)
writedlm(filename,step_supply_cost_coordinates_w_header , '\t')

println("\n...creating a saving plot...")
filename = ("$path\\$country Step Supply-Cost Curve $scenario.pdf")
plt = plot(x,y, title = "Step SC SPLINE $country", label="Raw data", xlabel="PJ",size=(1900,1080))
plot!(prod_coord, lcoh_coord, label="Clustered", ylabel="LCOH", size=(1900,1080))
savefig(filename)
