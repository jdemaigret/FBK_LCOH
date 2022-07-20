using DelimitedFiles, Dierckx

nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "POR" "ROA" "ROE" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR"]
#                 1     2     3     4     5     6     7     8     9     10    11    12   13    14    15    16     17    18   19     20    21    22    23    24    25   26    27    28     29   30    31
scenario = "(come PLEXOS)"

for reg in nomi_stringa

    println("reg")
    println(reg)

    local country = reg

    local supply_cost_curve_data = readdlm("C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Da mostrare\\V13\\File di testo\\$country results_matrix (come PLEXOS).txt", skipstart=1);

    if reg == "ARG"

        global supply_cost_curve_data_joined = zeros(size(supply_cost_curve_data)[1],2);
        local LCOH    = supply_cost_curve_data[:,1];
        local prod_kg = supply_cost_curve_data[:,2];

        local N = length(LCOH);

        local deleted_index = typemax(Int32);

        for i in 1:N
            if LCOH[i] == Inf
                LCOH[i]            = deleted_index;
                prod_kg[i]         = deleted_index;
            end
        end

        local LCOH_clean            = zeros(length(LCOH));
        local prod_kg_clean         = zeros(length(prod_kg));

        local LCOH_clean    = filter(!isequal(deleted_index), LCOH);
        local prod_kg_clean = filter(!isequal(deleted_index), prod_kg);

        local supply_cost_curve_data_clean = [LCOH_clean prod_kg_clean];

        println("\nsize(supply_cost_curve_data_clean)")
        println(size(supply_cost_curve_data_clean))

        global supply_cost_curve_data_joined = supply_cost_curve_data_clean;

        local prod_EJ = supply_cost_curve_data_joined[:,2] .* 120 ./ 10^12;
        println("sum(prod_EJ)")
        println(sum(prod_EJ))

    else

        local LCOH    = supply_cost_curve_data[:,1];
        local prod_kg = supply_cost_curve_data[:,2];

        local N = length(LCOH);

        local deleted_index = typemax(Int32);

        for i in 1:N
            if LCOH[i] == Inf
                LCOH[i]            = deleted_index;
                prod_kg[i]         = deleted_index;
            end
        end

        local LCOH_clean            = zeros(length(LCOH));
        local prod_kg_clean         = zeros(length(prod_kg));

        local LCOH_clean    = filter(!isequal(deleted_index), LCOH);
        local prod_kg_clean = filter(!isequal(deleted_index), prod_kg);

        local supply_cost_curve_data_clean = [LCOH_clean prod_kg_clean];

        println("\nsize(supply_cost_curve_data_clean)")
        println(size(supply_cost_curve_data_clean))

        local prod_EJ = supply_cost_curve_data_clean[:,2] .* 120 ./ 10^12;
        println("sum(prod_EJ)")
        println(sum(prod_EJ))

        global supply_cost_curve_data_joined = [supply_cost_curve_data_joined ; supply_cost_curve_data_clean]
        #vcat(supply_cost_curve_data_joined, supply_cost_curve_data_clean);

    end

    println("\nsize(supply_cost_curve_data_joined) fine ciclo for arrivati fino a $country")
    println(size(supply_cost_curve_data_joined))
    prod_EJ = supply_cost_curve_data_joined[:,2] .* 120 ./ 10^12;
    println("sum(prod_EJ)  fine ciclo for arrivati fino a $country")
    println(sum(prod_EJ))


end

println("\nsize(supply_cost_curve_data_joined) fuori dal for")
println(size(supply_cost_curve_data_joined))

supply_cost_curve_data_joined_sorted = sortslices(supply_cost_curve_data_joined, dims=1);

println("\nsize(supply_cost_curve_data_joined_sorted) fuori dal for")
println(size(supply_cost_curve_data_joined_sorted))

prod_EJ = supply_cost_curve_data_joined_sorted[:,2] .* 120 ./ 10^12;
println("sum(prod_EJ) fuori dal for")
println(sum(prod_EJ))

# println("\nsize(supply_cost_curve_data_clean_sorted)")
# println(size(supply_cost_curve_data_clean_sorted))

N_clean_sorted = size(supply_cost_curve_data_joined_sorted)[1];

prod_EJ         = supply_cost_curve_data_joined_sorted[:,2] .* 120 ./ 10^12;
# println("sum(prod_EJ)")
# println(sum(prod_EJ))
prog_prod_EJ    = zeros(length(prod_EJ));
prog_prod_EJ[1] = prod_EJ[1];

for i in 2:N_clean_sorted
    prog_prod_EJ[i] = prog_prod_EJ[i-1] + prod_EJ[i];
end

final_join_supply_cost_curve = hcat(supply_cost_curve_data_joined_sorted,prog_prod_EJ);

println("\nsize(final_join_supply_cost_curve) fuori dal for")
println(size(final_join_supply_cost_curve))

println("sum(final_join_supply_cost_curve[:,2]) fuori dal for")
println(sum(final_join_supply_cost_curve[:,2]))

println("last(final_join_supply_cost_curve[:,3]) fuori dal for")
println(last(final_join_supply_cost_curve[:,3]))

# println("\nsize(final_join_supply_cost_curve)")
# println(size(final_join_supply_cost_curve))

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

x = final_join_supply_cost_curve[:,3];
y = final_join_supply_cost_curve[:,1];

fit1 = Spline1D(x,y,w=ones(length(x)), k=1, bc="nearest", s=0.0); #LCOH in funzione di PJ

xs = range(0; length= 2000, stop=last(x))
ys = fit1(xs)

path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Da mostrare\\V13\\Supply Cost Curves Interpolate"
filename = ("$path\\Joined Supply-Cost Curve Interpolate $scenario.txt")
isfile(filename) && rm(filename)
writedlm(filename, [xs  ys], '\t')
