using DelimitedFiles

nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
nomi_stringa = ["ARG"]
scenario = "2050 Pess ws"
global read_path = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\Interpolate"

global supply_cost_matrix = zeros(34000,2)

global LCOH_vec = Float64[]
global potential_vec = Float32[]

@time for reg in nomi_stringa

    println("reg")
    println(reg)
    global data = readdlm("$read_path\\$reg Supply-Cost Curve Interpolate $scenario.txt")
    local LCOH_data = data[:,1]
    local potential_data = data[:,2]

    push!(LCOH_vec,LCOH_data)
    push!(potential_vec,potential_data)

end
#
# header = ["ARG" "" "AUS" "" "BRA" "" "CAN" "" "CHI" "" "CHL" "" "COL" "" "EAS" "" "EUR" "" "FRA" "" "GER" "" "IDN" "" "IND" "" "ITA" "" "JAP" "" "KOR" "" "LAM" "" "MEN" "" "MEX" "" "MOR" "" "OCE" "" "POR" "" "ROA" "" "ROE" "" "RUS" "" "SAU" "" "SEA" "" "SOU" "" "SPA" "" "SSA" "" "TUR" "" "UKG" "" "UKR" "" "USA" ""]
# path = "C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\Da mostrare\\V13\\Supply Cost Curves Interpolate"
# filename = ("$path\\Complete Supply-Cost Curve Interpolate $scenario.txt")
# isfile(filename) && rm(filename)
# writedlm(filename, [header ; supply_cost_matrix], '\t')
