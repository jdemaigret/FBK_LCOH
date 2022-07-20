using DelimitedFiles

nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
# nomi_stringa = ["ARG"]
global sc_mat = zeros(1000,68)
global scenario = "2030 Pess ws"
global i=1
global j=2
@time for country in nomi_stringa
    global sc_data = readdlm("E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\Interpolate\\$country Supply-Cost Curve Interpolate $scenario.txt")
    sc_mat[:,i]=sc_data[:,1]
    sc_mat[:,j]=sc_data[:,2]
    global i+=2
    global j+=2
end

filename_sc_mat = "SC Mat\\$scenario SC.txt"
header = ["ARG" "" "AUS" "" "BRA" "" "CAN" "" "CHI" "" "CHL" "" "COL" "" "EAS" "" "EUR" "" "FRA" "" "GER" "" "IDN" "" "IND" "" "ITA" "" "JAP" "" "KOR" "" "LAM" "" "MEN" "" "MEX" "" "MOR" "" "OCE" "" "POR" "" "ROA" "" "ROE" "" "RUS" "" "SAU" "" "SEA" "" "SOU" "" "SPA" "" "SSA" "" "TUR" "" "UKG" "" "UKR" "" "USA" ""];
isfile(filename_sc_mat) && rm(filename_sc_mat);
writedlm(filename_sc_mat, [header; sc_mat], '\t');
