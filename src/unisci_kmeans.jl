using DelimitedFiles

nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
# nomi_stringa = ["ARG"]

global kmeans_mat = zeros(10,68)
global scenario = "2030 Pess ws"
global s=1
global t=2
@time for country in nomi_stringa
    global kmeans_data = readdlm("E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\K-Means Cluster\\$country Step Supply-Cost Curve Coordinate $scenario.txt", skipstart=1)
    kmeans_mat[:,s]=kmeans_data[:,1]
    kmeans_mat[:,t]=kmeans_data[:,2]
    global s+=2
    global t+=2
end

filename_kmeans_mat = "SC Mat\\$scenario kmeans.txt"
header = ["ARG" "" "AUS" "" "BRA" "" "CAN" "" "CHI" "" "CHL" "" "COL" "" "EAS" "" "EUR" "" "FRA" "" "GER" "" "IDN" "" "IND" "" "ITA" "" "JAP" "" "KOR" "" "LAM" "" "MEN" "" "MEX" "" "MOR" "" "OCE" "" "POR" "" "ROA" "" "ROE" "" "RUS" "" "SAU" "" "SEA" "" "SOU" "" "SPA" "" "SSA" "" "TUR" "" "UKG" "" "UKR" "" "USA" ""];
header1 = ["Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH"]
isfile(filename_kmeans_mat) && rm(filename_kmeans_mat);
writedlm(filename_kmeans_mat,[ header;[header1; kmeans_mat]], '\t');
