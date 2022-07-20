using DelimitedFiles

nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
# nomi_stringa = ["ARG"]

global spline_mat = zeros(20,68)
global scenario = "2030 Pess ws"
global a=1
global b=2
@time for country in nomi_stringa
    global spline_data = readdlm("E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\Spline Cluster\\$country Step SC Coord SPLINE Coordinate $scenario.txt", skipstart=1)
    spline_mat[1:size(spline_data)[1],a]=spline_data[:,1]
    spline_mat[1:size(spline_data)[1],b]=spline_data[:,2]
    global a+=2
    global b+=2
end

filename_spline_mat = "SC Mat\\$scenario spline.txt"
header = ["ARG" "" "AUS" "" "BRA" "" "CAN" "" "CHI" "" "CHL" "" "COL" "" "EAS" "" "EUR" "" "FRA" "" "GER" "" "IDN" "" "IND" "" "ITA" "" "JAP" "" "KOR" "" "LAM" "" "MEN" "" "MEX" "" "MOR" "" "OCE" "" "POR" "" "ROA" "" "ROE" "" "RUS" "" "SAU" "" "SEA" "" "SOU" "" "SPA" "" "SSA" "" "TUR" "" "UKG" "" "UKR" "" "USA" ""];
header1 = ["Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH" "Potential" "LCOH"]
isfile(filename_spline_mat) && rm(filename_spline_mat);
writedlm(filename_spline_mat,[ header;[header1; spline_mat]], '\t');
