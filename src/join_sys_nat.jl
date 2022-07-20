using DelimitedFiles
nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
scenario = "2050 V ws"
read_path = "C:\\Users\\jdemaigret\\.julia\\dev\\GlobalEnergyGIS_lcoh\\src\\System Nature\\$scenario"
global sys_nat_4_excel = zeros(6,34)
global i=1
@time for country in nomi_stringa
    local sys_nat_data = readdlm("$read_path\\$country system nature areas $scenario.txt",'\t')
    global sys_nat_4_excel[:,i] = sys_nat_data[:,2]
    global i+=1
end
global write_path = "C:\\Users\\jdemaigret\\.julia\\dev\\GlobalEnergyGIS_lcoh\\src\\System Nature\\Complete matrices"
global filename_sys_nat_4_excel = "$write_path\\System nature all regions $scenario.txt";
isfile(filename_sys_nat_4_excel) && rm(filename_sys_nat_4_excel);
writedlm(filename_sys_nat_4_excel, sys_nat_4_excel, '\t');
