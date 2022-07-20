nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
nomi_stringa_ridotto = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "SAU" "SEA" "SOU" "SPA" "TUR" "UKG" "UKR"]

nomi_stringa = ["SAU" "SEA" "SOU" "SPA" "TUR" "UKG" "UKR" "USA"]

nomi = [ARG, AUS, BRA, CAN, CHI, CHL, COL, EAS, EUR, FRA,
        GER, IDN, IND, ITA, JAP, KOR, LAM, MEN, MEX, MOR,
        OCE, POR, ROA, ROE, RUS, SAU, SEA, SOU, SPA, SSA,
        TUR, UKG, UKR, USA]

nomi_string_big = ["RUS" "SSA"]

for i in 1:length(nomi)
    saveregions(nomi_stringa[i], nomi[i])
end

@time for i in 1:length(nomi_stringa)
    GISlcoh(gisregion=nomi_stringa[i])
end

for i in 30:34
    GISlcoh(gisregion=nomi_stringa[i])
end

for i in 1:length(nomi_stringa_ridotto)
    GISlcoh(gisregion=nomi_stringa_ridotto[i])
end

global LCOE_matrix = zeros(34,15)
for i in 1:length(nomi_stringa)
    gisregion = nomi_stringa[i]
    # LCOE_region = readdlm("C:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Characetristic CF Curves\\Version 8\\Text Files\\LCOEs V8\\$gisregion LCOEs V8.txt",',', skipstart=1)
    LCOE_region = readdlm("LCOEs V11a\\$gisregion LCOEs last.txt",',', skipstart=1)
    for j in 1:15
        LCOE_matrix[i,j] = LCOE_region[1,j]
    end
end

header = ["LCOE [USD/MWh]" "PV 1" "PV 2" "PV 3" "PV 4" "PV 5" "Onshore 1" "Onshore 2" "Onshore 3" "Onshore 4" "Onshore 5" "Offshore 1" "Offshore 2" "Offshore 3" "Offshore 4" "Offshore 5"]
nomi_stringa_colonna = ["ARG", "AUS", "BRA", "CAN", "CHI", "CHL", "COL", "EAS", "EUR", "FRA", "GER", "IDN", "IND", "ITA", "JAP", "KOR", "LAM", "MEN", "MEX", "MOR", "OCE", "POR", "ROA", "ROE", "RUS", "SAU", "SEA", "SOU", "SPA", "SSA", "TUR", "UKG", "UKR", "USA"]

LCOE_con_region = [nomi_stringa_colonna LCOE_matrix]

filename = GlobalEnergyGIS.in_datafolder("output", "LCOEs V11a.txt")
isfile(filename) && rm(filename);
writedlm(filename, [header ; LCOE_con_region], ',')




# SOLO INDONESIA

nomi_stringa_indonesia = ["SM_1" "SM_2" "SM_3" "Kal_1" "Kal_2" "Kal_3" "Kal_4" "JV_1" "JV_2" "JV_3" "LSI_1" "LSI_2" "LSI_3" "SW_1" "SW_2" "SW_3" "WNG" "MLK"]
nomi_indonesia = [SM_1, SM_2, SM_3, Kal_1, Kal_2, Kal_3, Kal_4, JV_1, JV_2, JV_3, LSI_1, LSI_2, LSI_3, SW_1, SW_2, SW_3, WNG, MLK]

for i in 1:length(nomi_indonesia)
    saveregions(nomi_stringa_indonesia[i], nomi_indonesia[i])
end

@time for i in 1:length(nomi_stringa_indonesia)
    GISlcoh(gisregion=nomi_stringa_indonesia[i])
end

nomi_stringa = ["ARG" "AUS" "BRA" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "SAU" "SEA" "SOU" "SPA" "TUR" "UKG" "UKR"]


nomi_stringa_grandi = ["CAN" "RUS" "SSA"]
nomi_stringa=[ "TUR" "UKG" "UKR"]
for i in 1:length(nomi_stringa)
    GISlcoh(gisregion=nomi_stringa[i])
end


nomi_stringa=[ "FRA" "GER" "IDN"]
for i in 1:length(nomi_stringa)
    GISlcoh(gisregion=nomi_stringa[i])
end


# Malaysia
nomi_stringa_malesia = ["PEN_1", "PEN_2", "PEN_3", "PEN_4", "PEN_5", "PEN_6", "PEN_7", "PEN_8", "PEN_9", "PEN_10", "PEN_11", "PEN_12", "PEN_13", "SAB_1", "SAB_2", "SAR_1", "SAR_2", "SAR_3"]
nomi_malesia = [PEN_1, PEN_2, PEN_3, PEN_4, PEN_5, PEN_6, PEN_7, PEN_8, PEN_9, PEN_10, PEN_11, PEN_12, PEN_13, SAB_1, SAB_2, SAR_1, SAR_2, SAR_3]


for i in 1:length(nomi_malesia)
    saveregions(nomi_stringa_malesia[i], nomi_malesia[i])
end

for i in 1:length(nomi_stringa_malesia)
    GISlcoh(gisregion=nomi_stringa_malesia[i])
end




for i in 1:length(nomi_stringa_africa)
    GISlcoh(gisregion=nomi_stringa_africa[i])
end


nomi_stringa = ["COD" "EGY" "ETH" "KEN" "MRT" "NAM"]

nomi = [COD EGY ETH KEN MRT NAM]

for i in 1:length(nomi)
    saveregions(nomi_stringa[i], nomi[i])
end

for i in 1:length(nomi_stringa)
    GISlcoh(gisregion=nomi_stringa[i])
end


nomi_stringa = ["Brunei" "Cambodia" "Laos" "Malaysia" "Myanmar" "Philippines" "Singapore" "Thailand" "Vietnam"]
nomi = [brunei, cambodia, laos, malaysia, myanmar, philippines, singapore, thailand, vietnam]
for i in 1:length(nomi)
    saveregions(nomi_stringa[i], nomi[i])
end
for i in 1:length(nomi_stringa)
    GISlcoh(gisregion=nomi_stringa[i])
end
