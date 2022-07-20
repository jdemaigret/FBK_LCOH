path = homedir()
excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
region_name = excel_file["Region Name"]
inputs = excel_file["Inputs"]
PV_power_dens = inputs[:][1,2]
PV_plant_area = inputs[:][2,2]
Onshore_power_dens = inputs[:][3,2]
Onshore_plant_area = inputs[:][4,2]
Offshore_power_dens = inputs[:][5,2]
Offshore_plant_area = inputs[:][6,2]
exclusion_zones = excel_file["Exclusion Zones"]
max_popdens_PV = exclusion_zones[:][2,2]
max_popdens_Onshore = exclusion_zones[:][3,2]
landtype_exclusion_PV_input = exclusion_zones[:][6:22,2]
landtype_exclusion_Onshore_input = exclusion_zones[:][6:22,3]
landtype_exclusion_PV = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
landtype_exclusion_Onshore = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
for i in 1:length(landtype_exclusion_PV)
    if  isequal(landtype_exclusion_PV_input[i],missing)
        landtype_exclusion_PV[i]=typemax(Int8)
    end
    if  isequal(landtype_exclusion_Onshore_input[i],missing)
        landtype_exclusion_Onshore[i]=typemax(Int8)
    end
end
protected_exclusion_PV_input = exclusion_zones[:][25:34,2]
protected_exclusion_Onshore_input = exclusion_zones[:][25:34,3]
protected_exclusion_PV = [1,2,3,4,5,6,7,8,9,10]
protected_exclusion_Onshore = [1,2,3,4,5,6,7,8,9,10]
for i in 1:length(protected_exclusion_PV)
    if  isequal(protected_exclusion_PV_input[i],missing)
        protected_exclusion_PV[i]=typemax(Int8)
    end
    if  isequal(protected_exclusion_Onshore_input[i],missing)
        protected_exclusion_Onshore[i]=typemax(Int8)
    end
end
filter(!isequal(typemax(Int8)),landtype_exclusion_PV)
filter(!isequal(typemax(Int8)),landtype_exclusion_Onshore)
filter(!isequal(typemax(Int8)),protected_exclusion_PV)
filter(!isequal(typemax(Int8)),protected_exclusion_Onshore)

max_depth_offshore = exclusion_zones[:][44,2]
min_shore_distance_offshore = exclusion_zones[:][45,2]

path = homedir()
excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
exclusion_zones = excel_file["Exclusion Zones"]
max_slope_pv = exclusion_zones[:][37,2]
max_slope_onshore = exclusion_zones[:][37,3]


path = homedir()
excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
exclusion_zones = excel_file["Exclusion Zones"]
BWS_threshold = exclusion_zones[:][41,2]

path = homedir()
excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
exclusion_zones = excel_file["Exclusion Zones"]
min_desalination_distance = exclusion_zones[:][48,2]

path = homedir()
excel_file = XLSX.readxlsx("$path\\save_regions.xlsx")
exclusion_zones = excel_file["Exclusion Zones"]
include_desal = exclusion_zones[:][40,2]
