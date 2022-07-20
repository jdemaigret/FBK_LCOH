using DelimitedFiles, GeoArrays

global scenario = "2050 V ws"
# global scenario = "2030 Opt ws"
mat2read_names = ["onshore_matrix" "pv_matrix"]
global read_path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\sys_nature\\Cap_matrices\\$scenario"

global region_names = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
global abs_matrix = zeros(36000,18000);

for mat2read in mat2read_names
    @time for region in region_names
        println("Region: $region")
        println("mat2read: $mat2read")
        println("$scenario")
        local rel_mat = readdlm("$read_path\\$region $mat2read ($scenario).txt" ,',')
        local regions,_,_,lonrange,latrange=loadregions(region)
        # Removing Offshore regions from rel_map, might only be confusing in final world map
        # This has to be excluded in case Offshore potential map is processed
        for i in 1:length(lonrange)
            for j in 1:length(latrange)
                if rel_mat[i,j] != 0
                    local abs_matrix[lonrange[i], latrange[j]]=rel_mat[i,j]
                end
            end
        end
    end

    local save_path = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario";
    # save_path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\Report\\Africa BOX\\$scenario";
    local ga = GeoArray(abs_matrix)
    bbox!(ga, (min_x=-180, min_y=90, max_x=180, max_y=-90))
    epsg!(ga, 4326)
    println("Saving GeoTiff...")
    GeoArrays.write!("$save_path\\$scenario $mat2read World Matrix.tif", ga)
end
