using DelimitedFiles, GeoArrays

global scenario = "2050 V ws"
# global scenario = "2030 Opt ws"

# global mat2read = "lcoh"

# global mat2read = "Ori PV"
# global mat2read = "Ori Onsh"
# global mat2read = "Hy to PV"
# global mat2read = "Hy PV Prev"
# global mat2read = "Hy to Onsh"
# global mat2read = "Hy Onsh Prev"

global mat2read_names = ["Ori PV" "Ori Onsh" "Hy to PV" "Hy PV Prev" "Hy to Onsh" "Hy Onsh Prev"]
# global mat2read_names = ["Ori PV" "Ori Onsh" "Hy to PV" "Hy PV Prev" "Hy Onsh Prev"]
# global mat2read_names = ["Hy to Onsh"]

global lcoh_threshold = 5
global read_path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\sys_nature\\Sys_nat_matrices\\$scenario"
# global read_path = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario\\$mat2read"
# global region_names = ["Sardegna"]
# global region_names = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
global region_names = ["CAN"]
# global region_names = ["COD" "EGY" "ETH" "KEN" "MOR" "MRT" "NAM"]
# global region_names = ["KEN"]
global abs_matrix = zeros(36000,18000);

# global regions_frame,_,_,lonrange_frame,latrange_frame=loadregions("africa_box_region")

# for i in 1:size(regions_frame)[1]
#     for j in 1:size(regions_frame)[2]
#         if (regions_frame[i,j]!=0) && (regions_frame[i,j]!=32767)
#             regions_frame[i,j]=1
#         end
#     end
# end
#
# for i in 1:length(lonrange_frame)
#     for j in 1:length(latrange_frame)
#         if regions_frame[i,j] == 1
#             global abs_matrix[lonrange_frame[i], latrange_frame[j]]=0
#         end
#     end
# end

# regions,_,_,lonrange,latrange=loadregions(region)

# println("size(regions)")
# println(size(regions))
#
# println("typeof(regions)")
# println(typeof(regions))
#
# println("length(unique(regions))")
# println(length(unique(regions)))
#
# if length(unique(regions))<10
#     println("unique(regions)")
#     println(unique(regions))
# end
for mat2read in mat2read_names
    @time for region in region_names
        println("Region: $region")
        println("mat2read: $mat2read")
        println("$scenario")

        # local rel_mat = readdlm("$read_path\\$region $(mat2read)_matrix ($scenario).txt" ,',')
        # global rel_mat = readdlm("C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\Report\\Africa BOX\\$scenario\\$region\\$region lcoh_matrix ($scenario).txt" ,',')
        local rel_mat = readdlm("$read_path\\$region $mat2read $scenario.txt" ,',')
        local regions,_,_,lonrange,latrange=loadregions(region)
        println("maximum(rel_mat)")
        println(maximum(rel_mat))
        for i in 1:size(rel_mat)[1]
            for j in 1:size(rel_mat)[2]
                if (rel_mat[i,j] > lcoh_threshold) || (rel_mat[i,j] == NaN)
                    local rel_mat[i,j] = 100
                end
            end
        end
        # Removing Offshore regions from rel_map, might only be confusing in final world map
        # This has to be excluded in case Offshore potential map is processed
        for i in size(regions)[1]
            for j in size(regions)[2]
                if regions[i,j] == 0
                    local rel_mat[i,j] = 0
                end
                # if regions[i,j] == typemax(Int32)
                #     local rel_mat[i,j] = 200
                # end
            end
        end
        for i in 1:length(lonrange)
            for j in 1:length(latrange)
                if rel_mat[i,j] != 0
                    global abs_matrix[lonrange[i], latrange[j]] = rel_mat[i,j]
                end
            end
        end

    end
    #
    # for i in 1:length(lonrange_frame)
    #     for j in 1:length(latrange_frame)
    #         if regions_frame[i,j] == 1
    #             global abs_matrix[lonrange_frame[i], latrange_frame[j]]=0
    #         end
    #     end
    # end
    local save_path = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario";
    # save_path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\Report\\Africa BOX\\$scenario";
    local ga = GeoArray(abs_matrix)
    bbox!(ga, (min_x=-180, min_y=90, max_x=180, max_y=-90))
    epsg!(ga, 4326)
    println("Saving GeoTiff...")
    GeoArrays.write!("$save_path\\$scenario $mat2read World Matrix.tif", ga)
    # @time GeoArrays.write!("$save_path\\$scenario $mat2read Africa Matrix.tif", ga)
end
