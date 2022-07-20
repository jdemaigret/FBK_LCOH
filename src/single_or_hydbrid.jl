using DelimitedFiles, Plots


nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
# nomi_stringa = ["ARG"]
@time for country in nomi_stringa

    println("Region: $country")

    local scenario = "2050 V ws"
    # local scenario = "2030 Opt ws"
    local supply_cost_curve_data = readdlm("E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\LCOH analysis (supply curves)\\Risultati\\Supply Cost Curves\\$scenario\\File di testo\\$country results_matrix ($scenario).txt", skipstart=1)
    local hybrid_threshold_PV = 0.01
    local hybrid_threshold_Onshore = 0.01
    local lcoh_threshold = 5

    local PV_density = 45000
    local Onshore_density = 5000

    local deleted_index = typemax(Int32)

    local lcoh_vec = supply_cost_curve_data[:,1]
    global sys_nature_mat = supply_cost_curve_data[:,3:4]

    local N = 1:size(sys_nature_mat)[1]

    for i in N
        if ((sys_nature_mat[i,1]==0) && (sys_nature_mat[i,2]==0)) || (lcoh_vec[i]>lcoh_threshold)
            global sys_nature_mat[i,1] = deleted_index
            global sys_nature_mat[i,2] = deleted_index
        end
    end

    local PV_vec = filter(!isequal(deleted_index), sys_nature_mat[:,1]);
    local Onshore_vec = filter(!isequal(deleted_index), sys_nature_mat[:,2]);

    global sys_nature_mat = [PV_vec Onshore_vec]
    local area_mat = zeros(size(sys_nature_mat))

    local N = 1:size(sys_nature_mat)[1]

    for i in N
        local area_mat[i,1] = sys_nature_mat[i,1] /  PV_density
        local area_mat[i,2] = sys_nature_mat[i,2] /  Onshore_density
    end

    local sys_nature_vec = zeros(size(sys_nature_mat)[1])


    # 1= Originally PV only
    # 2= Originally Onshore only
    # 3= PV prevalent hybrid
    # 4= Onshore prevalent hybrid
    # 5= Hybrid became PV only
    # 6= Hybrid became Onshore only

    global PV_area=0             # Originally PV only
    global Onshore_area=0        # Originally Onshore only
    global Hybrid_PV_area=0      # PV prevalent hybrid
    global Hybrid_Onshore_area=0 # Onshore prevalent hybrid
    global hyb2pv_area=0         # Hybrid became PV only
    global hyb2onsh_area=0       # Hybrid became Onshore only

    for i in N
        if sys_nature_mat[i,2] == 0 # Originally PV only
            local sys_nature_vec[i] = 1
            global PV_area+=area_mat[i,1]
        elseif sys_nature_mat[i,1] == 0 # Originally Onshore only
            local sys_nature_vec[i] = 2
            global Onshore_area+=area_mat[i,2]
        elseif (sys_nature_mat[i,1] > sys_nature_mat[i,2]) && (sys_nature_mat[i,2]>0)
            if (sys_nature_mat[i,2]/sys_nature_mat[i,1] <= hybrid_threshold_PV) # Hybrid became PV only
                local sys_nature_vec[i] = 5
                # global PV_area+=area_mat[i,1]
                global hyb2pv_area+=area_mat[i,1]
            else                                        # PV prevalent hybrid
                local sys_nature_vec[i] = 3
                global Hybrid_PV_area+=area_mat[i,1]
            end
        elseif (sys_nature_mat[i,2] > sys_nature_mat[i,1]) && (sys_nature_mat[i,1]>0)
            if (sys_nature_mat[i,1]/sys_nature_mat[i,2] <= hybrid_threshold_Onshore) # Hybrid became Onshore only
                local sys_nature_vec[i] = 6
                global hyb2onsh_area+=area_mat[i,2]
                # global Onshore_area+=area_mat[i,2]
            else                                        # Onshore prevalent hybrid
                local sys_nature_vec[i] = 4
                global Hybrid_Onshore_area+=area_mat[i,2]
            end
        end
    end

    # Quello di prima
    # for i in N
    #     if sys_nature_mat[i,2] == 0 # Originally PV only
    #         local sys_nature_vec[i] = 1
    #         global PV_area+=area_mat[i,1]
    #         continue
    #     elseif sys_nature_mat[i,1] == 0 # Originally Onshore only
    #         local sys_nature_vec[i] = 2
    #         global Onshore_area+=area_mat[i,2]
    #         continue
    #     end
    #     if sys_nature_mat[i,1] > sys_nature_mat[i,2]
    #         if (sys_nature_mat[i,2]/sys_nature_mat[i,1] <= hybrid_threshold_PV) # Hybrid became PV only
    #             local sys_nature_vec[i] = 5
    #             # global PV_area+=area_mat[i,1]
    #             global hyb2pv_area+=area_mat[i,1]
    #         else                                        # PV prevalent hybrid
    #             local sys_nature_vec[i] = 3
    #             global Hybrid_PV_area+=area_mat[i,1]
    #         end
    #     end
    #     if sys_nature_mat[i,2] > sys_nature_mat[i,1]
    #         if (sys_nature_mat[i,1]/sys_nature_mat[i,2] <= hybrid_threshold_Onshore) # Hybrid became Onshore only
    #             local sys_nature_vec[i] = 6
    #             global hyb2onsh_area+=area_mat[i,2]
    #             # global Onshore_area+=area_mat[i,2]
    #         else                                        # Onshore prevalent hybrid
    #             local sys_nature_vec[i] = 4
    #             global Hybrid_Onshore_area+=area_mat[i,2]
    #         end
    #     end
    # end

    global PV=0             # Originally PV only
    global Onshore=0        # Originally Onshore only
    global Hybrid_PV=0      # PV prevalent hybrid
    global Hybrid_Onshore=0 # Onshore prevalent hybrid
    global hyb2pv=0         # Hybrid became PV only
    global hyb2onsh=0       # Hybrid became Onshore only


    for i in 1:length(sys_nature_vec)
        if sys_nature_vec[i]==1
            global PV+=1
        elseif sys_nature_vec[i]==2
            global Onshore+=1
        elseif sys_nature_vec[i]==3
            global Hybrid_PV+=1
        elseif sys_nature_vec[i]==4
            global Hybrid_Onshore+=1
        elseif sys_nature_vec[i]==5
            global hyb2pv+=1
        elseif sys_nature_vec[i]==6
            global hyb2onsh+=1
        end
    end

    # if (Hybrid_PV!=0) && (Hybrid_Onshore!=0)
    #     local xs = ["PV" ,"Onshore", "Hybrid PV" ,"Hybrid Onshore" ]
    #     local ys = [PV, Onshore, Hybrid_PV,  Hybrid_Onshore]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # elseif (Hybrid_PV!=0)
    #     local xs = ["PV" ,"Onshore", "Hybrid PV"]
    #     local ys = [PV, Onshore, Hybrid_PV]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # elseif (Hybrid_Onshore!=0)
    #     local xs = ["PV" ,"Onshore", "Hybrid Onshore"]
    #     local ys = [PV, Onshore, Hybrid_Onshore]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # else
    #     local xs = ["PV" ,"Onshore"]
    #     local ys = [PV, Onshore]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # end
    #
    # if (Hybrid_PV!=0) && (Hybrid_Onshore!=0)
    #     local xs = ["PV" ,"Onshore", "Hybrid PV" ,"Hybrid Onshore" ]
    #     local ys = [PV_area, Onshore_area, Hybrid_PV_area,  Hybrid_Onshore_area]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # elseif (Hybrid_PV!=0)
    #     local xs = ["PV" ,"Onshore", "Hybrid PV"]
    #     local ys = [PV_area, Onshore_area, Hybrid_PV_area]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # elseif (Hybrid_Onshore!=0)
    #     local xs = ["PV" ,"Onshore", "Hybrid Onshore"]
    #     local ys = [PV_area, Onshore_area, Hybrid_Onshore_area]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # else
    #     local xs = ["PV" ,"Onshore"]
    #     local ys = [PV_area, Onshore_area]
    #     pie(xs, ys, title="$country Green hydrogen Generation Composition",l = 2)
    # end

    local filename_sys_type = "System Nature\\$country system nature $scenario.txt"



    local sys_type_names = ["Originally PV only", "Originally Onshore only", "PV prevalent hybrid", "Onshore prevalent hybrid", "Hybrid became PV only", "Hybrid became Onshore only"];
    local sys_type_nums = [PV, Onshore, Hybrid_PV, Hybrid_Onshore, hyb2pv, hyb2onsh];
    isfile(filename_sys_type) && rm(filename_sys_type);
    # writedlm(filename_sys_type, [sys_type_names sys_type_nums], '\t');

    local filename_sys_type_area = "System Nature\\$country system nature areas $scenario.txt"

    local sys_type_names_area = ["Originally PV only", "Originally Onshore only", "PV prevalent hybrid", "Onshore prevalent hybrid", "Hybrid became PV only", "Hybrid became Onshore only"];
    local sys_type_nums_area = [PV_area, Onshore_area, Hybrid_PV_area, Hybrid_Onshore_area, hyb2pv_area, hyb2onsh_area];
    isfile(filename_sys_type_area) && rm(filename_sys_type_area);
    writedlm(filename_sys_type_area, [sys_type_names_area sys_type_nums_area], '\t');

end
