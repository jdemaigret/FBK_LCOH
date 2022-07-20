using DelimitedFiles
nomi_stringa = ["ARG" "AUS" "BRA" "CAN" "CHI" "CHL" "COL" "EAS" "EUR" "FRA" "GER" "IDN" "IND" "ITA" "JAP" "KOR" "LAM" "MEN" "MEX" "MOR" "OCE" "POR" "ROA" "ROE" "RUS" "SAU" "SEA" "SOU" "SPA" "SSA" "TUR" "UKG" "UKR" "USA"]
# nomi_stringa = ["UKG"]
# global scenario = "2030 Opt ws"
global scenario = "2050 V ws"
global hybrid_threshold_PV = 0.01
global hybrid_threshold_Onshore = 0.01
global read_path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\sys_nature\\Cap_matrices"
global write_path = "C:\\Users\\jdemaigret\\Documents\\EF-IRENA\\sys_nature\\Sys_nat_matrices"

@time for country in nomi_stringa
    println("Country: $country")
    local pv_cap_mat = readdlm("$read_path\\$scenario\\$country pv_matrix ($scenario).txt",',')
    local onsh_cap_mat = readdlm("$read_path\\$scenario\\$country onshore_matrix ($scenario).txt",',')
    local lcoh_mat = readdlm("$read_path\\$scenario\\$country lcoh_matrix ($scenario).txt",',')
    local mat_rows = size(pv_cap_mat)[1]
    local mat_cols = size(pv_cap_mat)[2]
    global ori_pv_mat       = zeros(size(pv_cap_mat)) # Originally PV only
    global ori_onsh_mat     = zeros(size(pv_cap_mat)) # Originally Onshore only
    global hy_2_pv_mat      = zeros(size(pv_cap_mat)) # Hybrid became PV only
    global hy_pv_prev_mat   = zeros(size(pv_cap_mat)) # PV prevalent hybrid
    global hy_2_onsh_mat    = zeros(size(pv_cap_mat)) # Hybrid became Onshore only
    global hy_onsh_prev_mat = zeros(size(pv_cap_mat)) # Onshore prevalent hybrid

    for i in 1:mat_rows
        for j in 1:mat_cols
            if (pv_cap_mat[i,j]>0) && (onsh_cap_mat[i,j]==0)                         # Originally PV only
                local ori_pv_mat[i,j] = lcoh_mat[i,j]
            elseif (pv_cap_mat[i,j]==0) && (onsh_cap_mat[i,j]>0)                     # Originally Onshore only
                local ori_onsh_mat[i,j] = lcoh_mat[i,j]
            elseif (pv_cap_mat[i,j] > onsh_cap_mat[i,j]) && (onsh_cap_mat[i,j]>0)
                if (onsh_cap_mat[i,j] / pv_cap_mat[i,j]) <= hybrid_threshold_PV      # Hybrid became PV only
                    local hy_2_pv_mat[i,j] = lcoh_mat[i,j]
                else                                                                 # PV prevalent hybrid
                    local hy_pv_prev_mat[i,j] = lcoh_mat[i,j]
                end
            elseif (onsh_cap_mat[i,j] > pv_cap_mat[i,j]) && (pv_cap_mat[i,j]>0)
                if (pv_cap_mat[i,j] / onsh_cap_mat[i,j]) <= hybrid_threshold_Onshore # Hybrid became Onshore only
                    local hy_2_onsh_mat[i,j] = lcoh_mat[i,j]
                else                                                                 # Onshore prevalent hybrid
                    local hy_onsh_prev_mat[i,j] = lcoh_mat[i,j]
                end
            end
        end
    end

    local filename_ori_pv       = "$write_path\\$scenario\\$country Ori PV $scenario.txt";
    local filename_ori_onsh     = "$write_path\\$scenario\\$country Ori Onsh $scenario.txt";
    local filename_hy_2_pv      = "$write_path\\$scenario\\$country Hy to PV $scenario.txt";
    local filename_hy_pv_prev   = "$write_path\\$scenario\\$country Hy PV Prev $scenario.txt";
    local filename_hy_2_onsh    = "$write_path\\$scenario\\$country Hy to Onsh $scenario.txt";
    local filename_hy_onsh_prev = "$write_path\\$scenario\\$country Hy Onsh Prev $scenario.txt";
    isfile(filename_ori_pv) && rm(filename_ori_pv);
    isfile(filename_ori_onsh) && rm(filename_ori_onsh);
    isfile(filename_hy_2_pv) && rm(filename_hy_2_pv);
    isfile(filename_hy_pv_prev) && rm(filename_hy_pv_prev);
    isfile(filename_hy_2_onsh) && rm(filename_hy_2_onsh);
    isfile(filename_hy_onsh_prev) && rm(filename_hy_onsh_prev);
    writedlm(filename_ori_pv, ori_pv_mat, ',');
    writedlm(filename_ori_onsh, ori_onsh_mat, ',');
    writedlm(filename_hy_2_pv, hy_2_pv_mat, ',');
    writedlm(filename_hy_pv_prev, hy_pv_prev_mat, ',');
    writedlm(filename_hy_2_onsh, hy_2_onsh_mat, ',');
    writedlm(filename_hy_onsh_prev, hy_onsh_prev_mat, ',');

end
