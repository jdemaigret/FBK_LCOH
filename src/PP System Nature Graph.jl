using DelimitedFiles

region = "AUS"
time_horizon = 2050
scenario = "Optimistic"

# Country    Best  PV_Best  Onshore_Best
# AUS    ->   1-2    1-4     2-2
# CHL    ->   1-1    1-3     3-1
# CHI    ->   1-2    1-4     2-2
# GER    ->   3-2    1-4     4-2
# JAP    ->   2-2    2-4     4-2
# SAU    ->   1-3    1-5     3-3
# USA    ->   1-1    1-3     3-1

PV_quality = 1
Onshore_quality = 2

r_path = "C:\\Users\\jdemaigret\\.julia\\dev\\GlobalEnergyGIS_lcoh\\src\\HEATMAPS\\RAW"

r_data = readdlm("$r_path\\$region Gen Tech Ratio Matrix $time_horizon $scenario (pv$PV_quality, w$Onshore_quality).txt")

rows = 1:size(r_data)[1]
columns = 1:size(r_data)[2]

PV_bin_min      = 0.98
PV_bin_max      = 1

Onshore_bin_min = 0
Onshore_bin_max = 0.001

hybrid_bin_increment = (PV_bin_min-Onshore_bin_max)/3
hybrid_bin3_min = Onshore_bin_max
hybrid_bin3_max = hybrid_bin3_min + hybrid_bin_increment
hybrid_bin2_min = hybrid_bin3_max
hybrid_bin2_max = hybrid_bin2_min + hybrid_bin_increment
hybrid_bin1_min = hybrid_bin3_max
hybrid_bin1_max = hybrid_bin1_min + hybrid_bin_increment

for r in rows
    for c in columns
        if (PV_bin_min <= r_data[r,c]) && (r_data[r,c] <= PV_bin_max)  # 0.98<=x<=1
            r_data[r,c] = 1
        elseif (Onshore_bin_min <= r_data[r,c]) && (r_data[r,c] <= Onshore_bin_max) # 0<=x<=0.001
            r_data[r,c] = 0
        elseif (Onshore_bin_max < r_data[r,c]) && (r_data[r,c] <= hybrid_bin3_max) # 0.001<x<=0.32
            r_data[r,c] = 0.25
        elseif (hybrid_bin3_max < r_data[r,c]) && (r_data[r,c] <= hybrid_bin2_max) # 0.32<x<=0.64
            r_data[r,c] = 0.5
        elseif (hybrid_bin2_max < r_data[r,c]) && (r_data[r,c] < PV_bin_min) # 0.64<x<0.98
            r_data[r,c] = 0.75
        end
    end
end


w_path = "C:\\Users\\jdemaigret\\.julia\\dev\\GlobalEnergyGIS_lcoh\\src\\HEATMAPS\\PP"
filename_PP = "PP $region Gen Tech Ratio Matrix $time_horizon $scenario (pv$PV_quality, w$Onshore_quality).txt"
isfile(filename_PP) && rm(filename_PP);
writedlm("$w_path\\$filename_PP", r_data, '\t');
