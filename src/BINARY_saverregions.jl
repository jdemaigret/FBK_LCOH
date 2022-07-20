using XLSX

function read_save_region()
    excel_file = XLSX.readxlsx("save_regions.xlsx") ##Leggo Excel
    composed_region = excel_file["Region Definition"] ##Leggo Provincie
    n_entries = size(composed_region[:])[1]
    A = Array{Any, 2}(missing, n_entries, 2)
    for i in 1:n_entries
        if composed_region[i,3]=="missing"
            if composed_region[i,2]=="missing"
                #Case 1: State only
                A[i,1] = composed_region[i,1]
                A[i,2] = GADM(composed_region[i,1])
            else
                #Case 2: State and region
                A[i,1] = composed_region[i,2]
                A[i,2] = GADM([composed_region[i,1]], composed_region[i,2])
            end
        else
            #Case 3: State, region and province
            A[i,1] = composed_region[i,3]
            A[i,2] = GADM([composed_region[i,1] , composed_region[i,2]], composed_region[i,3])
        end
    end
    excel_file = XLSX.readxlsx("save_regions.xlsx") ##Leggo Excel
    region_name = excel_file["Region Name"]
    saveregions("$(region_name[:][1])",read_region)
end
