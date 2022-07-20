using DelimitedFiles, GeoArrays


scenario = "2050 V ws"
multiplier = 100

path_read_csv = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\2050 V ws\\big_raster"

filename_csv = "$path_read_csv\\big_raster $scenario (x$multiplier).txt"
data = readdlm(filename_csv, ',')

ga = GeoArrays.GeoArray(data)
GeoArrays.bbox!(ga, (min_x=-180, min_y=90, max_x=180, max_y=-90))
GeoArrays.epsg!(ga, 4326)
println("Saving multiplied GeoTiff...")
path_write_multiplied = path_read_csv
GeoArrays.write!("$path_write_multiplied\\$scenario big_raster (x$multiplier).tif", ga)
