using ArchGDAL, GDAL

scenario = "2050 V ws"

# Read the original geotiff

    path_read_tiff = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario"
    filename_read_tiff = "$scenario lcoh World Matrix.tif"
    path_filename_tiff = "$path_read_tiff\\$filename_read_tiff"
    println("Reading original GeoTiff...")
    raster = ArchGDAL.readraster(path_filename_tiff)

gt = ArchGDAL.getgeotransform(raster)
p = ArchGDAL.getproj(raster)

path_write_shapefile = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario\\Shapefile"

dataset = ArchGDAL.create(
    "$path_write_shapefile\\temporary.tif",
    driver = ArchGDAL.getdriver("GTiff"),
    width=36000,
    height=18000,
    nbands=1,
    dtype=Float64
)

ArchGDAL.write!(dataset, raster[:,:,1], 1)
ArchGDAL.setgeotransform!(dataset, gt)
ArchGDAL.setproj!(dataset, p)


dataset_V1 = ArchGDAL.create(
               "temporary.shp",
               driver = ArchGDAL.getdriver("ESRI Shapefile")
           )

           do ds
           ArchGDAL.createlayer(geom = GDAL.wkbPoint) do layer
               for (lon,lat) in zip(Lon, Lat)
                   ArchGDAL.createfeature(layer) do f
                       ArchGDAL.setgeom!(f, ArchGDAL.createpoint(lat, lon))
                   end
               end
               ArchGDAL.copy(layer, dataset = ds)
           end
       end;
