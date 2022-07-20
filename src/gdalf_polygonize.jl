using GDAL

scenario = "2050 V ws"


    path_read_tiff = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario"
    filename_read_tiff = "$scenario lcoh World Matrix.tif"
    path_filename_tiff = "$path_read_tiff\\$filename_read_tiff"

    path_write_shapefile = "E:\\Users\\jdemaigret\\Documents\\Hydrogen Project - EF_IRENA\\Global hydrogen and power model\\Visualization\\$scenario\\Shapefile"
    filename_write_shapefile = "$path_write_shapefile\\$scenario lcoh World Matrix.shp"


raster_file = path_filename_tiff
band = 1
out_file = "$filename_write_shapefile"
fieldname_1 = "LCOH"
nodata = NaN

    layername = "out"
    ds_raster = GDAL.gdalopen(raster_file, GDAL.GA_ReadOnly)
    band = GDAL.gdalgetrasterband(ds_raster, band)

    # if !isnan(nodata); GDAL.gdalsetrasternodatavalue(band, nodata); end
    # dst_ds = GDAL.ogr_dr_createdatasource(drv, dst_filename, C_NULL)
    drive = GDAL.gdalgetdriverbyname("ESRI Shapefile")
    ds_shp = GDAL.gdalcreate(drive, out_file, 0, 0, 0, GDAL.GDT_Unknown, C_NULL)

    REF = GDAL.osrnewspatialreference(C_NULL)
    GDAL.osrimportfromepsg(REF, 4326)
    # REF = C_NULL
    # REF = "WGS_1984"
    layer = GDAL.gdaldatasetcreatelayer(ds_shp, layername, REF, GDAL.wkbPolygon, C_NULL)

    fielddefn = GDAL.ogr_fld_create(fieldname_1, GDAL.OFTInteger)
    field = GDAL.ogr_l_createfield(layer, fielddefn, GDAL.TRUE)

    GDAL.gdalfpolygonize(
        band,   # band
        C_NULL, # mask
        layer, field,
        C_NULL, C_NULL, C_NULL)
    GDAL.gdalclose(ds_shp)
    GDAL.gdalclose(ds_raster)
