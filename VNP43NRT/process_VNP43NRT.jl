using Glob
using Dates
using Rasters
using DimensionalData.Dimensions.LookupArrays
import ArchGDAL
using VNP43NRT
using Logging

struct CustomLogger <: AbstractLogger
    stream::IO
    min_level::LogLevel
end

Logging.min_enabled_level(logger::CustomLogger) = logger.min_level

function Logging.shouldlog(logger::CustomLogger, level, _module, group, id)
    return level >= logger.min_level
end

function Logging.handle_message(logger::CustomLogger, level, message, _module, group, id, file, line; kwargs...)
    t = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(logger.stream, "[$t $(uppercase(string(level)))] $message")
end

global_logger(CustomLogger(stdout, Logging.Info))

SINUSOIDAL_CRS = WellKnownText("PROJCS[\"unknown\",GEOGCS[\"unknown\",DATUM[\"unknown\",SPHEROID[\"unknown\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"longitude_of_center\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")

band = ARGS[1]
@info "band: $(band)"
h = parse(Int64, ARGS[2])
v = parse(Int64, ARGS[3])
@info "h: $(h) v: $(v)"
tile_width_cells = parse(Int64, ARGS[4])
@info "tile width cells: $(tile_width_cells)"
start_date = Date(ARGS[5])
@info "start date: $(start_date)"
end_date = Date(ARGS[6])
@info "end date: $(end_date)"
reflectance_directory = ARGS[7]
@info "reflectance directory: $(reflectance_directory)"
solar_zenith_directory = ARGS[8]
@info "solar zenith directory: $(solar_zenith_directory)"
sensor_zenith_directory = ARGS[9]
@info "sensor zenith directory: $(sensor_zenith_directory)"
relative_azimuth_directory = ARGS[10]
@info "relative azimuth directory: $(relative_azimuth_directory)"
SZA_filename = ARGS[11]
@info "solar zenith noon file: $(SZA_filename)"
output_directory = ARGS[12]
@info "output directory: $(output_directory)"

reflectance_image_filenames = sort(glob("*.tif", reflectance_directory))

dates = [start_date + Day(d - 1) for d in 1:((end_date - start_date).value + 1)]

t = Ti(dates)

x_dim, y_dim = sinusoidal_tile_dims(h, v, tile_width_cells)

function load_timeseries(directory::String, variable::String, start_date::Date, end_date::Date, x_dim, y_dim)
    @info "searching directory: $(directory)"
    filenames = sort(glob("*.tif", directory))
    images = []

    dates = [start_date + Day(d - 1) for d in 1:((end_date - start_date).value + 1)]

    for date in dates
        date = Dates.format(date, dateformat"yyyy-mm-dd")
        match = findfirst(x -> occursin(date, x), filenames)
    
        if match === nothing
            @info "$(variable) image is not available on $(date)"
            image = Raster(fill(NaN, tile_width_cells, tile_width_cells, 1), dims=(x_dim, y_dim, Band(1:1)), missingval=NaN)
        else
            filename = filenames[match]
            @info "ingesting $(variable) image on $(date): $(filename)"
            image = Raster(filename, dims=(x_dim, y_dim, Band(1:1)))
        end
    
        push!(images, image)
    end
    
    return images
end

reflectance_images = load_timeseries(reflectance_directory, band, start_date, end_date, x_dim, y_dim)

function stack_timeseries(timeseries)
    permutedims(hcat([vec(image) for image in timeseries]...), [1,2])
end

reflectance_stack = stack_timeseries(reflectance_images)

solar_zenith_images = load_timeseries(solar_zenith_directory, band, start_date, end_date, x_dim, y_dim)
solar_zenith_stack = stack_timeseries(solar_zenith_images)

sensor_zenith_images = load_timeseries(sensor_zenith_directory, band, start_date, end_date, x_dim, y_dim)
sensor_zenith_stack = stack_timeseries(sensor_zenith_images)

relative_azimuth_images = load_timeseries(relative_azimuth_directory, band, start_date, end_date, x_dim, y_dim)
relative_azimuth_stack = stack_timeseries(relative_azimuth_images)

SZA = Raster(SZA_filename)
SZA_flat = vec(SZA)

results = NRT_BRDF_all(reflectance_stack, solar_zenith_stack, sensor_zenith_stack, relative_azimuth_stack, SZA_flat)

mkpath(output_directory)

date_format = dateformat"yyyy-mm-dd"
date_stamp = Dates.format(end_date, date_format)

WSA = Raster(reshape(results[:,1], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
WSA_filename = joinpath(output_directory, "$(date_stamp)_WSA.tif")
@info "writing WSA: $(WSA_filename)"
write(WSA_filename, WSA)

BSA = Raster(reshape(results[:,2], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
BSA_filename = joinpath(output_directory, "$(date_stamp)_BSA.tif")
@info "writing BSA: $(BSA_filename)"
write(BSA_filename, BSA)

NBAR = Raster(reshape(results[:,3], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
NBAR_filename = joinpath(output_directory, "$(date_stamp)_NBAR.tif")
@info "writing NBAR: $(NBAR_filename)"
write(NBAR_filename, NBAR)

WSA_SE = Raster(reshape(results[:,4], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
WSA_SE_filename = joinpath(output_directory, "$(date_stamp)_WSA_SE.tif")
@info "writing WSA_SE: $(WSA_SE_filename)"
write(WSA_SE_filename, WSA_SE)

BSA_SE = Raster(reshape(results[:,5], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
BSA_SE_filename = joinpath(output_directory, "$(date_stamp)_BSA_SE.tif")
@info "writing BSA_SE: $(BSA_SE_filename)"
write(BSA_SE_filename, BSA_SE)

NBAR_SE = Raster(reshape(results[:,6], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
NBAR_SE_filename = joinpath(output_directory, "$(date_stamp)_NBAR_SE.tif")
@info "writing NBAR_SE: $(NBAR_SE_filename)"
write(NBAR_SE_filename, NBAR_SE)

BRDF_SE = Raster(reshape(results[:,7], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
BRDF_SE_filename = joinpath(output_directory, "$(date_stamp)_BRDF_SE.tif")
@info "writing BRDF_SE: $(BRDF_SE_filename)"
write(BRDF_SE_filename, BRDF_SE)

BRDF_R2 = Raster(reshape(results[:,8], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
BRDF_R2_filename = joinpath(output_directory, "$(date_stamp)_BRDF_R2.tif")
@info "writing BRDF_R2: $(BRDF_R2_filename)"
write(BRDF_R2_filename, BRDF_R2)

count = Raster(reshape(results[:,9], (tile_width_cells, tile_width_cells)), dims=(x_dim, y_dim), missingval=NaN)
count_filename = joinpath(output_directory, "$(date_stamp)_count.tif")
@info "writing count: $(count_filename)"
write(count_filename, count)
