using Glob
using Dates
using Rasters
using Plots
using LinearAlgebra
using STARS
using STARS.BBoxes
using STARS.sentinel_tiles
using STARS.HLS
using STARS.VNP43
using STARS.STARS
using Logging

using Pkg

Pkg.add("OpenSSL")

using HTTP

BLAS.set_num_threads(1)

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

# command = f'cd "{STARS_source_directory}" && julia --project=@. "{julia_script_filename}" "{tile}" "{coarse_cell_size}" "{fine_cell_size}" "{VIIRS_start_date}" "{VIIRS_end_date}" "{HLS_start_date}" "{HLS_end_date}" "{coarse_directory}" "{fine_directory}" "{posterior_filename}" "{posterior_UQ_filename}" "{posterior_bias_filename}" "{posterior_bias_UQ_filename}" "{prior_filename}" "{prior_UQ_filename}" "{prior_bias_filename}" "{prior_bias_UQ_filename}"'

@info "processing STARS data fusion"
tile = ARGS[1]
@info "tile: $(tile)"
coarse_cell_size = parse(Int64, ARGS[2])
@info "coarse cell size: $(coarse_cell_size)"
fine_cell_size = parse(Int64, ARGS[3])
@info "fine cell size: $(fine_cell_size)"
VIIRS_start_date = Date(ARGS[4])
@info "VIIRS start date: $(VIIRS_start_date)"
VIIRS_end_date = Date(ARGS[5])
@info "VIIRS end date: $(VIIRS_end_date)"
HLS_start_date = Date(ARGS[6])
@info "HLS start date: $(HLS_start_date)"
HLS_end_date = Date(ARGS[7])
@info "HLS end date: $(HLS_end_date)"
coarse_directory = ARGS[8]
@info "coarse inputs directory: $(coarse_directory)"
fine_directory = ARGS[9]
@info "fine inputs directory: $(fine_directory)"
posterior_filename = ARGS[10]
@info "posterior filename: $(posterior_filename)"
posterior_UQ_filename = ARGS[11]
@info "posterior UQ filename: $(posterior_UQ_filename)"
posterior_bias_filename = ARGS[12]
@info "posterior bias filename: $(posterior_bias_filename)"
posterior_bias_UQ_filename = ARGS[13]
@info "posterior bias UQ filename: $(posterior_bias_UQ_filename)"

if size(ARGS)[1] >= 17
    prior_filename = ARGS[14]
    @info "prior filename: $(prior_filename)"
    prior_mean = Raster(prior_filename)
    prior_UQ_filename = ARGS[15]
    @info "prior UQ filename: $(prior_UQ_filename)"
    prior_sd = Raster(prior_UQ_filename)
    prior_bias_filename = ARGS[16]
    @info "prior bias filename: $(prior_bias_filename)"
    prior_bias_mean = Raster(prior_bias_filename)
    prior_bias_UQ_filename = ARGS[17]
    @info "prior bias UQ filename: $(prior_bias_UQ_filename)"
    prior_bias_sd = Raster(prior_bias_UQ_filename)
    prior = DataFusionState(prior_mean, prior_sd, prior_bias_mean, prior_bias_sd, nothing)
else
    prior = nothing
end

x_coarse, y_coarse = sentinel_tile_dims(tile, coarse_cell_size)
x_coarse_size = size(x_coarse)[1]
y_coarse_size = size(y_coarse)[1]
@info "coarse x size: $(x_coarse_size)"
@info "coarse y size: $(y_coarse_size)"
x_fine, y_fine = sentinel_tile_dims(tile, fine_cell_size)
x_fine_size = size(x_fine)[1]
y_fine_size = size(y_fine)[1]
@info "fine x size: $(x_fine_size)"
@info "fine y size: $(y_fine_size)"

coarse_image_filenames = sort(glob("*.tif", coarse_directory))
coarse_dates_found = [Date(split(basename(filename), "_")[3]) for filename in coarse_image_filenames]

fine_image_filenames = sort(glob("*.tif", fine_directory))
fine_dates_found = [Date(split(basename(filename), "_")[3]) for filename in fine_image_filenames]

coarse_start_date = VIIRS_start_date
coarse_end_date = VIIRS_end_date

fine_start_date = HLS_start_date
fine_end_date = HLS_end_date

dates = [fine_start_date + Day(d - 1) for d in 1:((fine_end_date - fine_start_date).value + 1)]
t = Ti(dates)
coarse_dims = (x_coarse, y_coarse, t)
fine_dims = (x_fine, y_fine, t)

covariance_dates = [coarse_start_date + Day(d - 1) for d in 1:((coarse_end_date - coarse_start_date).value + 1)]
t_covariance = Ti(covariance_dates)
covariance_dims = (x_coarse, y_coarse, t_covariance)

covariance_images = []

for (i, date) in enumerate(covariance_dates)
    date = Dates.format(date, dateformat"yyyy-mm-dd")
    match = findfirst(x -> occursin(date, x), coarse_image_filenames)
    timestep_index = Band(i:i)
    timestep_dims = (x_coarse, y_coarse, timestep_index)

    if match === nothing
        @info "coarse image is not available on $(date)"
        covariance_image = Raster(fill(NaN, x_coarse_size, y_coarse_size, 1), dims=timestep_dims, missingval=NaN)
        @info size(covariance_image)
    else
        filename = coarse_image_filenames[match]
        @info "ingesting coarse image on $(date): $(filename)"
        covariance_image = Raster(reshape(Raster(filename), x_coarse_size, y_coarse_size, 1), dims=timestep_dims, missingval=NaN)
        @info size(covariance_image)
    end

    push!(covariance_images, covariance_image)
end

@info "concatenating coarse images for covariance calculation"
covariance_images = Raster(cat(covariance_images..., dims=3), dims=covariance_dims, missingval=NaN)

# estimate spatial var parameter
n_eff = compute_n_eff(Int(round(coarse_cell_size / fine_cell_size)), 2, smoothness=1.5) ## Matern: range = 200m, smoothness = 1.5
sp_var = fast_var_est(covariance_images, n_eff_agg = n_eff)
cov_pars_raster = Raster(fill(NaN, size(covariance_images)[1], size(covariance_images)[2], 4), dims=(covariance_images.dims[1:2]...,Band(1:4)), missingval=covariance_images.missingval)
cov_pars_raster[:,:,1] = sp_var
cov_pars_raster[:,:,2] .= 200
cov_pars_raster[:,:,3] .= 1e-10
cov_pars_raster[:,:,4] .= 1.5

coarse_images = []

for (i, date) in enumerate(dates)
    date = Dates.format(date, dateformat"yyyy-mm-dd")
    match = findfirst(x -> occursin(date, x), coarse_image_filenames)
    timestep_index = Band(i:i)
    timestep_dims = (x_coarse, y_coarse, timestep_index)

    if match === nothing
        @info "coarse image is not available on $(date)"
        coarse_image = Raster(fill(NaN, x_coarse_size, y_coarse_size, 1), dims=timestep_dims, missingval=NaN)
        @info size(coarse_image)
    else
        filename = coarse_image_filenames[match]
        @info "ingesting coarse image on $(date): $(filename)"
        coarse_image = Raster(reshape(Raster(filename), x_coarse_size, y_coarse_size, 1), dims=timestep_dims, missingval=NaN)
        @info size(coarse_image)
    end

    push!(coarse_images, coarse_image)
end

@info "concatenating coarse image inputs"
coarse_images = Raster(cat(coarse_images..., dims=3), dims=coarse_dims, missingval=NaN)

fine_images = []

for (i, date) in enumerate(dates)
    date = Dates.format(date, dateformat"yyyy-mm-dd")
    match = findfirst(x -> occursin(date, x), fine_image_filenames)
    timestep_index = Band(i:i)
    timestep_dims = (x_fine, y_fine, timestep_index)

    if match === nothing
        @info "fine image is not available on $(date)"
        fine_image = Raster(fill(NaN, x_fine_size, y_fine_size, 1), dims=timestep_dims, missingval=NaN)
        @info size(fine_image)
    else
        filename = fine_image_filenames[match]
        @info "ingesting fine image on $(date): $(filename)"
        fine_image = Raster(reshape(Raster(filename), x_fine_size, y_fine_size, 1), dims=timestep_dims, missingval=NaN)
        @info size(fine_image)
    end

    push!(fine_images, fine_image)
end

@info "concatenating fine image inputs"
fine_images = Raster(cat(fine_images..., dims=3), dims=fine_dims, missingval=NaN)

target_date = dates[end]

@info "running data fusion"
fusion_results = coarse_fine_data_fusion(
    coarse_images,
    fine_images,
    cov_pars = cov_pars_raster,
    prior = prior,
    target_times = [target_date],
    buffer_distance = 100.,
    smooth = false,
    offset_ar = [1, 0.0],
    offset_var = [1e-5, 1e-5],
)

@info "writing fused mean: $(posterior_filename)"
write(posterior_filename, Raster(fusion_results.mean, dims=(x_fine, y_fine, Band(1:1)), missingval=NaN))
@info "writing fused SD: $(posterior_UQ_filename)"
write(posterior_UQ_filename, Raster(fusion_results.SD, dims=(x_fine, y_fine, Band(1:1)), missingval=NaN))

@info "writing bias mean: $(posterior_bias_filename)"
write(posterior_bias_filename, Raster(fusion_results.mean_bias, dims=(x_coarse, y_coarse, Band(1:1)), missingval=NaN))
@info "writing bias SD: $(posterior_bias_UQ_filename)"
write(posterior_bias_UQ_filename, Raster(fusion_results.SD_bias, dims=(x_coarse, y_coarse, Band(1:1)), missingval=NaN))

