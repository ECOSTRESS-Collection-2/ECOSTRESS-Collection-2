using Glob
using Dates
using Rasters
using Plots
using STARS
using STARS.BBoxes
using STARS.sentinel_tiles
using STARS.HLS
using STARS.VNP43
using STARS.STARS

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

x_coarse, y_coarse = sentinel_tile_dims(tile, coarse_cell_size)
x_fine, y_fine = sentinel_tile_dims(tile, fine_cell_size)

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

for date in covariance_dates
    date = Dates.format(date, dateformat"yyyy-mm-dd")
    match = findfirst(x -> occursin(date, x), coarse_image_filenames)

    if match === nothing
        @info "coarse image is not available on $(date)"
        covariance_image = Raster(fill(NaN, size(x_coarse)[1], size(y_coarse)[1], 1), dims=(x_coarse, y_coarse, Band(1:1)), missingval=NaN)
    else
        filename = coarse_image_filenames[match]
        @info "ingesting coarse image on $(date): $(filename)"
        covariance_image = Raster(filename)
    end

    push!(covariance_images, covariance_image)
end

covariance_images = Raster(cat(covariance_images..., dims=3), dims=covariance_dims, missingval=NaN)

# estimate spatial var parameter
n_eff = compute_n_eff(7,2,smoothness=1.5) ## Matern: range = 200m, smoothness = 1.5
sp_var = fast_var_est(covariance_images, n_eff_agg = n_eff)
cov_pars_raster = Raster(fill(NaN, size(covariance_images)[1], size(covariance_images)[2], 4), dims=(covariance_images.dims[1:2]...,Band(1:4)), missingval=covariance_images.missingval)
cov_pars_raster[:,:,1] = sp_var
cov_pars_raster[:,:,2] .= 200
cov_pars_raster[:,:,3] .= 1e-10
cov_pars_raster[:,:,4] .= 1.5

coarse_images = []

for date in dates
    date = Dates.format(date, dateformat"yyyy-mm-dd")
    match = findfirst(x -> occursin(date, x), coarse_image_filenames)

    if match === nothing
        @info "coarse image is not available on $(date)"
        coarse_image = Raster(fill(NaN, size(x_coarse)[1], size(y_coarse)[1], 1), dims=(x_coarse, y_coarse, Band(1:1)), missingval=NaN)
    else
        filename = coarse_image_filenames[match]
        @info "ingesting coarse image on $(date): $(filename)"
        coarse_image = Raster(filename)
    end

    push!(coarse_images, coarse_image)
end

coarse_images = Raster(cat(coarse_images..., dims=3), dims=coarse_dims, missingval=NaN)

fine_images = []

for date in dates
    date = Dates.format(date, dateformat"yyyy-mm-dd")
    match = findfirst(x -> occursin(date, x), fine_image_filenames)

    if match === nothing
        @info "fine image is not available on $(date)"
        fine_image = Raster(fill(NaN, size(x_fine)[1], size(y_fine)[1], 1), dims=(x_fine, y_fine, Band(1:1)), missingval=NaN)
    else
        filename = fine_image_filenames[match]
        @info "ingesting fine image on $(date): $(filename)"
        fine_image = Raster(filename)
    end

    push!(fine_images, fine_image)
end

fine_images = Raster(cat(fine_images..., dims=3), dims=fine_dims, missingval=NaN)

target_date = dates[end]

@info "running data fusion"
fusion_results = coarse_fine_data_fusion(
    coarse_images, 
    fine_images, 
    cov_pars = cov_pars_raster,
    target_times = [target_date],
    buffer_distance = 100.,
    smooth = false
)

@info "writing fused mean: $(posterior_filename)"
write(posterior_filename, Raster(fusion_results.mean, dims=(x_fine, y_fine, Band(1:1)), missingval=NaN))
@info "writing fused SD: $(posterior_UQ_filename)"
write(posterior_UQ_filename, Raster(fusion_results.SD, dims=(x_fine, y_fine, Band(1:1)), missingval=NaN))

@info "writing bias mean: $(posterior_bias_filename)"
write(posterior_bias_filename, Raster(fusion_results.mean_bias, dims=(x_coarse, y_coarse, Band(1:1)), missingval=NaN))
@info "writing bias SD: $(posterior_bias_UQ_filename)"
write(posterior_bias_UQ_filename, Raster(fusion_results.SD_bias, dims=(x_coarse, y_coarse, Band(1:1)), missingval=NaN))
