module VNP43NRT

using Dates
using LinearAlgebra
using Statistics
using Rasters
using DimensionalData.Dimensions.LookupArrays
using ProgressMeter

SINUSOIDAL_CRS = WellKnownText("PROJCS[\"unknown\",GEOGCS[\"unknown\",DATUM[\"unknown\",SPHEROID[\"unknown\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"longitude_of_center\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")

function sinusoidal_tile_dims(h::Int, v::Int, tile_width_cells::Int)::Tuple{X,Y}
    # boundaries of sinusodial projection
    GLOBE_UPPER_LEFT_X = -20015109.355798
    GLOBE_UPPER_LEFT_Y = 10007554.677899
    GLOBE_LOWER_RIGHT_X = 20015109.355798
    GLOBE_LOWER_RIGHT_Y = -10007554.677899

    # size across (width or height) of any equal-area sinusoidal tile
    TILE_SIZE = 1111950.5197665554

    # rows and columns of sinusoidal tile grid
    TOTAL_ROWS = 18
    TOTAL_COLUMNS = 36

    cell_size = TILE_SIZE / tile_width_cells
    tile_left_x = GLOBE_UPPER_LEFT_X + h * TILE_SIZE
    tile_right_x = GLOBE_UPPER_LEFT_X + (h + 1) * TILE_SIZE - cell_size
    tile_upper_y = GLOBE_LOWER_RIGHT_Y + (TOTAL_ROWS - v) * TILE_SIZE - cell_size
    tile_lower_y = GLOBE_LOWER_RIGHT_Y + (TOTAL_ROWS - 1 - v) * TILE_SIZE
    sampling = Intervals(Start())
    x_dim = X(Projected(LinRange(tile_left_x, tile_right_x, tile_width_cells), order=ForwardOrdered(), span=Regular(cell_size), sampling=sampling, crs=SINUSOIDAL_CRS))
    y_dim = Y(Projected(LinRange(tile_upper_y, tile_lower_y, tile_width_cells), order=ReverseOrdered(), span=Regular(-cell_size), sampling=sampling, crs=SINUSOIDAL_CRS))
    dims = (x_dim, y_dim)

    return dims
end

export sinusoidal_tile_dims

struct BRDFParameters
    brdf::AbstractMatrix
    se::AbstractVector
    R2::AbstractVector
    uncert::AbstractArray{Float64,3}
end

function zenith_from_solarnoon_time(tm::Float64, lat::Float64, lon::Float64)
    rad = π/180.0

    Jd = tm/86400.0 + 2440587.5
    Jc = (Jd - 2451545.0)/36525.0
    L0 = mod((280.46646 + Jc * (36000.76983 + 0.0003032 * Jc)),360)
    M = 357.52911 + Jc * (35999.05029 - 0.0001537 * Jc)
    e = 0.016708634 - Jc * (4.2037e-05 + 1.267e-07 * Jc)
    eqctr = sin(rad * M) * (1.914602 - Jc * (0.004817 + 1.4e-05 *
        Jc)) + sin(rad * 2.0 * M) * (0.019993 - 0.000101 * Jc) +
        sin(rad * 3.0 * M) * 0.000289
    lambda0 = L0 + eqctr
    omega = 125.04 - 1934.136 * Jc
    lambda = lambda0 - 0.00569 - 0.00478 * sin(rad * omega)
    seconds = 21.448 - Jc * (46.815 + Jc * (0.00059 - Jc * (0.001813)))
    obliq0 = 23.0 + (26.0 + (seconds/60.0))/60.0
    obliq = obliq0 + 0.00256 * cos(rad * omega)
    y = tan(rad * obliq/2)^2
    eqnTime = 4/rad * (y*sin(rad*2*L0) - 2*e*sin(rad*M) +
        4*e*y*sin(rad*M)*cos(rad*2*L0) - y^2/2*sin(rad*4*L0) -
        e^2*1.25*sin(rad*2*M))
    solarDec = asin(sin(rad*obliq)*sin(rad*lambda))
    sinSolarDec = sin(solarDec)
    cosSolarDec = cos(solarDec)
    solarTime = (mod(Jd-1/2,1)*1440+eqnTime)/4
    hourAngle = solarTime+lon-180
    cosZenith = sin(rad*lat)*sinSolarDec+cos(rad*lat)*cosSolarDec*cos(rad*hourAngle)

    if cosZenith < -1
        cosZenith=-1
    end

    if cosZenith > 1
        cosZenith=1
    end

    return acos(cosZenith)/rad
end

function zenith_from_solarnoon_vec(tm::Vector{Float64}, lat::Vector{Float64}, lon::Vector{Float64})
    n = length(tm)
    szn = Vector{Float64}(undef, n)
    
    for i in 1:n
        szn[i] = zenith_from_solarnoon_time(tm[i], lat[i], lon[i])
    end

    return szn
end

function Kvol(sz::Matrix, vz::Matrix, rz::Matrix)::Matrix
    sc = π / 180
    eps = acos.(cos.(sc .* sz) .* cos.(sc .* vz) .+ sin.(sc .* sz) .* sin.(sc .* vz) .* cos.(sc .* rz))
    Kvol = 1 ./ (cos.(sc .* sz) .+ cos.(sc .* vz)) .* ((π / 2 .- eps) .* cos.(eps) .+ sin.(eps)) .- π / 4

    return Kvol
end

function Kvol_vec(sz::AbstractVector, vz::AbstractVector, rz::AbstractVector)
    sc = pi/180
    eps = acos.(cos.(sc .* sz) .* cos.(sc .* vz) .+ sin.(sc .* sz) .* sin.(sc .* vz) .* cos.(sc .* rz))
    Kvol_vec = 1 ./ (cos.(sc .* sz) .+ cos.(sc .* vz)) .* ((pi/2 .- eps) .* cos.(eps) .+ sin.(eps)) .- pi/4

    return Kvol_vec
end

function Kvol_sc(sz::Real, vz::Real, rz::Real)
    sc = pi/180
    eps = acos(cos(sc*sz)*cos(sc*vz) + sin(sc*sz)*sin(sc*vz)*cos(sc*rz))
    Kvol_sc = 1/(cos(sc*sz) + cos(sc*vz))*((pi/2-eps)*cos(eps) + sin(eps)) - pi/4

    return Kvol_sc
end

function Kgeo_sc(sz::Real, vz::Real, rz::Real)
    sc = pi/180
    eps = acos(cos(sc*sz)*cos(sc*vz) + sin(sc*sz)*sin(sc*vz)*cos(sc*rz))
    D = sqrt(tan(sc*sz)^2 + tan(sc*vz)^2 - 2*tan(sc*sz)*tan(sc*vz)*cos(sc*rz))
    cost = 2*sqrt(D^2 + (tan(sc*sz)*tan(sc*vz)*sin(sc*rz))^2)/(1/cos(sc*sz) + 1/cos(sc*vz))
    cost = clamp(cost, -1.0, 1.0)
    t = acos(cost)
    O = 1/pi*(t-sin(t)*cos(t))*(1/cos(sc*sz) + 1/cos(sc*vz))
    O = max(O, 0)
    Kgeo_sc = O - 1/cos(sc*sz) - 1/cos(sc*vz) + 0.5*(1+cos(eps))/cos(sc*sz)/cos(sc*vz)

    return Kgeo_sc
end

function Kgeo(sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix)
    n, t = size(sz)
    Kg = zeros(n, t)
    for i in 1:n
        for j in 1:t
            Kg[i,j] = Kgeo_sc(sz[i,j], vz[i,j], rz[i,j])
        end
    end
    return Kg
end

function Kgeo_vec(sz::AbstractVector, vz::AbstractVector, rz::AbstractVector)
    t = length(vz)
    Kg = zeros(t)
    
    for j in 1:t
        Kg[j] = Kgeo_sc(sz[j], vz[j], rz[j])
    end

    return Kg
end

function NRT_BRDF(Y::AbstractMatrix, kv::AbstractMatrix, kg::AbstractMatrix, weighted::Bool, scale::Real)
    n, p = size(Y)
    brdf = fill(NaN, n, 3)
    R2 = fill(NaN, n)
    se = fill(NaN, n)
    uncert = fill(NaN, 3, 3, n)

    x = ones(3, p)
    for i in 1:n
        yt = Y[i,:]
        non_missing = findall(isfinite.(yt))
        nt = length(non_missing)

        if nt < 7
            continue
        else
            x[2,:] = kv[i,:]
            x[3,:] = kg[i,:]

            xt = x[:,non_missing]
            ytt = yt[non_missing]

            Si = inv(xt * xt')
            brdf[i,:] = (Si * xt * ytt)'
            yp = (brdf[i,:] * xt)'
            uncert[:,:,i] = Si
            se[i] = sqrt(sum((ytt - yp).^2)/(nt-3))
            R2[i] = 1 - se[i]^2 * (nt-3) / var(ytt; corrected=true) / nt
        end
    end

    return BRDFParameters(brdf, se, R2, uncert)
end

function NRT_BRDF_albedo(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool, scale::Real)
    # RossThick constant
    g0vol = -0.007574
    g1vol = -0.070987
    g2vol = 0.307588

    # LiSparseR constant
    g0geo = -1.284909
    g1geo = -0.166314
    g2geo = 0.041840

    gwsa = [1.0, 0.189184, -1.377622]
    gbsa = [1.0, 0.0, 0.0]

    n, p = size(Y)
    @info "reflectance n: $(n) p: $(p)"
    results = fill(NaN, n, 7) # (wsa, bsa, wsa_se, bsa_se, rmse, R2, nt)
    @info "results shape rows: $(size(results)[1]) cols: $(size(results)[2])"

    x = ones(3, p)
    for i in 1:n
        yt = Y[i,:]
        non_missing = findall(isfinite.(yt))
        nt = length(non_missing)

        if nt < 7
            continue
        else
            sznrad = soz_noon[i] * pi/180
            gbsa[2] = g0vol + g1vol * sznrad^2 + g2vol * sznrad^3
            gbsa[3] = g0geo + g1geo * sznrad^2 + g2geo * sznrad^3

            x[2,:] = Kvol_vec(sz[i,:], vz[i,:], rz[i,:])
            x[3,:] = Kgeo_vec(sz[i,:], vz[i,:], rz[i,:])

            xt = x[:,non_missing]
            ytt = yt[non_missing]

            if weighted
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing]

                SSi = Diagonal(xxt)
                Si = inv(xt * SSi * xt')
                brdf = Si * xt * SSi * ytt
                yp = (brdf' * xt)'
                se = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3))
            else
                Si = inv(xt * xt')
                brdf = Si * xt * ytt
                yp = (brdf' * xt)'
                se = sqrt(sum((ytt - yp).^2)/(nt-3))
            end

            results[i,1] = dot(gwsa, brdf) # wsa
            results[i,2] = dot(gbsa, brdf) # bsa

            results[i,3] = se * sqrt(dot(gwsa, Si * gwsa)) # wsa se
            results[i,4] = se * sqrt(dot(gbsa, Si * gbsa)) # bsa se

            results[i,5] = se # brdf rmse
            results[i,6] = 1 - se^2 * (nt-3) / var(ytt; corrected=true) / nt # brdf R2
            results[i,7] = nt # number of obs for brdf estimation
        end
    end

    return results
end

function NRT_BRDF_nadir(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool, scale::Real)
    n, p = size(Y)
    results = fill(NaN, n, 5) # (nadir, nadir_se, rmse, R2, nt)
    x = ones(3, p)
    Knadir = [1.0, 0.0, 0.0]

    for i in 1:n
        yt = Y[i,:]
        non_missing = findall(isfinite.(yt))
        nt = length(non_missing)

        if nt < 7
            continue
        else
            szn = soz_noon[i]
            x[2,:] = Kvol_vec(sz[i,:], vz[i,:], rz[i,:])
            x[3,:] = Kgeo_vec(sz[i,:], vz[i,:], rz[i,:])

            xt = x[:,non_missing]
            ytt = yt[non_missing]

            if weighted
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing]
                SSi = Diagonal(xxt)
                Si = inv(xt * SSi * xt')
                brdf = Si * xt * SSi * ytt
                yp = (brdf' * xt)'
                se = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3))
            else
                Si = inv(xt * xt')
                brdf = Si * xt * ytt
                yp = (brdf' * xt)'
                se = sqrt(sum((ytt - yp).^2)/(nt-3))
            end

            Knadir[2] = Kvol_sc(szn, 0.0, 0.0)
            Knadir[3] = Kgeo_sc(szn, 0.0, 0.0)
            results[i,1] = dot(Knadir, brdf) # nadir
            results[i,2] = se * sqrt(dot(Knadir, Si * Knadir)) # nadir se
            results[i,3] = se # brdf rmse
            results[i,4] = 1 - se^2 * (nt-3) / var(ytt; corrected=true) / nt # brdf R2
            results[i,5] = nt # number of obs for brdf estimation
        end
    end

    return results
end

function NRT_BRDF_all(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool = true, scale::Real = 1.87)
    # the rows are separate locations and the columns are separate times
    @info "processing BRDF"
    @info "reflectance rows: $(size(Y)[1]) cols: $(size(Y)[2])"
    @info "solar zenith rows: $(size(sz)[1]) cols: $(size(sz)[2])"
    @info "sensor zenith rows: $(size(vz)[1]) cols: $(size(vz)[2])"
    @info "relative azimuth rows: $(size(rz)[1]) cols: $(size(rz)[2])"
    @info "solar zenith noon size: $(size(soz_noon)[1])"

    # RossThick constants
    g0vol = -0.007574
    g1vol = -0.070987
    g2vol = 0.307588

    # LiSparseR constants
    g0geo = -1.284909
    g1geo = -0.166314
    g2geo = 0.041840

    gwsa = [1.0, 0.189184, -1.377622]
    gbsa = [1.0, 0.0, 0.0]
    gnbar = [1.0, 0.0, 0.0]

    n, p = size(Y)
    @info "n: $(n) p: $(p)"
    results = fill(NaN, n, 9) # (wsa, bsa, nadir, wsa_se, bsa_se, nadir_se, rmse, R2, nt)
    @info "results rows: $(size(results)[1]) cols: $(size(results)[2])"

    x = ones(3, p)

    @showprogress for i in 1:n
        yt = Y[i,:]
        non_missing = findall(isfinite.(yt))
        nt = length(non_missing)

        if nt < 7
            continue
        else
            sznrad = soz_noon[i] * pi/180
            gbsa[2] = g0vol + g1vol * sznrad^2 + g2vol * sznrad^3
            gbsa[3] = g0geo + g1geo * sznrad^2 + g2geo * sznrad^3

            gnbar[2] = Kvol_sc(soz_noon[i], 0.0, 0.0)
            gnbar[3] = Kgeo_sc(soz_noon[i], 0.0, 0.0)

            x[2,:] = Kvol_vec(sz[i,:], vz[i,:], rz[i,:])
            x[3,:] = Kgeo_vec(sz[i,:], vz[i,:], rz[i,:])

            xt = x[:,non_missing]
            ytt = yt[non_missing]

            if weighted
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing]

                SSi = Diagonal(xxt)
                Si = inv(xt * SSi * xt')
                brdf = Si * xt * SSi * ytt
                yp = (brdf' * xt)'
                se = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3))
            else
                Si = inv(xt * xt')
                brdf = Si * xt * ytt
                yp = (brdf' * xt)'
                se = sqrt(sum((ytt - yp).^2)/(nt-3))
            end

            results[i,1] = dot(gwsa, brdf) # wsa
            results[i,2] = dot(gbsa, brdf) # bsa
            results[i,3] = dot(gnbar, brdf) # nadir

            results[i,4] = se * sqrt(dot(gwsa, Si * gwsa)) # wsa se
            results[i,5] = se * sqrt(dot(gbsa, Si * gbsa)) # bsa se
            results[i,6] = se * sqrt(dot(gnbar, Si * gnbar)) # nadir se

            results[i,7] = se # brdf rmse
            results[i,8] = 1 - se^2 * (nt-3) / var(ytt; corrected=true) / nt # brdf R2
            results[i,9] = nt # number of obs for brdf estimation

            replace!(results, missing=>NaN)
        end
    end

    return results
end

export NRT_BRDF_all

function date_range(start::Date, stop::Date)::Vector{Date}
    return collect(start:Day(1):stop)
end

function calculate_SZA(day_of_year::Array{Int64,1}, hour_of_day::Union{Array{Float64,1}, Float64}, lat::Array{Float64,1})::Array{Float64,1}
    day_angle = (2 * pi * (day_of_year .- 1)) / 365
    lat = deg2rad.(lat)
    dec = deg2rad.((0.006918 .- 0.399912 .* cos.(day_angle) .+ 0.070257 .* sin.(day_angle) .- 0.006758 .* cos.(2 .* day_angle) .+ 0.000907 .* sin.(2 .* day_angle) .- 0.002697 .* cos.(3 .* day_angle) .+ 0.00148 .* sin.(3 .* day_angle)) .* (180 / pi))
    hour_angle = deg2rad.(hour_of_day * 15.0 .- 180.0)
    SZA = rad2deg.(acos.(sin.(lat) .* sin.(dec) .+ cos.(lat) .* cos.(dec) .* cos.(hour_angle)))
    return SZA
end

end
