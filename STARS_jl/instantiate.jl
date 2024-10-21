using Pkg
Pkg.activate(".")

if haskey(ENV, "CONDA_PREFIX")
    try
        run(pipeline(`curl "https://raw.githubusercontent.com/JuliaLang/MbedTLS.jl/master/src/cacert.pem"`, stdout="$CONDA_PREFIX/share/julia/cert.pem"))
    catch
    end
end

try
    run(pipeline(`curl "https://raw.githubusercontent.com/JuliaLang/MbedTLS.jl/master/src/cacert.pem"`, stdout="/opt/conda/share/julia/cert.pem"))
catch
end

Pkg.add("IJulia")
Pkg.build("IJulia")
Pkg.add("HDF5")
Pkg.build("HDF5")
Pkg.add("CoordinateTransformations")
Pkg.add(Pkg.PackageSpec(;name="DimensionalData",version="0.24.7"))
Pkg.add(Pkg.PackageSpec(;name="Rasters",version="0.5.3"))
Pkg.add(Pkg.PackageSpec(;name="GDAL",version="1.5.1"))
Pkg.instantiate()
Pkg.precompile()
