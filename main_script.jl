#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# Structural connectivity index #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

# for speed, use multithreading by launching julia with multiple threads (adjust according to your hardware)
# julia --threads 10

# 0. Preamble  -------------------
nothing

using Pkg
Pkg.activate(@__DIR__) 
Pkg.instantiate()

using Rasters
using CairoMakie

include("functions.jl")

kommuner = ["moss", "bodo", "holtalen", "sandefjord", "ullensaker", "nordreisa"]
okosystemer = ["naturskog", "broadleaved_forest", "type5_2"]

# 1. Demo: naturskog in Moss -------------------
kommune = kommuner[1]
ecosystem = okosystemer[1]

# A. Compute inputs: land use index and habitat quality ================
lui, qual = create_input_raster_layers(; kommune, ecosystem)
plot(lui)
plot(qual)

write("input/" * ecosystem * "_" * kommune * "_lui.tif", lui, force=true)
write("input/" * ecosystem * "_" * kommune * "_qualityBuffered.tif", qual, force=true)

msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(qual), buff=100.0) 
qual = crop(qual, to = msk)
qual = mask(qual, with = msk, missingval=NaN)
write("input/" * ecosystem * "_" * kommune * "_quality.tif", qual, force=true)

# B. Compute connected habitat and corridors ================
maps = compute_maps(; kommune, ecosystem, ref_value=false, n_threads=Threads.nthreads())
maps = mask_maps(; maps, kommune)
plot(maps[:ch])
plot(maps[:betk])

write("output/" * ecosystem * "_" * kommune * "_funksjonelthabitat.tif", maps[:ch], force=true)
write("output/" * ecosystem * "_" * kommune * "_korridorer.tif", maps[:betk], force=true)

# C. Compute reference connected habitat ================
maps = compute_maps(; kommune, ecosystem, ref_value=true, n_threads=Threads.nthreads())
maps = mask_maps(; maps, kommune)
plot(maps[:ch])

write("output/" * ecosystem * "_" * kommune * "_funksjonelthabitat_ref.tif", maps[:ch], force=true)


# 2. Loop over municipalities -------------------
kommuner = ["bodo", "holtalen", "sandefjord", "ullensaker", "nordreisa"]
ecosystem = okosystemer[1]

# compute inputs
for kommune in kommuner
    lui, qual = create_input_raster_layers(; kommune, ecosystem)
    write("input/" * ecosystem * "_" * kommune * "_lui.tif", lui, force=true)
    write("input/" * ecosystem * "_" * kommune * "_qualityBuffered.tif", qual, force=true)

    msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(qual), buff=100.0) 
    qual = crop(qual, to = msk)
    qual = mask(qual, with = msk, missingval=NaN)
    write("input/" * ecosystem * "_" * kommune * "_quality.tif", qual, force=true)
end

# compute connected habitat and corridors
for kommune in kommuner
    maps = compute_maps(; kommune, ecosystem, ref_value=false, n_threads=Threads.nthreads())
    maps = mask_maps(; maps, kommune)
    # plot(maps[:ch])
    write("output/" * ecosystem * "_" * kommune * "_funksjonelthabitat.tif", maps[:ch], force=true)
    write("output/" * ecosystem * "_" * kommune * "_korridorer.tif", maps[:betk], force=true)

    maps = compute_maps(; kommune, ecosystem, ref_value=true, n_threads=Threads.nthreads())
    maps = mask_maps(; maps, kommune)
    # plot(maps[:ch])
    write("output/" * ecosystem * "_" * kommune * "_funksjonelthabitat_ref.tif", maps[:ch], force=true)
end


# 3. Loop over municipalities and ecosystems -------------------
# note the other ecosystems were not computed for all municipalities
kommuner = ["moss", "bodo", "holtalen"]
okosystemer = ["broadleaved_forest", "type5_2"]

# compute inputs
for kommune in kommuner
    for ecosystem in okosystemer
        lui, qual = create_input_raster_layers(; kommune, ecosystem)
        write("input/" * ecosystem * "_" * kommune * "_lui.tif", lui, force=true)
        write("input/" * ecosystem * "_" * kommune * "_qualityBuffered.tif", qual, force=true)

        msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(qual), buff=100.0) 
        qual = crop(qual, to = msk)
        qual = mask(qual, with = msk, missingval=NaN)
        write("input/" * ecosystem * "_" * kommune * "_quality.tif", qual, force=true)
    end
end

# compute connected habitat and corridors
for kommune in kommuner
    for ecosystem in okosystemer
        maps = compute_maps(; kommune, ecosystem, ref_value=false, n_threads=Threads.nthreads())
        maps = mask_maps(; maps, kommune)
        # plot(maps[:ch])
        write("output/" * ecosystem * "_" * kommune * "_funksjonelthabitat.tif", maps[:ch], force=true)
        write("output/" * ecosystem * "_" * kommune * "_korridorer.tif", maps[:betk], force=true)

        maps = compute_maps(; kommune, ecosystem, ref_value=true, n_threads=Threads.nthreads())
        maps = mask_maps(; maps, kommune)
        # plot(maps[:ch])
        write("output/" * ecosystem * "_" * kommune * "_funksjonelthabitat_ref.tif", maps[:ch], force=true)
    end
end
