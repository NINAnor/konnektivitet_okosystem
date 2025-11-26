using Pkg
#Pkg.add(url = "https://github.com/ConScape/ConScape.jl.git", rev = "dev")

using Rasters, Shapefile, LibGEOS, ArchGDAL
using GeoInterface, GeometryBasics#, CoordinateTransformations
using CairoMakie
import GeoDataFrames as GDF
using GeometryOps
using Statistics
using Rasters: Center

using ConScape#dev install from dev branch
using SparseArrays, LinearAlgebra
using StatsBase, Optim

function asym_transf(x, mid, scl)
    y = (1 ./ (1 .+ exp.(.- ((x .- mid) .* scl))))
    mn = (1/(1 + exp(mid*scl)))
    y = (y .- mn) ./ (1-mn)
    return y
end

function municipal_mask(; kommune = "moss", tgt_crs, buff=0.0)
   # * Municipal boundary ==============
    if kommune == "moss"
        kommunenummer = "3002" # Moss
    elseif kommune == "bodo"
        kommunenummer = "1804" # Bodo
    elseif kommune == "holtalen"
        kommunenummer = "5026" # Holtalen
    elseif kommune == "sandefjord"
        kommunenummer = "3804" # Sandefjord
    elseif kommune == "ullensaker"
        kommunenummer = "3033" # Ullensaker
    else
        kommunenummer = "5428" # Nordreisa
    end


    df = GDF.read(joinpath(@__DIR__, "data","Norge_Kommuner_2023.shp"))
    kommune_shp = df[df.kommunenum .== kommunenummer,:]

    kommune_shp = GeometryOps.reproject(kommune_shp, crs(kommune_shp), tgt_crs)
    if buff > 0
        msk = buffer.(kommune_shp.geometry, buff)
    else
        msk = kommune_shp.geometry
    end
    return msk
end

# kommunenummer = "3002" # Moss
# ecosystem = "broadleaved_forest"

solver   = ConScape.VectorSolver()

# see report for description
theta=1.0
optimal_r = 88 # the distance conversion rate: converts the cost distance from ConScape into meter equivalent -- which, for ease of interpretation, allows alpha to be expressed in Euclidean distance
alpha = 1/2500 # distance decay rate for the Euclidean distance
Alpha = alpha * optimal_r # the actual distance decay used in ConScape for the cost (not Euclidean) distance.
res = 50
radius = Int(10_000 / res)
centersize = 20

params = Dict(
    :solver    => solver,
    :res       => res,
    :theta     => theta,
    :radius    => radius,
    :alpha     => alpha,
    :optimal_r => optimal_r,
    :Alpha     => Alpha
)

movement_mode = RandomisedShortestPath(ExpectedCost();
    distance_transformation=x -> exp(- Alpha * x),
    theta=theta
)

function create_input_raster_layers(; kommune = "moss", ecosystem = "broadleaved_forest")
    # * Land use intensity index ===============
    lui = Raster(joinpath(@__DIR__, "data", "ABI25x100.tif"), missingval=NaN, lazy=true) ./ 100
    
    msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(lui), buff=10_000.0)
    lui = crop(lui, to=msk)
    lui = mask(lui, with=msk, missingval=NaN)

    lui = aggregate(Center(), lui, 2)
    lui = map(x -> asym_transf(x, 4, 0.75), lui)
    #plot(lui)

    # * Ecosystem =============
    if ecosystem == "broadleaved_forest"
        raster_path = "data/Municipality_" * kommune * "_type_4_1_binary.tif"
    elseif ecosystem == "naturskog"
        raster_path = joinpath(@__DIR__, "data", "naturskog_v1_naturskognaerhet.tif")
    else
        raster_path = "data/Municipality_" * kommune * "_type_5_2_binary.tif"
    end   
    es = Raster(raster_path, lazy=true, missingval=NaN)

    msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(es), buff=11_000.0)
    es = crop(es, to=msk)

    if ecosystem == "naturskog"
        es = map(x -> isnan(x) ? 0 : x, es)
        max_es = 7.0

        es = resample(es, to=lui, method=:average)
        es = es ./ max_es 
    else
        es = resample(es, to=lui, method=:near)
    end

    es = crop(es, to=lui)
    es = mask(es, with=lui, missingval=NaN)
    #plot(es)

    qualities = es
    qualities = map(x -> isnan(x) ? 0 : x, qualities)
    qualities = Float64.(qualities)

    return lui, qualities
end

function create_input_stack(; kommune = "moss", ecosystem = "broadleaved_forest", return_sparse=false, ref_value=false)
    # * Land use intensity index ===============
    lui = Raster("input/" * ecosystem * "_" * kommune * "_lui.tif", missingval=NaN, lazy=false)

    # * Ecosystem =============
    qualities = Raster("input/" * ecosystem * "_" * kommune * "_qualityBuffered.tif", missingval=NaN, lazy=false)

    if ref_value && (ecosystem == "naturskog")
        map!(x -> x>0 ? 1 : 0, qualities, qualities)
    end

    # * Maps ==============

    # ideally, but I don't think that it works with a source-based functionality -- code needs updating 
    # msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(qualities), buff=100.0) 
    # So, until then we do a bit too much computing as to avoid edge effects:
    # msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(qualities), buff=10_000.0) 
    # tgts = mask(qualities, with = msk, missingval=NaN)
    tgts = copy(qualities)

    step_likl = (0.5 .+ 0.49 .* qualities)

    if !ref_value
        step_likl = (1 .- lui) .* step_likl
    end

    if return_sparse
        rast = RasterStack((; affinities=step_likl, source_qualities=qualities, target_qualities=modify(sparse, tgts)))
    else
        rast = RasterStack((; affinities=step_likl, source_qualities=qualities, target_qualities=tgts))
    end
    return rast
end

function compute_maps(; kommune = "moss", ecosystem = "broadleaved_forest", ref_value=false, n_threads=1)
    if ref_value
        measures = (; ch=ConnectedHabitat())
    else
        measures = (;
            ch=ConnectedHabitat(),
            betk=Betweenness(QualityAndProximityWeighted()),
        )
    end

    if n_threads==1
        rast = create_input_stack(; kommune, ecosystem, return_sparse=true, ref_value)

        problem = ConScape.Problem(; measures, movement_mode, solver=ConScape.VectorSolver())
        maps = ConScape.solve(problem, rast)
    else
        rast = create_input_stack(; kommune, ecosystem, return_sparse=false, ref_value)

        problem = ConScape.Problem(; measures, movement_mode, solver=ConScape.VectorSolver())
        windowed_problem = ConScape.WindowedProblem(problem; buffer=radius, centersize=20, threaded = true)
        maps = ConScape.solve(windowed_problem, rast)
    end

    return maps     
end

function mask_maps(; maps, kommune= "moss")
    msk = municipal_mask(; kommune, tgt_crs=Rasters.crs(maps), buff=100.0)
    maps = crop(maps, to = msk)
    maps = mask(maps, with = msk, missingval=NaN)

    return maps     
end

