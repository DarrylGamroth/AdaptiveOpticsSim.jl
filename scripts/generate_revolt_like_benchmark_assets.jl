using DelimitedFiles

const OUTPUT_DIR = joinpath(dirname(@__DIR__), "benchmarks", "assets", "revolt_like")

function revolt_like_mask(n::Int=19)
    radius = 9.35
    center = (n + 1) / 2
    mask = falses(n, n)
    @inbounds for j in 1:n, i in 1:n
        x = i - center
        y = j - center
        mask[i, j] = x^2 + y^2 <= radius^2
    end
    count(mask) == 277 || error("expected 277 active actuators, got $(count(mask))")
    return mask
end

function revolt_like_map(mask::BitMatrix)
    n_row, n_col = size(mask)
    grid = zeros(Int, n_row, n_col)
    idx = 1
    @inbounds for j in 1:n_col, i in 1:n_row
        if mask[i, j]
            grid[i, j] = idx
            idx += 1
        end
    end
    return grid
end

function nearest_interior_neighbor(mask::BitMatrix, i::Int, j::Int)
    n_row, n_col = size(mask)
    best = nothing
    best_score = typemax(Int)
    @inbounds for dj in -1:1, di in -1:1
        (di == 0 && dj == 0) && continue
        ii = i + di
        jj = j + dj
        1 <= ii <= n_row || continue
        1 <= jj <= n_col || continue
        mask[ii, jj] || continue
        neighbor_count = 0
        for odj in -1:1, odi in -1:1
            (odi == 0 && odj == 0) && continue
            iii = ii + odi
            jjj = jj + odj
            if 1 <= iii <= n_row && 1 <= jjj <= n_col && mask[iii, jjj]
                neighbor_count += 1
            end
        end
        score = -neighbor_count
        if score < best_score
            best_score = score
            best = (ii, jj)
        end
    end
    isnothing(best) && error("no interior neighbor found for edge actuator at ($i, $j)")
    return best
end

function revolt_like_extrapolation(map_grid::Matrix{Int}, mask::BitMatrix)
    n = maximum(map_grid)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    @inbounds for j in axes(map_grid, 2), i in axes(map_grid, 1)
        label = map_grid[i, j]
        label == 0 && continue
        neighbor_count = 0
        for dj in -1:1, di in -1:1
            (abs(di) + abs(dj) != 1) && continue
            ii = i + di
            jj = j + dj
            if 1 <= ii <= size(mask, 1) && 1 <= jj <= size(mask, 2) && mask[ii, jj]
                neighbor_count += 1
            end
        end
        if neighbor_count >= 4
            push!(rows, label - 1)
            push!(cols, label - 1)
            push!(vals, 1.0)
        else
            ni, nj = nearest_interior_neighbor(mask, i, j)
            neighbor_label = map_grid[ni, nj]
            push!(rows, label - 1)
            push!(cols, label - 1)
            push!(vals, 0.5)
            push!(rows, label - 1)
            push!(cols, neighbor_label - 1)
            push!(vals, 0.5)
        end
    end
    length(rows) >= n || error("expected at least one entry per row")
    return hcat(rows, cols, vals)
end

function write_map(path::AbstractString, grid::Matrix{Int})
    open(path, "w") do io
        println(io, "0 0 $(size(grid, 1)) $(size(grid, 2))")
        for i in axes(grid, 1)
            println(io, join(grid[i, :], ","))
        end
    end
end

function write_extrapolation(path::AbstractString, triplets::Matrix{Float64}, n::Int)
    open(path, "w") do io
        println(io, "0 0 $n $n")
        for row in eachrow(triplets)
            println(io, "$(Int(row[1])),$(Int(row[2])),$(row[3])")
        end
    end
end

function main()
    mkpath(OUTPUT_DIR)
    mask = revolt_like_mask()
    grid = revolt_like_map(mask)
    triplets = revolt_like_extrapolation(grid, mask)
    write_map(joinpath(OUTPUT_DIR, "revolt_like_dmActuatorMap_277.csv"), grid)
    write_extrapolation(joinpath(OUTPUT_DIR, "revolt_like_dmExtrapolation.csv"), triplets, maximum(grid))
    println("generated revolt-like assets in ", OUTPUT_DIR)
end

main()
