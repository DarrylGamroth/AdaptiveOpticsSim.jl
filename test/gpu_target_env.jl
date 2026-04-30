let root = dirname(@__DIR__)
    if !any(path -> abspath(path) == root, LOAD_PATH)
        push!(LOAD_PATH, root)
    end
end
