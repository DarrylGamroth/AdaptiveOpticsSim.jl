#!/usr/bin/env julia

# Prepare the repository's private, patched CUDACore source without loading CUDA.
using Pkg
using TOML
using UUIDs

const CUDACORE_UUID = UUID("bd0ed864-bdfe-4181-a5ed-ce625a5fdea2")
const CUDACORE_VERSION = v"6.5.0"
const CUDACORE_TREE = "9d504f23d39a56f827d1da4bfd1580990a28adbb"

function verify_source(path, expected_tree)
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    get(project, "uuid", nothing) == string(CUDACORE_UUID) ||
        error("Unexpected CUDACore UUID in $path")
    get(project, "version", nothing) == string(CUDACORE_VERSION) ||
        error("Unexpected CUDACore version in $path")
    actual_tree = bytes2hex(Pkg.GitTools.tree_hash(path))
    actual_tree == expected_tree ||
        error("CUDACore source tree mismatch in $path: expected $expected_tree, got $actual_tree")
    return nothing
end

function prepare_cuda_dependency()
    VERSION >= v"1.13" || error("Preparing CUDACore requires Julia 1.13 or newer")
    root = dirname(@__DIR__)
    deps = joinpath(root, "deps")
    patch = joinpath(deps, "patches", "cudacore-6.5.0.patch")
    metadata = TOML.parsefile(joinpath(deps, "patches", "cudacore-6.5.0.toml"))
    patched_tree = metadata["patched_tree"]
    patched_tree isa String && occursin(r"^[0-9a-f]{40}$", patched_tree) ||
        error("Invalid patched_tree in CUDACore patch metadata")
    destination = joinpath(deps, "cudacore-6.5.0")
    if ispath(destination) || islink(destination)
        verify_source(destination, patched_tree)
        println("Verified patched CUDACore: $destination")
        return nothing
    end

    # Pkg runs in a child process so the caller's active project is untouched.
    mktempdir() do environment
        source_file = joinpath(environment, "source.txt")
        bootstrap = """
        using Pkg, UUIDs
        uuid = UUID($(repr(string(CUDACORE_UUID))))
        Pkg.add(PackageSpec(name="CUDACore", uuid=uuid, version=v"6.5.0"))
        write($(repr(source_file)), Pkg.dependencies()[uuid].source)
        """
        command = `$(Base.julia_cmd()) --startup-file=no --project=$environment -e $bootstrap`
        run(addenv(command, "JULIA_PKG_PRECOMPILE_AUTO" => "0"))
        source = read(source_file, String)
        verify_source(source, CUDACORE_TREE)

        # Stage beside the destination so publishing is a same-filesystem rename.
        mktempdir(deps; prefix=".cudacore-") do staging
            staged_source = joinpath(staging, "source")
            cp(source, staged_source)
            for (directory, _, files) in walkdir(staged_source)
                chmod(directory, filemode(directory) | 0o700)
                for file in files
                    path = joinpath(directory, file)
                    islink(path) || chmod(path, filemode(path) | 0o200)
                end
            end
            # Apply outside the enclosing AOS Git worktree's path prefix.
            git = addenv(`git -C $staged_source`,
                "GIT_CEILING_DIRECTORIES" => staging,
                "GIT_DIR" => nothing, "GIT_WORK_TREE" => nothing)
            run(`$git apply --check $patch`)
            run(`$git apply $patch`)
            verify_source(staged_source, patched_tree)
            mv(staged_source, destination)
        end
    end
    println("Prepared patched CUDACore: $destination")
    return nothing
end

prepare_cuda_dependency()
