"""
An ordered finite collection of leaf sources.

Construction freezes every source and seals the collection at its final
cardinality. Heterogeneous source families use a concrete union element type
rather than abstract source storage.
"""
struct Asterism{V<:FixedSizeVector{<:AbstractSource}} <: AbstractSource
    sources::V
    function Asterism{V}(sources::V) where {
        V<:FixedSizeVector{<:AbstractSource},
    }
        _union_members_are_concrete(eltype(V)) || throw(
            InvalidConfiguration(
                "asterism storage must have concrete source type members",
            ),
        )
        return new{V}(sources)
    end
end

@inline function _freeze_asterism_source(source::AbstractSource)
    require_leaf_source(source, "Asterism child")
    return freeze_source(source)
end

function Asterism(sources::AbstractVector{<:AbstractSource})
    frozen = Tuple(_freeze_asterism_source(source) for source in sources)
    storage = _fixed_size_union_vector(frozen)
    return Asterism{typeof(storage)}(storage)
end

freeze_source(ast::Asterism) = Asterism(ast.sources)

@inline source_composition_style(::Asterism) = ExpandedSourceComposition()

Base.length(ast::Asterism) = length(ast.sources)

function wavelength(ast::Asterism)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    w0 = wavelength(ast.sources[1])
    @inbounds for i in 2:length(ast.sources)
        src = ast.sources[i]
        if wavelength(src) != w0
            throw(InvalidConfiguration("asterism sources must share a common wavelength"))
        end
    end
    return w0
end
