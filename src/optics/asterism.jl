struct Asterism{S<:AbstractSource,V<:AbstractVector{S}} <: AbstractSource
    sources::V
    function Asterism(sources::AbstractVector{S}) where {S<:AbstractSource}
        frozen = Vector{S}(undef, length(sources))
        @inbounds for i in eachindex(sources)
            require_leaf_source(sources[i], "Asterism child")
            frozen[i] = freeze_source(sources[i])
        end
        return new{S,typeof(frozen)}(frozen)
    end
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
