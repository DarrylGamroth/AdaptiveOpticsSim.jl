struct PETypeDeclaration
    path::String
    name::String
    declaration::String
    mutable::Bool
    field_count::Int
    line::Int
end

struct PEDebtSite
    path::String
    category::String
    scope::String
    binding::String
    line::Int
end

const PE_VOCABULARY = Set((
    "definition",
    "definitions",
    "params",
    "parameters",
    "plan",
    "prepared",
    "state",
    "workspace",
    "product",
    "products",
    "strategy",
))

const PE_NONSTORAGE_SIGNATURE_ROOTS = Set((
    :Tuple,
    :NamedTuple,
    :Vararg,
    :Type,
    :Val,
    :Union,
))

const PE_STORAGE_TYPE_ROOTS = Set((
    :Array,
    :Vector,
    :Matrix,
    :Memory,
    :GenericMemory,
    :AtomicMemory,
    :AbstractArray,
    :AbstractVector,
    :AbstractMatrix,
    :DenseArray,
    :DenseVector,
    :DenseMatrix,
    :FixedSizeArray,
    :FixedSizeVector,
    :FixedSizeMatrix,
    :Ref,
    :RefValue,
    :Dict,
    :IdDict,
    :WeakKeyDict,
    :Set,
))

const PE_IMPLICIT_ANY_EMPTY_CONSTRUCTORS = Set((
    :Dict,
    :IdDict,
    :WeakKeyDict,
    :Set,
))

function pe_source_paths(root::AbstractString)
    paths = String[]
    for relative_root in ("src", "ext")
        directory = joinpath(root, relative_root)
        isdir(directory) || continue
        for (parent, _, files) in walkdir(directory)
            append!(paths, joinpath(parent, file) for file in files
                if endswith(file, ".jl"))
        end
    end
    return sort!(paths)
end

pe_normalize_repository_path(path::AbstractString) =
    replace(String(path), '\\' => '/')

function pe_camel_tokens(name::AbstractString)
    tokens = String[]
    for component in split(strip(name, '_'), '_')
        for matched in eachmatch(
            r"[A-Z]+(?=[A-Z][a-z]|\d|$)|[A-Z]?[a-z]+|\d+",
            component,
        )
            push!(tokens, lowercase(matched.match))
        end
    end
    return tokens
end

pe_has_vocabulary(name::AbstractString) =
    !isdisjoint(Set(pe_camel_tokens(name)), PE_VOCABULARY)

function pe_type_name(expression)
    expression isa Symbol && return expression
    expression isa QuoteNode && expression.value isa Symbol &&
        return expression.value
    expression isa Expr || return nothing
    if expression.head in (:curly, :(<:), :(>:), :(::))
        return pe_type_name(first(expression.args))
    elseif expression.head == :.
        return pe_type_name(last(expression.args))
    end
    return nothing
end

function pe_declaration_name(expression)
    expression isa Symbol && return expression
    expression isa Expr || return nothing
    if expression.head in (:curly, :(<:), :(::))
        return pe_declaration_name(first(expression.args))
    end
    return nothing
end

function pe_is_type_alias_rhs(expression)
    expression isa Symbol && return true
    expression isa Expr || return false
    return expression.head in (:curly, :where, :., :(<:), :(>:))
end

function pe_struct_field(expression)
    value = expression
    while value isa Expr && value.head == :const
        length(value.args) == 1 || return nothing
        value = only(value.args)
    end

    has_default = value isa Expr && value.head == :(=)
    default = has_default ? last(value.args) : nothing
    declaration = has_default ? first(value.args) : value
    if declaration isa Symbol
        return (
            binding=declaration,
            annotation=nothing,
            has_default,
            default,
        )
    elseif declaration isa Expr && declaration.head == :(::) &&
            first(declaration.args) isa Symbol
        return (
            binding=first(declaration.args),
            annotation=last(declaration.args),
            has_default,
            default,
        )
    end
    return nothing
end

function _pe_collect_declarations!(
    declarations::Vector{PETypeDeclaration},
    expression,
    relative_path::String,
    line::Int,
)
    expression isa LineNumberNode && return expression.line
    expression isa Expr || return line

    if expression.head == :struct
        mutable = expression.args[1]::Bool
        name = pe_declaration_name(expression.args[2])
        body = expression.args[3]
        fields = body isa Expr && body.head == :block ?
            count(value -> pe_struct_field(value) !== nothing,
                body.args) : 0
        if name !== nothing && pe_has_vocabulary(String(name))
            push!(declarations, PETypeDeclaration(
                relative_path,
                String(name),
                mutable ? "mutable_struct" : "struct",
                mutable,
                fields,
                line,
            ))
        end
    elseif expression.head in (:abstract, :primitive)
        name = pe_declaration_name(first(expression.args))
        if name !== nothing && pe_has_vocabulary(String(name))
            push!(declarations, PETypeDeclaration(
                relative_path,
                String(name),
                expression.head == :abstract ?
                    "abstract_type" : "primitive_type",
                false,
                -1,
                line,
            ))
        end
    elseif expression.head == :macrocall &&
            pe_type_name(first(expression.args)) === Symbol("@enum")
        name = length(expression.args) >= 3 ?
            pe_declaration_name(expression.args[3]) : nothing
        if name !== nothing && pe_has_vocabulary(String(name))
            push!(declarations, PETypeDeclaration(
                relative_path,
                String(name),
                "enum",
                false,
                -1,
                line,
            ))
        end
    elseif expression.head == :const && length(expression.args) == 1
        assignment = first(expression.args)
        if assignment isa Expr && assignment.head == :(=)
            name = first(assignment.args)
            target = last(assignment.args)
            if name isa Symbol && pe_has_vocabulary(String(name)) &&
                    pe_is_type_alias_rhs(target)
                push!(declarations, PETypeDeclaration(
                    relative_path,
                    String(name),
                    "type_alias",
                    false,
                    -1,
                    line,
                ))
            end
        end
    end

    current_line = line
    for argument in expression.args
        current_line = _pe_collect_declarations!(
            declarations, argument, relative_path, current_line)
    end
    return current_line
end

function pe_type_declarations(root::AbstractString)
    declarations = PETypeDeclaration[]
    for path in pe_source_paths(root)
        relative_path = pe_normalize_repository_path(relpath(path, root))
        expression = Meta.parseall(read(path, String); filename=path)
        _pe_collect_declarations!(declarations, expression, relative_path, 1)
    end
    sort!(declarations; by=value -> (value.path, value.name))
    return declarations
end

function pe_contains_any(expression)
    expression === :Any && return true
    expression isa QuoteNode && return false
    expression isa Expr || return false
    return any(pe_contains_any, expression.args)
end

function pe_is_memory_type(expression)
    expression === :Memory && return true
    expression isa Expr || return false
    expression.head in (:curly, :.) || return false
    return pe_type_name(expression) === :Memory
end

function pe_storage_type_contains_any(expression)
    expression isa Expr && expression.head == :curly || return false
    root = pe_type_name(first(expression.args))
    root in PE_NONSTORAGE_SIGNATURE_ROOTS && return false
    root in PE_STORAGE_TYPE_ROOTS || return false
    return pe_contains_any(expression)
end

function pe_is_implicit_any_constructor(expression)
    expression isa Expr && expression.head == :call || return false
    callee = first(expression.args)
    (callee isa Symbol ||
        (callee isa Expr && callee.head == :.)) || return false
    root = pe_type_name(callee)
    positional_count = 0
    first_positional = nothing
    for argument in Iterators.drop(expression.args, 1)
        argument isa Expr && argument.head == :parameters && continue
        positional_count += 1
        positional_count == 1 && (first_positional = argument)
    end
    root in PE_IMPLICIT_ANY_EMPTY_CONSTRUCTORS &&
        return positional_count == 0
    root === :Vector &&
        return positional_count == 0 || first_positional === :undef
    root === :Matrix && return first_positional === :undef
    return false
end

function pe_expression_binding(expression)
    expression isa Symbol && return String(expression)
    expression isa Expr || return "<expression>"
    if expression.head == :(::)
        return pe_expression_binding(first(expression.args))
    elseif expression.head in (:ref, :.)
        return pe_expression_binding(first(expression.args))
    elseif expression.head == :tuple
        return "<destructured>"
    end
    return "<expression>"
end

function pe_function_name(signature)
    signature isa Symbol && return String(signature)
    signature isa Expr || return "<anonymous>"
    if signature.head in (:where, :(::))
        return pe_function_name(first(signature.args))
    elseif signature.head == :call
        name = first(signature.args)
        symbol = pe_type_name(name)
        return symbol === nothing ? string(name) : String(symbol)
    end
    return "<anonymous>"
end

function _pe_push_site!(
    sites::Vector{PEDebtSite},
    path::String,
    category::String,
    scope::String,
    binding::String,
    line::Int,
)
    push!(sites, PEDebtSite(path, category, scope, binding, line))
    return nothing
end

function _pe_memory_type_sites!(
    sites::Vector{PEDebtSite},
    expression,
    path::String,
    category::String,
    scope::String,
    binding::String,
    line::Int,
)
    expression === :Memory && begin
        _pe_push_site!(sites, path, category, scope, binding, line)
        return
    end
    expression isa Expr || return
    if pe_is_memory_type(expression)
        _pe_push_site!(sites, path, category, scope, binding, line)
        for argument in Iterators.drop(expression.args, 1)
            _pe_memory_type_sites!(sites, argument, path, category,
                scope, binding, line)
        end
        return
    end
    for argument in expression.args
        _pe_memory_type_sites!(sites, argument, path, category,
            scope, binding, line)
    end
    return
end

function _pe_scan_signature!(
    memory_sites::Vector{PEDebtSite},
    any_sites::Vector{PEDebtSite},
    signature,
    path::String,
    line::Int,
)
    function_name = pe_function_name(signature)
    callable = signature
    while callable isa Expr && callable.head == :where
        for constraint in Iterators.drop(callable.args, 1)
            _pe_memory_type_sites!(memory_sites, constraint, path,
                "signature_type", function_name, "<constraint>", line)
        end
        callable = first(callable.args)
    end
    if callable isa Expr && callable.head == :(::)
        return_type = last(callable.args)
        _pe_memory_type_sites!(memory_sites, return_type, path,
            "signature_type", function_name, "<return>", line)
        if pe_storage_type_contains_any(return_type)
            _pe_push_site!(any_sites, path, "return_erased_storage",
                function_name, "<return>", line)
        elseif occursin("prepare", lowercase(function_name)) &&
                pe_contains_any(return_type)
            _pe_push_site!(any_sites, path, "prepared_return_any",
                function_name, "<return>", line)
        end
        callable = first(callable.args)
    end
    callable isa Expr && callable.head == :call || return
    for argument in Iterators.drop(callable.args, 1)
        argument isa Expr || continue
        if argument.head == :parameters
            for keyword in argument.args
                _pe_scan_signature_argument!(memory_sites, any_sites,
                    keyword, path, function_name, line)
            end
        else
            _pe_scan_signature_argument!(memory_sites, any_sites,
                argument, path, function_name, line)
        end
    end
    return
end

function _pe_scan_signature_argument!(
    memory_sites::Vector{PEDebtSite},
    any_sites::Vector{PEDebtSite},
    argument,
    path::String,
    function_name::String,
    line::Int,
)
    value = argument isa Expr && argument.head in (:kw, :(=)) ?
        first(argument.args) : argument
    value isa Expr && value.head == :(... ) &&
        (value = first(value.args))
    value isa Expr && value.head == :(::) || return
    binding = pe_expression_binding(first(value.args))
    annotation = last(value.args)
    _pe_memory_type_sites!(memory_sites, annotation, path,
        "signature_type", function_name, binding, line)
    if pe_storage_type_contains_any(annotation)
        _pe_push_site!(any_sites, path, "signature_erased_storage",
            function_name, binding, line)
    end
    return
end

function _pe_scan_struct_body!(
    memory_sites::Vector{PEDebtSite},
    any_sites::Vector{PEDebtSite},
    body,
    path::String,
    type_name::String,
    line::Int,
)
    body isa Expr && body.head == :block ||
        return _pe_scan_expression!(memory_sites, any_sites, body,
            path, line, type_name, "<expression>")
    current_line = line
    for field in body.args
        if field isa LineNumberNode
            current_line = field.line
            continue
        end
        field_parts = pe_struct_field(field)
        if field_parts !== nothing && field_parts.annotation !== nothing
            binding = String(field_parts.binding)
            annotation = field_parts.annotation
            _pe_memory_type_sites!(memory_sites, annotation, path,
                "field_type", type_name, binding, current_line)
            if pe_contains_any(annotation)
                _pe_push_site!(any_sites, path, "field_stored_any",
                    type_name, binding, current_line)
            end
            if field_parts.has_default
                current_line = _pe_scan_expression!(memory_sites,
                    any_sites, field_parts.default, path, current_line,
                    type_name, binding)
            end
        elseif field_parts !== nothing
            binding = String(field_parts.binding)
            _pe_push_site!(any_sites, path, "field_stored_any",
                type_name, binding, current_line)
            if field_parts.has_default
                current_line = _pe_scan_expression!(memory_sites,
                    any_sites, field_parts.default, path, current_line,
                    type_name, binding)
            end
        else
            current_line = _pe_scan_expression!(memory_sites, any_sites,
                field, path, current_line, type_name, "<expression>")
        end
    end
    return current_line
end


function _pe_scan_call_arguments!(
    memory_sites::Vector{PEDebtSite},
    any_sites::Vector{PEDebtSite},
    expression::Expr,
    path::String,
    line::Int,
    scope::String,
)
    current_line = line
    for argument in Iterators.drop(expression.args, 1)
        current_line = _pe_scan_expression!(memory_sites, any_sites,
            argument, path, current_line, scope, "<argument>")
    end
    return current_line
end

function _pe_scan_expression!(
    memory_sites::Vector{PEDebtSite},
    any_sites::Vector{PEDebtSite},
    expression,
    path::String,
    line::Int,
    scope::String,
    binding::String,
)
    expression isa LineNumberNode && return expression.line
    expression isa Expr || return line

    if expression.head == :struct
        name = something(pe_declaration_name(expression.args[2]),
            Symbol("<anonymous>"))
        _pe_memory_type_sites!(memory_sites, expression.args[2], path,
            "type_constraint", String(name), "<type-parameter>", line)
        return _pe_scan_struct_body!(memory_sites, any_sites,
            expression.args[3], path, String(name), line)
    elseif expression.head in (:abstract, :primitive)
        _pe_memory_type_sites!(memory_sites, first(expression.args), path,
            "type_constraint", string(pe_declaration_name(
                first(expression.args))), "<type-parameter>", line)
        return line
    elseif expression.head == :function
        _pe_scan_signature!(memory_sites, any_sites,
            expression.args[1], path, line)
        length(expression.args) == 1 && return line
        return _pe_scan_expression!(memory_sites, any_sites,
            expression.args[2], path, line,
            pe_function_name(expression.args[1]), "<expression>")
    elseif expression.head == :(=) &&
            (first(expression.args) isa Expr) &&
            first(expression.args).head in (:call, :where, :(::))
        _pe_scan_signature!(memory_sites, any_sites,
            first(expression.args), path, line)
        return _pe_scan_expression!(memory_sites, any_sites,
            last(expression.args), path, line,
            pe_function_name(first(expression.args)), "<expression>")
    elseif expression.head == :(=)
        destination = pe_expression_binding(first(expression.args))
        current_line = _pe_scan_expression!(memory_sites, any_sites,
            first(expression.args), path, line, scope, destination)
        return _pe_scan_expression!(memory_sites, any_sites,
            last(expression.args), path, current_line, scope, destination)
    elseif expression.head == :call &&
            pe_is_memory_type(first(expression.args))
        _pe_push_site!(memory_sites, path, "constructor", scope,
            binding, line)
        pe_storage_type_contains_any(first(expression.args)) &&
            _pe_push_site!(any_sites, path,
                "erased_storage_constructor", scope, binding, line)
        return _pe_scan_call_arguments!(memory_sites, any_sites,
            expression, path, line, scope)
    elseif expression.head == :call &&
            pe_storage_type_contains_any(first(expression.args))
        _pe_push_site!(any_sites, path, "erased_storage_constructor",
            scope, binding, line)
        return _pe_scan_call_arguments!(memory_sites, any_sites,
            expression, path, line, scope)
    elseif pe_is_implicit_any_constructor(expression)
        _pe_push_site!(any_sites, path,
            "implicit_any_storage_constructor", scope, binding, line)
        return _pe_scan_call_arguments!(memory_sites, any_sites,
            expression, path, line, scope)
    elseif expression.head == :ref && !isempty(expression.args) &&
            first(expression.args) === :Any
        _pe_push_site!(any_sites, path, "any_array_literal", scope,
            binding, line)
    elseif expression.head == :vect && isempty(expression.args)
        _pe_push_site!(any_sites, path, "implicit_any_array_literal",
            scope, binding, line)
    elseif expression.head == :(::)
        annotation = last(expression.args)
        _pe_memory_type_sites!(memory_sites, annotation, path,
            "local_type", scope,
            pe_expression_binding(first(expression.args)), line)
        if pe_storage_type_contains_any(annotation)
            _pe_push_site!(any_sites, path, "local_erased_storage",
                scope, pe_expression_binding(first(expression.args)), line)
        end
        return line
    elseif pe_is_memory_type(expression)
        _pe_push_site!(memory_sites, path, "type_reference", scope,
            binding, line)
        return line
    end

    current_line = line
    for argument in expression.args
        current_line = _pe_scan_expression!(memory_sites, any_sites,
            argument, path, current_line, scope, binding)
    end
    return current_line
end

function pe_debt_sites(root::AbstractString)
    memory_sites = PEDebtSite[]
    any_sites = PEDebtSite[]
    for path in pe_source_paths(root)
        relative_path = pe_normalize_repository_path(relpath(path, root))
        expression = Meta.parseall(read(path, String); filename=path)
        _pe_scan_expression!(memory_sites, any_sites, expression,
            relative_path, 1, "<top-level>", "<expression>")
    end
    order(site) = (site.path, site.category, site.scope,
        site.binding, site.line)
    sort!(memory_sites; by=order)
    sort!(any_sites; by=order)
    return (memory=memory_sites, stored_any=any_sites)
end

function pe_debt_counts(sites)
    counts = Dict{Tuple{String,String,String,String},Int}()
    for site in sites
        key = (site.path, site.category, site.scope, site.binding)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end
