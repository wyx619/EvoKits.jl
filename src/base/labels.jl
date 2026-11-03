# Shared identifier and unique-label validation for tree and community inputs.
function _label_identifier(value, name)
    (value isa AbstractString || (value isa Integer && !(value isa Bool))) ||
        throw(ArgumentError("$name must be nonempty strings or integer identifiers"))
    label = String(string(value))
    isempty(strip(label)) && throw(ArgumentError("$name must not be empty"))
    return label
end

function _validated_labels(labels, n, prefix)
    labels === nothing && return ["$(prefix)_$i" for i in 1:n]
    labels isa AbstractVector || throw(ArgumentError("$(prefix)_labels must be a vector"))
    length(labels) == n || throw(ArgumentError("$(prefix)_labels length must match the matrix dimension"))
    result = [_label_identifier(label, "$(prefix)_labels") for label in labels]
    length(unique(result)) == n || throw(ArgumentError("$(prefix)_labels must not contain duplicates"))
    return result
end
