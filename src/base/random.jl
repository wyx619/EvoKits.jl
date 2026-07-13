function stable_grid_seed(base_seed::Integer, grid_id::AbstractString)::UInt64
    value = UInt64(base_seed) ⊻ UInt64(0xcbf29ce484222325)
    for byte in codeunits(grid_id)
        value = (value ⊻ UInt64(byte)) * UInt64(0x100000001b3)
    end
    return value
end
