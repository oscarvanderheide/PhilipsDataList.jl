"""
    to_kspace(path_to_data_or_list; drop_dims=true)

Read in the .{data/list} files and store the measured "STD" samples in a k-space.

## Note
- It is assumed that each readout has the same number of samples.
- The .data and .list file should have the same name (except for the extension). The `path` is not required to have an extension since this function will append .data and .list to the path to read the respective files.
- The k-space is stored as an OffsetArray s.t. we can use the ranges as indices. For example, `kspace[0, 0, 0, ...]` will give the k-space value(s) at the center of k-space.
- The k-space also has named dimensions (it is a `NamedDimsArray`), which in fact allows us to use the dimension names as indices. For example, `kspace[kx=0, ky=0, kz=0, ...]` will give the k-space value(s) at the center of k-space. 
- The ordering of the dimensions is the same as in the `DIMENSIONS_STD` constant. If `drop_dims=true`, dimensions of size 1 are dropped.
"""
function to_kspace(path_to_data_or_list; drop_dims=true)

    # Read in the data and list files
    samples_per_type, attributes_per_type, _ = read_data_list(path_to_data_or_list)

    # Get the samples and attributes for the "STD" type
    samples = samples_per_type.STD
    attributes = attributes_per_type.STD

    # Sort the data based on the attributes into a k-space
    kspace = _to_kspace(samples, attributes)

    # Drop dimensions of size 1
    if drop_dims
        kspace = squeeze(kspace)
    end
    return kspace
end

"""
    _to_kspace(samples::Vector{ComplexF32}, attributes::DataFrame)

Assemble a k-space as an OffsetArray with named dimensions (with names `DIMENSIONS_STD`) from the measured samples and attributes.

## Note
- The `attributes` DataFrame should only contain attributes of type "STD".
- Each readout is assumed to have the same number of samples per readout.

## Implementation details
- The k-space will be a high-dimensional array with dimensions `DIMENSIONS_STD`. 
- The ranges of the k-space dimensions are calculated from the attributes. 
- From the ranges, the sizes are calculated and then the actual k-space array is allocated. 
- We use an `OffsetArray` s.t. we can use the ranges as indices.
- We then loop over each readout (e.g. row of the `attributes` DataFrame) and extract the k-space location of the readout from the attributes. We then store the readout in the k-space at the correct location.
- The `kx` dimension is a special case. The `kx` dimension is not stored in the attributes but can be calculated from the size of the readout.
- The k-space is wrapped in a `NamedDimsArray` to add dimension names.
"""
function _to_kspace(samples::Vector{ComplexF32}, attributes::DataFrame)

    # Check that all attributes are of type "STD"
    @assert all(attributes[!, :typ] .== "STD")

    # Check that all readouts have the same number of samples
    @assert _all_readouts_same_size(attributes)

    # Get ranges of all the different attributes (dimensions) in the data
    ranges = _calculate_kspace_dimension_ranges(attributes)

    # Preallocate k-space
    kspace = _allocate_kspace(ranges)

    # Fill the k-space with the measured STD samples
    _fill_kspace!(kspace, samples, attributes)

    return kspace
end

"""
    _all_readouts_same_size(attributes::DataFrame)

Check that all readouts have the same number of samples
"""
function _all_readouts_same_size(attributes::DataFrame)
    return length(unique(attributes[!, :size])) == 1
end

"""
    _calculate_kspace_dimension_ranges(attributes::DataFrame)

Calculate range for each k-space dimension and store in a NamedTuple
"""
function _calculate_kspace_dimension_ranges(attributes::DataFrame)
    NamedTuple{DIMENSIONS_STD}(
        _calculate_range(dim, attributes) for dim in DIMENSIONS_STD
    )
end

"""
    _calculate_range(dimension::Symbol, attributes::DataFrame) 

Calculate range for a single k-space dimension
"""
function _calculate_range(dimension::Symbol, attributes::DataFrame)

    # Special case for kx
    dimension == :kx && (return _calculate_range_kx(attributes))

    # Get the range of the specified dimension
    minimum(attributes[!, dimension]):maximum(attributes[!, dimension])

end

"""
    _calculate_range_kx(attributes::DataFrame)

Calculate range for the kx dimension which is assumed to be the same for all readouts.
"""
function _calculate_range_kx(attributes::DataFrame)

    # Get the size in bytes of a readout
    size_readout = attributes[begin, :size]

    # Convert to number of samples per readout
    samples_per_readout = _num_bytes_to_num_samples(size_readout)

    return -(samples_per_readout ÷ 2):((samples_per_readout÷2)-1)
end

"""
    _allocate_kspace(ranges::NamedTuple{DIMENSIONS_STD})

Allocate k-space based on the `ranges` for each dimension in k-space
"""
function _allocate_kspace(ranges::NamedTuple{DIMENSIONS_STD})

    # From ranges, calculate sizes
    sizes = [length(range) for range in values(ranges)]

    @show sizes
    # Allocate k-space
    kspace = zeros(ComplexF32, sizes...)

    # Wrap in an OffsetArray using the ranges
    kspace = OffsetArray(kspace, ranges...)

    # Add dimension names
    kspace = NamedDimsArray(kspace, propertynames(ranges))

    return kspace
end

"""
    _fill_kspace!(kspace, samples::Vector{ComplexF32}, attributes::DataFrame)

Fill `k-space` with measured `samples` with the k-space locations extracted from the `attributes`.
"""
function _fill_kspace!(kspace, samples::Vector{ComplexF32}, attributes::DataFrame)

    # Reshape samples into a matrix
    samples = reshape(samples, size(kspace, :kx), :)

    # Validate that num readouts in data is the same as num readouts in attributes
    @assert nrow(attributes) == size(samples, 2)

    # Loop over each row of the attributes DataFrame to extract k-space location of each readout
    kx_range = axes(kspace, :kx)

    for i in ProgressBar(1:nrow(attributes))

        # Extract current readout from measured samples
        current_readout = OffsetVector((@view samples[:, i]), kx_range)
        # Extract k-space location
        idx = _get_kspace_idx(attributes, i)
        # Store current readout in k-space
        kspace[:, idx] .= current_readout
    end

    # TODO: Multiply by :sign value

    return nothing
end

"""
    _get_kspace_idx(attributes::DataFrame, i::Int)

Get index in k-space for a given readout
"""
function _get_kspace_idx(attributes::DataFrame, i::Int)
    CartesianIndex(Tuple(attributes[i, dim] for dim in DIMENSIONS_STD if dim != :kx))
end

"""
    squeeze(x::AbstractArray)

Drop dimensions of size 1 from an array `x`.
"""
squeeze(x::AbstractArray) = dropdims(x; dims=tuple(findall(size(x) .== 1)...))