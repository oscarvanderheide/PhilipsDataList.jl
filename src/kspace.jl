"""
    data_list_to_kspace(path_to_data_or_list; drop_dims=true, remove_readout_oversampling=false, offset_array=false)

Read in the .{data/list} files and store the measured "STD" samples in a k-space.

## Warning
This function should only be used on data from Cartesian acquisitions. I don't know what exactly happens for non-Cartesian data.

## Note
- It is assumed that each readout has the same number of samples.
- The .data and .list file should have the same name (except for the extension). The `path` is not required to have an extension since this function will append .data and .list to the path to read the respective files.
- The k-space has named dimensions (it is a `NamedDimsArray`) which allows, for example, to 
retrieve all samples from a channel with index `i` by indexing with `kspace[chan=i]`.
- The ordering of the dimensions is the same as in the `DIMENSIONS_STD` constant. However, if `drop_dims=true`, dimensions of size 1 are dropped.
- If `offsetarray` is `true`, the k-space is stored as an OffsetArray s.t. we can use the ranges contained in the list file as indices. For example, `kspace[kx=0, ky=0, kz=0, ...]` will give the k-space value(s) at the center of k-space.
- If `remove_readout_oversampling` is `true`, the readout oversampling is removed by cropping the k-space and applying an ifft along the readout direction. 
"""
function data_list_to_kspace(path_to_data_or_list; drop_dims=true,
    remove_readout_oversampling=false, offset_array=false)

    # Read in the data and list files
    samples_per_type, attributes_per_type, general_info = read_data_list(path_to_data_or_list)

    # Considered only the measured sampled of type "STD"
    samples = samples_per_type.STD

    # Extract the attributes for the "STD" type samples
    attributes = attributes_per_type.STD

    # Sort the data based on the attributes into a k-space
    kspace = _samples_to_kspace(samples, attributes)

    # Drop dimensions of size 1
    drop_dims && (kspace = squeeze(kspace))

    # Remove readout oversampling
    if remove_readout_oversampling

        # Get the oversampling factor
        oversampling_factor = _extract_readout_oversampling_factor(general_info)

        # Remove readout oversampling
        kspace = _remove_readout_oversampling(kspace, oversampling_factor)
    end

    # Remove custom indices (e.g. axes(kspace, :kx) = -128:127) if `offset_array` is false
    if !offset_array
        kspace = _remove_custom_indices_keep_nameddims(kspace)
    end

    return kspace
end

"""
    _samples_to_kspace(samples::Vector{ComplexF32}, attributes::DataFrame)

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
function _samples_to_kspace(samples::Vector{ComplexF32}, attributes::DataFrame)

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

Calculate range for each k-space dimension, including handling the special case for the kx dimension, which is assumed to be the same for all readouts, and store in a NamedTuple.
"""
function _calculate_kspace_dimension_ranges(attributes::DataFrame)

    # Calculate range for the special case kx first

    # Get the size in bytes of a readout
    size_readout = attributes[begin, :size]

    # Convert to number of samples per readout
    samples_per_readout = _num_bytes_to_num_samples(size_readout)

    # Calculate range for kx
    range_kx = -(samples_per_readout ÷ 2):((samples_per_readout÷2)-1)

    # Helper function to calculate range for any dimension
    _calculate_range(dim) = begin
        if dim == :kx
            return range_kx
        else
            range_min = minimum(attributes[!, dim])
            range_max = maximum(attributes[!, dim])
            return range_min:range_max
        end
    end

    # Calculate range for all dimensions and store in a NamedTuple
    NamedTuple{DIMENSIONS_STD}(
        _calculate_range(dim) for dim in DIMENSIONS_STD
    )
end

"""
    _allocate_kspace(ranges::NamedTuple{DIMENSIONS_STD})

Allocate k-space based on the `ranges` for each dimension in k-space. The k-space is stored as an `OffsetArray` s.t. we can use the ranges as indices. The k-space is also wrapped in a `NamedDimsArray` to add dimension names.
"""
function _allocate_kspace(ranges::NamedTuple{DIMENSIONS_STD})

    # From ranges, calculate sizes
    sizes = [length(range) for range in values(ranges)]

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

Fill `k-space` with measured `samples` with the k-space locations extracted from the `attributes`. Because `kspace` is an `OffsetArray`, we can use the ranges obtained from the list file as indices.
"""
function _fill_kspace!(kspace, samples::Vector{ComplexF32}, attributes::DataFrame)

    # Reshape samples into a matrix with dimensions (kx, readouts)
    samples = reshape(samples, size(kspace, :kx), :)

    # Validate that num readouts in data is the same as num readouts in attributes
    @assert nrow(attributes) == size(samples, 2)

    # Loop over each row of the attributes DataFrame to extract k-space location of each readout
    kx_range = axes(kspace, :kx)

    @info "Sorting data into k-space"

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

Get index in k-space for the `i`-th readout.
"""
function _get_kspace_idx(attributes::DataFrame, i::Int)
    CartesianIndex(Tuple(attributes[i, dim] for dim in DIMENSIONS_STD if dim != :kx))
end

"""
    squeeze(x::AbstractArray)

Drop dimensions of size 1 from an array `x`.
"""
squeeze(x::AbstractArray) = dropdims(x; dims=tuple(findall(size(x) .== 1)...))

"""
    _extract_readout_oversampling_factor(general_info::Vector{String})

Retrieve the readout oversampling factor from the part of the list file that looks like this:

"# mix  echo n.a.  k-space oversample factors           value"
"# ---- ---- ----  ----------------------------------   ---------"
".    0    0    0  kx_oversample_factor               :    2.0000"
"""
function _extract_readout_oversampling_factor(general_info::Vector{String})

    # Iterate over the lines in the general_info
    for line in general_info
        # Check if the string contains "kx_oversample_factor"
        if occursin("kx_oversample_factor", line)
            # Use a regex to find the floating point number in the string
            m = match(r"(\d+\.\d+)", line)
            # If a match was found, convert it to a float and then to an integer
            if m !== nothing
                # Read as floating point number
                oversampling_factor = parse(Float64, only(m))
                # Check that it is actually an integer
                @assert oversampling_factor ≈ round(oversampling_factor)
                # Convert to integer
                oversampling_factor = round(Int, oversampling_factor)

                return oversampling_factor
            else
                error("Could not find a floating point number in the string: $str")
            end
            # Exit the loop once the number is found
            break
        end
    end

    error("Could not find the kx_oversample_factor in the general_info")
end

"""
    _remove_readout_oversampling(kspace::NamedDimsArray, oversampling_factor::Int)

Remove readout oversampling from data by fft-ing and cropping.

## Note:
- The `kx` dimension is assumed to be the first dimension. 
- It is also assumed no partial Fourier technique is applied (in the readout direction). 
- The `kspace` is assumed to be an `OffsetArray` with named dimensions. 
- Before applying the fft, the offsets are actually removed because the offsets will cause unwanted phases when applying the fft. After cropping and fft-ing back, offsets are added back.
"""
function _remove_readout_oversampling(kspace::NamedDimsArray, oversampling_factor::Int)

    # Get current number of samples in readout direction
    n = size(kspace, :kx)

    # New kx_max index afte croppin
    kx_max = n ÷ (2 * oversampling_factor)

    # Store dimnames to add them later again
    kspace_dimnames = dimnames(kspace)

    # Store axes to add them later again
    kspace_axes = axes(kspace)

    # Get parent data array without dimension names and offsets
    kspace = _remove_custom_indices_keep_nameddims(kspace)

    # fft to data along readout direction
    tmp = ifftshift(ifft(kspace, :kx), :kx∿) # Made possible by NamedDims

    # Rename fft turns :kx into :kx∿, rename to :x here
    tmp = NamedDims.rename(tmp, :kx∿ => :x)

    # Determine central region that is kept (no offsets here)
    central_region = ((n÷2)-(kx_max)+1):((n÷2)+(kx_max))

    # Crop the data in the `x` domain
    # tmp = tmp[central_region, fill(:, ndims(kspace)-1)...]
    tmp = tmp[x=central_region]

    # ifft back to `kx` domain
    kspace = fft(fftshift(tmp, :x), :x)

    # Rename :x∿ to :kx
    kspace = NamedDims.rename(kspace, :x∿ => :kx)

    # Calculate new kx range (taking into account offset)
    kx_range = -kx_max:(kx_max-1)

    # Add offsets and dimension names
    kspace = OffsetArray(kspace, kx_range, kspace_axes[2:end]...)

    return kspace
end

"""
    _remove_custom_indices_keep_nameddims(kspace)

Remove the custom indices from the k-space array while preserving the named dimensions.
"""
function _remove_custom_indices_keep_nameddims(kspace)

    # Check if kspace is a NamedDimsArray wrapped around an OffsetArray or vice versa
    if kspace isa NamedDimsArray && parent(kspace) isa OffsetArray
        # NamedDimsArray wrapping an OffsetArray
        kspace_dimnames = dimnames(kspace)
        kspace = parent(parent(kspace))  # Get the raw array without offsets and names
    elseif kspace isa OffsetArray && parent(kspace) isa NamedDimsArray
        # OffsetArray wrapping a NamedDimsArray
        kspace_dimnames = dimnames(parent(kspace))
        kspace = parent(parent(kspace))  # Get the raw array without offsets and names
    else
        throw(ArgumentError("kspace must be a NamedDimsArray wrapping an OffsetArray or vice versa"))
    end

    # Reconstruct the NamedDimsArray with original dimension names
    kspace = NamedDimsArray(kspace, kspace_dimnames)

    return kspace
end