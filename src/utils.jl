"""
    _validate_path(path_to_file::String, expected_extension)

Check that the extension of `path_to_file` is as expected and that the file at that location is not empty. Throw an error if either of these conditions is not met.
"""
function _validate_path(path_to_file::String, expected_extension::String)

    # Get the extension of the file
    extension = splitext(path_to_file) |> last

    # If expected_extension starts with a dot, remove it
    if startswith(expected_extension, ".")
        expected_extension = expected_extension[2:end]
    end

    # Check that extension is as expected
    if extension != ".$expected_extension"
        throw(DomainError("Invalid file extension: expected '.data', got '.$expected_extension'"))
    end

    # Check that the file is not empty
    if (filesize(path_to_file) == 0)
        throw(DomainError("File is empty: $path_to_file"))
    end

    return nothing
end

"""
    _extract_general_info(list_lines::Vector{String})

Extract the lines that contain "general information" (Philips' words).

These lines always start with either '#' or '.'
"""
function _extract_general_info(list_lines::Vector{String})

    # Extract lines that start with "#" or "."
    is_info_line = line -> startswith(line, "#") || startswith(line, ".")
    general_info = filter(is_info_line, list_lines)

    # Throw an error if no general information is found
    if isempty(general_info)
        error("No general information found in .list file")
    end

    return general_info
end

"""
_extract_attributes(list_lines::Vector{String})

Extract the attributes for each of the _complex data vectors_ from the .list file.

The attributes are found in the lines that do *not* start with "#" or ".". The attributes are first read in as a `Vector` of `Strings`. Then, the attributes are converted into a `DataFrame`.

## Implementation details
Using DataFrame() to read the attributes into a DataFrame was not successful because it required parsing all the values to their appropriate types. Instead, we use CSV.read(),
which takes care of the parsing for us, but it requires a file-like object. We can use IOBuffer to create a file-like object from the attribute lines. The lines need to be cleaned up first through.
"""
function _extract_attributes(list_lines)

    # These are lines that do *not* start with "#" or "."
    is_attributes_line = line -> !(startswith(line, "#") || startswith(line, "."))
    attributes_lines = filter(is_attributes_line, list_lines)

    # Throw an error if no attributes are found
    if isempty(attributes_lines)
        error("No attributes found in .list file")
    end

    # The attributes are stored as a Vector{String}.
    # Now we convert them into a DataFrame.

    # Extract the names of the attributes from the .list
    # file to be used as header in the DataFrame
    attribute_names = _get_attributes_header(list_lines)

    # Remove leading whitespace from each line
    attributes_lines = lstrip.(attributes_lines, ' ')

    # Replace all remaining whitespace with commas
    attributes_lines = replace.(attributes_lines, r"\s+" => ",")

    # Join the attributes into a single string
    attributes_str = join(attributes_lines, "\n")

    # Create a file-like IOBuffer from the string
    attributes_io = IOBuffer(attributes_str)

    # Finally, create the DataFrame
    attributes_df = CSV.read(
        attributes_io, # source
        DataFrame, # destination type
        header=attribute_names, # column names
        normalizenames=true # normalize column names (i.e. remove dots and spaces, etc.)
    )

    return attributes_df
end

"""
    _get_attributes_header(list_lines::Vector{String})

Extracts the names of the "attributes" from the .list file.

It searches for the line containing "START OF DATA VECTOR INDEX" and extracts the attribute names from the line two lines below that. The attribute names are stored as a vector of strings.

The attribute names are intended to be used to create column names of a DataFrame holding the attributes.
"""
function _get_attributes_header(list_lines::Vector{String})

    # Look for the line containing "START OF DATA VECTOR INDEX"
    i = findfirst(contains("START OF DATA VECTOR INDEX"), list_lines)

    if isnothing(i)
        error("Could not find 'START OF DATA VECTOR INDEX' in .list file")
    end

    # Two lines below the "START OF DATA VECTOR INDEX" line
    # is the line containing the names of the attributes.
    ATTRIBUTES_HEADER_OFFSET = 2
    attributes_header = list_lines[i+ATTRIBUTES_HEADER_OFFSET]

    # Remove the "# " from the beginning of the line
    attributes_header = replace(attributes_header, r"#\s+" => "")
    # Remove remaining whitespace and replace with commas
    attributes_header = replace(attributes_header, r"\s+" => ",")
    # Split into a vector of strings
    attributes_header = split(attributes_header, ',')
    # Convert SubString{String} to String
    attributes_header = String.(attributes_header)

    return attributes_header
end

"""
    _split_attributes_per_type(attributes::DataFrame)

Split the attributes DataFrame into a NamedTuple of DataFrames for each type.

# Arguments
- `attributes::DataFrame`: The DataFrame containing the attributes data.

# Returns
- `NamedTuple`: A NamedTuple where each field corresponds to a type and contains a DataFrame with the attributes of that type.
"""
function _split_attributes_per_type(attributes)

    # Helper function to select rows from `attributes` based on the type
    _filter_type(type) = filter(row -> row.typ .== type, attributes)

    # Initialize empty dict
    attributes_per_type = Dict()

    # Create dict entry for each type and fill with corresponding attributes
    for type in COMPLEX_DATA_VECTOR_TYPES
        attributes_per_type[Symbol(type)] = _filter_type(type)
    end

    # Convert to NamedTuple
    attributes_per_type = NamedTuple(attributes_per_type)

    return attributes_per_type
end

"""
    _num_bytes_to_num_samples(bytes::Int)

Converts the number of `bytes` to the number of samples based on the size of the raw data type (`ComplexF32`).
"""
function _num_bytes_to_num_samples(bytes::Int)
    # Check that the number of bytes is divisible by the size of the raw data type
    if bytes % sizeof(COMPLEX_ELTYPE) != 0
        error("bytes is not divisible by sizeof(COMPLEX_ELTYPE), result may not be an integer")
    end
    return bytes ÷ sizeof(COMPLEX_ELTYPE)
end

"""
    total_num_samples(type::Union{Symbol,String}, attributes::DataFrame)

Calculate the total number of samples for a given `type` from the `:size` column of the `attributes` DataFrame.
"""
function _total_num_samples(type::Union{Symbol,String}, attributes::DataFrame)
    # Calculate total number of bytes from the size column of attributes
    total_bytes = sum(attributes[attributes.typ.==String(type), :size])
    # Convert bytes to samples
    _num_bytes_to_num_samples(total_bytes)
end

"""
    _preallocate_samples(attributes::DataFrame)

For each of the different types of "complex data vectors", we will store the samples in a separate array. We will also "preallocate" the arrays with `sizehint!` to avoid resizing them each time we `append!` new samples.

# Returns
- `samples_per_type::NamedTuple`: A NamedTuple where each field corresponds to a type and contains a preallocated array for the samples of that type.
"""
function _preallocate_samples(attributes::DataFrame)

    # Initialize an empty dictionary
    samples_per_type = Dict{Symbol,Any}()

    # Fill the dictionary with empty arrays for each type
    for type in COMPLEX_DATA_VECTOR_TYPES
        samples_per_type[Symbol(type)] = ComplexF32[]
    end

    # Convert the dictionary to a NamedTuple
    samples_per_type = NamedTuple(samples_per_type)

    # Preallocate the arrays with sizehint!
    for (type, samples_array) in pairs(samples_per_type)
        sizehint!(samples_array, _total_num_samples(type, attributes))
    end

    return samples_per_type
end

"""
    _read_and_store_samples_per_type!(
        samples_per_type::NamedTuple,
        path_to_datafile::String,
        attributes::DataFrame)

Reads "complex data vectors" from a .data file and stores them in the corresponding array in the NamedTuple.

# Arguments
- `samples_per_type::NamedTuple`: A NamedTuple where each field corresponds to a type and contains a preallocated array for the samples of that type.
- `path_to_datafile::String`: The path to the .data file.
- `attributes::DataFrame`: The DataFrame containing the attributes data.

# Returns
- This function does not return anything. It modifies the `samples_per_type` NamedTuple in-place.
"""
function _read_and_store_samples_per_type!(
    samples_per_type::NamedTuple,
    path_to_datafile::String,
    attributes::DataFrame)

    @info "Reading in .data file"
    
    open(path_to_datafile, "r") do datafile

        for row in ProgressBar(eachrow(attributes))
            # Determine size of current "complex data vector" in bytes
            num_bytes_to_read = row.size
            # Read in raw bytes (read uses UInt8 by default)
            raw_bytes = read(datafile, num_bytes_to_read, all=false)
            # Reinterpret as ComplexF32
            complex_samples = reinterpret(COMPLEX_ELTYPE, raw_bytes)
            # Determine type of current "complex data vector"
            type = Symbol(row.typ)
            # Append the samples to the array of the corresponding type
            append!(samples_per_type[type], complex_samples)
        end

        # Make sure we reached the end of file. Otherwise we might have
        # missed some samples.
        if !eof(datafile)
            error("Did not reach end of file")
        end
    end
end

"""
    _offset_and_size_to_range(offset::Int, size::Int)

The .list file contains `offset` and `size` of each complex data vector in the .data file (in bytes).
This function converts the `offset` and `size` to a range to extract the complex data vector
from the `samples` vector that is returned by the `read_data_list` function.
"""
function _offset_and_size_to_range(offset::Int, size::Int)
    start = _num_bytes_to_num_samples(offset) + 1 # correct for 0-based indexing
    stop = start + _num_bytes_to_num_samples(size)
    return start:stop
end

