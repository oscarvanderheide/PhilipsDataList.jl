"""
    _read_data_file(path::String)

Reads a .data file from the specified path and returns all the samples
as a vector of ComplexF32. Note that it contains both noise samples as well as
data samples for each coil and requires additional processing
based on information in the corresponding .list file.

# Arguments
- `path::String`: The path to the data file.

# Returns
- `samples::Vector{ComplexF32}`: The noise and data samples read from the data file.

# Errors
- Throws an error if the file extension is not '.data'.
- Throws an error if the file cannot be opened or read.
"""
function _read_data_file(path::String)

    # Check that extension is .data
    extension = splitext(path) |> last

    if extension != ".data"
        error("Invalid file extension: expected '.data', got '$extension'")
    end

    # Check that the file is not empty
    if (filesize(path) == 0)
        error("File is empty: $path")
    end

    # Allocate a Vector{COMPLEX_ELTYPE} (= Vector{ComplexF32})
    # to hold all the acquired samples during the scan
    num_samples = _num_bytes_to_num_samples(filesize(path))
    samples = Vector{COMPLEX_ELTYPE}(undef, num_samples)

    # Open the .data file and read in all the samples
    try
        println("Reading in samples from .data file")
        read!(path, samples)
    catch e
        error("Failed to read in complex data samples from .data file: $path")
        rethrow(e)
    end

    return samples
end

"""
    _read_list_file(path::String)

First, read a .list file into a vector of strings, each representing a line in the .list file.

The .list file contains (in Philips' words):

1. "general information"

    Examples of general are the maximum kx, ky and kz indices, oversampling factors and the number of channels.

    The general information is found in lines that start with "#" or ".".

2. "attributes" of the "complex data vectors" in the corresponding .data file.

    For example, each "complex data vector" could be a readout and the attributes tell which channel and ky index the readout corresponds to, how many bytes of data it contains and at what offset in the .data file it starts.

    The attributes are found in lines that do *not* start with "#" or ".".

# Arguments
- `path::String`: The path to the .list file.

# Returns
- `attributes::DataFrame`: A DataFrame containing the attributes of the "complex data vectors" in the .data file.
- `general_info::Vector{String}`: A vector of strings, each representing a line with general information contained in the .list file.


# Errors
- Throws an error if the file extension is not .list.
- Throws an error if the file cannot be opened.
- Throws an error if general information lines are found.
- Throws an error if no attributes lines are found.

# Example
```julia
attributes, general_info = _read_list_file("/path/to/file.list")
```
"""
function _read_list_file(path::String)

    # Check that extension is .list
    extension = splitext(path) |> last

    if extension != ".list"
        error("Invalid file extension: expected '.list', got '$extension'")
    end

    # Open the .list file and store each line as a String into a Vector
    list_lines = try
        readlines(path)
    catch e
        error("Failed to open .list file:\n$path")
        rethrow(e)
    end

    # Extract the lines that contain "general information" (Philips' words).
    # These lines always start with either '#' or '.'
    is_info_line = line -> startswith(line, "#") || startswith(line, ".")
    general_info = filter(is_info_line, list_lines)

    if isempty(general_info)
        error("No general information found in .list file")
    end

    # Now extract the "attributes" from the .list file
    # These are lines that do *not* start with "#" or "."
    is_attributes_line = line -> !(startswith(line, "#") || startswith(line, "."))
    attributes = filter(is_attributes_line, list_lines)

    if isempty(attributes)
        error("No attributes found in .list file")
    end

    # At this point, the attributes are stored as a vector of strings.
    # Now convert it to a DataFrame to make it easier to work with the attributes later on

    # Extract the names of the attributes from the .list
    # file to be used as header in the DataFrame
    header = _get_attributes_header(general_info)

    attributes = _attributes_lines_to_df(attributes, header)

    return attributes, general_info
end

"""
    _attributes_lines_to_df(attributes_lines)

Extracts and stores the attribute lines from a .list file in a DataFrame.

# Arguments
- `attributes_lines`: An array of strings representing the lines of a .list file containing attributes of the "complex data vectors".

# Returns
- `attributes_df`: A DataFrame containing the attribute for all of the "complex data vectors" in a .data file.

The function first extracts the attribute names from the list file to be used as the header in the DataFrame.
Then, it filters out the lines containing the attribute values.
Next, it cleans up the lines by removing leading whitespace and replacing remaining whitespace with commas.
After that, it joins the cleaned attribute lines into a single string.
Finally, it creates a file-like IOBuffer from the string and uses CSV.read() to parse the attributes into a DataFrame.
"""
function _attributes_lines_to_df(attributes_lines, attribute_names)

    # Using DataFrame() to read the attributes into a DataFrame
    # was not successful because it required parsing all the values
    # to their appropriate types. Instead, we use CSV.read(),
    # which takes care of the parsing for us, but it requires a
    # file-like object. We can use IOBuffer to create a file-like object
    # from the attribute lines. The lines need to be cleaned up first through.

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

Extracts the names of the "attributes" from the .list file. It searches for the line
containing "START OF DATA VECTOR INDEX" and extracts the attribute names from
the line two lines below that. The attribute names are stored as a vector of strings.

The attribute names are intended to be used to create column names of a DataFrame holding the attributes.

# Arguments
- `list_lines::Vector{String}`: A vector of strings representing the lines in a .list file.

# Returns
- `attribute_names::Vector{String}`: A vector of strings representing the attribute names.

# Errors
- Throws an error if the 'START OF DATA VECTOR INDEX' is not found in the list file.

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
    split_attributes_per_type(attributes::DataFrame)

Split the attributes DataFrame into a NamedTuple of DataFrames for each type.

# Arguments
- `attributes::DataFrame`: The DataFrame containing the attributes data.

# Returns
- `NamedTuple`: A NamedTuple where each field corresponds to a type and contains a DataFrame with the attributes of that type.
"""
function split_attributes_per_type(attributes)

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

Converts the number of bytes to the number of samples based on the size of the raw data type (COMPLEX_ELTYPE = ComplexF32).

# Arguments
- `bytes::Int`: The number of bytes.

# Returns
The number of samples.

# Errors
- Throws an error if `bytes` is not divisible by the size of the raw data type.

"""
function _num_bytes_to_num_samples(bytes::Int)
    # check that the number of bytes is divisible by the size of the raw data type
    if bytes % sizeof(COMPLEX_ELTYPE) != 0
        error("bytes is not divisible by sizeof(COMPLEX_ELTYPE), result may not be an integer")
    end
    return bytes ÷ sizeof(COMPLEX_ELTYPE)
end

"""
    total_num_samples(type, attributes)

Calculate the total number of samples for a given type from the `:size` column of the `attributes` DataFrame.

# Arguments
- `type::Union{Symbol,String}`: The type for which to calculate the total number of samples.
- `attributes::DataFrame`: The DataFrame containing the attributes data.

# Returns
- `Int`: The total number of samples for the given type.
"""
function _total_num_samples(type::Union{Symbol,String}, attributes::DataFrame)
    # Calculate total number of bytes from the size column of attributes
    total_bytes = sum(attributes[attributes.typ.==String(type), :size])
    # Convert bytes to samples
    _num_bytes_to_num_samples(total_bytes)
end

"""
    preallocate_samples(attributes::DataFrame)

For each of the different types of "complex data vectors",
we will store the samples in a separate array. We will also
"preallocate" the arrays with `sizehint!` to avoid resizing them
each time we `append!` new samples. 

# Arguments
- `attributes::DataFrame`: The DataFrame containing the attributes data.

# Returns
- `NamedTuple`: A NamedTuple where each field corresponds to a type and contains a preallocated array for the samples of that type.
"""
function preallocate_samples(attributes)

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
    read_and_store_samples_per_type!(
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
function read_and_store_samples_per_type!(
    samples_per_type::NamedTuple,
    path_to_datafile::String,
    attributes::DataFrame)

    open(path_to_datafile, "r") do datafile

        for row in eachrow(attributes)
            # Determine size of current "complex data vector" in bytes
            num_bytes_to_read = row.size
            # Read in raw bytes (read uses UInt8 by default)
            raw_bytes = read(datafile, num_bytes_to_read)
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

The .list file contains `offset` and `size` of each complex data vector in the .data file.
This function converts the `offset` and `size` to a range to extract the complex data vector
from the `samples` vector that is returned by the `read_data_list` function.

# Arguments
- `offset::Int`: The offset in bytes.
- `size::Int`: The size in bytes.

# Returns
A `UnitRange` representing the range of indices corresponding to the offset and size.
"""
function _offset_and_size_to_range(offset::Int, size::Int)
    start = _num_bytes_to_num_samples(offset) + 1 # correct for 0-based indexing
    stop = start + _num_bytes_to_num_samples(size)
    return start:stop
end

