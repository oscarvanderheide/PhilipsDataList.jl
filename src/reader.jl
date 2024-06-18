"""
    read_data_list(path_to_data_or_list::String)

The main entry point for reading *.{data,list} files from Philips' MR systems.

This function only really _reads_ in the data but does not sort it into, for example, k-spaces.

# Arguments
- `path_to_data_or_list::String`: The path to the data and list files (without extension).

# Returns
- `samples_per_type::NamedTuple`: A NamedTuple with each field being an array of samples corresponding to one of the the types of "complex data vectors" (see `COMPLEX_DATA_VECTOR_TYPES`).
- `attributes_per_type::NamedTuple`: A NamedTuple with each field being a DataFrame containing attributes corresponding to one of the the types of "complex data vectors" (see `COMPLEX_DATA_VECTOR_TYPES`).
- `info::Vector{String}`: General information about the scan.

# Notes
The .data and .list file should have the same name (except for the extension). The `path` is not required to have an extension since this function will append .data and .list to the path to read the respective files.
"""
function read_data_list(path_to_data_or_list::String)

    # Remove extension from path if any
    path, _ = splitext(path_to_data_or_list)

    # Read in the .list file to get a DataFrame with attributes and general info
    attributes, info = _read_list_file("$path.list")

    # Preallocate arrays for samples of different types and store them in a NamedTuple
    samples_per_type = _preallocate_samples(attributes)

    # Read in samples from .data file and store them in the pre-allocated arrays
    _read_and_store_samples_per_type!(samples_per_type, "$path.data", attributes)

    # Split the attributes DataFrame into a DataFrame for each type
    attributes_per_type = _split_attributes_per_type(attributes)

    return samples_per_type, attributes_per_type, info
end

"""
    _read_list_file(path_to_list_file::String)

First, read a .list file into a vector of strings, each representing a line in the .list file.

The .list file contains (in Philips' words):

1. "general information"

    Examples of general are the maximum kx, ky and kz indices, oversampling factors and the number of channels.

    The general information is found in lines that start with "#" or ".".

2. "attributes" of the "complex data vectors" in the corresponding .data file.

    For example, each "complex data vector" could be a readout and the attributes tell which channel and ky index the readout corresponds to, how many bytes of data it contains and at what offset in the .data file it starts.

    The attributes are found in lines that do *not* start with "#" or ".".

# Returns
- `attributes::DataFrame`: A DataFrame containing the attributes of the "complex data vectors" in the .data file.
- `general_info::Vector{String}`: A vector of strings, each representing a line with general information contained in the .list file.
"""
function _read_list_file(path_to_list_file::String)

    # Validate that the file extension is .list and that the file is not empty
    _validate_path(path_to_list_file, ".list")

    # Open the .list file and store each line as a String into a Vector
    list_lines = readlines(path_to_list_file)

    # Extract the lines with general scan information at the start of the list file
    general_info = _extract_general_info(list_lines)

    # Extract attributes of the "complex data vectors" and store as DataFrame
    attributes = _extract_attributes(list_lines)

    return attributes, general_info
end