"""
    read_data_list(path_to_data_or_list::String)

This function is the main entry point for reading *.{data,list} files from Philips' MR systems. 

It reads measured samples (ComplexF32) from a .data file and stores them in a NamedTuple where each key corresponds to a type of "complex data vector" and the value is an array of samples of that type. 

General scan information and attributes of all the "complex data vectors" are read from the list file. The attributes are stored in a DataFrame and also separated by type into a Dict of DataFrames.

# Important
This function only really reads in the data but does not necessarily put it into a useful format for further processing. For that, users should implement their own routines that make use of the data, attributes and general info returned by this function.

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
    samples_per_type = preallocate_samples(attributes)

    # Read in samples from .data file and store them in the pre-allocated arrays
    read_and_store_samples_per_type!(samples_per_type, "$path.data", attributes)

    # Split the attributes DataFrame into a DataFrame for each type
    attributes_per_type = split_attributes_per_type(attributes)

    return samples_per_type, attributes_per_type, info
end
