"""
    _validate_path(path_to_file::String, expected_extension)

Check that the extension of `path_to_file` is as expected and that the file at that location
is not empty. Throw an error if either of these conditions is not met.
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
        throw(DomainError("Invalid file extension: expected '.data',
        got '.$expected_extension'"))
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

The attributes are found in the lines that do *not* start with "#" or ".". The attributes
are first read in as a `Vector` of `Strings`. Then, the attributes are converted into a
`DataFrame`.

## Implementation details
Using DataFrame() to read the attributes into a DataFrame was not successful because it
required parsing all the values to their appropriate types. Instead, we use CSV.read(),
which takes care of the parsing for us, but it requires a file-like object. We can use
IOBuffer to create a file-like object from the attribute lines. The lines need to be cleaned
up first through.
"""
function _extract_attributes(list_lines)

    # These are lines that do *not* start with "#" or "."
    is_attributes_line = line -> !(startswith(line, "#") || startswith(line, "."))
    attributes_lines = filter(is_attributes_line, list_lines)

    # Throw an error if no attributes are found
    if isempty(attributes_lines)
        error("No attributes found in .list file")
    end

    # The attributes are stored as a Vector{String}. Now we convert them into a DataFrame.

    # Extract the names of the attributes from the .list file to be used as header in the
    # DataFrame
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

It searches for the line containing "START OF DATA VECTOR INDEX" and extracts the attribute
names from the line two lines below that. The attribute names are stored as a vector of
strings.

The attribute names are intended to be used to create column names of a DataFrame holding
the attributes.
"""
function _get_attributes_header(list_lines::Vector{String})

    # Look for the line containing "START OF DATA VECTOR INDEX"
    i = findfirst(contains("START OF DATA VECTOR INDEX"), list_lines)

    if isnothing(i)
        error("Could not find 'START OF DATA VECTOR INDEX' in .list file")
    end

    # Two lines below the "START OF DATA VECTOR INDEX" line is the line containing the names
    # of the attributes.
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
- `NamedTuple`: A NamedTuple where each field corresponds to a type and contains a DataFrame
  with the attributes of that type.
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

Converts the number of `bytes` to the number of samples based on the size of the raw data
type (`ComplexF32`).
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

Calculate the total number of samples for a given `type` from the `:size` column of the
`attributes` DataFrame.
"""
function _total_num_samples(type::Union{Symbol,String}, attributes::DataFrame)
    # Calculate total number of bytes from the size column of attributes
    total_bytes = sum(attributes[attributes.typ.==String(type), :size])
    # Convert bytes to samples
    _num_bytes_to_num_samples(total_bytes)
end

"""
    _preallocate_samples(attributes::DataFrame)

For each of the different types of "complex data vectors", we will store the samples in a
separate array. We will also "preallocate" the arrays with `sizehint!` to avoid resizing
them each time we `append!` new samples.

# Returns
- `samples_per_type::NamedTuple`: A NamedTuple where each field corresponds to a type and
  contains a preallocated array for the samples of that type.
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

Reads "complex data vectors" from a .data file and stores them in the corresponding array in
the NamedTuple.

# Arguments
- `samples_per_type::NamedTuple`: A NamedTuple where each field corresponds to a type and
  contains a preallocated array for the samples of that type.
- `path_to_datafile::String`: The path to the .data file.
- `attributes::DataFrame`: The DataFrame containing the attributes data.

# Returns
- This function does not return anything. It modifies the `samples_per_type` NamedTuple
  in-place.
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

        # Make sure we reached the end of file. Otherwise we might have missed some samples.
        if !eof(datafile)
            @warn "Did not reach end of file"
        end
    end
end

"""
    _offset_and_size_to_range(offset::Int, size::Int)

The .list file contains `offset` and `size` of each complex data vector in the .data file
(in bytes). This function converts the `offset` and `size` to a range to extract the complex
data vector from the `samples` vector that is returned by the `read_data_list` function.
"""
function _offset_and_size_to_range(offset::Int, size::Int)
    start = _num_bytes_to_num_samples(offset) + 1 # correct for 0-based indexing
    stop = start + _num_bytes_to_num_samples(size)
    return start:stop
end

"""
    _remove_empty_fields(samples_per_type::NamedTuple, attributes_per_type::NamedTuple)

Remove the fields of both `samples_per_type` and `attributes_per_type` for which there are
no samples/rows in the DataFrame.
"""
function _remove_empty_fields(samples_per_type::NamedTuple, attributes_per_type::NamedTuple)

    # Check that both NamedTuples have the same fields before filtering
    @assert propertynames(samples_per_type) == propertynames(attributes_per_type)

    # Create Dicts that will hold the non-empty fields only
    samples_dict = Dict{Symbol,Any}()
    attributes_dict = Dict{Symbol,Any}()

    for (key, value) in pairs(samples_per_type)
        if !isempty(value)
            samples_dict[key] = value
        end
    end

    for (key, value) in pairs(attributes_per_type)
        if !isempty(value)
            attributes_dict[key] = value
        end
    end

    # Convert dictionaries back to NamedTuples
    samples_per_type_filtered = NamedTuple(samples_dict)
    attributes_per_type_filtered = NamedTuple(attributes_dict)

    # Check that both NamedTuples have the same fields after filtering
    @assert propertynames(samples_per_type_filtered) ==
            propertynames(attributes_per_type_filtered)

    return samples_per_type_filtered, attributes_per_type_filtered
end

function _repair_list_lines(list_lines::Vector{String})

    nchan = _get_num_coil_channels(list_lines)
    start_idx, end_idx = _get_data_vector_index_bounds(list_lines)

    prefix_lines = list_lines[1:start_idx]
    index_lines = list_lines[start_idx+1:end_idx-1]
    suffix_lines = list_lines[end_idx:end]

    repaired_index_lines, did_repair = _repair_data_vector_index_lines(index_lines, nchan)

    if !did_repair
        return list_lines, false
    end

    return vcat(prefix_lines, repaired_index_lines, suffix_lines), true
end

function _write_list_lines(path::String, lines::Vector{String})
    open(path, "w") do io
        for (index, line) in enumerate(lines)
            if index > 1
                write(io, "\n")
            end
            write(io, line)
        end
    end
end

function _repair_list_file_if_needed(path_to_list_file::String, list_lines::Vector{String})
    _needs_list_repair(list_lines) || return list_lines

    repaired_lines, did_repair = _repair_list_lines(list_lines)
    did_repair || return list_lines

    backup_path = "$path_to_list_file.corrupt"
    if isfile(backup_path)
        error("Backup file already exists: $backup_path")
    end

    temp_path, temp_io = mktemp(dirname(path_to_list_file))
    close(temp_io)
    moved_original = false

    try
        _write_list_lines(temp_path, repaired_lines)
        mv(path_to_list_file, backup_path)
        moved_original = true
        mv(temp_path, path_to_list_file; force=true)
    catch
        isfile(temp_path) && rm(temp_path; force=true)
        if moved_original && !isfile(path_to_list_file) && isfile(backup_path)
            mv(backup_path, path_to_list_file; force=true)
        end
        rethrow()
    end

    return repaired_lines
end

function _needs_list_repair(list_lines::Vector{String})
    release_major = _try_get_gyroscan_release_major(list_lines)
    return !isnothing(release_major) && release_major >= 12
end

function _try_get_gyroscan_release_major(list_lines::Vector{String})
    for line in list_lines
        if occursin("Gyroscan SW release", line)
            match_result = match(r":\s*(\d+)", line)
            if isnothing(match_result)
                error("Could not parse Gyroscan SW release from .list file")
            end
            return parse(Int, match_result.captures[1])
        end
    end

    return nothing
end

function _get_gyroscan_release_major(list_lines::Vector{String})
    release_major = _try_get_gyroscan_release_major(list_lines)
    isnothing(release_major) && error("Could not find Gyroscan SW release in .list file")
    return release_major
end

function _get_num_coil_channels(list_lines::Vector{String})
    for line in list_lines
        if occursin("number of coil channels", line)
            match_result = match(r":\s*(\d+)\s*$", line)
            if isnothing(match_result)
                error("Could not parse number of coil channels from .list file")
            end
            return parse(Int, match_result.captures[1])
        end
    end
    error("Could not find number of coil channels in .list file")
end

function _get_data_vector_index_bounds(list_lines::Vector{String})
    start_idx = findfirst(contains("START OF DATA VECTOR INDEX"), list_lines)
    end_idx = findfirst(contains("END OF DATA VECTOR INDEX"), list_lines)

    if isnothing(start_idx) || isnothing(end_idx) || end_idx <= start_idx
        error("Could not determine data vector index bounds in .list file")
    end

    return start_idx, end_idx
end

function _repair_data_vector_index_lines(index_lines::Vector{String}, nchan::Int)

    entries = [_parse_list_index_line(line) for line in index_lines]
    repaired_entries = Any[]
    inserted_bytes = 0
    did_repair = false

    for i in eachindex(entries)
        entry = entries[i]
        entry isa String && push!(repaired_entries, entry)
        entry isa NamedTuple || continue

        corrected_entry = _shift_offset(entry, inserted_bytes)
        push!(repaired_entries, corrected_entry)

        next_index = _find_next_entry(entries, i + 1)
        if isnothing(next_index)
            continue
        end

        next_entry = entries[next_index]
        @assert next_entry isa NamedTuple

        if _should_insert_missing_std_row(corrected_entry, next_entry, nchan)
            missing_entry = _make_missing_std_row(corrected_entry, next_entry, nchan)
            push!(repaired_entries, missing_entry)
            inserted_bytes += missing_entry.size
            did_repair = true
        end
    end

    if !did_repair
        return index_lines, false
    end

    return [_format_list_index_entry(entry) for entry in repaired_entries], true
end

function _parse_list_index_line(line::String)
    stripped = strip(line)
    if isempty(stripped) || startswith(stripped, "#") || startswith(stripped, ".")
        return line
    end

    fields = split(stripped)
    if length(fields) != 21
        error("Unexpected number of columns ($(length(fields))) in data vector index line: $line")
    end

    numeric = parse.(Int, fields[2:end])

    return (
        typ=fields[1],
        mix=numeric[1],
        dyn=numeric[2],
        card=numeric[3],
        echo=numeric[4],
        loca=numeric[5],
        chan=numeric[6],
        extr1=numeric[7],
        extr2=numeric[8],
        ky=numeric[9],
        kz=numeric[10],
        na=numeric[11],
        aver=numeric[12],
        sign=numeric[13],
        rf=numeric[14],
        grad=numeric[15],
        enc=numeric[16],
        rtop=numeric[17],
        rr=numeric[18],
        size=numeric[19],
        offset=numeric[20],
    )
end

function _format_list_index_entry(entry::String)
    return entry
end

function _format_list_index_entry(entry::NamedTuple)
    return lpad(entry.typ, 5) *
           lpad(string(entry.mix), 6) *
           lpad(string(entry.dyn), 6) *
           lpad(string(entry.card), 6) *
           lpad(string(entry.echo), 6) *
           lpad(string(entry.loca), 6) *
           lpad(string(entry.chan), 6) *
           lpad(string(entry.extr1), 6) *
           lpad(string(entry.extr2), 6) *
           lpad(string(entry.ky), 6) *
           lpad(string(entry.kz), 6) *
           lpad(string(entry.na), 6) *
           lpad(string(entry.aver), 6) *
           lpad(string(entry.sign), 6) *
           lpad(string(entry.rf), 6) *
           lpad(string(entry.grad), 6) *
           lpad(string(entry.enc), 6) *
           lpad(string(entry.rtop), 6) *
           lpad(string(entry.rr), 6) *
           lpad(string(entry.size), 7) *
           " " * string(entry.offset)
end

function _find_next_entry(entries::Vector, start_index::Int)
    for i in start_index:length(entries)
        entries[i] isa NamedTuple && return i
    end
    return nothing
end

function _shift_offset(entry::NamedTuple, inserted_bytes::Int)
    return merge(entry, (offset=entry.offset + inserted_bytes,))
end

function _should_insert_missing_std_row(current::NamedTuple, next::NamedTuple, nchan::Int)
    current.typ == "STD" || return false
    next.typ == "STD" || return false
    expected_next = mod(current.chan + 1, nchan)
    missing_next = mod(current.chan + 2, nchan)

    if _same_std_group_except_channel(current, next)
        return next.chan == missing_next && next.chan != expected_next
    end

    return _is_missing_end_of_std_group(current, next, nchan) ||
           _is_missing_start_of_std_group(current, next, nchan)
end

function _same_std_group_except_channel(a::NamedTuple, b::NamedTuple)
    return a.mix == b.mix &&
           a.dyn == b.dyn &&
           a.card == b.card &&
           a.echo == b.echo &&
           a.loca == b.loca &&
           a.extr1 == b.extr1 &&
           a.extr2 == b.extr2 &&
           a.ky == b.ky &&
           a.kz == b.kz &&
           a.na == b.na &&
           a.aver == b.aver &&
           a.sign == b.sign &&
           a.rf == b.rf &&
           a.grad == b.grad &&
           a.enc == b.enc &&
           a.rtop == b.rtop &&
           a.rr == b.rr &&
           a.size == b.size
end

function _make_missing_std_row(current::NamedTuple, next::NamedTuple, nchan::Int)
    if _same_std_group_except_channel(current, next)
        missing_chan = current.chan + 1
        current.chan < next.chan || (missing_chan = 0)
        return merge(current, (chan=missing_chan, offset=current.offset + current.size,))
    elseif _is_missing_end_of_std_group(current, next, nchan)
        return merge(current, (chan=current.chan + 1, offset=current.offset + current.size,))
    elseif _is_missing_start_of_std_group(current, next, nchan)
        return merge(next, (chan=0, offset=current.offset + current.size,))
    end

    error("Could not reconstruct missing STD row")
end

function _is_missing_end_of_std_group(current::NamedTuple, next::NamedTuple, nchan::Int)
    return current.chan == nchan - 2 && next.chan == 0
end

function _is_missing_start_of_std_group(current::NamedTuple, next::NamedTuple, nchan::Int)
    return current.chan == nchan - 1 && next.chan == 1
end
