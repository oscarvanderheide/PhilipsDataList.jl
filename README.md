# PhilipsDataList

[![Build Status](https://github.com/Oscar/PhilipsDataList.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Oscar/PhilipsDataList.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package contains basic functionally to read in *.{data,list} files exported using Philips' MR systems.

## Installation

Activate the environment where you want to use this package, enter Pkg mode and add the package with the following command:
```julia
add git@gitlab.op.umcutrecht.nl:computational-imaging-lab/philipsdatalist.git
```

## Package Functionality

#### `read_data_list(path::String)`
- This function reads samples from the .data file. There are a limited number of _complex data vector types_ (e.g. `STD`, `NOI`). For each type, the samples are stored separately in a `Vector{ComplexF32}`. The vectors with samples are stored together in a `NamedTuple` with fieldnames corresponding to the different types (e.g. `STD`, `NOI`).
- The general scan information is extracted from the .list file and stored as `info::Vector{String}`.
- The `attributes` of each of the _complex data vectors_ is extracted from the .list file. For each _complex data vector type_ the attributes are stored in a `DataFrame` and the `DataFrames` are stored in a `NamedTuple`.
- The function returns `samples_per_type`, `attributes_per_type` and `info`.

This package only really _reads_ the .{data,list} files and it does not _process_ (e.g. sort) the data. 

#### `data_list_to_kspace(path::String)`

- This function reads in the .{data,list} files using `read_data_list` and then it sorts the samples of type `STD` into a k-space. 
- The k-space has named dimensions. The dimension names are found in `DIMENSIONS_STD`. By default, dimensions of size 1 are dropped. 
- The `kspace` array is also an `OffsetArray`.  This allows, for example, the k-space center(s) to be extracted as `kspace[kx=0, ky=0, kz=0]`. After FFT-ing, it make sense to get rid of the offsets.
- For highly undersampled scans, this approach is memory-inefficient because non-sampled readouts will be stored as zeros.


## Background

The .data file is a binary file that contain measured samples in (complex) single precision floating point format. 

For further processing of the data (e.g. sorting), note that the samples are grouped as _complex data vectors_ (Philips' phrasing). There are several complex data vector types:  
```
STD = Standard data vector (image data or spectroscopy data)
REJ = Rejected standard data vector
    (only for scans with arrhythmia rejection)
PHX = Correction data vector for EPI/GraSE phase correction
FRX = Correction data vector for frequency spectrum correction
NOI = Preparation phase data vector for noise determination
NAV = Phase navigator data vector
DNA = Dynamic phase navigator vector
```
Besides the type, there are several other attributes that _complex data vectors_ have:
```
mix      = mixed sequence number
dyn      = dynamic scan number
card     = cardiac phase number
echo     = echo number
loca     = location number
chan     = synco channel number
extr1    = extra attribute 1 (semantics depend on type of scan)
extr2    = extra attribute 2 (   ''       ''   ''  ''  ''  '' )
ky,kz    = k-space coordinates in 1st and 2nd preparation direction (image data)
kx,ky,kz = k-space coordinates in 1st, 2nd and 3rd preparation direction (spectroscopy data)
aver     = sequence number of this signal average
sign     = sign of measurement gradient used for this data vector (1 = positive, -1 = negative)
rf       = sequence number of this rf echo (only for TSE, TFE, GraSE)
grad     = sequence number of this gradient echo (only for EPI/GraSE)
enc      = encoding time (only for EPI/GraSE)
rtop     = R-top offset in ms
rr       = RR-interval length in ms
size   = data vector size   in number of bytes (1 complex element = 2 floats = 8 bytes)
offset = data vector offset in number of bytes (first data vector starts at offset 0)
```

The .list file contains the attributes for each of the complex data vectors and which samples belong to which _complex data vector_ (using `size` and `offset`). It also contains some _general information_ of the scan such as the k-space and image space coordinate ranges, the number of locations, dynamics and averages.

## Todo
- Store `info` in a `DataFrame` as well.
- Make a faster version in case all the `STD` samples are contiguous.

