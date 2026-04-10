using Test
using PhilipsDataList

function format_fixture_row(typ, mix, dyn, card, echo, loca, chan, extr1, extr2, ky, kz,
    na, aver, sign, rf, grad, enc, rtop, rr, size, offset)
    return join((typ, mix, dyn, card, echo, loca, chan, extr1, extr2, ky, kz, na, aver,
        sign, rf, grad, enc, rtop, rr, size, offset), " ")
end

function write_fixture_list(path, rows; release="12.1-1")
    lines = String[
        "# Synthetic Philips list fixture",
        "# Gyroscan SW release                  :  $release",
        ".    0    0    0  number of coil channels            :    4",
        "# === START OF DATA VECTOR INDEX ===============================================",
        "#",
        "# typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  enc   rtop  rr    size   offset",
        "# --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ------ ------",
        "#",
    ]

    append!(lines, rows)

    append!(lines, [
        "#",
        "# === END OF DATA VECTOR INDEX ================================================",
        "#",
        "# === END OF DATA DESCRIPTION FILE ============================================",
    ])

    open(path, "w") do io
        for (index, line) in enumerate(lines)
            if index > 1
                write(io, "\n")
            end
            write(io, line)
        end
    end
end

@testset "repair_list_file skips older releases" begin
    mktempdir() do dir
        list_path = joinpath(dir, "synthetic.list")
        original_rows = [
            format_fixture_row("STD", 0, 0, 0, 0, 0, 0, 0, 0, 10, 1, 0, 0, 1, 0, 0, 0, 0, 0, 8, 0),
            format_fixture_row("STD", 0, 0, 0, 0, 0, 2, 0, 0, 10, 1, 0, 0, 1, 0, 0, 0, 0, 0, 8, 8),
        ]

        write_fixture_list(list_path, original_rows; release="5.7")
        original_contents = read(list_path, String)

        repaired_path = repair_list_file(list_path)

        @test repaired_path == list_path
        @test read(list_path, String) == original_contents
        @test !isfile("$list_path.corrupt")
    end
end

@testset "repair_list_file skips healthy release 12 files" begin
    mktempdir() do dir
        list_path = joinpath(dir, "synthetic.list")

        rows = String[]
        offset = 0

        for ky in 10:11
            for chan in 0:3
                push!(rows, format_fixture_row("STD", 0, 0, 0, 0, 0, chan, 0, 0, ky, 1, 0, 0, 1, 0, 0, 0, 0, 0, 8, offset))
                offset += 8
            end
        end

        push!(rows, format_fixture_row("NOI", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 16, offset))

        write_fixture_list(list_path, rows)
        original_contents = read(list_path, String)

        repaired_path = repair_list_file(list_path)

        @test repaired_path == list_path
        @test read(list_path, String) == original_contents
        @test !isfile("$list_path.corrupt")
    end
end

@testset "read_data_list repairs corrupted list automatically" begin
    mktempdir() do dir
        list_path = joinpath(dir, "synthetic.list")
        data_path = joinpath(dir, "synthetic.data")

        rows = String[]
        offset = 0
        size_std = 8

        for chan in (0, 1, 3)
            push!(rows, format_fixture_row("STD", 0, 0, 0, 0, 0, chan, 0, 0, 10, 1, 0, 0, 1, 0, 0, 0, 0, 0, size_std, offset))
            offset += size_std
        end

        for chan in 0:3
            push!(rows, format_fixture_row("STD", 0, 0, 0, 0, 0, chan, 0, 0, 11, 1, 0, 0, 1, 0, 0, 0, 0, 0, size_std, offset))
            offset += size_std
        end

        for chan in 0:2
            push!(rows, format_fixture_row("STD", 0, 0, 0, 0, 0, chan, 0, 0, 12, 1, 0, 0, 1, 0, 0, 0, 0, 0, size_std, offset))
            offset += size_std
        end

        for chan in 0:3
            push!(rows, format_fixture_row("STD", 0, 0, 0, 0, 0, chan, 0, 0, 13, 1, 0, 0, 1, 0, 0, 0, 0, 0, size_std, offset))
            offset += size_std
        end

        for chan in 1:3
            push!(rows, format_fixture_row("STD", 0, 0, 0, 0, 0, chan, 0, 0, 14, 1, 0, 0, 1, 0, 0, 0, 0, 0, size_std, offset))
            offset += size_std
        end

        push!(rows, format_fixture_row("NOI", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 16, offset))

        write_fixture_list(list_path, rows)

        open(data_path, "w") do io
            write(io, zeros(UInt8, 176))
        end

        corrupt_path = "$list_path.corrupt"
        samples_per_type, attributes, _ = read_data_list(list_path)

        @test isfile(list_path)
        @test isfile(corrupt_path)
        @test length(samples_per_type.STD) == 20

        @test attributes.STD.chan == repeat(collect(0:3), 5)

        repaired_lines = readlines(list_path)
        start_idx = findfirst(contains("START OF DATA VECTOR INDEX"), repaired_lines)
        end_idx = findfirst(contains("END OF DATA VECTOR INDEX"), repaired_lines)

        data_rows = filter(line -> begin
            stripped = strip(line)
            !isempty(stripped) && !startswith(stripped, "#") && !startswith(stripped, ".")
        end, repaired_lines[start_idx+1:end_idx-1])

        parsed_rows = split.(strip.(data_rows))
        parsed_std = filter(fields -> fields[1] == "STD", parsed_rows)
        parsed_noi = filter(fields -> fields[1] == "NOI", parsed_rows)

        @test length(parsed_std) == 20
        @test parse.(Int, getindex.(parsed_std, 7)) == repeat(collect(0:3), 5)
        @test parse(Int, last(last(parsed_std))) == 152
        @test parse(Int, last(only(parsed_noi))) == 160
    end
end

@testset "repair_list_file refuses existing backup" begin
    mktempdir() do dir
        list_path = joinpath(dir, "synthetic.list")

        write_fixture_list(list_path, [
            format_fixture_row("STD", 0, 0, 0, 0, 0, 0, 0, 0, 10, 1, 0, 0, 1, 0, 0, 0, 0, 0, 8, 0),
            format_fixture_row("STD", 0, 0, 0, 0, 0, 2, 0, 0, 10, 1, 0, 0, 1, 0, 0, 0, 0, 0, 8, 8),
        ])

        open("$list_path.corrupt", "w") do io
            write(io, "already backed up")
        end

        error_text = ""
        try
            repair_list_file(list_path)
            @test false
        catch error
            error_text = sprint(showerror, error)
        end

        @test occursin("Backup file already exists", error_text)
        @test read(list_path, String) == join(readlines(list_path), "\n")
        @test read("$list_path.corrupt", String) == "already backed up"
    end
end