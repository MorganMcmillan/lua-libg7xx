local band, bor, brshift, blshift = bit.band, bit.bor, bit.brshift or bit.rshift, bit.blshift or bit.lshift
local pack = string.pack

--[[
 * decode.c
 *
 * CCITT ADPCM decoder
 *
 * Usage : decode [-3|4|5] [-a|u|l] < infile > outfile
]]

local g72x = require("g72x")
local g721 = require("g721")
local g723_24 = require("g723_24")
local g723_40 = require("g723_40")

local AUDIO_ENCODING_ULAW = 1
local AUDIO_ENCODING_ALAW = 2
local AUDIO_ENCODING_LINEAR = 3

--- prints usage info to stdout and exits
local function print_usage()
    print("CCITT ADPCM Decoder -- usage:")
    print("\tdecode [-3|4|5] [-a|u|l] [infile] [outfile]")
    print("where:")
    print("\t-3\tProcess G.723 24kbps (3-bit) input data")
    print("\t-4\tProcess G.721 32kbps (4-bit) input data [default]")
    print("\t-5\tProcess G.723 40kbps (5-bit) input data")
    print("\t-a\tGenerate 8-bit A-law data")
    print("\t-u\tGenerate 8-bit u-law data [default]")
    print("\t-l\tGenerate 16-bit linear PCM data")
    error()
end

--- Unpack input codes and pass them back as bytes.
--- Returns 1 if there is residual input, returns -1 if eof, else returns 0.
---@param code integer[] pointer to a single integer
---@param bits integer
---@param in_file file*
---@return boolean
local function unpack_input(code, bits, in_file)
    local in_buffer
    in_buffer = in_buffer or 0
    local in_bits
    in_bits = in_bits or 0
    local in_byte

    if in_bits < bits then
        in_byte = in_file:read(1)
        if not in_byte then
            code[1] = 0
            return false
        end
        in_buffer = bor(in_buffer, blshift(in_byte, in_bits))
        in_bits = in_bits + 8
    end
    code[1] = band(in_buffer, (blshift(1, bits) - 1))
    in_buffer = brshift(in_buffer ,brshift(in_buffer, bits))
    in_bits = in_bits - bits
    return in_bits > 0
end

--- Usage : decode [-3|4|5] [-a|u|l] [infile] [outfile]
--- @param args string[]
local function main(args)
    local code = {}
    local state = {}

    g72x.init_state(state)
    local out_coding = AUDIO_ENCODING_ULAW
    local out_size = 1 -- sizeof(char)
    local dec_routine = g721.decoder
    local dec_bits = 4

    -- Process encoding argument, if any
    local max_i
    for i=1, #args do
        max_i = i
        local arg = args[i]
        if arg:sub(1,1) ~= '-' then break end
        local option = arg:sub(2,2)
        
        if option == '3' then
            dec_routine = g723_24.decoder
            dec_bits = 3
        elseif option == '4' then
            dec_routine = g721.decoder
            dec_bits = 4
        elseif option == '5' then
            dec_routine = g723_40.decoder
            dec_bits = 5
        elseif option == 'u' then
            out_coding = AUDIO_ENCODING_ULAW
            out_size = 1
        elseif option == 'a' then
            out_coding = AUDIO_ENCODING_ALAW
            out_size = 1
        elseif option == 'l' then
            out_coding = AUDIO_ENCODING_LINEAR
            out_size = 2
        else
            print_usage()
        end
    end

    if #args - max_i < 2 then
        print_usage()
    end

    local in_file, out_file = assert(io.open(args[#args - 1], "rb")), assert(io.open(args[#args], "wb"))

    -- Read and unpack input codes and process them
    while unpack_input(code, dec_bits, in_file) do
        local sample = dec_routine(code, out_coding, state)
        if out_size == 2 then
            out_file:write(pack("h", sample))
        else
            out_file:write(pack("b", sample))
        end
    end

    in_file:close()
    out_file:close()
end

main({...})