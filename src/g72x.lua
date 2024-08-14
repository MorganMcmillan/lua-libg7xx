local abs ,band, bxor = math.abs ,bit.band, bit.bxor

local bitutil = require("bitutil")
local blshift, brshift = bitutil.blshift, bitutil.brshift

--[[
 * This source code is a product of Sun Microsystems, Inc. and is provided
 * for unrestricted use.  Users may copy or modify this source code without
 * charge.
 *
 * SUN SOURCE CODE IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING
 * THE WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
 *
 * Sun source code is provided with no support and without any obligation on
 * the part of Sun Microsystems, Inc. to assist in its use, correction,
 * modification or enhancement.
 *
 * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
 * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY THIS SOFTWARE
 * OR ANY PART THEREOF.
 *
 * In no event will Sun Microsystems, Inc. be liable for any lost revenue
 * or profits or other special, indirect and consequential damages, even if
 * Sun has been advised of the possibility of such damages.
 *
 * Sun Microsystems, Inc.
 * 2550 Garcia Avenue
 * Mountain View, California  94043
 */
#include <stdlib.h>
/*
 * g72x.c
 *
 * Common routines for G.721 and G.723 conversions.
]]
local g72x = {}

local g711 = require("g711")
local linear2alaw, alaw2linear, linear2ulaw, ulaw2linear = g711.linear2alaw, g711.alaw2linear, g711.linear2ulaw, g711.ulaw2linear

local power2 = {1, 2, 4, 8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000}

--- quantizes the input val against the table of size short integers.
--- It returns i if table[i - 1] <= val < table[i].
---
--- Using linear search for simple coding.
---@param val integer
---@param table integer[]
---@param size integer the size of the table
---@return integer index the quantization index, or size + 1 if not found
local function quan(val, table, size)
    for i = 1, size do
        if val < table[i] then return i end
    end
    return size + 1
end

--- returns the integer product of the 14-bit integer "an" and
--- "floating point" representation (4-bit exponent, 6-bit mantessa) "srn".
---@param an integer
---@param srn integer
---@return integer
local function fmult(an, srn)
    local anmag = an > 0 and an or band(-an, 0x1FFF)
    local anexp = quan(anmag, power2, 15) - 6
    local anmant = (anmag == 0) and 32 or (anexp >= 0 and brshift(anmag, anexp) or blshift(anmag, -anexp))
    local wanexp = anexp + band(brshift(srn, 6), 0xF) - 13

    local wanmant = brshift(anmant * band(srn, 077) + 0x30, 4)
    local retval = (wanexp >= 0) and band(blshift(wanmant, wanexp), 0x7FFF) or brshift(wanmant, -wanexp)

    return bxor(an, srn) < 0 and -retval or retval
end

--- The following is the definition of the state structure
--- used by the G.721/G.723 encoder and decoder to preserve their internal
--- state between successive calls.  The meanings of the majority
--- of the state structure fields are explained in detail in the
--- CCITT Recommendation G.721.  The field names are essentially identical
--- to variable names in the bit level description of the coding algorithm
--- included in this Recommendation.
---@class g72x_state
---@field yl integer Locked or steady state step size multiplier.
---@field yu integer Unlocked or non-steady state step size multiplier.
---@field dms integer Short term energy estimate.
---@field dml integer Long term energy estimate.
---@field ap integer Linear weighting coefficient of 'yl' and 'yu'.
---@field a integer[] Coefficients of pole portion of prediction filter. Is 2 elements long
---@field b integer[] Coefficients of zero portion of prediction filter. Is 6 elements long
---@field pk integer[] Signs of previous two samples of a partially reconstructed signal. Is 2 elements long
---@field dq integer[] Previous 6 samples of the quantized difference signal represented in an internal floating point format. Is 6 elements long
---@field sr integer[] Previous 2 samples of the quantized difference signal represented in an internal floating point format. Is 2 elements long
---@field td integer delayed tone detect, new in 1988 version

--- This routine initializes and/or resets the g72x_state structure
--- pointed to by 'state'.
--- All the initial state values are specified in the CCITT G.721 document.
---@param state g72x_state
function g72x.init_state(state)
    state.yl = 34816
    state.yu = 544
    state.dms = 0
    state.dml = 0
    state.ap = 0
    state.a = {0, 0}
    state.pk = {0, 0}
    state.sr = {0, 0}
    state.b = {0, 0, 0, 0, 0, 0}
    state.dq = {0, 0, 0, 0, 0, 0}
    state.td = 0
end

--- computes the estimated signal from 6-zero predictor.
---@param state g72x_state
---@return integer sezi the estimated signal
function g72x.predictor_zero(state)
    local sezi = fmult(brshift(state.b[1], 2), state.dq[1] )
    for i=2,6 do -- ACCUM
        sezi = sezi + fmult(brshift(state.b[i], 2), state.dq[i])
    end
    return sezi
end

--- computes the estimated signal from 2-pole predictor.
---@param state g72x_state
---@return integer sezi the estimated signal
function g72x.predictor_pole(state)
    return fmult(brshift(state.a[2], 2), state.sr[2] ) +
            fmult(brshift(state.a[1], 2), state.sr[1] )
end

--- computes the quantization step size of the adaptive quantizer.
---@param state g72x_state
---@return integer step_size the quantization step size
function g72x.step_size(state)
    if state.ap >= 256 then
        return state.yu
    else
        local y = brshift(state.yl, 6)
        local dif = state.yu - y
        local al = brshift(state.ap, 2)
        if dif > 0 then
            y = y + brshift(dif * al, 6)
        elseif dif < 0 then
            y = y + brshift(dif * al + 0x3F, 6)
        end
        return y
    end
end

--- Given a raw sample, 'd', of the difference signal and a
--- quantization step size scale factor, 'y', this routine returns the
--- ADPCM codeword to which that sample gets quantized.  The step
--- size scale factor division operation is done in the log base 2 domain
--- as a subtraction.
---@param d integer Raw difference signal sample
---@param y integer Step size multiplier
---@param table integer[] quantization table
---@param size integer DEPRECATED: table size of short integers
---@return integer codeword the ADPCM codeword
local function quantize(d, y, table, size)
    local dqm  -- Magnitude of 'd'
    local exp  -- Integer part of base 2 log of 'd'
    local mant -- Fractional part of base 2 log
    local dl   -- Log of magnitude of 'd'
    local dln  -- Step size scale factor normalized log
    local i

    -- LOG
    --
    -- Compute base 2 log of 'd', and store in 'dl'.

    dqm = abs(d)
    exp = quan(brshift(dqm, 1), power2, 15)
    mant = band(brshift(brshift(dqm, 7), exp), 0x7F) -- Fractional portion.
    dl = brshift(exp, 7) + mant

    -- SUBTB
    --
    -- "Divide" by step size multiplier.
    dln = dl - brshift(y, 2)

    -- QUAN
    --
    -- Obtain codword i for 'd'.

    i = quan(dln, table, size)
    if d < 0 then -- take 1's complement of i
        return brshift(size, 1) + 1 - i
    elseif i == 0 then              -- take 1's complement of 0
        return brshift(size, 1) + 1 -- new in 1988
    else
        return i
    end
end

g72x.quantize = quantize

--- Returns reconstructed difference signal 'dq' obtained from
--- codeword 'i' and quantization step size scale factor 'y'.
--- Multiplication is performed in log base 2 domain as addition.
---@param sign integer 0 for non-negative value
---@param dqln integer G.72x codeword
---@param y integer Step size multiplier
---@return integer dq the reconstructed difference signal
function g72x.reconstruct(sign, dqln, y)
    local dql -- Log of 'dq' magnitude
    local dex -- Integer part of log
    local dqt
    local dq -- Reconstructed difference signal sample

    dql = dqln + brshift(y, 2) -- ADDA

    if dql < 0 then
        return sign and -0x8000 or 0
    else -- ANTILOG
        dex = band(brshift(dql, 7), 15)
        dqt = 128 + band(dql, 127)
        dq = brshift(brshift(dqt, 7), 14 - dex)
        return sign and (dq - 0x8000) or dq
    end
end

--- updates the state variables for each output code
---@param code_size integer distinguish 723_40 with others
---@param y integer quantizer step size
---@param wi integer scale factor multiplier
---@param fi integer for long/short term energies
---@param dq integer quantized prediction difference
---@param sr integer reconstructed signal
---@param dqsez integer difference from 2-pole predictor
---@param state g72x_state coder state pointer
function g72x.update(code_size, y, wi, fi, dq, sr, dqsez, state)
    local mag, exp -- Adaptive predictor, FLOAT A
    local a2p = 0  -- LIMC
    local a1ul     -- UPA1
    local pks1     -- UPA2
    local fa1
    local tr -- tone/transition detector
    local ylint, thr2, dqthr
    local ylfrac, thr1
    local pk0

    pk0 = (dqsez < 0) and 1 or 0 -- needed in updating predictor poles

    mag = band(dq, 0x7FFF) -- prediction difference magnitude
    -- TRANS
    ylint = brshift(state.yl, 15)           -- exponent part of yl
    ylfrac = band(brshift(state.yl, 10), 0x1F) -- fractional part of yl
    thr1 = blshift(32 + ylfrac, ylint)         -- threshold
    thr2 = (ylint > 9) and brshift(31, 10) or thr1  -- limit thr2 to 31 << 10
    dqthr = brshift(thr2 + brshift(thr2, 1), 1)     -- dqthr = 0.75 * thr2
    if state.td == 0 then                -- signal supposed voice
        tr = 0
    elseif mag <= dqthr then -- supposed data, but small mag
        tr = 0            -- treated as voice
    else                   -- signal is data (modem)
        tr = 1
    end

    -- Quantizer scale factor adaptation.

    -- FUNCTW & FILTD & DELAY
    -- update non-steady state step size multiplier
    state.yu = y + brshift(wi - y, 5)

    -- LIMB
    if state.yu < 544 then -- 544 <= yu <= 5120
        state.yu = 544
    elseif state.yu > 5120 then
        state.yu = 5120
    end
    -- FILTE & DELAY
    -- update steady state step size multiplier
    state.yl = state.yl + state.yu + brshift(-state.yl, 6)

    -- Adaptive predictor coefficients.

    if tr == 1 then -- reset a's and b's for modem signal
        for i=1,2 do state.a[i] = 0 end
        for i=1,6 do state.b[i] = 0 end
    else                             -- update a's and b's
        pks1 = bxor(pk0, state.pk[1])  -- UPA2

        -- update predictor pole a[2]
        a2p = state.a[2]  - brshift(state.a[2], 7)
        if dqsez ~= 0 then
            fa1 = (pks1) and state.a[1]  or -state.a[1]
            if fa1 < -8191 then -- a2p = function of fa1
                a2p = a2p - 0x100
            elseif fa1 > 8191 then
                a2p = a2p + 0xFF
            else
                a2p = a2p + brshift(fa1, 5)
            end

            if bxor(pk0, state.pk[2] ) ~= 0 then
                -- LIMC
                if a2p <= -12160 then
                    a2p = -12288
                elseif a2p >= 12416 then
                    a2p = 12288
                else
                    a2p = a2p - 0x80
                end
            elseif a2p <= -12416 then
                a2p = -12288
            elseif a2p >= 12160 then
                a2p = 12288
            else
                a2p = a2p + 0x80
            end
        end

        -- TRIGB & DELAY
        state.a[2]  = a2p

        -- UPA1
        -- update predictor pole a[1]
        state.a[1] = state.a[1] - brshift(state.a[1], 8)
        if dqsez ~= 0 then
            if pks1 == 0 then
                state.a[1] = state.a[1] + 192
            else
                state.a[1] = state.a[1] - 192
            end
        end

        -- LIMD
        a1ul = 15360 - a2p
        if state.a[1]  < -a1ul then
            state.a[1]  = -a1ul
        elseif state.a[1]  > a1ul then
            state.a[1]  = a1ul
        end

        -- UPB : update predictor zeros b[6]
        for cnt=1,6 do
            if code_size == 5 then -- for 40Kbps G.723
                state.b[cnt] = state.b[cnt] - brshift(state.b[cnt], 9)
            else -- for G.721 and 24Kbps G.723
                state.b[cnt] = state.b[cnt] - brshift(state.b[cnt], 8)
            end
            if band(dq, 0x7FFF) ~= 0 then -- XOR
                if bxor(dq, state.dq[cnt]) >= 0 then
                    state.b[cnt] = state.b[cnt] + 128
                else
                    state.b[cnt] = state.b[cnt] - 128
                end
            end
        end
    end

    for cnt=6,2,-1 do
        state.dq[cnt] = state.dq[cnt - 1]
    end
    -- FLOAT A : convert dq[1]  to 4-bit exp, 6-bit mantissa f.p.
    if mag == 0 then
        state.dq[1]  = (dq >= 0) and 0x20 or 0xFC20
    else
        exp = quan(mag, power2, 15)
        state.dq[1]  = (dq >= 0) and brshift(exp, 6) + brshift(brshift(mag, 6), exp) or brshift(exp, 6) + brshift(brshift(mag, 6), exp) - 0x400
    end

    state.sr[2]  = state.sr[1]
    -- FLOAT B : convert sr to 4-bit exp., 6-bit mantissa f.p.
    if sr == 0 then
        state.sr[1]  = 0x20
    elseif sr > 0 then
        exp = quan(sr, power2, 15)
        state.sr[1]  = brshift(exp, 6) + brshift(brshift(sr, 6), exp)
    elseif sr > -32768 then
        mag = -sr
        exp = quan(mag, power2, 15)
        state.sr[1]  = brshift(exp, 6) + brshift(brshift(mag, 6), exp) - 0x400
    else
        state.sr[1]  = 0xFC20
    end

    -- DELAY A
    state.pk[2]  = state.pk[1]
    state.pk[1]  = pk0

    -- TONE
    if tr == 1 then           -- this sample has been treated as data
        state.td = 0 -- next one will be treated as voice
    elseif a2p < -11776 then -- small sample-to-sample correlation
        state.td = 1 -- signal may be data
    else                   -- signal is voice
        state.td = 0
    end

    -- Adaptation speed control.
    state.dms = state.dms + brshift(fi - state.dms, 5)          -- FILTA
    state.dml = state.dml + brshift(brshift(fi, 2) - state.dml, 7) -- FILTB

    if tr == 1 then
        state.ap = 256
    elseif y < 1536 then -- SUBTC
        state.ap = state.ap + brshift(0x200 - state.ap, 4)
    elseif state.td == 1 then
        state.ap = state.ap + brshift(0x200 - state.ap, 4)
    elseif abs(brshift(state.dms, 2) - state.dml) >= brshift(state.dml, 3) then
        state.ap = state.ap + brshift(0x200 - state.ap, 4)
    else
        state.ap = state.ap + brshift(-state.ap, 4)
    end
end



--- tandem_adjust(sr, se, y, i, sign)
---
--- At the end of ADPCM decoding, it simulates an encoder which may be receiving
--- the output of this decoder as a tandem process. If the output of the
--- simulated encoder differs from the input to this decoder, the decoder output
--- is adjusted by one level of A-law or u-law codes.
---@param sr integer decoder output linear PCM sample
---@param se integer predictor estimate sample
---@param y integer quantizer step size
---@param i integer decoder input code
---@param sign integer sign bit of code i
---@param qtab integer[]
---@return integer sp adjusted A-law or u-law compressed sample
function g72x.tandem_adjust_alaw(sr, se, y, i, sign, qtab)
    local sp -- A-law compressed 8-bit code
    local dx -- prediction error
    local id -- quantized prediction error
    local sd -- adjusted A-law decoded sample value
    local im -- biased magnitude of i
    local imx -- biased magnitude of id

    if sr <= -32768 then
        sr = -1
    end
    sp = linear2alaw(blshift(brshift(sr, 1), 3)) -- short to A-law compression
    dx = brshift(alaw2linear(sp), 2) - se -- 16-bit prediction error
    id = quantize(dx, y, qtab, sign - 1)

    if id == i then -- no adjustment on sp
        return sp
    else -- sp adjustment needed
        -- ADPCM codes : 8, 9, ... F, 0, 1, ... , 6, 7
        im = bxor(i, sign) -- 2's complement to biased unsigned
        imx = bxor(id, sign)

        if imx > im then -- sp adjusted to next lower value
            if band(sp, 0x80) ~= 0 then
                sd = (sp == 0xD5) and 0x55 or bxor(bxor(sp, 0x55) - 1, 0x55)
            else
                sd = (sp == 0x2A) and 0x2A or bxor(bxor(sp, 0x55) + 1, 0x55)
            end
        else -- sp adjusted to next higher value
            if band(sp, 0x80) ~= 0 then
                sd = (sp == 0xAA) and 0xAA or bxor(bxor(sp, 0x55) + 1, 0x55)
            else
                sd = (sp == 0x55) and 0xD5 or bxor(bxor(sp, 0x55) - 1, 0x55)
            end
        end
        return sd
    end
end

---@param sr integer decoder output linear PCM sample
---@param se integer predictor estimate sample
---@param y integer quantizer step size
---@param i integer decoder input code
---@param sign integer
---@param qtab integer[]
---@return integer sp adjusted u-law compressed sample
function g72x.tandem_adjust_ulaw(sr, se, y, i, sign, qtab)
    local sp -- u-law compressed 8-bit code
    local dx -- prediction error
    local id -- quantized prediction error
    local sd -- adjusted u-law decoded sample value
    local im -- biased magnitude of i
    local imx -- biased magnitude of id

    if sr <= -32768 then
        sr = 0
    end
    sp = linear2ulaw(brshift(sr, 2))        -- short to u-law compression
    dx = brshift(ulaw2linear(sp), 2) - se -- 16-bit prediction error
    id = quantize(dx, y, qtab, sign - 1)
    if id == i then
        return sp
    else
        -- ADPCM codes : 8, 9, ... F, 0, 1, ... , 6, 7
        im = bxor(i, sign) -- 2's complement to biased unsigned
        imx = bxor(id, sign)
        if imx > im then -- sp adjusted to next lower value
            if band(sp, 0x80) ~= 0 then
                sd = (sp == 0xFF) and 0x7E or sp + 1
            else
                sd = (sp == 0) and 0 or sp - 1
            end

        else -- sp adjusted to next higher value
            if band(sp, 0x80) ~= 0 then
                sd = (sp == 0x80) and 0x80 or sp - 1
            else
                sd = (sp == 0x7F) and 0xFE or sp + 1
            end
        end
        return sd
    end
end

return g72x