local bitutil = require("bitutil")
local band, brshift, blshift = bit.band, bitutil.brshift, bitutil.blshift

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

 * g723_40.c
 *
 * Description:
 *
 * g723_40_encoder(), g723_40_decoder()
 *
 * These routines comprise an implementation of the CCITT G.723 40Kbps
 * ADPCM coding algorithm.  Essentially, this implementation is identical to
 * the bit level description except for a few deviations which
 * take advantage of workstation attributes, such as hardware 2's
 * complement arithmetic.
 *
 * The deviation from the bit level specification (lookup tables),
 * preserves the bit level performance specifications.
 *
 * As outlined in the G.723 Recommendation, the algorithm is broken
 * down into modules.  Each section of code below is preceded by
 * the name of the module which it is implementing.
 *
]]
local g723_40 = {}

local g72x = require("g72x.lua")
local predictor_zero, predictor_pole, step_size, reconstruct, quantize, update, tandem_adjust_alaw, tandem_adjust_ulaw  = g72x.predictor_zero, g72x.predictor_pole, g72x.step_size, g72x.reconstruct, g72x.quantize, g72x.update, g72x.tandem_adjust_alaw, g72x.tandem_adjust_ulaw

local g711 = require("g711")
local alaw2linear, ulaw2linear = g711.alaw2linear, g711.ulaw2linear

local AUDIO_ENCODING_ULAW = 1
local AUDIO_ENCODING_ALAW = 2
local AUDIO_ENCODING_LINEAR = 3

-- Maps G.723_40 code word to ructeconstructed scale factor normalized log
-- magnitude values.
local _dqlntab = { -2048, -66, 28, 104, 169, 224, 274, 318, 358, 395, 429, 459, 488, 514, 539, 566, 566, 539, 514, 488, 459, 429, 395, 358, 318, 274, 224, 169, 104, 28, -66, -2048 }

-- Maps G.723_40 code word to log of scale factor multiplier.
local _witab = { 448, 448, 768, 1248, 1280, 1312, 1856, 3200, 4512, 5728, 7008, 8960, 11456, 14080, 16928, 22272,22272, 16928, 14080, 11456, 8960, 7008, 5728, 4512,3200, 1856, 1312, 1280, 1248, 768, 448, 448 }

-- Maps G.723_40 code words to a set of values whose long and short
-- term averages are computed and then compared to give an indication
-- how stationary (steady state) the signal is.
local _fitab = { 0, 0, 0, 0, 0, 0x200, 0x200, 0x200, 0x200, 0x200, 0x400, 0x600, 0x800, 0xA00, 0xC00, 0xC00, 0xC00, 0xC00, 0xA00, 0x800, 0x600, 0x400, 0x200, 0x200, 0x200, 0x200, 0x200, 0, 0, 0, 0, 0 }

local qtab_723_40 = { -122, -16, 68, 139, 198, 250, 298, 339, 378, 413, 445, 475, 502, 528, 553 }

--- g723_40_encoder()
---
--- Encodes a 16-bit linear PCM, A-law or u-law input sample and retuens
--- the resulting 5-bit CCITT G.723 40Kbps code.
--- Returns -1 if the input coding value is invalid.
---@param sl integer input sample
---@param in_coding 1|2|3 the encoding to use
---@param state g72x_state
function g723_40.encoder(sl, in_coding, state)
    local sei, sezi, se, sez -- ACCUM 
    local d                  -- SUBTA 
    local y                  -- MIX 
    local sr                 -- ADDB 
    local dqsez              -- ADDC 
    local dq, i

    -- linearize input sample to 14-bit PCM 
    if in_coding == AUDIO_ENCODING_ALAW then
        sl = brshift(alaw2linear(sl), 2)
    elseif in_coding == AUDIO_ENCODING_ULAW then
        sl = brshift(ulaw2linear(sl), 2)
    elseif in_coding == AUDIO_ENCODING_LINEAR then
        sl = brshift(sl, 2) -- sl of 14-bit dynamic range 
    else
        return -1
    end

    sezi = predictor_zero(state)
    sez = brshift(sezi, 1)
    sei = sezi + predictor_pole(state)
    se = brshift(sei, 1) -- se = estimated signal 

    d = sl - se -- d = estimation difference 

    -- quantize prediction difference 
    y = step_size(state)            -- adaptive quantizer step size 
    i = quantize(d, y, qtab_723_40, 15) -- i = ADPCM code 

    dq = reconstruct(band(i, 0x10), _dqlntab[i], y) -- quantized diff 

    sr = dq < 0 and se - band(dq, 0x7FFF) or se + dq -- reconstructed signal 

    dqsez = sr + sez - se -- dqsez = pole prediction diff. 

    update(5, y, _witab[i], _fitab[i], dq, sr, dqsez, state)

    return i
end

--- g723_40_decoder()
---
--- Decodes a 5-bit CCITT G.723 40Kbps code and returns
--- the resulting 16-bit linear PCM, A-law or u-law sample value.
--- -1 is returned if the output coding is unknown.
---@param i integer the encoded sample
---@param out_coding 1|2|3 the encoding to use
---@param state g72x_state
---@return integer the decoded sample
function g723_40.decoder(i, out_coding, state)
    local sezi, sei, sez, se -- ACCUM 
    local y                  -- MIX 
    local sr                 -- ADDB 
    local dq
    local dqsez

    i = band(i, 0x1f) -- mask to get proper bits 
    sezi = predictor_zero(state)
    sez = brshift(sezi, 1)
    sei = sezi + predictor_pole(state)
    se = brshift(sei, 1) -- se = estimated signal 

    y = step_size(state)                   -- adaptive quantizer step size 
    dq = reconstruct(band(i, 0x10), _dqlntab[i], y) -- estimation diff. 

    sr = dq < 0 and (se - band(dq, 0x7FFF)) or (se + dq) -- reconst. signal 

    dqsez = sr - se + sez -- pole prediction diff. 

    update(5, y, _witab[i], _fitab[i], dq, sr, dqsez, state)

    if out_coding == AUDIO_ENCODING_ALAW then
        return (tandem_adjust_alaw(sr, se, y, i, 0x10, qtab_723_40))
    elseif out_coding == AUDIO_ENCODING_ULAW then
        return (tandem_adjust_ulaw(sr, se, y, i, 0x10, qtab_723_40))
    elseif out_coding == AUDIO_ENCODING_LINEAR then
        return blshift(sr, 2) -- sr was of 14-bit dynamic range 
    else
        return -1
    end
end

return g723_40