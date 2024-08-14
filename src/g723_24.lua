local bitutil = require("bitutil")
local band, brshift = bit.band, bitutil.brshift

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

 * g723_24.c
 *
 * Description:
 *
 * g723_24_encoder(), g723_24_decoder()
 *
 * These routines comprise an implementation of the CCITT G.723 24 Kbps
 * ADPCM coding algorithm.  Essentially, this implementation is identical to
 * the bit level description except for a few deviations which take advantage
 * of workstation attributes, such as hardware 2's complement arithmetic.
 *
]]
local g723_24 = {}

local g72x = require("g72x.lua")
local predictor_zero, predictor_pole, step_size, reconstruct, quantize, update, tandem_adjust_alaw, tandem_adjust_ulaw  = g72x.predictor_zero, g72x.predictor_pole, g72x.step_size, g72x.reconstruct, g72x.quantize, g72x.update, g72x.tandem_adjust_alaw, g72x.tandem_adjust_ulaw

local g711 = require("g711")
local alaw2linear, ulaw2linear = g711.alaw2linear, g711.ulaw2linear

local AUDIO_ENCODING_ULAW = 1
local AUDIO_ENCODING_ALAW = 2
local AUDIO_ENCODING_LINEAR = 3

-- Maps G.723_24 code word to reconstructed scale factor normalized log
-- magnitude values.
local _dqlntab = { -2048, 135, 273, 373, 373, 273, 135, -2048 }

-- Maps G.723_24 code word to log of scale factor multiplier. 
local _witab = { -128, 960, 4384, 18624, 18624, 4384, 960, -128 }

-- Maps G.723_24 code words to a set of values whose long and short
-- term averages are computed and then compared to give an indication
-- how stationary (steady state) the signal is.
local _fitab = { 0, 0x200, 0x400, 0xE00, 0xE00, 0x400, 0x200, 0 }

local qtab_723_24 = { 8, 218, 331 }

--- g723_24_encoder()
---
--- Encodes a linear PCM, A-law or u-law input sample and returns its 3-bit code.
--- Returns -1 if invalid input coding value.
---@param sl integer input sample
---@param in_coding 1|2|3 the encoding to use
---@param state g72x_state
function g723_24.encoder(sl, in_coding, state)
    local sei, sezi, se, sez; -- ACCUM 
    local d;                  -- SUBTA 
    local y;                  -- MIX 
    local sr;                 -- ADDB 
    local dqsez;              -- ADDC 
    local dq, i;

    -- linearize input sample to 14-bit PCM 
    if in_coding == AUDIO_ENCODING_ALAW then
        sl = brshift(alaw2linear(sl), 2);
    elseif in_coding == AUDIO_ENCODING_ULAW then
        sl = brshift(ulaw2linear(sl), 2);
    elseif in_coding == AUDIO_ENCODING_LINEAR then
        sl = brshift(sl, 2); -- sl of 14-bit dynamic range 
    else
        return -1;
    end

    sezi = predictor_zero(state);
    sez = brshift(sezi, 1);
    sei = sezi + predictor_pole(state);
    se = brshift(sei, 1); -- se = estimated signal 

    d = sl - se; -- d = estimation diff. 

    -- quantize prediction difference d 
    y = step_size(state);                -- quantizer step size 
    i = quantize(d, y, qtab_723_24, 3);      -- i = ADPCM code 
    dq = reconstruct(i & 4, _dqlntab[i], y); -- quantized diff. 

    sr = dq < 0 and se - band(dq, 0x3FFF) or se + dq; -- reconstructed signal 

    dqsez = sr + sez - se; -- pole prediction diff. 

    update(3, y, _witab[i], _fitab[i], dq, sr, dqsez, state);

    return i;
end

--- g723_24_decoder()
---
--- Decodes a 3-bit CCITT G.723_24 ADPCM code and returns
--- the resulting 16-bit linear PCM, A-law or u-law sample value.
--- -1 is returned if the output coding is unknown.
---@param i integer
---@param out_coding 1|2|3 the encoding to use
---@param state g72x_state
---@return integer the decoded sample
function g723_24.decoder(i, out_coding, state)
    local sezi, sei, sez, se; -- ACCUM 
    local y;                  -- MIX 
    local sr;                 -- ADDB 
    local dq;
    local dqsez;

    i = band(i, 0x07); -- mask to get proper bits 
    sezi = predictor_zero(state);
    sez = brshift(sezi, 1);
    sei = sezi + predictor_pole(state);
    se = brshift(sei, 1); -- se = estimated signal 

    y = step_size(state);                   -- adaptive quantizer step size 
    dq = reconstruct(band(i, 0x04), _dqlntab[i], y); -- unquantize pred diff 

    sr = dq < 0 and se - band(dq, 0x3FFF) or se + dq; -- reconst. signal 

    dqsez = sr - se + sez; -- pole prediction diff. 

    update(3, y, _witab[i], _fitab[i], dq, sr, dqsez, state);

    if out_coding == AUDIO_ENCODING_ALAW then
        return (tandem_adjust_alaw(sr, se, y, i, 4, qtab_723_24));
    elseif out_coding == AUDIO_ENCODING_ULAW then
        return (tandem_adjust_ulaw(sr, se, y, i, 4, qtab_723_24));
    elseif out_coding == AUDIO_ENCODING_ULAW then
        return brshift(sr, 2); -- sr was of 14-bit dynamic range 
    else
        return -1;
    end
end

return g723_24