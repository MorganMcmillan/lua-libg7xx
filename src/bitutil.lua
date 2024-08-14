local bit_blshift, bit_brshift = bit.blshift, bit.brshift

---@param n integer
---@param bits integer
---@return integer y
local function blshift(n, bits)
  return n > 0 and bit_blshift(n, bits) or -bit_blshift(-n, bits)
end

---@param n integer
---@param bits integer
---@return integer y
local function brshift(n, bits)
  return n > 0 and bit_brshift(n, bits) or -bit_brshift(-n, bits)
end

return {
  blshift = blshift,
  brshift = brshift
}