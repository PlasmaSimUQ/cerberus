function binary_search(list, value)
  local low = 1
  local high = #list

  local mid_old = 0
  local mid = 1

  while (low < high) and (mid_old ~= mid) do
    mid_old = mid
    mid = math.floor((low + high) / 2)
    if (list[mid] <= value) and (list[mid + 1] >= value) then
      return mid
    elseif list[mid] > value then
      high = mid
    elseif list[mid] < value then
      low = mid
    end
  end
  return false
end

function linear_interp(x, xlist, ylist)
  local idx = binary_search(xlist, x)
  if idx then
    local y0 = ylist[idx]
    local y1 = ylist[idx + 1]
    local x0 = xlist[idx]
    local x1 = xlist[idx + 1]
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)
  else
    error(string.format('data out of bounds x = %g not in [%g, %g]', x, xlist[1], xlist[#xlist]))
  end
end
