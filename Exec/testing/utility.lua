function make_circles(x, y, collection)
  local d, dd

  for i, v in ipairs(collection) do
    dd = v[1] ^ 2 - ((x - v[2]) ^ 2 + (y - v[3]) ^ 2)
    if i == 1 then
      d = dd
    end
    d = math.max(d, dd)
  end

  return -d
end

-- collection = {{{x_lo, x_hi},{y_lo,y_hi}}, ... }
function make_rectangles(x, y, collection)
  local d, d1, dx, dy, dx_, dy_

  local coords = { x, y }

  for i, v in ipairs(collection) do
    dx = math.max(coords[1] - v[1][2], v[1][1] - coords[1])
    dy = math.max(coords[2] - v[2][2], v[2][1] - coords[2])

    dx_ = math.max(dx, 0.0)
    dy_ = math.max(dy, 0.0)

    d1 = math.sqrt(dx_ * dx_ + dy_ * dy_) + math.min(0.0, math.max(dx, dy))

    if i == 1 then
      d = d1
    else
      d = math.min(d, d1)
    end
  end

  return d
end

-- vertices = {{x,y},{xy},...}
function polygon_distance(x, y, vertices)
  local dd = math.huge
  local wn = 0 -- winding number (inside/outside)

  local pax, pay, pbx, pby, ax, ay, bx, by, vx, vy, nsegs, t, dn, l

  -- iterate over the line segments
  nsegs = #vertices - 1

  pbx = x - vertices[1][1]
  pby = y - vertices[1][2]

  for i = 1, nsegs do
    ax = vertices[i][1]
    ay = vertices[i][2]

    bx = vertices[i + 1][1]
    by = vertices[i + 1][2]

    -- vector a -> b
    vx = bx - ax
    vy = by - ay

    -- length of segment
    l = math.sqrt(vx * vx + vy * vy)

    -- normalise vector
    vx = vx / l
    vy = vy / l

    pax = pbx
    pay = pby

    pbx = x - bx
    pby = y - by

    t = pax * vx + pay * vy -- t-parameter of projection onto line
    dn = pax * vy - pay * vx -- normal distance from p to line

    -- Distance to line segment
    if t < 0 then
      dd = math.min(dd, pax * pax + pay * pay) -- distance to vertex[0] of line
    elseif t > l then
      dd = math.min(dd, pbx * pbx + pby * pby) -- distance to vertex[1] of line
    else
      dd = math.min(dd, dn * dn) -- normal distance to line
    end

    -- Is the point in the polygon?
    -- See: http://geomalgorithms.com/a03-_inclusion.html
    if ay <= y then
      if (by > y) and (dn < 0) then
        wn = wn + 1
      end
    else
      if (by <= y) and (dn > 0) then
        wn = wn - 1
      end
    end
  end

  -- normalise d*d to d
  d = math.sqrt(dd)
  if wn ~= 0.0 then
    -- p is inside the polygon
    return -d
  end

  return d
end

-- square = {{-1,-1}, {1,-1}, {1,1}, {-1,1}, {-1,-1}}
-- print(poly_distance(0, 2, square))
