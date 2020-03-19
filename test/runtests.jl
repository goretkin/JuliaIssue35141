using GridVAMP
using Test
const GGeo = GridVAMP.GridGeometry

import GeometryTypes
import GeometryTypes: Vec, Point
import LinearAlgebra: dot
# import OffsetArrays: OffsetArray

function project_linesegment(onto, point)
  (p2, p1) = onto
  d = (p2 - p1)
  vp1 = (d * dot(point - p1, d)) / (dot(d, d))
  vp1 + p1
end

"""from `onto` to `point`"""
function perp_vec(onto, point)
  from = project_linesegment(onto, point)
  return point - from
end

"""conservative rasterization
Pixels _is_ a little 1-by-1 squares, centered at interger coordinates.
"""
function BROKEN_rasterize_segment!(out::AbstractArray, (p1, p2))
  for ci in CartesianIndices(out)
    i = Point(Tuple(ci)...)
    v = perp_vec((p1, p2), i)
    out[ci] = all(-0.5 .< v .< 0.5)
  end
  return out
end

@testset "integer line" begin
    for i = 1:100
        r = -100:100
        p1 = (rand(r), rand(r))
        p2 = (rand(r), rand(r))
        int_line1 = GGeo.cells_overlapped_nd(p1, p2)
        int_line2 = GGeo.cells_overlapped_2d(p1..., p2...)
        int_line3 = collect(GGeo.LineSegmentGridIterator(GGeo.LineSegment(p1, p2)))

        #=
        corners = (first(int_line1), last(int_line1))
        ranges = map(x->UnitRange(minmax(x...)...), zip(corners...))
        bmp = OffsetArray(zeros(Bool, length.(ranges)...), ranges...)
        rasterize_segment!(bmp, (Point(p1...), Point(p2...)))
        int_line3 = Tuple.(findall(bmp))
        =#

        @test Set(int_line1) == Set(int_line2)
        @test Set(int_line1) == Set(int_line3)
    end
    p1 = (1, 1)
    p2 = (10, 4)

    ref = [(1, 1), (2, 1), (2, 1), (3, 2), (4, 2), (5, 2), (5, 2), (6, 3), (7, 3), (8, 3), (8, 3), (9, 4), (10, 4)]
    # test that (5, 2) -> (6, 3) transition happens
    @test ref == GridVAMP.GridGeometry.cells_overlapped_nd(p1, p2)
end

@testset "coo to im" begin
  b = rand(Bool, 10, 10)
  b2 = GridVAMP.GridGeometry.coo_to_bool_array(Tuple.(findall(b)))
  @test sort(findall(b)) == sort(findall(b2))
end

@testset "rotgrid" begin
  a = reshape(1:9, 3, 3)
  for n = 0:3, i=1:3, j=1:3
    # set 0,0 to each of the 9 positions
    b = GridVAMP.GridGeometry.OffsetArray(a, ((1:3) .- i, (1:3) .- j))
    br = GridVAMP.GridGeometry.rotgrid(b, n)
    @test br[0, 0] == b[0, 0]
  end
end

#=
p1 = (1, 1)
p2 = (10, 4)
int_line1 = GridVAMP.GridGeometry.cells_overlapped_nd(p1, p2)
corners = (first(int_line1), last(int_line1))
ranges = map(x->UnitRange(minmax(x...)...), zip(corners...))
bmp = OffsetArray(zeros(Bool, length.(ranges)...), ranges...)
rasterize_segment!(bmp, (Point(p1...), Point(p2...)))
int_line3 = Tuple.(findall(bmp))

bmp = bmp .* 2
bmp[CartesianIndex.(setdiff(Set(int_line1), Set(int_line3)))] .= 3


f = GRUtils.Figure()
GRUtils.aspectratio!(f, 1)
#GRUtils.grid!(f, true)
GRUtils.hold!(f, true)
GRUtils.heatmap!(f, (map(r->UnitRange(first(r) - 0.5, last(r) + 0.5), axes(bmp)))..., bmp'; grid=true)
GRUtils.plot!(f, [p1[1], p2[1]], [p1[2], p2[2]]; grid=true)
#GRUtils.GR.grid(1, 1, 0, 0, 0, 0)
display(f)
=#

include("collision_motion_planning.jl")
include("visibility_geometry.jl")
include("visibility.jl")
include("local_vamp.jl")
