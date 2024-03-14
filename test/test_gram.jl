using Test

using CompScienceMeshes
using BEAST
using StaticArrays
using LinearAlgebra

ω = 1.0

l1 = meshsegment(1.0,1/2)
l2 = meshsegment(1.0,1/4)

idcs = BEAST.interior_and_junction_vertices(l1, boundary(l1))
@test length(idcs) == 3

# lag1 = lagrangec0d1(l1, boundary(l1))
lag1 = lagrangec0d1(l1, skeleton(l1,0))
lag2 = lagrangecxd0(l2)

# @show numfunctions(lag1)
@test numfunctions(lag1) == 3
@test numfunctions(lag2) == 4

id = Identity()
G = assemble(id, lag1, lag2)
H = assemble(id, lag2, lag1)

Gt = [
    3 1 0 0
    1 3 3 1
    0 0 1 3] // 16

@test norm(G-H') ≈ 0.0  atol = sqrt(eps()) #G == H'
@test norm(G-Gt) ≈ 0.0  atol = sqrt(eps())

# Test whether localoperator and localoperator 2 give the same result
# for test and basis functions defined on the same mesh.
id = Identity()
nc = NCross()

#fn = Pkg.dir("BEAST","test","sphere2.in")
fn = joinpath(dirname(@__FILE__),"assets","sphere316.in")
m = readmesh(fn)
rt = raviartthomas(m)

G1 = zeros(ComplexF64, numfunctions(rt), numfunctions(rt))
BEAST.assemble_local_matched!(id, rt, rt, (v,m,n)->(G1[m,n]+=v))

G2 = zeros(ComplexF64, numfunctions(rt), numfunctions(rt))
BEAST.assemble_local_mixed!(id, rt, rt, (v,m,n)->(G2[m,n]+=v))

# Test wether the gram matrix assembled works for
# different but overlapping meshes
m1 = meshrectangle(1.0,1.0,1.0)
m2 = translate(m1, point(0.5,0.5,0.0))
Id = Identity()
x1 = lagrangecxd0(m1)
x2 = lagrangecxd0(m2)
G = assemble(Id,x1,x2)
@test sum(G) ≈ 1/4
