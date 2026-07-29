abstract type ImpedanceOperator3D{T,K} <: IntegralOperator end

scalartype(op::ImpedanceOperator3D{T,K}) where {T, K <: Val{0}} = T
scalartype(op::ImpedanceOperator3D{T,K}) where {T, K} = promote_type(T, K)

gamma(op::ImpedanceOperator3D{T,Val{0}}) where {T} = zero(T)
gamma(op::ImpedanceOperator3D{T,K}) where {T, K} = op.gamma


struct IBCOperator3D{T,U,V} <: ImpedanceOperator3D{T,U}
  gamma::T
  α::U
  β::U
  z::V
end

struct PCOperator3D{T,U,V,W} <: ImpedanceOperator3D{T,U}
  gamma::T
  α::U
  β::U
  z_alpha::V
  z_beta::W
end

struct ImpedancePCOperator3D{T,U,V} <: ImpedanceOperator3D{T,U}
  gamma::T
  α::U
  β::U
  z::V
end

defaultquadstrat(op::ImpedanceOperator3D, tfs::RTRefSpace, bfs::RTRefSpace) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)


"""
    ImpedanceOperator3D <: LocalOperator

The impedance operator to model spatially varying impedance.
"""
struct LIBCOperator3D{T} <: LocalOperator
  z::T
end

kernelvals(biop::LIBCOperator3D, x) = nothing
function integrand(op::LIBCOperator3D, kernel, x, g, f)
  x = cartesian(x)
  dot(f[1], g[1])*op.z(x)
end
scalartype(op::LIBCOperator3D) = ComplexF64

################################################################################
#
#  Kernel definitions
#
################################################################################


function (igd::Integrand{<:IBCOperator3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green
    βG = β * green
    z = igd.operator.z(cartesian(x))

    _integrands(f,g) do fi,gj
        z * αG * dot(fi.value, gj.value) + z * βG * dot(fi.divergence, gj.divergence)
    end
end

function (igd::Integrand{<:PCOperator3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green
    βG = β * green
    z_alpha = igd.operator.z_alpha(cartesian(x))
    z_beta = igd.operator.z_beta(cartesian(x))

    _integrands(f,g) do fi,gj
        z_alpha * αG * dot(fi.value, gj.value) + z_beta * βG * dot(fi.divergence, gj.divergence)
    end
end

function (igd::Integrand{<:ImpedancePCOperator3D})(x,y,f,g)
    α = igd.operator.α
    β = igd.operator.β
    γ = igd.operator.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(i4pi*iR)

    αG = α * green
    βG = β * green
    z = igd.operator.z(cartesian(x))

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + z * βG * dot(fi.divergence, gj.divergence)
    end
end

mutable struct GaussianInc{T,P,U}
  direction::P
  polarisation::P
  gamma::T
  amplitude::T
  origin::U
  variance::U
end


scalartype(x::GaussianInc{T,P,U}) where {T,P,U} = promote_type(T, eltype(P)) 


function (e::GaussianInc)(x)
  γ = e.gamma
  d = e.direction
  u = e.polarisation
  a = e.amplitude
  o = e.origin
  s = e.variance
  #@show x
  gfactor = 0.5*sum(i -> ((x[i]-o[i])/s[i])^2, eachindex(o))
  #@show gfactor
  a * exp(-γ * dot(d, x)) * u * gfactor
end

function curl(field::PlaneWaveMW)
  γ = field.gamma
  d = field.direction
  u = field.polarisation
  a = field.amplitude
  o = field.origin
  s = field.variance
  v = d × u
  b = -a * γ
  GaussianInc(d, v, γ, b, o, s)
end

*(a::Number, e::GaussianInc) = GaussianInc(e.direction, e.polarisation, e.gamma, a*e.amplitude, e.origin, e.variance)