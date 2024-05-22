


"""
    StagedTimeStep{T,N,NN}

T: the value type of the basis function.
N: the number of stages.
NN: the number of stages squared
Each time step has intermediary stages given by the vertor c in a Butcher tableau (A,b,c)
"""
struct StagedTimeStep{T, N, NN, I}
	Δt :: T
	Nt :: I
	c  :: SVector{N,T}
	A  :: SArray{Tuple{N,N},T,2,NN}
	b  :: SVector{N,T}
	zTransformedTermCount :: I
	contourRadius         :: T
end

"""
	FiniteDiffTimeStep{T,N,NN}

To be used in single step finite difference method schemes such as 
the backward euler or backward difference formula 2

"""

struct FiniteDiffTimeStep{T, U, V}
	Δt::T
	Nt::U
	zTransformedTermCount ::U
	contourRadius::T
	method::V
end

scalartype(sts :: StagedTimeStep{T, N, NN, I}) where {T, N, NN, I} = T
temporalbasis(sts :: StagedTimeStep{T, N, NN, I}) where {T, N, NN, I} = timebasisdelta(sts.Δt, sts.Nt)

numfunctions(s::StagedTimeStep) = s.Nt
numfunctions(s::FiniteDiffTimeStep) = s.Nt

numstages(s) = 1
numstages(s::StagedTimeStep) = size(s.c,1)

scalartype(fdts :: FiniteDiffTimeStep{T,I}) where {T,I} = T
temporalbasis(fdts :: FiniteDiffTimeStep{T,I}) where {T,I} = timebasisdelta(fdts.Δt, fdts.Nt)