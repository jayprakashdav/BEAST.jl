struct FiniteDiffConvolutionQuadrature{}
	timedomainkernel :: AbstractSpaceTimeOperator
end

abstract  type FiniteDiffMethod end

struct BDF2{T} <:FiniteDiffMethod where T
	p::BEAST.Polynomial
	dt::T
end

struct BE{T} <: FiniteDiffMethod where T
	p::BEAST.Polynomial
	dt::T
end

BDF2(dt) = BDF2(BEAST.Polynomial(1.5/dt, -2.0/dt, 0.5/dt), dt) #2-step backward difference formula

BE(dt) = BE(BEAST.Polynomial(1.0/dt,-1.0/dt), dt) #Backward Euler or 1-step backward difference formula

BDF1(dt) = BE(dt) #1-step backward difference formula


scalartype(etcq::FiniteDiffConvolutionQuadrature) = ComplexF64;

function laplacekernel(op::FiniteDiffConvolutionQuadrature, s, T)
	ws_weight = op.timedomainkernel.ws_weight
	hs_weight = op.timedomainkernel.hs_weight
	c = op.timedomainkernel.speed_of_light
	ws_diffs = op.timedomainkernel.ws_diffs
	hs_diffs = op.timedomainkernel.hs_diffs
	if isa(op.timedomainkernel, MWSingleLayerTDIO)
		@info "converting time domain single layer operator to laplace domain"
		return MWSingleLayer3D(s/c, ws_weight*s^ws_diffs/c, hs_weight*s^hs_diffs*T(c))
	elseif isa(op.timedomainkernel, AMWSingleLayerTDIO)
		@info "converting scalar potential part of the time domain single layer to laplace domain"
		return AugmentedMaxwellOperator3D(s/c, ws_weight*s^ws_diffs*T(1.0), hs_weight*s^hs_diffs*T(c))
	end
end

function assemble(cqop :: FiniteDiffConvolutionQuadrature,
                  testfns :: SpaceTimeBasis,
                  trialfns :: SpaceTimeBasis, threading::Type{Threading{:multi}}=Threading{:multi})

	LaplaceOp(s::T) where {T} = laplacekernel(cqop, s, T)
	method = testfns.time.method
	Δt = testfns.time.Δt
	Q = testfns.time.zTransformedTermCount
	rho = testfns.time.contourRadius

	test_spatial_basis  = testfns.space
	trial_spatial_basis = trialfns.space

	# Compute the Z transformed sequence.
	# Assume that the operator applied on the conjugate of s is the same as the
	# conjugate of the operator applied on s,
	# so that only half of the values are computed
	Qmax = Q>>1+1
	M = numfunctions(test_spatial_basis)
	N = numfunctions(trial_spatial_basis)
	Tz = promote_type(scalartype(cqop), scalartype(testfns), scalartype(trialfns))
	Zz = Vector{Array{Tz,2}}(undef,Qmax)

	if threading.parameters[1]==:single
		@info "single threaded assembly"
		for q = 0:Qmax-1
			# Build a temporary matrix for each eigenvalue
			s = laplace_to_z(rho, q, Q, Δt, method)
			Zz[q+1] = assemble(LaplaceOp(s), test_spatial_basis, trial_spatial_basis)
		end
	elseif threading.parameters[1]==:multi
		@info "multi threading assembly"
		T = Threads.nthreads()
		Qsplits = [round(Int, s) for s in range(0, stop=Qmax, length=T+1)]
		#function that correctly stores the computed values of Z for each Q
		Threads.@threads for idx in 1:T
			for q in Qsplits[idx]:Qsplits[idx+1]-1
				s = laplace_to_z(rho, q, Q, Δt, method)
				Zz[q+1] = assemble(LaplaceOp(s), test_spatial_basis, trial_spatial_basis)
			end
		end
	end

	# return the inverse Z transform
	kmax = Q
	T = real(Tz)
	Z = zeros(T, M, N, kmax)
	#return Zz
	print("Inverse z transform dots out of 10:")
	for q = 0:kmax-1
		Z[:,:,q+1] = real_inverse_z_transform(q, rho, Q, Zz)
		isinteger((q+1)*10/kmax) ? print(".") : nothing
	end
	return ConvolutionOperators.DenseConvOp(Z)

end