using BEAST
using BEAST.StaticArrays

mutable struct SurfaceFourierTransform{T}
    xi::T
end

"""
    SurfaceFourierTransform(;gamma)

Create the surface fourier transform operator where gamma is the fourier domain parameter.
"""
function SurfaceFourierTransform(;xi=nothing)
    SurfaceFourierTransform(gamma)
end

defaultquadstrat(op::SurfaceFourierTransform, basis) = BEAST.SingleNumQStrat(8)
quaddata(op::SurfaceFourierTransform,rs,els,qs::BEAST.SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::SurfaceFourierTransform,refspace,p,y,q,el,qdata,qs::BEAST.SingleNumQStrat) = qdata[1,q]

function kernelvals(op::SurfaceFourierTransform,y)

        expn = exp(-im*dot(op.xi,y.cart))
        krn = (phase=expn)
end

function integrand(op::SurfaceFourierTransform, krn, fp)

    j = fp.value
    ρ = fp.divergence

    j*krn
end

function surfaceft(op::SurfaceFourierTransform, points, coeffs, basis;
	type=SVector{3,ComplexF64},
	quadstrat=defaultquadstrat(op, basis))

	ff = zeros(type, size(points))
	store(v,m,n) = (ff[m] += v*coeffs[n])
	surfaceft!(store, op, points, basis; type, quadstrat)
	return ff
end


function surfaceft!(store, op::SurfaceFourierTransform, points, basis;
	type=SVector{3,ComplexF64},
	quadstrat=defaultquadstrat(op, basis))

	z = zeros(type,length(points))

	els, ad = assemblydata(basis)
	rs = refspace(basis)

	geo = geometry(basis)
	nf = numfunctions(rs, domain(chart(geo, first(geo))))
	zlocal = Array{type}(undef, nf)
	qdata = quaddata(op,rs,els,quadstrat)

	print("dots out of 10: ")

	todo, done, pctg = length(points), 0, 0

	for (p,y) in enumerate(points)
        op.xi = y
		for (q,el) in enumerate(els)

			fill!(zlocal,zero(type))
			qr = quadrule(op,rs,p,y,q,el,qdata,quadstrat)
			fouriertransformlocal!(zlocal,op,rs,y,el,qr)

			# assemble from local contributions
			for (r,z) in enumerate(zlocal)
                for (n,b) in ad[q,r]
					store(z*b,p,n)
				end
			end
		end

		done += 1
		new_pctg = round(Int, done / todo * 100)
		if new_pctg > pctg + 9
				#println(todo," ",done," ",new_pctg)
				print(".")
				pctg = new_pctg
		end
	end

	println("")

end

function fouriertransformlocal!(zlocal,op,refspace,y,el,qr)

    for q in qr
        x = q.point
        F = q.value
        dx = q.weight

        krn = kernelvals(op, x)
        for r in eachindex(zlocal)
            zlocal[r] += integrand(op,krn,F[r]) * dx
        end

    end

end