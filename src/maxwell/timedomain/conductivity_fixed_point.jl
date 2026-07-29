mutable struct ConductivityTDFunc_fp{L,M,N,P,T} <: ConductivityTD_Functionaltype{T}
    chr::L
    dchr::M
    numdiffs::N
    efield::P
    jflux::P
end

mutable struct ConductivityTDOp_fp <: ConductivityTD_Operatortype
    op::ConductivityTDFunc_fp
end

ConductivityTDFunc_fp(chr::L, dchr::M, numdiffs::N, efield::P, jflux::P) where {L,M, N, P} = ConductivityTDFunc_fp{L,M,N,P,eltype(efield[1])}(chr,dchr,numdiffs,efield,jflux)

scalartype(p::ConductivityTDFunc_fp) = eltype(p.efield[1])
scalartype(p::ConductivityTDOp_fp) = eltype(p.op.efield[1])

function conductivityfunc_fp(chr::BEAST.Polynomial;numdiffs=0)
    dchr = derive(chr)
    efield = fill(SVector(0.0,0.0,0.0), (2,2))
    charge = fill(SVector(0.0), (2,2))
    ConductivityTDFunc_fp(chr, dchr, numdiffs, efield, efield)
end

function (f::ConductivityTDFunc_fp)(cell, cqdpt, mp)
    ei = f.efield[cell, cqdpt]
    return f.chr(norm(ei))*ei
end

#= function (f::ConductivityTDFunc)(cell, cqdpt, mp)
    y = mp.cart[2]/10
    ei = [0.0,1.0*sqrt(y),0.0]
    #return ei
    return ei*f.chr.(norm(ei))
end =#

function kernelvals(f::ConductivityTDOp_fp, mp, cell, cqdpt)
    ei = f.op.efield[cell, cqdpt]
    if norm(ei)==0
        ei = 1e-9.+ei
    end
    #dsigma = f.op.dchr(norm(ei))*kron(ei, ei')/norm(ei)+(f.op.chr(norm(ei)))*I(3)
    dsigma = (f.op.chr(norm(ei)))
    return dsigma
end

function assemble_local_matched!(biop::ConductivityTDOp_fp, tfs::BEAST.Space, bfs::BEAST.Space, store;
    quadstrat=BEAST.defaultquadstrat(biop, tfs, bfs), kwargs...)

    tels, tad, ta2g = BEAST.assemblydata(tfs)
    bels, bad, ba2g = BEAST.assemblydata(bfs)

    bg2a = zeros(Int, length(BEAST.geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    trefs = BEAST.refspace(tfs)
    brefs = BEAST.refspace(bfs)
    tgeo = geometry(tfs)
    tdom = domain(chart(tgeo, first(tgeo)))

    qd = BEAST.quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    verbose = length(tels) > 10_000
    verbose && print("dots out of 20: ")
    todo, done, pctg = length(tels), 0, 0
    locmat = zeros(BEAST.scalartype(biop, trefs, brefs), BEAST.numfunctions(trefs, tdom), numfunctions(brefs, tdom))
    for (p,cell) in enumerate(tels)
        P = ta2g[p]
        q = bg2a[P]
        q == 0 && continue

        qr = BEAST.quadrule(biop, trefs, brefs, cell, qd, quadstrat)
        fill!(locmat, 0)
        BEAST.cellinteractions_matched!(locmat, biop, trefs, brefs, cell, qr,p)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p,i], (n,b) in bad[q,j]
                store(a * locmat[i,j] * b, m, n)
        
        end end

        new_pctg = round(Int, (done += 1) / todo * 100)
        verbose && new_pctg > pctg + 4 && (print("."); pctg = new_pctg)
    end
end

function cellinteractions(biop::ConductivityTDOp_fp, trefs::U, brefs::V, cell, qr, p) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    for (i,q) in enumerate(qr)

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * BEAST.jacobian(mp)
        kernel = BEAST.kernelvals(biop, mp, p, i)

        for m in 1 : num_tshs
            tval = tvals[m]

            for n in 1 : num_bshs
                bval = bvals[n]

                igd = BEAST.integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end

function cellinteractions_matched!(zlocal, biop::ConductivityTDOp_fp, trefs, brefs, cell, qr, p)

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    # zlocal = zeros(Float64, num_tshs, num_bshs)
    for (i,q) in enumerate(qr)

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * BEAST.jacobian(mp)
        kernel = BEAST.kernelvals(biop, mp, p,i)
        
        for n in 1 : num_bshs
            bval = bvals[n]
            for m in 1 : num_tshs
                tval = tvals[m]

                igd = BEAST.integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd
            end
        end
    end

    return zlocal
end

function assemble!(field::ConductivityTDFunc_fp, tfs::BEAST.Space, store;
    quadstrat=BEAST.defaultquadstrat(field, tfs), kwargs...)

    tels, tad = BEAST.assemblydata(tfs)

    trefs = BEAST.refspace(tfs)
    tgeo = geometry(tfs)
    tdom = domain(chart(tgeo, first(tgeo)))

    qd = BEAST.quaddata(field, trefs, tels, quadstrat)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = BEAST.quadrule(field, trefs, t, tcell, qd, quadstrat)
        blocal = BEAST.celltestvalues(trefs, t, tcell, field, qr)

        for i in 1 : BEAST.numfunctions(trefs, tdom)
            for (m,a) in tad[t,i]
                store(a*blocal[i], m)
            end
        end

    end

end

function celltestvalues(tshs::BEAST.RefSpace{T}, t, tcell, field::ConductivityTDFunc_fp, qr) where {T}

    num_tshs = numfunctions(tshs, domain(tcell))
    interactions = zeros(BEAST.scalartype(field, tshs), num_tshs)

    num_oqp = length(qr)

    for p in 1 : num_oqp
        mp = qr[p].point

        dx =qr[p].weight

        fval = field(t,p,mp)
        tvals = qr[p].value

        for m in 1 : num_tshs
            tval = tvals[m]

            igd = BEAST.integrand(field, tval, fval)
            interactions[m] += igd * dx
        end
    end

    return interactions
end