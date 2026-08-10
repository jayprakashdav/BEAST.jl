## ConductivityTDFunc3 is named to fixed point iteration ##
mutable struct ConductivityTDFunc3{M,P,T} <: ConductivityTD_Functionaltype{T}
    e::M #ElementaryCharge
    k_B::M #BoltzmannConstant
    h_cross::M #ReducedPlanckConstant
    Z0::M #FreespaceImpedance
    delta::M #miniband width
    l_p::M #lattice period length
    en_DL_y::M #charge times doping density times thickness
    wb_red::M #Reduced bloch frequency (electric field to be multiplied)
    v_p::M #Peak drift velocity in Esaki-Tsu formula
    L_y::M #Thickness of superlattice
    tp::M #temperature
    nu_e::M #elastic collision frequency
    nu_p::M #inelastic collision frequency
    scl::M #SI to normalized scaling
    efield::P
    jflux::P
end

mutable struct ConductivityTDOp3 <: ConductivityTD_Operatortype
    op::ConductivityTDFunc3
end

ConductivityTDFunc3(e::M, k_B::M, h_cross::M, Z0::M, delta::M, l_p::M, en_DL_y::M, wb_red::M, v_p::M, L_y::M, tp::M, nu_e::M, nu_p::M, scl::M, efield::P, jflux::P) where {M, P} = ConductivityTDFunc3{M,P,eltype(efield[1])}(e, k_B, h_cross, Z0, delta, l_p, en_DL_y, wb_red, v_p, L_y, tp, nu_e, nu_p, scl, efield, jflux)

scalartype(p::ConductivityTDFunc3) = eltype(p.efield[1])
scalartype(p::ConductivityTDOp3) = eltype(p.op.efield[1])

function conductivityfunc3(delta, l_p, n_D, L_y, tp, nu_e, nu_p, scl, I_0, I_1)
    #constant definition
    e = 1.60217e-19 #ElementaryCharge
    k_B = 1.38065e-23 #BoltzmannConstant
    h_cross = 1.05457e-34 #ReducedPlanckConstant
    Z0 = 376.7303 #FreespaceImpedance

    delta = delta*e
    en_DL_y = e*n_D*L_y
    wb_red = e*l_p/h_cross
    v_p = (delta/4)*(l_p/h_cross)*(I_1/I_0)

    efield = fill(SVector(0.0,0.0,0.0), (2,2))
    charge = fill(SVector(0.0), (2,2))
    ConductivityTDFunc3(e,k_B,h_cross,Z0,delta, l_p, en_DL_y, wb_red,v_p, L_y, tp, nu_e, nu_p, scl, efield, efield)
end

function (f::ConductivityTDFunc3)(cell, cqdpt, mp)
    ei = f.efield[cell, cqdpt]
    # e -> 0 analytic limit: j -> 0 (see conductivity2.jl)
    if norm(ei) < 1e-12
        return zero(ei)
    end
    #return f.chr(norm(ei))*ei
    E = norm(ei)*f.scl #normalized to SI unit

    ν_e = f.nu_e
    ν_p = f.nu_p
    v_p = f.v_p

    ω_B = f.wb_red*E
    v_d = v_p*(2*ν_e*ω_B)./(ν_e*(ν_e+ν_p) .+ (ω_B).^2)
    v_d2 = 7.5*v_p*sqrt(ν_e/(ν_e+ν_p))*(E/1e7).^1.4
    j_s = f.en_DL_y*v_d
    j_s2 = f.en_DL_y*v_d2

    scaled_j_st = f.Z0*(j_s+j_s2)*(1/f.scl) #from SI to normalized
    return scaled_j_st*ei/norm(ei)
end

#= function (f::ConductivityTDFunc)(cell, cqdpt, mp)
    y = mp.cart[2]/10
    ei = [0.0,1.0*sqrt(y),0.0]
    #return ei
    return ei*f.chr.(norm(ei))
end =#

function kernelvals(f::ConductivityTDOp3, mp, cell, cqdpt)
    ei = f.op.efield[cell, cqdpt]
    # e -> 0 analytic limit: sigma -> ohmic sheet conductance sigma_s(0)
    # (see conductivity2.jl); threshold in normalized units.
    if norm(ei) < 1e-12
        return f.op.Z0 * f.op.en_DL_y * f.op.v_p * 2 * f.op.nu_e * f.op.wb_red /
               (f.op.nu_e * (f.op.nu_e + f.op.nu_p))
    end
    #dsigma = f.op.dchr(norm(ei))*kron(ei, ei')/norm(ei)+(f.op.chr(norm(ei)))*I(3)
    #return dsigma
    Ei = ei*f.op.scl
    E = norm(ei)*f.op.scl #normalized to SI unit

    ν_e = f.op.nu_e
    ν_p = f.op.nu_p
    v_p = f.op.v_p    

    ω_B = f.op.wb_red*E

    v_d = v_p*(2*ν_e*ω_B)/(ν_e*(ν_e+ν_p) + (ω_B)^2)
    v_d2 = 7.5*v_p*sqrt(ν_e/(ν_e+ν_p))*(E/1e7)^1.4
    j_s = f.op.en_DL_y*v_d
    j_s2 = f.op.en_DL_y*v_d2

    scaled_j_st = f.op.Z0*(j_s+j_s2)*(1/f.op.scl) #from SI to normalized
    #dj_s1 = f.op.en_DL_y*v_p*2*ν_e*f.op.wb_red*(ν_e*(ν_e+ν_p) - (ω_B)^2)/(ν_e*(ν_e+ν_p) + (ω_B)^2)^2
    #dj_s2 = f.op.en_DL_y*7.5*1.4*1e-7*v_p*sqrt(ν_e/(ν_e+ν_p))*(E/1e7)^(0.4)
    #dj_s1 = -2*f.op.en_DL_y*v_p*2*ν_e*(f.op.wb_red)^2*((ω_B))/(ν_e*(ν_e+ν_p) + (ω_B)^2)^2
    #dj_s2 = f.op.en_DL_y*7.5*0.4*(1e-7)^2*v_p*sqrt(ν_e/(ν_e+ν_p))*(E/1e7)^(-0.6)
    #dj_s1 = e*N_D*L_y*v_p*2*ν_e*(e*(L/h_cross))^2*((ω_B))/(ν_e*(ν_e+ν_p) + (ω_B)^2)^2*f.op.scl
    #dj_s2 = e*N_D*L_y*7.5*0.4*v_p*1e-7*(E/1e7)^(-0.6)
    #dj_s = f.op.Z0*(dj_s1*+dj_s2)*f.op.scl
    #dsigma1 = dj_s*kron(ei, ei')/norm(ei)
    sigma = (scaled_j_st/norm(ei))
    return sigma
end

function assemble_local_matched!(biop::ConductivityTDOp3, tfs::BEAST.Space, bfs::BEAST.Space, store;
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

function cellinteractions(biop::ConductivityTDOp3, trefs::U, brefs::V, cell, qr, p) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

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

function cellinteractions_matched!(zlocal, biop::ConductivityTDOp3, trefs, brefs, cell, qr, p)

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

function assemble!(field::ConductivityTDFunc3, tfs::BEAST.Space, store;
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

function celltestvalues(tshs::BEAST.RefSpace{T}, t, tcell, field::ConductivityTDFunc3, qr) where {T}

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