function td_solve_augmented(eq1::BEAST.DiscreteEquation, eq2::BEAST.DiscreteEquation, cq=true)
    X = _spacedict_to_directproductspace(eq1.test_space_dict)
    Y = _spacedict_to_directproductspace(eq1.trial_space_dict)
    #Z = assemble(eq1.equation.lhs, X, Y)
    if cq
        V = eq1.trial_space_dict[1]
        W = eq1.test_space_dict[1]
        Δt = V.time.Δt
        Nt = V.time.Nt
        LaplaceTs(s::T) where {T} = MWSingleLayer3D(s/c, -s*s/c, T(0.0))
        LaplaceTh(s::T) where {T} = AugmentedMaxwellOperator3D(s/c, T(0.0), -s*c)
        Div_cq = DivergenceOp(false, true)
        Id = Identity()
        Ts_cq = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
        Ts_cq.hs_weight=0.0
        Th_cq = BEAST.AMWSingleLayerTDIO(1.0,1.0,1.0,0,2)
        Z_Ts = BEAST.assemble(Ts_cq, V, V)
        Z_Th = BEAST.assemble(Th_cq, W, V)
        Z_Div = BEAST.assemble(Div_cq, W.space, V.space)
        Z_Diff = BEAST.assemble(Id, W.space, W.space)/Δt
        B = BEAST.assemble(eq1.equation.rhs, eq1.trial_space_dict)
        return motnlaug(Z_Ts, Z_Th, Z_Div, Z_Diff, B)#march on time nonlinear augmented
    end
end

function marchonintimenl(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    Z0 = zeros(eltype(Z), size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = zeros(eltype(Ġ), size(Ġ)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Ġ0,Ġ,1) 
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    V0 = zeros(N+Ne,N+Ne)
    V0[1:N, 1:N] = Z0
    V0[1:N, N+1:N+Ne] = -Ġ0
    V0[N+1:N+Ne, 1:N] = G_j0
    sol = zeros(2*N)
    xj_all = zeros(N)
    xe_all = zeros(N)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq2.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDOp(σ)
    bσ = zeros(T, N)
    iZ = BEAST.GMRESSolver(Z0, restart=0, reltol=1e-6)
    invZ = inv(Z0)
    C1 = G_j0*invZ
    C2 = C1*Ġ0
    #= try =#
    for i in 1:Nt
        println(i)
        R = inc[:,i]
        fill!(yj,0)
        fill!(ye,0)
        BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 100
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            V0[N+1:N+Ne,N+1:N+Ne] = -Q
            b = [rhs1; rhs2]
            sol = inv(Matrix(V0))*b
            xj[:,i] = sol[1:N]
            xe[:,i] = sol[N+1:end]
            #= rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev - C1*rhs1
            iw = BEAST.GMRESSolver(C2-Q, restart=0, reltol=1e-6)
            u, ch = solve(iw, rhs2)
            xe[:,i] .= u
            xj[:,i] .= invZ*(rhs1+Ġ0*xe[:,i]) =#
            #xe_all = hcat(xe_all, xe[:,i])
            #xj_all = hcat(xj_all, xj[:,i])
            println("norm xe ", norm(xe[:,i]-xeprev)/M)
            if norm(xe[:,i]-xeprev)/M < 1e-4
                break
            end
            #= println("L_inf norm of diff of xe ", maximum(abs.((xe[:,i]-xeprev))))
            if maximum(abs.((xe[:,i]-xeprev)))/M < 1e-10
                break
            end =#
        end
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xe, xj_all, xe_all
    #= catch e
        return xj, xe, xj_all, xe_all
    end =#
end
