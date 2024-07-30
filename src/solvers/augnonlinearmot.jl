function td_nl_solve_augmented(Ts, Th, E, s, tfs, bfs)
    Δt = bfs[1].time.Δt
    Nt = bfs[1].time.Nt
    nrows = cumsum(numfunctions.(spatialbasis.(tfs)))
    ncols = cumsum(numfunctions.(spatialbasis.(bfs)))

    Zs = BEAST.assemble(Ts, tfs[1], bfs[1])
    Zs0=BEAST.ConvolutionOperators.timeslice(Zs,1)
    Zh = BEAST.assemble(Th, tfs[1], bfs[2])
    Zh0=BEAST.ConvolutionOperators.timeslice(Zh,1)
    
    M = zeros(Float64, sum(nrows), sum(ncols))
    M[1:nrows[1], 1:ncols[1]] = Zs0 #M[1,1]
    M[1:nrows[1], ncols[1]+1:ncols[2]] = Zh0 #M[1,2]
    
    Div = DivergenceOp(false, true)
    Z_div = BEAST.assemble(Div, tfs[2].space, bfs[1].space)
    Z_diff = BEAST.assemble(Identity(), tfs[2].space, bfs[2].space)/Δt
    
    Ge = BEAST.assemble(Identity(), tfs[1].space, bfs[3].space)/Δt
    Gj = BEAST.assemble(Identity(), tfs[3].space, bfs[1].space)
    
    M[1:nrows[1], ncols[2]+1:ncols[3]] = Ge #M[1,3]
    M[nrows[1]+1:nrows[2], 1:ncols[1]] = Z_div #M[2,1]
    M[nrows[1]+1:nrows[2], ncols[1]+1:ncols[2]] = Z_diff #M[2,2]
    M[nrows[2]+1:nrows[3], 1:ncols[1]] = Gj #M[3,1]

    B = BEAST.assemble(E, tfs[1])
    return motnlaug(M, Zs, Zh, Z_diff, Ge, Gj, B, s, tfs, bfs) #march on time nonlinear augmented
end

function motnlaug(M, Zs, Zh, Z_diff, Ge, Gj, B, s,tfs, bfs)
    nrows = cumsum(numfunctions.(spatialbasis.(tfs)))
    ncols = cumsum(numfunctions.(spatialbasis.(bfs)))
    
    T = scalartype(bfs[1])
    sol = zeros(T,sum(ncols)-1)
    xj = zeros(T,numfunctions(bfs[1].space),Nt)
    xq = zeros(T,numfunctions(bfs[2].space),Nt)
	xe = zeros(T,numfunctions(bfs[3].space),Nt)
    
    xqprev = zeros(T,numfunctions(bfs[2].space))
	xeprev = zeros(T,numfunctions(bfs[3].space))
    
    yj = zeros(T,numfunctions(tfs[1].space))
    yq = zeros(T,numfunctions(tfs[1].space))
	ye = zeros(T,numfunctions(tfs[1].space))
    
    jk_start = 2
    qk_start = 2
    ek_start = 2
    jk_stop = Nt
    qk_stop = Nt
    ek_stop = BEAST.numfunctions(eq2.trial_space_dict[1].time)+1
    csxj = zeros(T,numfunctions(bfs[1].space),Nt)
    csxq = zeros(T,numfunctions(bfs[2].space),Nt)
    csxe = zeros(T,numfunctions(bfs[3].space),Nt)
    
    se = BEAST.ConductivityTDOp(s)
    sq = BEAST.ConductivityTDOpch(s)

    bs = zeros(T, N)
    #= iZ = BEAST.GMRESSolver(Z0, restart=0, reltol=1e-6)
    invZ = inv(Z0)
    C1 = G_j0*invZ
    C2 = C1*Ġ0 =#
    #= try =#
    for i in 1:Nt
        println(i)
        R = inc[:,i]
        fill!(yj,0)
        fill!(ye,0)
        BEAST.ConvolutionOperators.convolve!(yj,Zs,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(yq,Zh,xq,csxq,i,qk_start,qk_stop)
        if i>1 
            ye = Ge*xe[:,i-1]
        end
        iter_max = 100
        for l in 1:iter_max
            if (l==1) && (i>1)
                xeprev = xe[:,i-1]
                xqprev = xr[:,i-1]
                update!(s, xj, xq, xe, i-1, tfs, bfs)
            else
                xeprev = xe[:,i]
                xqprev = xq[:,i]
                update!(s, xj, xq, xe, i, tfs, bfs)
            end
            
            Qe = BEAST.assemble(se, tfs[3].space, bfs[3].space)
            Qq = BEAST.assemble(sq, tfs[3].space, bfs[2].space)
            M[nrows[2]+1:nrows[3], ncols[1]+1:ncols[2]] = Qq #M[3,2]
            M[nrows[2]+1:nrows[3], ncols[2]+1:ncols[3]] = Qe #M[3,3]

            bs = BEAST.assemble(s, tfs[3])
            rhs1 = R - yj - yq + ye
            rhs2 = Z_diff*xq[:,i-1]
            rhs3 = bs + Qe*xeprev + Qq*xqprev
            b = [rhs1; rhs2; rhs3]
            
            sol = inv(Matrix(V0))*b

            xj[:,i] = sol[1:ncols[1]] #currents
            xq[:,i] = sol[ncols[1]+1:ncols[2]] #charges
            xe[:,i] = sol[ncols[2]+1:ncols[3]] #electric fields

            #= rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev - C1*rhs1
            iw = BEAST.GMRESSolver(C2-Q, restart=0, reltol=1e-6)
            u, ch = solve(iw, rhs2)
            xe[:,i] .= u
            xj[:,i] .= invZ*(rhs1+Ġ0*xe[:,i]) =#

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
            csxq[:,i] .= csxq[:,i-1] .+ xq[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxq[:,i] .= xq[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xq, xe
    #= catch e
        return xj, xe, xj_all, xe_all
    end =#
end
