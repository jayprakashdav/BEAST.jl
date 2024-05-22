function marchonintimeaug(Zs, Zh, Zdiv, Zdiff, B)

    Zs0 = ConvolutionOperators.timeslice(Zs,1)
    Zh0 = ConvolutionOperators.timeslice(Zh,1)
    T = eltype(Zs0)
    M,N = size(Zs0)
    @assert M == size(B,1)

    Nj = size(Zs0,2)
    Nq = size(Zh0,2)
    Nt = size(B,2)

    xj = zeros(T,Nj,Nt)
    yj = zeros(T,Nj)
    csxj = zeros(T,Nj,Nt)

    xq = zeros(T,Nq,Nt)
    yq = zeros(T,Nj)
    csxq = zeros(T,Nq,Nt)

    invZsys = inv([Zs0 Zh0; Matrix(Zdiv) Matrix(Zdiff)])

    for i in 1:Nt
        R = B[:,i]
        k_start = 2
        k_stop = Nt

        fill!(yj,0)
        fill!(yq,0)
        ConvolutionOperators.convolve!(yj,Zs,xj,csxj,i,k_start,k_stop)
        ConvolutionOperators.convolve!(yq,Zh,xq,csxq,i,k_start,k_stop)
        b1 = R - yj - yq
        if i==1
            b2 = zeros(T,size(xq,1))
        else
            b2 = Zdiff*xq[:,i-1]
        end
        b = [b1;b2]
        x = invZsys*b
        #@assert nj+nq = length(x)
        xj[:,i] = x[1:Nj]
        xq[:,i] = x[Nj+1:Nj+Nq]
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxq[:,i] .= csxq[:,i-1] .+ xq[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxq[:,i] .= xq[:,i]
        end

        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end

    return xj, xq
end
