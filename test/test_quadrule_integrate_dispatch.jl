using Test
using LinearAlgebra
using CompScienceMeshes, BEAST

@testitem "integrate!: fusing quadrule avoids allocating a boxed quadrature rule" begin
    using CompScienceMeshes
    using LinearAlgebra

    fn = joinpath(dirname(pathof(BEAST)), "../examples/assets/sphere45.in")
    m = BEAST.readmesh(fn)
    X = raviartthomas(m)
    op = Maxwell3D.singlelayer(gamma=1.0)

    # This quadrature strategy's `integrate!` method picks between five distinct
    # rule types (CommonFace/CommonEdge/CommonVertex/WiltonSERule/DoubleQuadRule)
    # for a given pair of triangles -- enough concrete types to defeat Julia's
    # union-splitting optimization if the chosen rule were allowed to escape as a
    # plain return value, instead of being consumed within the branch that built it.
    qs = BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 4, 3)

    tels, tad = assemblydata(X)
    bels, bad = assemblydata(X)
    trefs = brefs = refspace(X)
    qd = BEAST.quaddata(op, trefs, brefs, tels, bels, qs)

    zlocal = zeros(scalartype(op, X, X),
        numfunctions(trefs, CompScienceMeshes.domain(tels[1])),
        numfunctions(brefs, CompScienceMeshes.domain(bels[1])))

    fused!(p, q) = begin
        tcell, bcell = tels[p], bels[q]
        fill!(zlocal, 0)
        BEAST.integrate!(op, trefs, brefs, p, tcell, q, bcell, qd, qs,
            zlocal, X, p, X, q; action=BEAST.ApplyIntegrate())
    end

    # Reproduces the pre-2.10 two-step call: build the rule, return it to the
    # caller, and only then dispatch on its (now boxed) runtime type to apply it.
    unfused!(p, q) = begin
        tcell, bcell = tels[p], bels[q]
        fill!(zlocal, 0)
        qrule = BEAST.integrate!(op, trefs, brefs, p, tcell, q, bcell, qd, qs;
            action=BEAST.ReturnQRule())
        BEAST.integrate!(zlocal, op, X, p, tcell, X, q, bcell, qrule)
    end

    # Exercise several branches (well-separated, and the closest few triangles,
    # which will hit the near/touching branches) rather than relying on any one
    # of them individually: the exact byte count of a given branch can include
    # allocations unrelated to this fix (e.g. from SauterSchwabQuadrature's own
    # internals), but fusing should never allocate *more* than the two-step
    # call it replaces, and for the union-boxing this test targets, strictly less.
    center(el) = sum(el.vertices) / length(el.vertices)
    p = 1
    distances = sortperm([norm(center(tels[p]) - center(bels[q])) for q in eachindex(bels)])
    test_qs = unique([distances[1], distances[2], distances[end]])

    for q in test_qs
        fused!(p, q);
        unfused!(p, q) # compile before measuring
    end

    fused_bytes = sum(q -> @allocated(fused!(p, q)), test_qs)
    unfused_bytes = sum(q -> @allocated(unfused!(p, q)), test_qs)

    @show fused_bytes, unfused_bytes

    @test fused_bytes < unfused_bytes
end
