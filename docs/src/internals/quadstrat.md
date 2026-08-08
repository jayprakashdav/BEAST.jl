# Quadrature strategies

There are many ways to approximately compute the singular integrals that appear in boundary element discretisations of surface and volume integral euqations.

BEAST.jl is configured to select reasonable defaults, but advanced users may want to select their own quadrature rules. This section provides information on how to do this.

!!! warning
    TODO: discuss directly available/implemented singularity treatment

## quaddata and integrate!

Numerical quadrature is governed by a pair of functions that need to be designed to work together:

- `quaddata`: this function is executed before the assembly loop is entered. It's job is to compute all data needed for quadrature that the developer wants to be cached. Typically this is all geometric information such as the parametric and cartesian coordinates of all quadratures rules for all elemenents. It makes sense to cache this data as it will be used many times over in the double for loop that governs assembly. Typically, near singular interactions require more careful treatment than far interactions. This means that multiple quadrature rules per elements can be required. In such cases, the developer may want to opt to sture quadrature points and weights for all these rules. The function returns a quaddata object that holds all the cached data.
- `integrate!`: for `IntegralOperator` assembly, the method of `integrate!` dispatching on the quadrature strategy is executed inside the assembly hotloop. It receives a pair of elements and the quaddata object as its arguments. Based on this, the relevant cached data is extracted and a quadrature-rule object is built, describing which numerical integration routine is appropriate for this pair of elements. Rather than returning that object, this method immediately hands it to the *other* method of `integrate!` (the one dispatching on `action` and the rule object) from within the same method/branch, so that the choice of which integration routine to run resolves statically rather than through a dynamic dispatch. Pass `action=ReturnQRule()` to get the rule object back unevaluated instead (the default is `action=ApplyIntegrate()`, which evaluates it). Before BEAST 2.10 this was two separate functions, `quadrule` (returning the rule object) and `momintegrals!` (consuming it); `quadrule` remains a distinct, still-existing function for a few operator families outside `IntegralOperator` (local operators, excitations, farfield/nearfield postprocessing).

  A third action, `ApplyIntegrateNonConforming()`, exists for `NonConformingOverlapQRule` and `NonConformingTouchQRule`. Their own `integrate!` methods are reached only after the generic dispatcher has already stripped the test/trial `Space` down to local refspaces, so there's no full `Space` left for the `ApplyIntegrate` path to work with; `ApplyIntegrateNonConforming` skips straight to the refspace-level `integrate!` instead.

## quadstrat

The pair of `quaddata`/`integrate!` methods that is used is determined by the type of the operator and finite elements, and a `quadstrat` object. This object is passed to the assembly routine and passed on to `quaddata` and `integrate!`, so it can be considered during dispatch.

Parameters, such as those that determine the accuracy of the numerical quadrature, are part of the runtime payload of the quadstrat object. This is usefull when the user is interested on the impact of these parameters on the performance and the accuracy of the solver without having to supply a new pair of `quadstrat`/`quaddata` methods for each possible value of these parameters.

Roughly this leads to the following (simplified) assembly routine:

```julia
function assemble(op, tfs, bfs, store; quadstrat=QS)

    tad, tels = assemblydata(op, tfs)
    bad, bels = assemblydata(op, bfs)
    
    qd = quadata(op,tels,bels,quadstrat)
    for tel in tels
        for bel in bels
            zlocal = zeros(...)
            integrate!(op,tel,bel,qd,quadstrat,zlocal,tfs,tel,bfs,bel; action=ApplyIntegrate())

            for i in axes(zlocal,1)
                for j in axes(zlocal,2)
                    m, a = tad[tel,i]
                    n, b = bad[bel,j]
                    store(a*zlocal[i,j]*b,m,n)
end end end end end
```

It is conceivable that the types and functions described above look like this:

```julia
struct DoubleNumQS
    test_precision
    trial_precision
end

function quaddata(op, tels, bels, quadstrat::DoubleNumQS)
    tqps = [quadpoints(tel,precision=quadstrat.test_precision) for tel in tels]
    bqps = [quadpoints(bel,precision=quadstrat.basis_precision) for bel in bels]
    return (test_quadpoints=tqps, basis_quadpoints=bqps)
end


struct DoubleNumQR
    test_quadpoints
    trial_quadpoints
end

struct HighPrecisionQR end

function integrate!(op, tel, bel, qd, quadstrat::DoubleNumQs,
        out=nothing, tfs=nothing, tptr=nothing, bfs=nothing, bptr=nothing;
        action::QuadRuleAction=ApplyIntegrate())
    if wellseparated(tel, bel)
        qr = DoubleNumQR(qd.test_quadpoints[tel], qd.basis_quadpoints[bel])
    else
        qr = HighPrecisionQR(tel, bel)
    end
    integrate!(action, out, op, tfs, tptr, tel, bfs, bptr, bel, qr)
end

function integrate!(out, op, tfs, tptr, tel, bfs, bptr, bel, qr::DoubleNumQR)
    ...
end

function integrate!(out, op, tfs, tptr, tel, bfs, bptr, bel, qr::HighPrecisionQR)
    ...
end
```

Building `qr` and calling `integrate!` on it from within the same method/branch (rather than returning `qr` for a separate caller to dispatch on) is what avoids the dynamic dispatch: `qr`'s type is known to the compiler at the point of that call, even though it varies across calls depending on the runtime geometry of `tel`/`bel`. Passing `action=ReturnQRule()` skips that call and just returns `qr`, matching the pre-2.10 behaviour of the old `quadrule` function.