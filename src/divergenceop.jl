struct DivergenceOp <: LocalOperator
    outer::Bool
    inner::Bool
end
struct CurlOp <: LocalOperator
    outer::Bool
    inner::Bool
end

kernelvals(biop::DivergenceOp, x) = nothing
kernelvals(biop::CurlOp, x) = nothing

function integrand(op::DivergenceOp, kernel, x, test_fn, trial_fn)
    if op.outer
        dot(test_fn.divergence, trial_fn.value)
    elseif op.inner
        dot(test_fn.value, trial_fn.divergence)
    end
end

function integrand(op::CurlOp, kernel, x, test_fn, trial_fn)
    if op.outer
        dot(test_fn.curl, trial_fn.value)
    elseif op.inner
        dot(test_fn.value, trial_fn.curl)
    end
end

scalartype(op::DivergenceOp) = Union{}
scalartype(op::CurlOp) = Union{}