using LinearAlgebra, Tables

"""
    resgrad!(μ, resid, Jac, df, β)

Update the mean response, `μ`, the residual, `resid` and the transpose of
the Jacobian, `Jac`, given `df`, an object like a `DataFrame` for which `Tables.rows`
returns a NamedTuple including names `dose`, `time`, and `conc`.

Returns the sum of squared residuals.
"""
function resgrad!(μ, resid, Jac, cov, φ, f)
    rss = zero(eltype(μ))
    for (i,r) in enumerate(Tables.rows(cov))
        μi, g = f(r, φ)
        μ[i] = μi
        resi = resid[i] = r.conc - μi
        rss += abs2(resi)
        Jac[i] = g
    end
    rss
end

function resgradre!(μ, resid, Jac, df, β, b)
    ϕ = similar(β)
    grad = similar(β)
    oldsubj = zero(eltype(df.subj))
    for (i,r) in enumerate(Tables.rows(df))
        if r.subj ≠ oldsubj
            oldsubj = r.subj
            for j in eachindex(β)
                ϕ[j] = β[j] + b[j, oldsubj]
            end
        end
        μ[i] = onecompartment!(grad, r.dose, r.time, ϕ)
        resid[i] = r.conc - μ[i]
        for j in eachindex(grad)
            Jac[i,j] = grad[j]
        end
    end
    resid
end

## Simple nonlinear least squares
function increment!(δ, resid, Jac)
    m, n = size(Jac)
    fac = qr!(Jac)
    lmul!(fac.Q', resid)
    for i in eachindex(δ)
        δ[i] = resid[i]
    end
    ldiv!(UpperTriangular(fac.R), δ)
    sum(abs2, view(resid, 1:n)) / sum(abs2, view(resid, (n+1):m))
end

function nls!(β::NamedTuple, df, f)
    n = size(df, 1)
    Jac = repeat([β], n)
    μ = zeros(typeof(β[1]), n)
    resid = similar(μ)
    rss = resgrad!(μ, resid, Jac, df, β, f)
    Jacobian = similar(μ, n, length(β))
    sch = Tables.schema(β)
    for (i, r) in enumerate(Jac)
        Tables.eachcolumn(sch, r) do val, j::Int, nm::symbol
            Jacobian[i,j] = val
        end
    end
    qrJac = qr!(Jacobian)
    rss
end
#=
    δ = similar(β)     # parameter increment
    b = copy(β)        # trial parameter value
    n = size(df, 1)
    μ = similar(β, n)
    resid = similar(μ)
    Jac = similar(β, (n, 3))
    oldrss = resgrad!(μ, resid, Jac, df, β)
    cvg = increment!(δ, resid, Jac)
    tol = 0.0001       # convergence criterion tolerance
    minstep = 0.001    # minimum step factor
    maxiter = 100      # maximum number of iterations
    iter = 1
    while cvg > tol && iter ≤ maxiter
        step = 1.0     # step factor
        b .= β .+ step .* δ
        rss = resgrad!(μ, resid, Jac, df, b)
        while rss > oldrss && step ≥ minstep  # step-halving to ensure reduction of rss
            step /= 2
            b .= β .+ step .* δ
            rss = resgrad!(μ, resid, Jac, df, b)
        end
        if step < minstep
            throw(ErrorException("Step factor reduced below minstep of $minstep"))
        end
        copy!(β, b)
        cvg = increment!(δ, resid, Jac)
        iter += 1
        oldrss = rss
    end
    if iter > maxiter
        throw(ErrorException("Maximum number of iterations, $maxiter, exceeded"))
    end
    (lK = b[1], lKa = b[2], lCl = b[3])
end
=#
function updateTerms!(m::LinearMixedModel, resid, Jac)
    copyto!(first(m.feterms).x, Jac)
    copyto!(last(m.feterms).x, resid)
    re = first(m.reterms)
    transpose!(re.z, Jac)
    copyto!(re.adjA.nzval, re.z)
    terms = vcat(m.reterms, m.feterms)
    k = length(terms)
    A = m.A
    for j in 1:k
        for i in j:k
            mul!(A[Block(i,j)], terms[i]', terms[j])
        end
    end
    updateL!(m)
    L = m.L
    prss = abs2(first(L[Block(3,3)]))
    prss, (sum(abs2, L[Block(3,1)]) + sum(abs2, L[Block(3,2)])) / prss
end

function pnls!(m::LinearMixedModel, β, b, df)
    δ = fill!(similar(β), 0)     # parameter increment
    β₀ = copy(β)       # trial parameter value
    b₀ = copy(b)
    δb = [fill!(similar(b), 0)]
    n = size(df, 1)
    μ = similar(β, n)
    resid = similar(μ)
    Jac = similar(β, (n, 3))
    resgradre!(μ, resid, Jac, df, β, b)
    oldprss, cvg = updateTerms!(m, resid, Jac)
    fixef!(δ, m)
    ranef!(δb, m, δ, false)
    tol = 0.0001       # convergence criterion tolerance
    minstep = 0.001    # minimum step factor
    maxiter = 100      # maximum number of iterations
    iter = 1
    while cvg > tol && iter ≤ maxiter
        step = 1.0                              # step factor
        β .= β₀ .+ step .* δ
        b .= b₀ .+ step .* first(δb)
        resgradre!(μ, resid, Jac, df, β, b)
        prss, cvg = updateTerms!(m, resid, Jac)
        while prss > oldprss && step ≥ minstep  # step-halving to ensure reduction of rss
            step /= 2
            β .= β₀ .+ step .* δ
            b .= b₀ .+ step .* first(δb)
            resgradre!(μ, resid, Jac, df, β, b)
            prss, cvg = updateTerms!(m, resid, Jac)
        end
        if step < minstep
            throw(ErrorException("Step factor reduced below minstep of $minstep"))
        end
        copyto!(β₀, β)
        copyto!(b₀, b)
        fixef!(δ, m)
        ranef!(δb, m, δ, false)
        iter += 1
        oldprss = prss
    end
    if iter > maxiter
        throw(ErrorException("Maximum number of iterations, $maxiter, exceeded"))
    end
    objective(lmm)
end
