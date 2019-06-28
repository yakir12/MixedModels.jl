"""
    onecompartment(covariates, pars)

Return the predicted concentration for a one-compartment model given `covariates`, with
properties `dose` and `time` and parameters `pars` with properties `K`, `Ka`, and `Cl`.
The gradient is also returned as a `NamedTuple`.

This function is based on the output of the R function call
    deriv(~ dose*K*Ka*(exp(-K*time)-exp(-Ka*time))/(Cl*(Ka - K)), c("K","Ka","Cl"))
which performs symbolic differentiation followed by common subexpression elimination.
"""
function onecompartment(covariates, pars)
    dose = covariates.dose
    t = covariates.time
    K = pars.K
    Ka = pars.Ka
    Cl = pars.Cl
    _expr1 = dose * K
    _expr2 = _expr1 * Ka
    _expr5 = exp(-K * t)
    _expr8 = exp(-Ka * t)
    _expr9 = _expr5 - _expr8
    _expr10 = _expr2 * _expr9
    _expr11 = Ka - K
    _expr12 = Cl * _expr11
    _expr21 = abs2(_expr12)
    _expr22 = _expr10 * Cl/_expr21
    _expr10/_expr12,
    (K = (dose * Ka * _expr9 - _expr2 * (_expr5 * t))/_expr12 + _expr22,
     Ka = (_expr1 * _expr9 + _expr2 * (_expr8 * t))/_expr12 - _expr22,
     Cl = -(_expr10 * _expr11/_expr21))
end

function onecompartmentlog(covariates, pars)
    K = exp(pars.lK)
    Ka = exp(pars.lKa)
    Cl = exp(pars.lCl)
    μ, g = onecompartment(covariates, (K=K, Ka=Ka, Cl=Cl))
    μ, (lK = g.K*K, lKa = g.Ka*Ka, lCl = g.Cl*Cl)
end
