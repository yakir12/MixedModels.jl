var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "MixedModels.jl Documentation",
    "title": "MixedModels.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MixedModels.jl-Documentation-1",
    "page": "MixedModels.jl Documentation",
    "title": "MixedModels.jl Documentation",
    "category": "section",
    "text": "CurrentModule = MixedModelsMixedModels.jl is a Julia package providing capabilities for fitting and examining linear and generalized linear mixed-effect models. It is similar in scope to the lme4 package for R.Pages = [\"constructors.md\",\n         \"extractors.md\",\n         \"bootstrap.md\"]\nDepth = 2"
},

{
    "location": "constructors.html#",
    "page": "Model constructors",
    "title": "Model constructors",
    "category": "page",
    "text": ""
},

{
    "location": "constructors.html#MixedModels.lmm",
    "page": "Model constructors",
    "title": "MixedModels.lmm",
    "category": "Function",
    "text": "lmm(f::DataFrames.Formula, fr::DataFrames.DataFrame; weights = [], contrasts = Dict())\n\nCreate a LinearMixedModel from f, a formula that contains both fixed-effects terms and random effects, and fr.\n\nThe return value is ready to be fit! but has not yet been fit.\n\n\n\n"
},

{
    "location": "constructors.html#Model-constructors-1",
    "page": "Model constructors",
    "title": "Model constructors",
    "category": "section",
    "text": "The lmm function creates a linear mixed-effects model representation from a Formula and an appropriate data type. At present a DataFrame is required but that is expected to change.lmm"
},

{
    "location": "constructors.html#Examples-of-linear-mixed-effects-model-fits-1",
    "page": "Model constructors",
    "title": "Examples of linear mixed-effects model fits",
    "category": "section",
    "text": "For illustration, several data sets from the lme4 package for R are made available in .rda format in this package. These include the Dyestuff and Dyestuff2 data sets.julia> using DataFrames, RData, MixedModels\n\njulia> const dat = convert(Dict{Symbol,DataFrame}, load(Pkg.dir(\"MixedModels\", \"test\", \"dat.rda\")));\n\njulia> # dat[:Dyestuff]\nThe columns in these data sets have been renamed for convenience in comparing models between examples. The response is always named Y. Potential grouping factors for random-effects terms are named G, H, etc."
},

{
    "location": "constructors.html#Models-with-simple,-scalar-random-effects-1",
    "page": "Model constructors",
    "title": "Models with simple, scalar random effects",
    "category": "section",
    "text": "The formula language in Julia is similar to that in R except that the formula must be enclosed in a call to the @formula macro. A basic model with simple, scalar random effects for the levels of G (the batch of an intermediate product, in this case) is declared and fit asjulia> fm1 = fit!(lmm(@formula(Y ~ 1 + (1|G)), dat[:Dyestuff]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G)\n   logLik   -2 logLik     AIC        BIC    \n -163.66353  327.32706  333.32706  337.53065\n\nVariance components:\n              Column    Variance  Std.Dev. \n G        (Intercept)  1388.3333 37.260345\n Residual              2451.2500 49.510100\n Number of obs: 30; levels of grouping factors: 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    1527.5   17.6946  86.326  <1e-99\n\n(If you are new to Julia you may find that this first fit takes an unexpectedly long time, due to Just-In-Time (JIT) compilation of the code. The second and subsequent calls to such functions are much faster.)julia> @time fit!(lmm(@formula(Y ~ 1 + (1|G)), dat[:Dyestuff2]))\n  0.001088 seconds (1.59 k allocations: 86.670 KiB)\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G)\n   logLik   -2 logLik     AIC        BIC    \n -81.436518 162.873037 168.873037 173.076629\n\nVariance components:\n              Column    Variance  Std.Dev. \n G        (Intercept)   0.000000 0.0000000\n Residual              13.346099 3.6532314\n Number of obs: 30; levels of grouping factors: 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    5.6656  0.666986 8.49433  <1e-16\n\n"
},

{
    "location": "constructors.html#Simple,-scalar-random-effects-1",
    "page": "Model constructors",
    "title": "Simple, scalar random effects",
    "category": "section",
    "text": "A simple, scalar random effects term in a mixed-effects model formula is of the form (1|G). All random effects terms end with |G where G is the grouping factor for the random effect. The name or, more generally, the expression G should evaluate to a categorical array that has a distinct set of levels. The random effects are associated with the levels of the grouping factor.A scalar random effect is, as the name implies, one scalar value for each level of the grouping factor. A simple, scalar random effects term is of the form, (1|G). It corresponds to a shift in the intercept for each level of the grouping factor."
},

{
    "location": "constructors.html#Models-with-vector-valued-random-effects-1",
    "page": "Model constructors",
    "title": "Models with vector-valued random effects",
    "category": "section",
    "text": "The sleepstudy data are observations of reaction time, Y, on several subjects, G, after 0 to 9 days of sleep deprivation, U. A model with random intercepts and random slopes for each subject, allowing for within-subject correlation of the slope and intercept, is fit asjulia> fm2 = fit!(lmm(@formula(Y ~ 1 + U + (1+U|G)), dat[:sleepstudy]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + U + ((1 + U) | G)\n   logLik   -2 logLik     AIC        BIC    \n -875.96967 1751.93934 1763.93934 1783.09709\n\nVariance components:\n              Column    Variance  Std.Dev.   Corr.\n G        (Intercept)  565.51067 23.780468\n          U             32.68212  5.716828  0.08\n Residual              654.94145 25.591824\n Number of obs: 180; levels of grouping factors: 18\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)   251.405   6.63226 37.9064  <1e-99\nU             10.4673   1.50224 6.96781  <1e-11\n\nA model with uncorrelated random effects for the intercept and slope by subject is fit asjulia> fm3 = fit!(lmm(@formula(Y ~ 1 + U + (1|G) + (0+U|G)), dat[:sleepstudy]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + U + (1 | G) + ((0 + U) | G)\n   logLik   -2 logLik     AIC        BIC    \n -876.00163 1752.00326 1762.00326 1777.96804\n\nVariance components:\n              Column    Variance  Std.Dev.   Corr.\n G        (Intercept)  584.258968 24.17145\n          U             33.632805  5.79938  0.00\n Residual              653.115782 25.55613\n Number of obs: 180; levels of grouping factors: 18\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)   251.405   6.70771   37.48  <1e-99\nU             10.4673   1.51931 6.88951  <1e-11\n\nAlthough technically there are two random-effects terms in the formula for fm3 both have the same grouping factor and, internally, are amalgamated into a single vector-valued term."
},

{
    "location": "constructors.html#Models-with-multiple,-scalar-random-effects-terms-1",
    "page": "Model constructors",
    "title": "Models with multiple, scalar random-effects terms",
    "category": "section",
    "text": "A model for the Penicillin data incorporates random effects for the plate, G, and for the sample, H. As every sample is used on every plate these two factors are crossed.julia> fm4 = fit!(lmm(@formula(Y ~ 1 + (1|G) + (1|H)), dat[:Penicillin]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G) + (1 | H)\n   logLik   -2 logLik     AIC        BIC    \n -166.09417  332.18835  340.18835  352.06760\n\nVariance components:\n              Column    Variance  Std.Dev. \n G        (Intercept)  0.7149795 0.8455646\n H        (Intercept)  3.1351920 1.7706474\n Residual              0.3024264 0.5499331\n Number of obs: 144; levels of grouping factors: 24, 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)   22.9722  0.744596 30.8519  <1e-99\n\nIn contrast the sample, G, grouping factor is nested within the batch, H, grouping factor in the Pastes data. That is, each level of G occurs in conjunction with only one level of H.julia> fm5 = fit!(lmm(@formula(Y ~ 1 + (1|G) + (1|H)), dat[:Pastes]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G) + (1 | H)\n   logLik   -2 logLik     AIC        BIC    \n -123.99723  247.99447  255.99447  264.37184\n\nVariance components:\n              Column    Variance  Std.Dev.  \n G        (Intercept)  8.4336166 2.90406897\n H        (Intercept)  1.1991794 1.09507048\n Residual              0.6780021 0.82340884\n Number of obs: 60; levels of grouping factors: 30, 10\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)   60.0533  0.642136 93.5212  <1e-99\n\nIn observational studies it is common to encounter partially crossed grouping factors. For example, the InstEval data are course evaluations by students, G, of instructors, H. Additional covariates include the academic department, H, in which the course was given and A, whether or not it was a service course.julia> fm6 = fit!(lmm(@formula(Y ~ 1 + A * I + (1|G) + (1|H)), dat[:InstEval]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + A * I + (1 | G) + (1 | H)\n     logLik        -2 logLik          AIC             BIC       \n -1.18792777×10⁵  2.37585553×10⁵  2.37647553×10⁵  2.37932876×10⁵\n\nVariance components:\n              Column     Variance   Std.Dev.  \n G        (Intercept)  0.105417976 0.32468135\n H        (Intercept)  0.258416368 0.50834670\n Residual              1.384727771 1.17674457\n Number of obs: 73421; levels of grouping factors: 2972, 1128\n\n  Fixed-effects parameters:\n                Estimate Std.Error   z value P(>|z|)\n(Intercept)      3.22961  0.064053   50.4209  <1e-99\nA: 1            0.252025 0.0686507   3.67112  0.0002\nI: 5            0.129536  0.101294   1.27882  0.2010\nI: 10          -0.176751 0.0881352  -2.00545  0.0449\nI: 12          0.0517102 0.0817524  0.632522  0.5270\nI: 6           0.0347319  0.085621  0.405647  0.6850\nI: 7             0.14594 0.0997984   1.46235  0.1436\nI: 4            0.151689 0.0816897   1.85689  0.0633\nI: 8            0.104206  0.118751  0.877517  0.3802\nI: 9           0.0440401 0.0962985  0.457329  0.6474\nI: 14          0.0517546 0.0986029  0.524879  0.5997\nI: 1           0.0466719  0.101942  0.457828  0.6471\nI: 3           0.0563461 0.0977925   0.57618  0.5645\nI: 11          0.0596536  0.100233   0.59515  0.5517\nI: 2          0.00556281  0.110867 0.0501756  0.9600\nA: 1 & I: 5    -0.180757  0.123179  -1.46744  0.1423\nA: 1 & I: 10   0.0186492  0.110017  0.169513  0.8654\nA: 1 & I: 12   -0.282269 0.0792937  -3.55979  0.0004\nA: 1 & I: 6    -0.494464 0.0790278  -6.25683   <1e-9\nA: 1 & I: 7    -0.392054  0.110313  -3.55403  0.0004\nA: 1 & I: 4    -0.278547 0.0823727  -3.38154  0.0007\nA: 1 & I: 8    -0.189526  0.111449  -1.70056  0.0890\nA: 1 & I: 9    -0.499868 0.0885423  -5.64553   <1e-7\nA: 1 & I: 14   -0.497162 0.0917162  -5.42065   <1e-7\nA: 1 & I: 1     -0.24042 0.0982071   -2.4481  0.0144\nA: 1 & I: 3    -0.223013 0.0890548  -2.50422  0.0123\nA: 1 & I: 11   -0.516997 0.0809077  -6.38997   <1e-9\nA: 1 & I: 2    -0.384773  0.091843  -4.18946   <1e-4\n\n"
},

{
    "location": "constructors.html#MixedModels.glmm",
    "page": "Model constructors",
    "title": "MixedModels.glmm",
    "category": "Function",
    "text": "glmm(f::Formula, fr::ModelFrame, d::Distribution[, l::GLM.Link])\n\nReturn a GeneralizedLinearMixedModel object.\n\nThe value is ready to be fit! but has not yet been fit.\n\n\n\n"
},

{
    "location": "constructors.html#Fitting-generalized-linear-mixed-models-1",
    "page": "Model constructors",
    "title": "Fitting generalized linear mixed models",
    "category": "section",
    "text": "To create a GLMM usingglmmthe distribution family for the response, given the random effects, must be specified.julia> gm1 = fit!(glmm(@formula(r2 ~ 1 + a + g + b + s + m + (1|id) + (1|item)), dat[:VerbAgg],\n    Bernoulli()))\nGeneralized Linear Mixed Model fit by minimizing the Laplace approximation to the deviance\n  Formula: r2 ~ 1 + a + g + b + s + m + (1 | id) + (1 | item)\n  Distribution: Distributions.Bernoulli{Float64}\n  Link: GLM.LogitLink()\n\n  Deviance (Laplace approximation): 8135.8329\n\nVariance components:\n          Column     Variance   Std.Dev. \n id   (Intercept)  1.793470989 1.3392054\n item (Intercept)  0.117151977 0.3422747\n\n Number of obs: 7584; levels of grouping factors: 316, 24\n\nFixed-effects parameters:\n              Estimate Std.Error  z value P(>|z|)\n(Intercept)   0.553345  0.385363  1.43591  0.1510\na            0.0574211 0.0167527  3.42757  0.0006\ng: M          0.320792  0.191206  1.67773  0.0934\nb: scold      -1.05975   0.18416 -5.75448   <1e-8\nb: shout       -2.1038  0.186519 -11.2793  <1e-28\ns: self       -1.05429  0.151196   -6.973  <1e-11\nm: do         -0.70698  0.151009 -4.68172   <1e-5\n\nThe canonical link, which is the GLM.LogitLink for the Bernoulli distribution, is used if no explicit link is specified.In the GLM package the appropriate distribution for a 0/1 response is the Bernoulli distribution. The Binomial distribution is only used when the response is the fraction of trials returning a positive, in which case the number of trials must be specified as the case weights."
},

{
    "location": "extractors.html#",
    "page": "Extractor functions",
    "title": "Extractor functions",
    "category": "page",
    "text": ""
},

{
    "location": "extractors.html#StatsBase.coef",
    "page": "Extractor functions",
    "title": "StatsBase.coef",
    "category": "Function",
    "text": "coef(obj::StatisticalModel)\n\nReturn the coefficients of the model.\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.coeftable",
    "page": "Extractor functions",
    "title": "StatsBase.coeftable",
    "category": "Function",
    "text": "coeftable(obj::StatisticalModel)\n\nReturn a table of class CoefTable with coefficients and related statistics.\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.dof",
    "page": "Extractor functions",
    "title": "StatsBase.dof",
    "category": "Function",
    "text": "dof(obj::StatisticalModel)\n\nReturn the number of degrees of freedom consumed in the model, including when applicable the intercept and the distribution's dispersion parameter.\n\n\n\ndof(d::UnivariateDistribution)\n\nGet the degrees of freedom.\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.deviance",
    "page": "Extractor functions",
    "title": "StatsBase.deviance",
    "category": "Function",
    "text": "deviance(obj::StatisticalModel)\n\nReturn the deviance of the model relative to a reference, which is usually when applicable the saturated model. It is equal, up to a constant, to -2 log L, with L the likelihood of the model.\n\n\n\ndeviance(obj::LinearModel)\n\nFor linear models, the deviance is equal to the residual sum of squares (RSS).\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.fitted",
    "page": "Extractor functions",
    "title": "StatsBase.fitted",
    "category": "Function",
    "text": "fitted(obj::RegressionModel)\n\nReturn the fitted values of the model.\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.loglikelihood",
    "page": "Extractor functions",
    "title": "StatsBase.loglikelihood",
    "category": "Function",
    "text": "loglikelihood(obj::StatisticalModel)\n\nReturn the log-likelihood of the model.\n\n\n\nloglikelihood(d::UnivariateDistribution, X::AbstractArray)\n\nThe log-likelihood of distribution d w.r.t. all samples contained in array x.\n\n\n\nloglikelihood(d::MultivariateDistribution, x::AbstractMatrix)\n\nThe log-likelihood of distribution d w.r.t. all columns contained in matrix x.\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.stderr",
    "page": "Extractor functions",
    "title": "StatsBase.stderr",
    "category": "Function",
    "text": "stderr(obj::StatisticalModel)\n\nReturn the standard errors for the coefficients of the model.\n\n\n\n"
},

{
    "location": "extractors.html#StatsBase.vcov",
    "page": "Extractor functions",
    "title": "StatsBase.vcov",
    "category": "Function",
    "text": "vcov(obj::StatisticalModel)\n\nReturn the variance-covariance matrix for the coefficients of the model.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.fixef",
    "page": "Extractor functions",
    "title": "MixedModels.fixef",
    "category": "Function",
    "text": "fixef(m::MixedModel)\n\nReturns the fixed-effects parameter vector estimate.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.fnames",
    "page": "Extractor functions",
    "title": "MixedModels.fnames",
    "category": "Function",
    "text": "fnames(m::MixedModel)\n\nReturn the names of the grouping factors for the random-effects terms.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.getΛ",
    "page": "Extractor functions",
    "title": "MixedModels.getΛ",
    "category": "Function",
    "text": "getΛ(m::MixedModel)\n\nReturn a vector of covariance template matrices for the random effects of m\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.getθ",
    "page": "Extractor functions",
    "title": "MixedModels.getθ",
    "category": "Function",
    "text": "getθ(A::FactorReTerm)\n\nReturn a vector of the elements of the lower triangle blocks in A.Λ (column-major ordering)\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.lowerbd",
    "page": "Extractor functions",
    "title": "MixedModels.lowerbd",
    "category": "Function",
    "text": "lowerbd{T}(A::FactorReTerm{T})\nlowerbd{T}(A::MatrixTerm{T})\nlowerbd{T}(v::Vector{AbstractTerm{T}})\n\nReturn the vector of lower bounds on the parameters, θ.\n\nThese are the elements in the lower triangle in column-major ordering. Diagonals have a lower bound of 0.  Off-diagonals have a lower-bound of -Inf.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.objective",
    "page": "Extractor functions",
    "title": "MixedModels.objective",
    "category": "Function",
    "text": "objective(m::LinearMixedModel)\n\nReturn negative twice the log-likelihood of model m\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.pwrss",
    "page": "Extractor functions",
    "title": "MixedModels.pwrss",
    "category": "Function",
    "text": "pwrss(m::LinearMixedModel)\n\nThe penalized residual sum-of-squares.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.ranef",
    "page": "Extractor functions",
    "title": "MixedModels.ranef",
    "category": "Function",
    "text": "ranef(m::MixedModel; uscale=false, named=true)\n\nReturn, as a Vector{Vector{T}} (Vector{NamedVector{T}} if named=true) the conditional modes of the random effects in model m.\n\nIf uscale is true the random effects are on the spherical (i.e. u) scale, otherwise on the original scale.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.sdest",
    "page": "Extractor functions",
    "title": "MixedModels.sdest",
    "category": "Function",
    "text": "sdest(m::LinearMixedModel)\n\nReturn the estimate of σ, the standard deviation of the per-observation noise.\n\n\n\n"
},

{
    "location": "extractors.html#MixedModels.varest",
    "page": "Extractor functions",
    "title": "MixedModels.varest",
    "category": "Function",
    "text": "varest(m::LinearMixedModel)\n\nReturns the estimate of σ², the variance of the conditional distribution of Y given B.\n\n\n\n"
},

{
    "location": "extractors.html#Extractor-functions-1",
    "page": "Extractor functions",
    "title": "Extractor functions",
    "category": "section",
    "text": "LinearMixedModel and GeneralizedLinearMixedModel are subtypes of StatsBase.RegressionModel Many of the generic extractors defined in the StatsBase package have methods for these models.StatsBase.coef\nStatsBase.coeftable\nStatsBase.dof\nStatsBase.deviance\nStatsBase.fitted\nStatsBase.loglikelihood\nStatsBase.stderr\nStatsBase.vcovOther extractors are defined in the MixedModels package itself.fixef\nfnames\ngetΛ\ngetθ\nlowerbd\nobjective\npwrss\nranef\nsdest\nvarestApplied to one of the models previously fit these yieldjulia> using DataFrames, RData, MixedModels\n\njulia> const dat = convert(Dict{Symbol,DataFrame}, load(Pkg.dir(\"MixedModels\", \"test\", \"dat.rda\")));\n\njulia> fm1 = fit!(lmm(@formula(Y ~ 1 + (1|G)), dat[:Dyestuff]))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G)\n   logLik   -2 logLik     AIC        BIC    \n -163.66353  327.32706  333.32706  337.53065\n\nVariance components:\n              Column    Variance  Std.Dev. \n G        (Intercept)  1388.3333 37.260345\n Residual              2451.2500 49.510100\n Number of obs: 30; levels of grouping factors: 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    1527.5   17.6946  86.326  <1e-99\n\n\njulia> fixef(fm1)\n1-element Array{Float64,1}:\n 1527.5\n\njulia> coef(fm1)\n1-element Array{Float64,1}:\n 1527.5\n\njulia> coeftable(fm1)\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    1527.5   17.6946  86.326  <1e-99\n\n\njulia> getΛ(fm1)\n1-element Array{Float64,1}:\n 0.752581\n\njulia> getθ(fm1)\n1-element Array{Float64,1}:\n 0.752581\n\njulia> loglikelihood(fm1)\n-163.6635299405672\n\njulia> pwrss(fm1)\n73537.50049200655\n\njulia> showall(ranef(fm1))\nArray{Float64,2}[[-16.6282 0.369516 26.9747 -21.8014 53.5798 -42.4943]]\njulia> showall(ranef(fm1, uscale=true))\nArray{Float64,2}[[-22.0949 0.490999 35.8429 -28.9689 71.1948 -56.4648]]\njulia> sdest(fm1)\n49.51010014532609\n\njulia> std(fm1)\n2-element Array{Array{Float64,1},1}:\n [37.2603]\n [49.5101]\n\njulia> stderr(fm1)\n1-element Array{Float64,1}:\n 17.6946\n\njulia> varest(fm1)\n2451.2500164002186\n\njulia> vcov(fm1)\n1×1 Array{Float64,2}:\n 313.097\n"
},

{
    "location": "bootstrap.html#",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "Parametric bootstrap for linear mixed-effects models",
    "category": "page",
    "text": ""
},

{
    "location": "bootstrap.html#Parametric-bootstrap-for-linear-mixed-effects-models-1",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "Parametric bootstrap for linear mixed-effects models",
    "category": "section",
    "text": "Julia is well-suited to implementing bootstrapping and other simulation-based methods for statistical models. The bootstrap! function in the MixedModels package provides an efficient parametric bootstrap for linear mixed-effects models, assuming that the results of interest from each simulated response vector can be incorporated into a vector of floating-point values."
},

{
    "location": "bootstrap.html#The-parametric-bootstrap-1",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "The parametric bootstrap",
    "category": "section",
    "text": "Bootstrapping is a family of procedures for generating sample values of a statistic, allowing for visualization of the distribution of the statistic or for inference from this sample of values.A _parametric bootstrap_ is used with a parametric model, m, that has been fitted to data. The procedure is to simulate n response vectors from m using the estimated parameter values and refit m to these responses in turn, accumulating the statistics of interest at each iteration.The parameters of a linear mixed-effects model as fit by the lmm function are the fixed-effects parameters, β, the standard deviation, σ, of the per-observation noise, and the covariance parameter, θ, that defines the variance-covariance matrices of the random effects.For example, a simple linear mixed-effects model for the Dyestuff data in the lme4 package for R is fit byjulia> using DataFrames, Gadfly, MixedModels, RData\njulia> ds = names!(dat[:Dyestuff], [:Batch, :Yield])\n30×2 DataFrames.DataFrame\n│ Row │ Batch │ Yield  │\n├─────┼───────┼────────┤\n│ 1   │ \"A\"   │ 1545.0 │\n│ 2   │ \"A\"   │ 1440.0 │\n│ 3   │ \"A\"   │ 1440.0 │\n│ 4   │ \"A\"   │ 1520.0 │\n│ 5   │ \"A\"   │ 1580.0 │\n│ 6   │ \"B\"   │ 1540.0 │\n│ 7   │ \"B\"   │ 1555.0 │\n│ 8   │ \"B\"   │ 1490.0 │\n⋮\n│ 22  │ \"E\"   │ 1630.0 │\n│ 23  │ \"E\"   │ 1515.0 │\n│ 24  │ \"E\"   │ 1635.0 │\n│ 25  │ \"E\"   │ 1625.0 │\n│ 26  │ \"F\"   │ 1520.0 │\n│ 27  │ \"F\"   │ 1455.0 │\n│ 28  │ \"F\"   │ 1450.0 │\n│ 29  │ \"F\"   │ 1480.0 │\n│ 30  │ \"F\"   │ 1445.0 │\n\njulia> m1 = fit!(lmm(@formula(Yield ~ 1 + (1 | Batch)), ds))\nLinear mixed model fit by maximum likelihood\n Formula: Yield ~ 1 + (1 | Batch)\n   logLik   -2 logLik     AIC        BIC    \n -163.66353  327.32706  333.32706  337.53065\n\nVariance components:\n              Column    Variance  Std.Dev. \n Batch    (Intercept)  1388.3333 37.260345\n Residual              2451.2500 49.510100\n Number of obs: 30; levels of grouping factors: 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    1527.5   17.6946  86.326  <1e-99\n\n"
},

{
    "location": "bootstrap.html#Using-the-bootstrap!-function-1",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "Using the bootstrap! function",
    "category": "section",
    "text": "This quick explanation is provided for those who only wish to use the bootstrap! method and do not need detailed explanations of how it works. The three arguments to bootstrap! are the matrix that will be overwritten with the results, the model to bootstrap, and a function that overwrites a vector with the results of interest from the model.Suppose the objective is to obtain 100,000 parametric bootstrap samples of the estimates of the \"variance components\", σ² and σ₁², in this model.  In many implementations of mixed-effects models the estimate of σ₁², the variance of the scalar random effects, is reported along with a standard error, as if the estimator could be assumed to have a Gaussian distribution. Is this a reasonable assumption?A suitable function to save the results isjulia> function saveresults!(v, m)\n    v[1] = varest(m)\n    v[2] = abs2(getθ(m)[1]) * v[1]\nend\nsaveresults! (generic function with 1 method)\nThe varest extractor function returns the estimate of σ².  As seen above, the estimate of the σ₁ is the product of Θ and the estimate of σ.  The expression abs2(getΘ(m)[1]) evaluates to Θ². The [1] is necessary because the value returned by getθ is a vector and a scalar is needed here.As with any simulation-based method, it is advisable to set the random number seed before calling bootstrap! for reproducibility.julia> srand(1234321);\njulia> results = bootstrap!(zeros(2, 100000), m1, saveresults!);\nThe results for each bootstrap replication are stored in the columns of the matrix passed in as the first argument.  A density plot of the first row using the Gadfly package is created as\nplot(x = view(results, 1, :), Geom.density(), Guide.xlabel(\"Parametric bootstrap estimates of σ²\"))(Image: Density of parametric bootstrap estimates of σ² from model m1)(Image: Density of parametric bootstrap estimates of σ₁² from model m1)The distribution of the bootstrap samples of σ² is a bit skewed but not terribly so.  However, the distribution of the bootstrap samples of the estimate of σ₁² is highly skewed and has a spike at zero."
},

]}
