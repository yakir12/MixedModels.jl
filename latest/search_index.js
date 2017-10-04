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
    "text": "CurrentModule = MixedModelsMixedModels.jl is a Julia package providing capabilities for fitting and examining linear and generalized linear mixed-effect models. It is similar in scope to the lme4 package for R.Pages = [\"SimpleLMM.md\",\n         \"constructors.md\",\n         \"extractors.md\",\n         \"bootstrap.md\"]\nDepth = 2"
},

{
    "location": "SimpleLMM.html#",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "A Simple, Linear, Mixed-effects Model",
    "category": "page",
    "text": ""
},

{
    "location": "SimpleLMM.html#A-Simple,-Linear,-Mixed-effects-Model-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "A Simple, Linear, Mixed-effects Model",
    "category": "section",
    "text": "In this book we describe the theory behind a type of statistical model called mixed-effects models and the practice of fitting and analyzing such models using the MixedModels package for Julia. These models are used in many different disciplines. Because the descriptions of the models can vary markedly between disciplines, we begin by describing what mixed-effects models are and by exploring a very simple example of one type of mixed model, the linear mixed model.This simple example allows us to illustrate the use of the lmm function in the MixedModels package for fitting such models and other functions  for analyzing the fitted model. We also describe methods of assessing the precision of the parameter estimates and of visualizing the conditional distribution of the random effects, given the observed data."
},

{
    "location": "SimpleLMM.html#Mixed-effects-Models-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Mixed-effects Models",
    "category": "section",
    "text": "Mixed-effects models, like many other types of statistical models, describe a relationship between a response variable and some of the covariates that have been measured or observed along with the response. In mixed-effects models at least one of the covariates is a categorical covariate representing experimental or observational “units” in the data set. In the example from the chemical industry discussed in this chapter, the observational unit is the batch of an intermediate product used in production of a dye. In medical and social sciences the observational units are often the human or animal subjects in the study. In agriculture the experimental units may be the plots of land or the specific plants being studied.In all of these cases the categorical covariate or covariates are observed at a set of discrete levels. We may use numbers, such as subject identifiers, to designate the particular levels that we observed but these numbers are simply labels. The important characteristic of a categorical covariate is that, at each observed value of the response, the covariate takes on the value of one of a set of distinct levels.Parameters associated with the particular levels of a covariate are sometimes called the “effects” of the levels. If the set of possible levels of the covariate is fixed and reproducible we model the covariate using fixed-effects parameters. If the levels that were observed represent a random sample from the set of all possible levels we incorporate random effects in the model.There are two things to notice about this distinction between fixed-effects parameters and random effects. First, the names are misleading because the distinction between fixed and random is more a property of the levels of the categorical covariate than a property of the effects associated with them. Secondly, we distinguish between “fixed-effects parameters”, which are indeed parameters in the statistical model, and “random effects”, which, strictly speaking, are not parameters. As we will see shortly, random effects are unobserved random variables.To make the distinction more concrete, suppose the objective is to model the annual reading test scores for students in a school district and that the covariates recorded with the score include a student identifier and the student’s gender. Both of these are categorical covariates. The levels of the gender covariate, male and female, are fixed. If we consider data from another school district or we incorporate scores from earlier tests, we will not change those levels. On the other hand, the students whose scores we observed would generally be regarded as a sample from the set of all possible students whom we could have observed. Adding more data, either from more school districts or from results on previous or subsequent tests, will increase the number of distinct levels of the student identifier.Mixed-effects models or, more simply, mixed models are statistical models that incorporate both fixed-effects parameters and random effects. Because of the way that we will define random effects, a model with random effects always includes at least one fixed-effects parameter. Thus, any model with random effects is a mixed model.We characterize the statistical model in terms of two random variables: a q-dimensional vector of random effects represented by the random variable mathcalB and an n-dimensional response vector represented by the random variable mathcalY. (We use upper-case “script” characters to denote random variables. The corresponding lower-case upright letter denotes a particular value of the random variable.) We observe the value, bfy, of mathcalY. We do not observe the value, bfb, of mathcalB.When formulating the model we describe the unconditional distribution of mathcalB and the conditional distribution, (mathcalYmathcalB=bfb). The descriptions of the distributions involve the form of the distribution and the values of certain parameters. We use the observed values of the response and the covariates to estimate these parameters and to make inferences about them.That’s the big picture. Now let’s make this more concrete by describing a particular, versatile class of mixed models called linear mixed models and by studying a simple example of such a model. First we describe the data in the example."
},

{
    "location": "SimpleLMM.html#The-Dyestuff-and-Dyestuff2-Data-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "The Dyestuff and Dyestuff2 Data",
    "category": "section",
    "text": "Models with random effects have been in use for a long time. The first edition of the classic book, Statistical Methods in Research and Production, edited by O.L. Davies, was published in 1947 and contained examples of the use of random effects to characterize batch-to-batch variability in chemical processes. The data from one of these examples are available as Dyestuff in the MixedModels package. In this section we describe and plot these data and introduce a second example, the Dyestuff2 data, described in Box and Tiao (1973)."
},

{
    "location": "SimpleLMM.html#The-Dyestuff-Data-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "The Dyestuff Data",
    "category": "section",
    "text": "The data are described in Davies (), the fourth edition of the book mentioned above, as coming froman investigation to find out how much the variation from batch to batch in the quality of an intermediate product (H-acid) contributes to the variation in the yield of the dyestuff (Naphthalene Black 12B) made from it. In the experiment six samples of the intermediate, representing different batches of works manufacture, were obtained, and five preparations of the dyestuff were made in the laboratory from each sample. The equivalent yield of each preparation as grams of standard colour was determined by dye-trial.First attach the packages to be usedjulia> using DataFrames, Distributions, Gadfly, GLM, MixedModels, RData\n\njulia> using Gadfly.Geom: point, line, histogram, density, vline\n\njulia> using Gadfly.Guide: xlabel, ylabel, yticks\nThe Dyestuff data are available in the lme4 package for R. This data frame and others have been stored in saved RData format in the test directory within the MixedModels package.Access the Dyestuff datajulia> const dat = convert(Dict{Symbol,DataFrame}, load(Pkg.dir(\"MixedModels\", \"test\", \"dat.rda\")));\n\njulia> dyestuff = dat[:Dyestuff];\nand plot itplot(dyestuff, x = \"Y\", y = \"G\", point, xlabel(\"Yield of dyestuff (g)\"), ylabel(\"Batch\"))(Image: Yield versus Batch for the Dyestuff data)In the dotplot we can see that there is considerable variability in yield, even for preparations from the same batch, but there is also noticeable batch-to-batch variability. For example, four of the five preparations from batch F provided lower yields than did any of the preparations from batches C and E.Recall that the labels for the batches are just labels and that their ordering is arbitrary.  In a plot, however, the order of the levels influences the perception of the pattern.  Rather than providing an arbitrary pattern it is best to order the levels according to some criterion for the plot.  In this case a good choice is to order the batches by increasing mean yield, which can be easily done in R.(Note: at present this plot fails because of the ongoing DataFrames conversion.)julia> #dyestuff = rcopy(\"within(Dyestuff, Batch <- reorder(Batch, Yield, mean))\");\n#plot(dyestuff, x=\"Y\", y=\"G\", point, xlabel(\"Yield of dyestuff (g)\"))\nIn Sect. [sec:DyestuffLMM] we will use mixed models to quantify the variability in yield between batches. For the time being let us just note that the particular batches used in this experiment are a selection or sample from the set of all batches that we wish to consider. Furthermore, the extent to which one particular batch tends to increase or decrease the mean yield of the process — in other words, the “effect” of that particular batch on the yield — is not as interesting to us as is the extent of the variability between batches. For the purposes of designing, monitoring and controlling a process we want to predict the yield from future batches, taking into account the batch-to-batch variability and the within-batch variability. Being able to estimate the extent to which a particular batch in the past increased or decreased the yield is not usually an important goal for us. We will model the effects of the batches as random effects rather than as fixed-effects parameters."
},

{
    "location": "SimpleLMM.html#The-Dyestuff2-Data-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "The Dyestuff2 Data",
    "category": "section",
    "text": "The data are simulated data presented in Box and Tiao (1973), where the authors stateThese data had to be constructed for although examples of this sort undoubtedly occur in practice they seem to be rarely published.The structure and summary are intentionally similar to those of the Dyestuff data. As can be seen in Fig. [fig:Dyestuff2dot]julia> dyestuff2 = dat[:Dyestuff2];\n(Image: )the batch-to-batch variability in these data is small compared to the within-batch variability. In some approaches to mixed models it can be difficult to fit models to such data. Paradoxically, small “variance components” can be more difficult to estimate than large variance components.The methods we will present are not compromised when estimating small variance components."
},

{
    "location": "SimpleLMM.html#Fitting-Linear-Mixed-Models-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Fitting Linear Mixed Models",
    "category": "section",
    "text": "Before we formally define a linear mixed model, let’s go ahead and fit models to these data sets using lmm which takes, as its first two arguments, a formula specifying the model and the data with which to evaluate the formula. The structure of the formula will be explained after showing the example."
},

{
    "location": "SimpleLMM.html#A-Model-For-the-Dyestuff-Data-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "A Model For the Dyestuff Data",
    "category": "section",
    "text": "A model allowing for an overall level of the Yield and for an additive random effect for each level of Batch can be fit asjulia> mm1 = fit!(lmm(@formula(Y ~ 1 + (1 | G)), dyestuff))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G)\n   logLik   -2 logLik     AIC        BIC    \n -163.66353  327.32706  333.32706  337.53065\n\nVariance components:\n              Column    Variance  Std.Dev. \n G        (Intercept)  1388.3333 37.260345\n Residual              2451.2500 49.510100\n Number of obs: 30; levels of grouping factors: 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    1527.5   17.6946  86.326  <1e-99\n\nAs shown in the summary of the model fit, the default estimation criterion is maximum likelihood.  The summary provides several other model-fit statistics such as Akaike’s Information Criterion (AIC), Schwarz’s Bayesian Information Criterion (BIC), the log-likelihood at the parameter estimates, and the objective function (negative twice the log-likelihood) at the parameter estimates. These are all statistics related to the model fit and are used to compare different models fit to the same data.The third section is the table of estimates of parameters associated with the random effects. There are two sources of variability in this model, a batch-to-batch variability in the level of the response and the residual or per-observation variability — also called the within-batch variability. The name “residual” is used in statistical modeling to denote the part of the variability that cannot be explained or modeled with the other terms. It is the variation in the observed data that is “left over” after determining the estimates of the parameters in the other parts of the model.Some of the variability in the response is associated with the fixed-effects terms. In this model there is only one such term, labeled the (Intercept). The name “intercept”, which is better suited to models based on straight lines written in a slope/intercept form, should be understood to represent an overall “typical” or mean level of the response in this case. (For those wondering about the parentheses around the name, they are included so that a user cannot accidentally name a variable in conflict with this name.) The line labeled Batch in the random effects table shows that the random effects added to the intercept term, one for each level of the factor, are modeled as random variables whose unconditional variance is estimated as 1388.33 g^2. The corresponding standard deviations is 37.26 g for the ML fit.Note that the last column in the random effects summary table is the estimate of the variability expressed as a standard deviation rather than as a variance. These are provided because it is usually easier to visualize the variability in standard deviations, which are on the scale of the response, than it is to visualize the magnitude of a variance. The values in this column are a simple re-expression (the square root) of the estimated variances. Do not confuse them with the standard errors of the variance estimators, which are not given here. As described in later sections, standard errors of variance estimates are generally not useful because the distribution of the estimator of a variance is skewed - often badly skewed.The line labeled Residual in this table gives the estimate, 2451.25 g^2, of the variance of the residuals and the corresponding standard deviation, 49.51 g. In written descriptions of the model the residual variance parameter is written sigma^2 and the variance of the random effects is sigma_1^2.  Their estimates are widehatsigma^2 and widehatsigma_1^2The last line in the random effects table states the number of observations to which the model was fit and the number of levels of any “grouping factors” for the random effects. In this case we have a single random effects term, (1 | Batch), in the model formula and the grouping factor for that term is Batch. There will be a total of six random effects, one for each level of Batch.The final part of the printed display gives the estimates and standard errors of any fixed-effects parameters in the model. The only fixed-effects term in the model formula is the (Intercept). The estimate of this parameter is 1527.5 g. The standard error of the intercept estimate is 17.69 g."
},

{
    "location": "SimpleLMM.html#A-Model-For-the-Dyestuff2-Data-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "A Model For the Dyestuff2 Data",
    "category": "section",
    "text": "Fitting a similar model to the dyestuff2 data produces an estimate widehatsigma_1^2=0.julia> mm2 = fit!(lmm(@formula(Y ~ 1 + (1 | G)), dyestuff2))\nLinear mixed model fit by maximum likelihood\n Formula: Y ~ 1 + (1 | G)\n   logLik   -2 logLik     AIC        BIC    \n -81.436518 162.873037 168.873037 173.076629\n\nVariance components:\n              Column    Variance  Std.Dev. \n G        (Intercept)   0.000000 0.0000000\n Residual              13.346099 3.6532314\n Number of obs: 30; levels of grouping factors: 6\n\n  Fixed-effects parameters:\n             Estimate Std.Error z value P(>|z|)\n(Intercept)    5.6656  0.666986 8.49433  <1e-16\n\nAn estimate of 0 for sigma_1 does not mean that there is no variation between the groups. Indeed Fig. [fig:Dyestuff2dot] shows that there is some small amount of variability between the groups. The estimate, widehatsigma_1, is a measure of the “between-group” variability that is in excess of the variability induced by the \"within-group\" or residual variability in the responses.  If 30 observations were simulated from a \"normal\" (also called \"Gaussian\") distribution and divided arbitrarily into 6 groups of 5, a plot of the data would look very much like Fig. [fig:Dyestuff2dot].  (In fact, it is likely that this is how that data set was generated.) It is only where there is excess variability between the groups that widehatsigma_10.  Obtaining widehatsigma_1=0 is not a mistake; it is simply a statement about the data and the model.The important point to take away from this example is the need to allow for the estimates of variance components that are zero. Such a model is said to be singular, in the sense that it corresponds to a linear model in which we have removed the random effects associated with Batch. Singular models can and do occur in practice. Even when the final fitted model is not singular, we must allow for such models to be expressed when determining the parameter estimates through numerical optimization.It happens that this model corresponds to the linear model (i.e. a model with fixed-effects only)julia> lm1 = lm(@formula(Y ~ 1), dyestuff2)\nDataFrames.DataFrameRegressionModel{GLM.LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,Base.LinAlg.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}\n\nFormula: Y ~ +1\n\nCoefficients:\n             Estimate Std.Error t value Pr(>|t|)\n(Intercept)    5.6656  0.678388 8.35156    <1e-8\n\nThe log-likelihood for this modeljulia> loglikelihood(lm1)\n-81.43651832691287\nis the same as that of fm2. The standard error of the intercept in lm1 is a bit larger than that of fm2 because the estimate of the residual variance is evaluated differently in the linear model."
},

{
    "location": "SimpleLMM.html#Further-Assessment-of-the-Fitted-Models-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Further Assessment of the Fitted Models",
    "category": "section",
    "text": "The parameter estimates in a statistical model represent our “best guess” at the unknown values of the model parameters and, as such, are important results in statistical modeling. However, they are not the whole story. Statistical models characterize the variability in the data and we must assess the effect of this variability on the parameter estimates and on the precision of predictions made from the model.In Sect. [sec:variability] we introduce a method of assessing variability in parameter estimates using the “profiled log-likelihood” and in Sect. [sec:assessRE] we show methods of characterizing the conditional distribution of the random effects given the data. Before we get to these sections, however, we should state in some detail the probability model for linear mixed-effects and establish some definitions and notation. In particular, before we can discuss profiling the log-likelihood, we should define the log-likelihood. We do that in the next section."
},

{
    "location": "SimpleLMM.html#The-Linear-Mixed-effects-Probability-Model-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "The Linear Mixed-effects Probability Model",
    "category": "section",
    "text": "In explaining some of parameter estimates related to the random effects we have used terms such as “unconditional distribution” from the theory of probability. Before proceeding further we clarify the linear mixed-effects probability model and define several terms and concepts that will be used throughout the book. Readers who are more interested in practical results than in the statistical theory should feel free to skip this section."
},

{
    "location": "SimpleLMM.html#Definitions-and-Results-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Definitions and Results",
    "category": "section",
    "text": "In this section we provide some definitions and formulas without derivation and with minimal explanation, so that we can use these terms in what follows. In Chap. [chap:computational] we revisit these definitions providing derivations and more explanation.As mentioned in Sect. [sec:memod], a mixed model incorporates two random variables: mathcalB, the q-dimensional vector of random effects, and mathcalY, the n-dimensional response vector. In a linear mixed model the unconditional distribution of mathcalB and the conditional distribution, (mathcalY  mathcalB=bfb), are both multivariate Gaussian distributions,\\begin{aligned}   (\\mathcal{Y} | \\mathcal{B}=\\bf{b}) &\\sim\\mathcal{N}(\\bf{ X\\beta + Z b},\\sigma^2\\bf{I})\\\n  \\mathcal{B}&\\sim\\mathcal{N}(\\bf{0},\\Sigma_\\theta) . \\end{aligned}The conditional mean of mathcal Y, given mathcal B=bf b, is the linear predictor, bf Xbfbeta+bf Zbf b, which depends on the p-dimensional fixed-effects parameter, bf beta, and on bf b. The model matrices, bf X and bf Z, of dimension ntimes p and ntimes q, respectively, are determined from the formula for the model and the values of covariates. Although the matrix bf Z can be large (i.e. both n and q can be large), it is sparse (i.e. most of the elements in the matrix are zero).The relative covariance factor, Lambda_theta, is a qtimes q lower-triangular matrix, depending on the variance-component parameter, bftheta, and generating the symmetric qtimes q variance-covariance matrix, Sigma_theta, as\\begin{equation}   \\Sigma_\\theta=\\sigma^2\\Lambda_\\theta\\Lambda_\\theta' \\end{equation}The spherical random effects, mathcalUsimmathcalN(bf 0sigma^2bf I_q), determine mathcal B according to\\begin{equation}   \\mathcal{B}=\\Lambda_\\theta\\mathcal{U}$. \\end{equation}The penalized residual sum of squares (PRSS),\\begin{equation}   r^2(\\bf\\theta,\\bf\\beta,\\bf u)=\\|\\bf y -\\bf X\\bf\\beta -\\bf Z\\Lambda_\\theta\\bf u\\|^2+\\|\\bf u\\|^2, \\end{equation}is the sum of the residual sum of squares, measuring fidelity of the model to the data, and a penalty on the size of bf u, measuring the complexity of the model. Minimizing r^2 with respect to bf u,\\begin{equation}   r^2_{\\beta,\\theta} =\\min_{\\bf u}\\left{\\|{\\bf y} -{\\bf X}{\\bf\\beta} -{\\bf Z}\\Lambda_\\theta{\\bf u}\\|^2+\\|{\\bf u}\\|^2\\right} \\end{equation}is a direct (i.e. non-iterative) computation.  The particular method used to solve this generates a blocked Choleksy factor, bfR_theta, which is an upper triangular qtimes q matrix satisfying\\begin{equation}   \\bf R_\\theta'\\bf R_\\theta=\\Lambda_\\theta'\\bf Z'\\bf Z\\Lambda_\\theta+\\bf I_q . \\end{equation}where bf I_q is the qtimes q identity matrix.Negative twice the log-likelihood of the parameters, given the data, bf y, is\\begin{equation} d({\\bf\\theta},{\\bf\\beta},\\sigma|{\\bf y}) =n\\log(2\\pi\\sigma^2)+\\log(|{\\bf R}_\\theta|^2)+\\frac{r^2_{\\beta,\\theta}}{\\sigma^2}. \\end{equation}where bf R_theta denotes the determinant of bf R_theta. Because bf R_theta is triangular, its determinant is the product of its diagonal elements.Negative twice the log-likelihood will be called the objective in what follows.  It is the value to be minimized by the parameter estimates.  It is, up to an additive factor, the deviance of the parameters.  Unfortunately, it is not clear what the additive factor should be in the case of linear mixed models.  In many applications, however, it is not the deviance of the model that is of interest as much the change in the deviance between two fitted models.  When calculating the change in the deviance the additive factor will cancel out so the change in the deviance when comparing models is the same as the change in this objective.Because the conditional mean, bfmu_mathcal Ymathcal B=bf b=bf Xbfbeta+bf ZLambda_thetabf u, is a linear function of both bfbeta and bf u, minimization of the PRSS with respect to both bfbeta and bf u to produce\\begin{equation} r^2_\\theta =\\min_{{\\bf\\beta},{\\bf u}}\\left{\\|{\\bf y} -{\\bf X}{\\bf\\beta} -{\\bf Z}\\Lambda_\\theta{\\bf u}\\|^2+\\|{\\bf u}\\|^2\\right} \\end{equation}is also a direct calculation. The values of bf u and bfbeta that provide this minimum are called, respectively, the conditional mode, tildebf u_theta, of the spherical random effects and the conditional estimate, widehatbfbeta_theta, of the fixed effects. At the conditional estimate of the fixed effects the objective is\\begin{equation} d({\\bf\\theta},\\widehat{\\beta}_\\theta,\\sigma|{\\bf y}) =n\\log(2\\pi\\sigma^2)+\\log(|{\\bf L}_\\theta|^2)+\\frac{r^2_\\theta}{\\sigma^2}. \\end{equation}Minimizing this expression with respect to sigma^2 produces the conditional estimate\\begin{equation} \\widehat{\\sigma^2}_\\theta=\\frac{r^2_\\theta}{n} \\end{equation}which provides the profiled log-likelihood on the deviance scale as\\begin{equation} \\tilde{d}(\\bf\\theta|\\bf y)=d(\\bf\\theta,\\widehat{\\beta}_\\theta,\\widehat{\\sigma}_\\theta|\\bf y) =\\log(|\\bf R_\\theta|^2)+n\\left[1+\\log\\left(\\frac{2\\pi r^2_\\theta}{n}\\right)\\right], \\end{equation}a function of bftheta alone.The MLE of bftheta, written widehatbftheta, is the value that minimizes this profiled objective. We determine this value by numerical optimization. In the process of evaluating tilded(widehatbfthetabf y) we determine widehatbfbeta=widehatbfbeta_widehattheta, tildebf u_widehattheta and r^2_widehattheta, from which we can evaluate widehatsigma=sqrtr^2_widehatthetan.The elements of the conditional mode of mathcal B, evaluated at the parameter estimates,\\begin{equation} \\tilde{\\bf b}_{\\widehat{\\theta}}= \\Lambda_{\\widehat{\\theta}}\\tilde{\\bf u}_{\\widehat{\\theta}} \\end{equation}are sometimes called the best linear unbiased predictors or BLUPs of the random effects. Although BLUPs an appealing acronym, I don’t find the term particularly instructive (what is a “linear unbiased predictor” and in what sense are these the “best”?) and prefer the term “conditional modes”, because these are the values of bf b that maximize the density of the conditional distribution (mathcalB  mathcalY = bf y.  For a linear mixed model, where all the conditional and unconditional distributions are Gaussian, these values are also the conditional means."
},

{
    "location": "SimpleLMM.html#Fields-of-a-LinearMixedModel-object-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Fields of a LinearMixedModel object",
    "category": "section",
    "text": "The optional second argument, verbose, in a call to fit! of a LinearMixedModel object produces output showing the progress of the iterative optimization of tilded(bfthetabf y).julia> mm1 = fit!(lmm(@formula(Y ~ 1 + (1 | G)), dyestuff), true);\nf_1: 327.76702 [1.0]\nf_2: 331.03619 [1.75]\nf_3: 330.64583 [0.25]\nf_4: 327.69511 [0.97619]\nf_5: 327.56631 [0.928569]\nf_6: 327.3826 [0.833327]\nf_7: 327.35315 [0.807188]\nf_8: 327.34663 [0.799688]\nf_9: 327.341 [0.792188]\nf_10: 327.33253 [0.777188]\nf_11: 327.32733 [0.747188]\nf_12: 327.32862 [0.739688]\nf_13: 327.32706 [0.752777]\nf_14: 327.32707 [0.753527]\nf_15: 327.32706 [0.752584]\nf_16: 327.32706 [0.752509]\nf_17: 327.32706 [0.752591]\nf_18: 327.32706 [0.752581]\nThe algorithm converges after 18 function evaluations to a profiled deviance of 327.32706 at theta=0752581. In this model the parameter theta is of length 1, the single element being the ratio sigma_1sigma.Whether or not verbose output is requested, the optsum field of a LinearMixedModel contains information on the optimization.  The various tolerances or the optimizer name can be changed between creating a LinearMixedModel and calling fit! on it to exert finer control on the optimization process.julia> mm1.optsum\nInitial parameter vector: [1.0]\nInitial objective value:  327.76702162461663\n\nOptimizer (from NLopt):   LN_BOBYQA\nLower bounds:             [0.0]\nftol_rel:                 1.0e-12\nftol_abs:                 1.0e-8\nxtol_rel:                 0.0\nxtol_abs:                 [1.0e-10]\ninitial_step:             [0.75]\nmaxfeval:                 -1\n\nFunction evaluations:     18\nFinal parameter vector:   [0.752581]\nFinal objective value:    327.3270598811344\nReturn code:              FTOL_REACHED\n\nThe full list of fields in a LinearMixedModel object isjulia> fieldnames(LinearMixedModel)\n6-element Array{Symbol,1}:\n :formula\n :trms   \n :sqrtwts\n :A      \n :L      \n :optsum \nThe formula field is a copy of the model formulajulia> mm1.formula\nFormula: Y ~ 1 + (1 | G)\nThe mf field is a ModelFrame constructed from the formula and data arguments to lmm. It contains a DataFrame formed by reducing the data argument to only those columns needed to evaluate the formula and only those rows that have no missing data. It also contains information on the terms in the model formula and on any \"contrasts\" associated with categorical variables in the fixed-effects terms.The trms field is a vector of numerical objects representing the terms in the model, including the response vector. As the names imply, the sqrtwts and wttrms fields are for incorporating case weights. These fields are not often used when fitting linear mixed models but are vital to the process of fitting a generalized linear mixed model, described in Chapter [sec:glmm].  When used, sqrtwts is a diagonal matrix of size n. A size of 0 indicates weights are not used.julia> mm1.sqrtwts\n0-element Array{Float64,1}\nThe trms field is a vector of length ge 3.julia> length(mm1.trms)\n3\nThe last two elements are bf X, the ntimes p model matrix for the fixed-effects parameters, bfbeta, and bf y, the response vector stored as a matrix of size ntimes 1.  In mm1, bf X consists of a single column of 1'sjulia> mm1.trms[end - 1]\nMixedModels.MatrixTerm{Float64,Array{Float64,2}}([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0], [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0], String[\"(Intercept)\"])\njulia> mm1.trms[end]\nMixedModels.MatrixTerm{Float64,Array{Float64,2}}([1545.0; 1440.0; 1440.0; 1520.0; 1580.0; 1540.0; 1555.0; 1490.0; 1560.0; 1495.0; 1595.0; 1550.0; 1605.0; 1510.0; 1560.0; 1445.0; 1440.0; 1595.0; 1465.0; 1545.0; 1595.0; 1630.0; 1515.0; 1635.0; 1625.0; 1520.0; 1455.0; 1450.0; 1480.0; 1445.0], [1545.0; 1440.0; 1440.0; 1520.0; 1580.0; 1540.0; 1555.0; 1490.0; 1560.0; 1495.0; 1595.0; 1550.0; 1605.0; 1510.0; 1560.0; 1445.0; 1440.0; 1595.0; 1465.0; 1545.0; 1595.0; 1630.0; 1515.0; 1635.0; 1625.0; 1520.0; 1455.0; 1450.0; 1480.0; 1445.0], String[\"\"])\nThe elements of trms before the last two represent vertical sections of bf Z associated with the random effects terms in the model.  In mm1 there is only one random effects term, (1 | Batch), and bf Z has only one section, the one generated by this term, of type ScalarReMat.julia> mm1.trms[1]\nMixedModels.ScalarFactorReTerm{Float64,String,UInt8}(String[\"A\", \"A\", \"A\", \"A\", \"A\", \"B\", \"B\", \"B\", \"B\", \"B\", \"C\", \"C\", \"C\", \"C\", \"C\", \"D\", \"D\", \"D\", \"D\", \"D\", \"E\", \"E\", \"E\", \"E\", \"E\", \"F\", \"F\", \"F\", \"F\", \"F\"], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], :G, String[\"(Intercept)\"], 0.7525806752871207)\nIn practice these matrices are stored in a highly condensed form because, in some models, they can be very large. In small examples the structure is more obvious when the ScalarReMat is converted to a sparse or a dense matrix.julia> sparse(mm1.trms[1])\n30×6 SparseMatrixCSC{Float64,Int64} with 30 stored entries:\n  [1 ,  1]  =  1.0\n  [2 ,  1]  =  1.0\n  [3 ,  1]  =  1.0\n  [4 ,  1]  =  1.0\n  [5 ,  1]  =  1.0\n  [6 ,  2]  =  1.0\n  [7 ,  2]  =  1.0\n  [8 ,  2]  =  1.0\n  [9 ,  2]  =  1.0\n  [10,  2]  =  1.0\n  [11,  3]  =  1.0\n  [12,  3]  =  1.0\n  [13,  3]  =  1.0\n  [14,  3]  =  1.0\n  [15,  3]  =  1.0\n  [16,  4]  =  1.0\n  [17,  4]  =  1.0\n  [18,  4]  =  1.0\n  [19,  4]  =  1.0\n  [20,  4]  =  1.0\n  [21,  5]  =  1.0\n  [22,  5]  =  1.0\n  [23,  5]  =  1.0\n  [24,  5]  =  1.0\n  [25,  5]  =  1.0\n  [26,  6]  =  1.0\n  [27,  6]  =  1.0\n  [28,  6]  =  1.0\n  [29,  6]  =  1.0\n  [30,  6]  =  1.0\njulia> full(mm1.trms[1])\n30×6 Array{Float64,2}:\n 1.0  0.0  0.0  0.0  0.0  0.0\n 1.0  0.0  0.0  0.0  0.0  0.0\n 1.0  0.0  0.0  0.0  0.0  0.0\n 1.0  0.0  0.0  0.0  0.0  0.0\n 1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  0.0  1.0  0.0  0.0\n 0.0  0.0  0.0  1.0  0.0  0.0\n 0.0  0.0  0.0  1.0  0.0  0.0\n 0.0  0.0  0.0  1.0  0.0  0.0\n 0.0  0.0  0.0  1.0  0.0  0.0\n 0.0  0.0  0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  0.0  0.0  1.0\n 0.0  0.0  0.0  0.0  0.0  1.0\n 0.0  0.0  0.0  0.0  0.0  1.0\n 0.0  0.0  0.0  0.0  0.0  1.0\n 0.0  0.0  0.0  0.0  0.0  1.0\nThe A field is a representation of the blocked, square, symmetric matrix bf A = Z  X  yZ  X  y. Only the upper triangle of A is stored. The number of blocks of the rows and columns of A is the number of vertical sections of bf Z (i.e. the number of random-effects terms) plus 2.julia> nblocks(mm1.A)\n(3, 3)\njulia> mm1.A[Block(1, 1)]\n6×6 Diagonal{Float64}:\n 5.0   ⋅    ⋅    ⋅    ⋅    ⋅ \n  ⋅   5.0   ⋅    ⋅    ⋅    ⋅ \n  ⋅    ⋅   5.0   ⋅    ⋅    ⋅ \n  ⋅    ⋅    ⋅   5.0   ⋅    ⋅ \n  ⋅    ⋅    ⋅    ⋅   5.0   ⋅ \n  ⋅    ⋅    ⋅    ⋅    ⋅   5.0\njulia> mm1.A[Block(2, 1)]\n1×6 Array{Float64,2}:\n 5.0  5.0  5.0  5.0  5.0  5.0\njulia> mm1.A[Block(2, 2)]\n1×1 Array{Float64,2}:\n 30.0\njulia> mm1.A[Block(3, 1)]\n1×6 Array{Float64,2}:\n 7525.0  7640.0  7820.0  7490.0  8000.0  7350.0\njulia> mm1.A[Block(3, 2)]\n1×1 Array{Float64,2}:\n 45825.0\njulia> mm1.A[Block(3, 3)]\n1×1 Array{Float64,2}:\n 7.01129e7\n"
},

{
    "location": "SimpleLMM.html#Fields-modified-during-the-optimization-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Fields modified during the optimization",
    "category": "section",
    "text": "Changing the value of theta changes the Lambda field. (Note: to input a symbol like Lambda in a Jupyter code cell or in the Julia read-eval-print loop (REPL), type \\Lambda followed by a tab character.  Such \"latex completions\" are available for many UTF-8 characters used in Julia.)  The matrix Lambda has a special structure.  It is a block diagonal matrix where the diagonal blocks are Kronecker products of an identity matrix and a (small) lower triangular matrix.  The diagonal blocks correspond to the random-effects terms.  For a scalar random-effects term, like (1 | Batch) the diagonal block is the Kronecker product of an identity matrix and a 1times 1 matrix.  This result in this case is just a multiple of the identity matrix.It is not necessary to store the full Lambda matrix.  Storing the small lower-triangular matrices is sufficient.julia> getΛ(mm1)\n1-element Array{Float64,1}:\n 0.752581\nThe L field is a blocked matrix like the A field containing the upper Cholesky factor of\\begin{bmatrix}   \\bf{\\Lambda'Z'Z\\Lambda + I} & \\bf{\\Lambda'Z'X} & \\bf{\\Lambda'Z'y} \\\n  \\bf{X'Z\\Lambda} & \\bf{X'X} & \\bf{X'y} \\\n  \\bf{y'Z\\Lambda} & \\bf{y'Z} & \\bf{y'y} \\end{bmatrix}   julia> mm1.L\n8×8 LowerTriangular{Float64,BlockArrays.BlockArray{Float64,2,AbstractArray{Float64,2}}}:\n    1.95752      ⋅           ⋅           ⋅           ⋅           ⋅           ⋅          ⋅   \n    0.0         1.95752      ⋅           ⋅           ⋅           ⋅           ⋅          ⋅   \n    0.0         0.0         1.95752      ⋅           ⋅           ⋅           ⋅          ⋅   \n    0.0         0.0         0.0         1.95752      ⋅           ⋅           ⋅          ⋅   \n    0.0         0.0         0.0         0.0         1.95752      ⋅           ⋅          ⋅   \n    0.0         0.0         0.0         0.0         0.0         1.95752      ⋅          ⋅   \n    1.92228     1.92228     1.92228     1.92228     1.92228     1.92228     2.79804     ⋅   \n 2893.03     2937.24     3006.45     2879.58     3075.65     2825.75     4274.01     271.178\njulia> mm1.L.data[Block(1, 1)]\n6×6 Diagonal{Float64}:\n 1.95752   ⋅        ⋅        ⋅        ⋅        ⋅     \n  ⋅       1.95752   ⋅        ⋅        ⋅        ⋅     \n  ⋅        ⋅       1.95752   ⋅        ⋅        ⋅     \n  ⋅        ⋅        ⋅       1.95752   ⋅        ⋅     \n  ⋅        ⋅        ⋅        ⋅       1.95752   ⋅     \n  ⋅        ⋅        ⋅        ⋅        ⋅       1.95752\njulia> mm1.L.data[Block(2, 1)]\n1×6 Array{Float64,2}:\n 1.92228  1.92228  1.92228  1.92228  1.92228  1.92228\njulia> mm1.L.data[Block(2, 2)]\n1×1 Array{Float64,2}:\n 2.79804\njulia> mm1.L.data[Block(3, 1)]\n1×6 Array{Float64,2}:\n 2893.03  2937.24  3006.45  2879.58  3075.65  2825.75\njulia> mm1.L.data[Block(3, 2)]\n1×1 Array{Float64,2}:\n 4274.01\njulia> mm1.L.data[Block(3, 3)]\n1×1 Array{Float64,2}:\n 271.178\nAll the information needed to evaluate the profiled log-likelihood is available in the R field; log(bf R_theta^2) isjulia> 2 * sum(log.(diag(mm1.L.data[Block(1,1)])))\n8.060146362820694\nIt can also be evaluated as logdet(mm1) or 2 * logdet(mm1.R[1, 1])julia> logdet(mm1) == (2*logdet(mm1.L.data[Block(1, 1)])) == (2*sum(log.(diag(mm1.L.data[Block(1, 1)]))))\ntrue\nThe penalized residual sum of squares is the square of the single element of the lower-right block, R[3, 3] in this casejulia> abs2(mm1.L.data[Block(3, 3)][1, 1])\n73537.50049200655\njulia> pwrss(mm1)\n73537.50049200655\nThe objective isjulia> logdet(mm1) + nobs(mm1) * (1 + log(2π * pwrss(mm1) / nobs(mm1)))\n327.3270598811344\n"
},

{
    "location": "SimpleLMM.html#Assessing-variability-of-parameter-estimates-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Assessing variability of parameter estimates",
    "category": "section",
    "text": ""
},

{
    "location": "SimpleLMM.html#Parametric-bootstrap-samples-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Parametric bootstrap samples",
    "category": "section",
    "text": "One way to assess the variability of the parameter estimates is to generate a parametric bootstrap sample from the model.  The technique is to simulate response vectors from the model at the estimated parameter values and refit the model to each of these simulated responses, recording the values of the parameters.  The bootstrap method for these models performs these simulations and returns 4 arrays: a vector of objective (negative twice the log-likelihood) values, a vector of estimates of sigma^2, a matrix of estimates of the fixed-effects parameters and a matrix of the estimates of the relative covariance parameters.  In this case there is only one fixed-effects parameter and one relative covariance parameter, which is the ratio of the standard deviation of the random effects to the standard deviation of the per-sample noise.First set the random number seed for reproducibility.julia> srand(1234321);\n\njulia> mm1bstp = bootstrap(10000, mm1);\n\njulia> size(mm1bstp)\n(10000, 5)\njulia> show(names(mm1bstp))\nSymbol[:obj, :σ, :β₁, :θ₁, :σ₁]"
},

{
    "location": "SimpleLMM.html#Histograms,-kernel-density-plots-and-quantile-quantile-plots-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Histograms, kernel density plots and quantile-quantile plots",
    "category": "section",
    "text": "I am a firm believer in the value of plotting results before summarizing them. Well chosen plots can provide insights not available from a simple numerical summary. It is common to visualize the distribution of a sample using a histogram, which approximates the shape of the probability density function. The density can also be approximated more smoothly using a kernel density plot. Finally, the extent to which the distribution of a sample can be approximated by a particular distribution or distribution family can be assessed by a quantile-quantile (qq) plot, the most common of which is the normal probability plot.The Gadfly package for Julia uses a \"grammar of graphics\" specification, similar to the ggplot2 package for R.  A histogram or a kernel density plot are describes as geometries and specified by Geom.histogram and Geom.density, respectively.julia> plot(mm1bstp, x = :β₁, Geom.histogram)\nPlot(...)\n(Image: )julia> plot(mm1bstp, x = :σ, Geom.histogram)\nPlot(...)\n(Image: )julia> plot(mm1bstp, x = :σ₁, Geom.histogram)\nPlot(...)\n(Image: )The last two histograms show that, even if the models are defined in terms of variances, the variance is usually not a good scale on which to assess the variability of the parameter estimates.  The standard deviation or, in some cases, the logarithm of the standard deviation is a more suitable scale.The histogram of sigma_1^2 has a \"spike\" at zero.  Because the value of sigma^2 is never zero, a value of sigma_1^2=0 must correspond to theta=0.  There are several ways to count the zeros in theta1.  The simplest is to use countnz, which counts the non-zeros, and subtrack that value from the total number of values in theta1.julia> length(mm1bstp[:θ₁]) - countnz(mm1bstp[:θ₁])\n941\nThat is, nearly 1/10 of the theta1 values are zeros.  Because such a spike or pulse will be spread out or diffused in a kernel density plot,julia> plot(mm1bstp, x = :θ₁, Geom.density)\nPlot(...)\n(Image: )such a plot is not suitable for a sample of a bounded parameter that includes values on the boundary.The density of the estimates of the other two parameters, beta_1 and sigma, are depicted well in kernel density plots.julia> plot(mm1bstp, x = :β₁, Geom.density)\nPlot(...)\n(Image: )julia> plot(mm1bstp, x = :σ, Geom.density)\nPlot(...)\n(Image: )The standard approach of summarizing a sample by its mean and standard deviation, or of constructing a confidence interval using the sample mean, the standard error of the mean and quantiles of a t or normal distribution, are based on the assumption that the sample is approximately normal (also called Gaussian) in shape.  A normal probability plot, which plots sample quantiles versus quantiles of the standard normal distribution, mathcalN(01), can be used to assess the validity of this assumption.  If the points fall approximately along a straight line, the assumption of normality should be valid.  Systematic departures from a straight line are cause for concern.In Gadfly a normal probability plot can be  constructed by specifying the statistic to be generated as Stat.qq and either x or y as the distribution Normal(). For the present purposes it is an advantage to put the theoretical quantiles on the x axis.This approach is suitable for small to moderate sample sizes, but not for sample sizes of 10,000.  To smooth the plot and to reduce the size of the plot files, we plot quantiles defined by a sequence of n \"probability points\".  These are constructed by partitioning the interval (0, 1) into n equal-width subintervals and returning the midpoint of each of the subintervals.julia> function ppoints(n)\n    if n ≤ 0\n        throw(ArgumentError(\"n = $n should be positive\"))\n    end\n    width = inv(n)\n    width / 2 : width : one(width)\nend\nppoints (generic function with 1 method)\n\njulia> const ppt250 = ppoints(250)\n0.002:0.004:0.998\nThe kernel density estimate of sigma is more symmetric(Image: )and the normal probability plot of sigma is also reasonably straight.(Image: )The normal probability plot of sigma_1 has a flat section at sigma_1 = 0.(Image: )In terms of the variances, sigma^2 and sigma_1^2, the normal probability plots are(Image: )(Image: )"
},

{
    "location": "SimpleLMM.html#Confidence-intervals-based-on-bootstrap-samples-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Confidence intervals based on bootstrap samples",
    "category": "section",
    "text": "When the distribution of a parameter estimator is close to normal or to a T distribution, symmetric confidence intervals are an appropriate representation of the uncertainty in the parameter estimates.  However, they are not appropriate for skewed and/or bounded estimator distributions, such as those for sigma^2 and sigma_2^1 shown above.The fact that a symmetric confidence interval is not appropriate for sigma^2 should not be surprising.  In an introductory statistics course the calculation of a confidence interval on sigma^2 from a random sample of a normal distribution using quantiles of a chi^2 distribution is often introduced. So a symmetric confidence interval on sigma^2 is inappropriate in the simplest case but is often expected to be appropriate in much more complicated cases, as shown by the fact that many statistical software packages include standard errors of variance component estimates in the output from a mixed model fitting procedure.  Creating confidence intervals in this way is optimistic at best. Completely nonsensical would be another way of characterizing this approach.A more reasonable alternative for constructing a 1 - alpha confidence interval from a bootstrap sample is to report a contiguous interval that contains a 1 - alpha proportion of the sample values.But there are many such intervals.  Suppose that a 95% confidence interval was to be constructed from one of the samples of size 10,000 of bootstrapped values.  To get a contigous interval the sample should be sorted. The sorted sample values, also called the order statistics of the sample, are denoted by a bracketed subscript.  That is, sigma_1 is the smallest value in the sample, sigma_2 is the second smallest, up to sigma_10000, which is the largest.One possible interval containing 95% of the sample is (sigma_1 sigma_9500).  Another is (sigma_2 sigma_9501) and so on up to (sigma_501sigma_10000).  There needs to be a method of choosing one of these intervals.  On approach would be to always choose the central 95% of the sample.  That is, cut off 2.5% of the sample on the left side and 2.5% on the right side.  julia> sigma95 = quantile(mm1bstp[:σ], [0.025, 0.975])\n2-element Array{Float64,1}:\n 35.3694\n 62.9763\nThis approach has the advantage that the endpoints of a 95% interval on sigma^2 are the squares of the endpoints of a 95% interval on sigma.julia> isapprox(abs2.(sigma95), quantile(abs2.(mm1bstp[:σ]), [0.025, 0.975]))\ntrue\nThe intervals are compared with isapprox rather than exact equality because, in floating point arithmetic, it is not always the case that left(sqrtxright)^2 = x.  This comparison can also be expressed in Julia asjulia> abs2.(sigma95) ≈ quantile(abs2.(mm1bstp[:σ]), [0.025, 0.975])\ntrue\nAn alternative approach is to evaluate all of the contiguous intervals containing, say, 95% of the sample and return the shortest shortest such interval.  This is the equivalent of a Highest Posterior Density (HPD) interval sometimes used in Bayesian analysis.  If the procedure is applied to a unimodal (i.e. one that has only one peak or mode) theoretical probability density the resulting interval has the property that the density at the left endpoint is equal to the density at the right endpoint and that the density at any point outside the interval is less than the density at any point inside the interval.  Establishing this equivalence is left as an exercise for the mathematically inclined reader.  (Hint: Start with the interval defined by the \"equal density at the endpoints\" property and consider what happens if you shift that interval while maintaining the same area under the density curve.  You will be replacing a region of higher density by one with a lower density and the interval must become wider to maintain the same area.)With large samples a brute-force enumeration approach works.julia> function hpdinterval(v, level=0.95)\n    n = length(v)\n    if !(0 < level < 1)\n        throw(ArgumentError(\"level = $level must be in (0, 1)\"))\n    end\n    if (lbd = floor(Int, (1 - level) * n)) < 2\n        throw(ArgumentError(\n            \"level = $level is too large from sample size $n\"))\n    end\n    ordstat = sort(v)\n    leftendpts = ordstat[1:lbd]\n    rtendpts = ordstat[(1 + n - lbd):n]\n    (w, ind) = findmin(rtendpts - leftendpts)\n    return [leftendpts[ind], rtendpts[ind]]\nend\nhpdinterval (generic function with 2 methods)\nFor example, the 95% HPD interval calculated from the sample of beta_1 values isjulia> hpdinterval(mm1bstp[:β₁])\n2-element Array{Float64,1}:\n 1492.49\n 1561.32\nwhich is very close to the central probability interval ofjulia> quantile(mm1bstp[:β₁], [0.025, 0.975])\n2-element Array{Float64,1}:\n 1492.45\n 1561.28\nbecause the empirical distribution of the beta_1 sample is very similar to a normal distribution.  In particular, it is more-or-less symmetric and also unimodal.The HPD interval on sigma^2 is julia> hpdinterval(abs2.(mm1bstp[:σ]))\n2-element Array{Float64,1}:\n 1068.03\n 3745.88\nwhich is shifted to the left relative to the central probability intervaljulia> quantile(abs2.(mm1bstp[:σ]), [0.025, 0.975])\n2-element Array{Float64,1}:\n 1250.99\n 3966.02\nbecause the distribution of the sigma^2 sample is skewed to the right.  The HPD interval will truncate the lower density, long, right tail and include more of the higher density, short, left tail.The HPD interval does not have the property that the endpoints of the interval on sigma^2 are the squares of the endpoints of the intervals on sigma, because \"shorter\" on the scale of sigma does not necessarily correspond to shorter on the scale of sigma^2.julia> sigma95hpd = hpdinterval(mm1bstp[:σ])\n2-element Array{Float64,1}:\n 35.4844\n 63.0209\njulia> abs2.(sigma95hpd)\n2-element Array{Float64,1}:\n 1259.14\n 3971.64\nFinally, a 95% HPD interval on sigma_1 includes the boundary value sigma_1=0.julia> hpdinterval(mm1bstp[:σ₁])\n2-element Array{Float64,1}:\n  0.0   \n 54.7193\nIn fact, the confidence level or coverage probability must be rather small before the boundary value is excludedjulia> hpdinterval(mm1bstp[:σ₁], 0.798)\n2-element Array{Float64,1}:\n  9.83921\n 52.2513 \njulia> hpdinterval(mm1bstp[:σ₁], 0.799)\n2-element Array{Float64,1}:\n  0.0  \n 42.525\n"
},

{
    "location": "SimpleLMM.html#Empirical-cumulative-distribution-function-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Empirical cumulative distribution function",
    "category": "section",
    "text": "The empirical cumulative distribution function (ecdf) of a sample maps the range of the sample onto [0,1] by x → proportion of sample ≤ x.  In general this is a \"step function\", which takes jumps of size 1/length(samp) at each observed sample value.  For large samples, we can plot it as a qq plot where the theoretical quantiles are the probability points and are on the vertical axis.(Image: )The orange lines added to the plot show the construction of the central probability 80% confidence interval on sigma_1 and the red lines show the 80% HPD interval.  Comparing the spacing of the left end points to that of the right end points shows that the HPD interval is shorter, because, in switching from the orange to the red lines, the right end point moves further to the left than does the left end point.The differences in the widths becomes more dramatic on the scale of sigma_1^2(Image: )ADD: similar analysis for mm2"
},

{
    "location": "SimpleLMM.html#Assessing-the-Random-Effects-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Assessing the Random Effects",
    "category": "section",
    "text": "In Sect. [sec:definitions] we mentioned that what are sometimes called the BLUPs (or best linear unbiased predictors) of the random effects, mathcal B, are the conditional modes evaluated at the parameter estimates, calculated as tildeb_widehattheta=Lambda_widehatthetatildeu_widehattheta.These values are often considered as some sort of “estimates” of the random effects. It can be helpful to think of them this way but it can also be misleading. As we have stated, the random effects are not, strictly speaking, parameters—they are unobserved random variables. We don’t estimate the random effects in the same sense that we estimate parameters. Instead, we consider the conditional distribution of mathcal B given the observed data, (mathcal Bmathcal Y=mathbf  y).Because the unconditional distribution, mathcal BsimmathcalN(mathbf  0Sigma_theta) is continuous, the conditional distribution, (mathcal Bmathcal Y=mathbf  y) will also be continuous. In general, the mode of a probability density is the point of maximum density, so the phrase “conditional mode” refers to the point at which this conditional density is maximized. Because this definition relates to the probability model, the values of the parameters are assumed to be known. In practice, of course, we don’t know the values of the parameters (if we did there would be no purpose in forming the parameter estimates), so we use the estimated values of the parameters to evaluate the conditional modes.Those who are familiar with the multivariate Gaussian distribution may recognize that, because both mathcal B and (mathcal Ymathcal B=mathbf  b) are multivariate Gaussian, (mathcal Bmathcal Y=mathbf  y) will also be multivariate Gaussian and the conditional mode will also be the conditional mean of mathcal B, given mathcal Y=mathbf  y. This is the case for a linear mixed model but it does not carry over to other forms of mixed models. In the general case all we can say about tildemathbf    u or tildemathbf  b is that they maximize a conditional density, which is why we use the term “conditional mode” to describe these values. We will only use the term “conditional mean” and the symbol, mathbf mu, in reference to mathrmE(mathcal Ymathcal B=mathbf  b), which is the conditional mean of mathcal Y given mathcal B, and an important part of the formulation of all types of mixed-effects models.The ranef extractor returns the conditional modes.julia> ranef(mm1)  # FIXME return an ordered dict\n1-element Array{Array{Float64,2},1}:\n [-16.6282 0.369516 26.9747 -21.8014 53.5798 -42.4943]\nThe result is an array of matrices, one for each random effects term in the model.  In this case there is only one matrix because there is only one random-effects term, (1 | Batch), in the model. There is only one row in this matrix because the random-effects term, (1 | Batch), is a simple, scalar term.To make this more explicit, random-effects terms in the model formula are those that contain the vertical bar () character. The variable is the grouping factor for the random effects generated by this term. An expression for the grouping factor, usually just the name of a variable, occurs to the right of the vertical bar. If the expression on the left of the vertical bar is , as it is here, we describe the term as a _simple, scalar, random-effects term_. The designation “scalar” means there will be exactly one random effect generated for each level of the grouping factor. A simple, scalar term generates a block of indicator columns — the indicators for the grouping factor — in mathbf Z. Because there is only one random-effects term in this model and because that term is a simple, scalar term, the model matrix, 𝐙, for this model is the indicator matrix for the levels of Batch.In the next chapter we fit models with multiple simple, scalar terms and, in subsequent chapters, we extend random-effects terms beyond simple, scalar terms. When we have only simple, scalar terms in the model, each term has a unique grouping factor and the elements of the list returned by can be considered as associated with terms or with grouping factors. In more complex models a particular grouping factor may occur in more than one term, in which case the elements of the list are associated with the grouping factors, not the terms.Given the data, 𝐲, and the parameter estimates, we can evaluate a measure of the dispersion of (mathcal Bmathcal Y=mathbf y). In the case of a linear mixed model, this is the conditional standard deviation, from which we can obtain a prediction interval. The extractor is named condVar.julia> condVar(mm1)\n1-element Array{Array{Float64,3},1}:\n [362.31]\n\n[362.31]\n\n[362.31]\n\n[362.31]\n\n[362.31]\n\n[362.31]\n"
},

{
    "location": "SimpleLMM.html#Chapter-Summary-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Chapter Summary",
    "category": "section",
    "text": "A considerable amount of material has been presented in this chapter, especially considering the word “simple” in its title (it’s the model that is simple, not the material). A summary may be in order.A mixed-effects model incorporates fixed-effects parameters and random effects, which are unobserved random variables, mathcal B. In a linear mixed model, both the unconditional distribution of mathcal B and the conditional distribution, (mathcal Ymathcal B=mathbf b), are multivariate Gaussian distributions. Furthermore, this conditional distribution is a spherical Gaussian with mean, mathbfmu, determined by the linear predictor, mathbf Zmathbf b+mathbf Xmathbfbeta. That is,\\begin{equation}(\\mathcal Y|\\mathcal B=\\mathbf b)\\sim   \\mathcal{N}(\\mathbf Z\\mathbf b+\\mathbf X\\mathbf\\beta, \\sigma^2\\mathbf I_n) .\\end{equation}The unconditional distribution of mathcal B has mean mathbf 0 and a parameterized qtimes q variance-covariance matrix, Sigma_theta.In the models we considered in this chapter, Sigma_theta, is a simple multiple of the identity matrix, mathbf I_6. This matrix is always a multiple of the identity in models with just one random-effects term that is a simple, scalar term. The reason for introducing all the machinery that we did is to allow for more general model specifications.The maximum likelihood estimates of the parameters are obtained by minimizing the deviance. For linear mixed models we can minimize the profiled deviance, which is a function of mathbftheta only, thereby considerably simplifying the optimization problem.To assess the precision of the parameter estimates, we profile the deviance function with respect to each parameter and apply a signed square root transformation to the likelihood ratio test statistic, producing a profile zeta function for each parameter. These functions provide likelihood-based confidence intervals for the parameters. Profile zeta plots allow us to visually assess the precision of individual parameters. Density plots derived from the profile zeta function provide another way of examining the distribution of the estimators of the parameters.Prediction intervals from the conditional distribution of the random effects, given the observed data, allow us to assess the precision of the random effects."
},

{
    "location": "SimpleLMM.html#Notation-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Notation",
    "category": "section",
    "text": ""
},

{
    "location": "SimpleLMM.html#Random-Variables-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Random Variables",
    "category": "section",
    "text": "mathcal Y\n: The responses (n-dimensional Gaussian)\nmathcal B\n: The random effects on the original scale (q-dimensional Gaussian with mean mathbf 0)\nmathcal U\n: The orthogonal random effects (q-dimensional spherical Gaussian)Values of these random variables are denoted by the corresponding bold-face, lower-case letters: mathbf y, mathbf b and mathbf u. We observe mathbf y. We do not observe mathbf b or mathbf u."
},

{
    "location": "SimpleLMM.html#Parameters-in-the-Probability-Model-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Parameters in the Probability Model",
    "category": "section",
    "text": "mathbfbeta\n: The p-dimension fixed-effects parameter vector.\nmathbftheta\n: The variance-component parameter vector. Its (unnamed) dimension is typically very small. Dimensions of 1, 2 or 3 are common in practice.\nsigma\n: The (scalar) common scale parameter, sigma0. It is called the common scale parameter because it is incorporated in the variance-covariance matrices of both mathcal Y and mathcal U.\nmathbftheta\nThe \"covariance\" parameter vector which determines the qtimes q lower triangular matrix Lambda_theta, called the relative covariance factor, which, in turn, determines the qtimes q sparse, symmetric semidefinite variance-covariance matrix Sigma_theta=sigma^2Lambda_thetaLambda_theta that defines the distribution of mathcal B."
},

{
    "location": "SimpleLMM.html#Model-Matrices-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Model Matrices",
    "category": "section",
    "text": "mathbf X\n: Fixed-effects model matrix of size ntimes p.\nmathbf Z\n: Random-effects model matrix of size ntimes q."
},

{
    "location": "SimpleLMM.html#Derived-Matrices-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Derived Matrices",
    "category": "section",
    "text": "mathbf L_theta\n: The sparse, lower triangular Cholesky factor of Lambda_thetamathbf Zmathbf ZLambda_theta+mathbf I_q"
},

{
    "location": "SimpleLMM.html#Vectors-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Vectors",
    "category": "section",
    "text": "In addition to the parameter vectors already mentioned, we definemathbf y\n: the n-dimensional observed response vector\nmathbfeta\n: the n-dimension linear predictor,\\begin{equation}   \\mathbf{\\eta=X\\beta+Zb=Z\\Lambda_\\theta u+X\\beta} \\end{equation}mathbfmu\n: the n-dimensional conditional mean of mathcal Y given mathcal B=mathbf b (or, equivalently, given mathcal U=mathbf u)\\begin{equation}   \\mathbf\\mu=\\mathrm{E}[\\mathcal Y|\\mathcal B=\\mathbf b]=\\mathrm{E}[\\mathcal Y|\\mathcal U=\\mathbf u] \\end{equation}tildeu_theta\n: the q-dimensional conditional mode (the value at which the conditional density is maximized) of mathcal U given mathcal Y=mathbf y."
},

{
    "location": "SimpleLMM.html#Exercises-1",
    "page": "A Simple, Linear, Mixed-effects Model",
    "title": "Exercises",
    "category": "section",
    "text": "These exercises and several others in this book use data sets from the package for . You will need to ensure that this package is installed before you can access the data sets.To load a particular data set,or load just the one data setCheck the documentation, the structure () and a summary of the Rail data (Fig. [fig:Raildot]).\nFit a model with as the response and a simple, scalar random-effects term for the variable Rail. Create a dotplot of the conditional modes of the random effects.\nCreate a bootstrap simulation from the model and construct 95% bootstrap-based confidence intervals on the parameters. Is the confidence interval on sigma_1 close to being symmetric about the estimate? Is the corresponding interval on log(sigma_1) close to being symmetric about its estimate?\nCreate the profile zeta plot for this model. For which parameters are there good normal approximations?\nPlot the prediction intervals on the random effects from this model. Do any of these prediction intervals contain zero? Consider the relative magnitudes of widehatsigma_1 and widehatsigma in this model compared to those in model for the data. Should these ratios of sigma_1sigma lead you to expect a different pattern of prediction intervals in this plot than those in Fig. [fig:fm01preddot]?"
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
