import pandas as pd
import numpy as np
import math
from functools import partial
from scipy.optimize import fmin_l_bfgs_b as minimize_func
from scipy.stats import chi2
from scipy.special import gammaln

######### matrix2AB ######################


# the matrices for beta-binomial calculation
KS_matrix = np.array([[1, 0, 1, 1, 0, 1, 0, 0, 0], [0, 1, -1, 0, 1, -1, 0, 0, 0]])
gamma_reduce = np.array([1, -1, -1, -1, 1, 1, 1, -1, -1])


def bb_loglikelihood(params, count_matrix):
    [a, b] = params
    ab_matrix = np.array([1, 1, 1, a + b, a, b, a + b, a, b])

    # perform matrix multiplication to get inputs to log-gamma
    input_matrix = np.matmul(count_matrix, KS_matrix) + ab_matrix
    # get corresponding log-gamma values and reduce over pon-values
    # make sure, count_df is 2d to safe the check
    # if is_1d:  # check whether gammatrix is 2-dim - otherwise sum aggregation over axis 0 is faulty
    #     gamma_matrix = gammaln(input_matrix)
    # else:
    gamma_matrix = np.sum(gammaln(input_matrix), axis=0)
    # add or subtract using gamma_reduce matrix and sum to loglikelihood (scalar)
    log_likelihood = np.sum(gamma_matrix * gamma_reduce)
    return log_likelihood


def bb_loglikelihood_fitting(params, count_matrix, penalty):
    """
    Fitting params [alpha, beta] to maximize loglikelihood
    """

    # Here, we apply the penalty term of alpha and beta (default 0.5 is slightly arbitray...)
    result = penalty * math.log(sum(params)) - bb_loglikelihood(
        params, count_matrix
    )  # matrix is dim2
    return result


def fit_bb(count_matrix, pen=0.5):
    """
    Obtaining maximum likelihood estimator of beta-binomial distribution
    count_df is the array of depth-mismatch (trials, success) pairs over the PoN list for either strand
    during minimization of fitting function (max for loglikelihood) penalty term is applied to constrain alpha and beta
        Ref for L-BFGS-B algorithm:
        A Limited Memory Algorithm for Bound Constrained Optimization
        R. H. Byrd, P. Lu and J. Nocedal. , (1995),
        SIAM Journal on Scientific and Statistical Computing, 16, 5, pp. 1190-1208.
    """

    # get the matrix from the PON_string
    #  "6|0|1|0=6|0|17|0" --> array([[6,6],[0,0], [17,1],[0,0]])

    # minimize loglikelihood using L-BFGS-B algorithm
    ab_params = minimize_func(
        bb_loglikelihood_fitting,
        [20, 20],
        args=(count_matrix, pen),
        approx_grad=True,
        bounds=[(0.1, 10000000), (1, 10000000)],
    )[0]
    return "|".join([str(round(param, 5)) for param in ab_params])


def bb_loglikelihood_1d(obs_row, params):
    """
    specialized 1-d version of bb_loglikelihood for p_value of targets
    copy of code is justified by omitting one if clause in the heavily used 2-d version
    """

    [a, b] = params
    ab_matrix = np.array([1, 1, 1, a + b, a, b, a + b, a, b])
    # convert df into matrix for np.array operations that change dims
    count_matrix = obs_row.values
    # perform matrix multiplication to get inputs to log-gamma
    input_matrix = np.matmul(count_matrix, KS_matrix) + ab_matrix
    # get corresponding log-gamma values and reduce over pon-values
    gamma_matrix = gammaln(input_matrix)
    # else:
    # gamma_matrix = np.sum(gammaln(input_matrix), axis=0)
    # add or subtract using gamma_reduce matrix and sum to loglikelihood (scalar)
    log_likelihood = np.sum(gamma_matrix * gamma_reduce)
    return log_likelihood


def bb_pvalue(obs_df, params):
    """
    get the sum of exponentials of loglikelihoods (densities) per observation
    params is strand-specific [A,B]
    """

    # get the loglikelihood per observation
    obs_df["p"] = obs_df.apply(bb_loglikelihood_1d, params=params, axis=1)

    # sum up the exponentials
    p_value = np.exp(obs_df["p"]).sum()

    return p_value


def fisher_combination(p_values):

    if 0 in p_values.values():
        return 0
    else:
        return 1 - chi2.cdf(
            sum([-2 * math.log(x) for x in p_values.values()]),
            2 * len(p_values.values()),
        )
