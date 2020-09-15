import pandas as pd
import math

from ebutils import get_count_df, retrieveABdata, get_obs_df
from ebcore import fit_bb, bb_pvalue, fisher_combination

default_config = {
    "fit_pen": 0.5,
    "count_dict": {
        0: "alt+",
        1: "alt-",
        2: "depth+",
        3: "depth-"}
}


def matrix2AB(row, config=default_config):
    count_df = get_count_df(row['PON:ALT=Depth'], config["count_dict"])
    return fit_bb(count_df, config['fit_pen'])


######### EB2AB ######################

def AB2EBscore(row, config=default_config):
    '''
    takes a df containing an AB column of shape "A+|A--B+|B-""
    '''
    # retrieve the data from the row
    target_s, AB_dict = retrieveABdata(row, config)

    # get the p_values for each strand
    p_values = {}
    p_values['+'] = bb_pvalue(get_obs_df(target_s,
                                         ['depth+', 'alt+']), AB_dict['+'])
    p_values['-'] = bb_pvalue(get_obs_df(target_s,
                                         ['depth-', 'alt-']), AB_dict['-'])
    # combine p_values with fisher combination
    EB_p = fisher_combination(p_values)
    if EB_p < 1e-60:
        return 60
    if EB_p > 1.0 - 1e-10:
        return 0
    return -round(math.log10(EB_p), 3)
