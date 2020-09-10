import pandas as pd
from ebutils import get_count_df, retrieveABdata, get_obs_df
from ebcore import fit_bb, bb_pvalue, fisher_combination


def matrix2AB(row, config=None):
    count_df = get_count_df(row['PON:ALT=Depth'], config["count_dict"])
    return fit_bb(count_df, config['fit_pen'])


######### EB2AB ######################

def AB2EBscore(row, EBconfig):
    '''
    takes a df containing an AB column of shape "A+|A--B+|B-""
    '''
    # retrieve the data from the row
    target_s, AB_list = retrieveABdata(row, EBconfig)
    
    # get the p_values for each strand
    p_values = {}
    p_values['+'] = bb_pvalue(get_obs_df(target_s, ['depth+', 'alt+']), params_split[:2])
    p_values['-'] = bb_pvalue(get_obs_df(target_s, ['depth-', 'alt-']), params_split[2:])
    # combine p_values with fisher combination
    EB_p = fisher_combination(p_values)
    if EB_p < 1e-60:
        return 60
    if EB_p > 1.0 - 1e-10:
        return 0
    return -round(math.log10(EB_p), 3)