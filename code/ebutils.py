import pandas as pd


def get_count_df(pon_string, count_dict):
    '''
    converts the base-wise read coverage to a matrix
    '''

    data_split = [s for ad in pon_string.split("=") for s in ad.split("-")]
    count_df = pd.DataFrame.from_dict(
        {count_dict[i]: v.split("|") for i, v in enumerate(data_split)})
    return count_df.astype(int)


######### AB2EB ######################

def retrieveABdata(row, EBconfig):
    # retrieve the data from the row
    # target_s:
    # turn string "0-5=12-42" into target_s:
    # pd.Series
    # alt+       0
    # alt-       5
    # depth+    12
    # depth-    42
    target = row['Tumor:Alt=Depth']
    count_dict = EBconfig["count_dict"]
    target_split = [s for ad in target.split("=") for s in ad.split("-")]
    target_s = pd.Series(
        {count_dict[i]: v for i, v in enumerate(target_split)}).astype(int)
    # params:
    # turn string A+|B+-A-|B- into AB dict {'+':[A+, B+], '-':[A-, B-]}
    params = row['ABparams']
    AB_list = [float(ab) for s in params.split("-") for ab in s.split("|")]
    AB_dict = {'+': AB_list[:2], '-': AB_list[2:]}
    return target_s, AB_dict


def get_obs_df(target_s, cols):
    '''
    turn the target_s into obs_df
    '''

    # cols is either ['depth+', 'alt+'] or ['depth-', 'alt-']
    alt_type = cols[1]
    # creates an observation df for each observation from depth-alt to depth-depth
    n_minus_k = target_s[cols[0]] - target_s[alt_type]
    # obs_df is instantiated from target_s dict with n_minus_k + 1 rows
    obs_df = pd.DataFrame(target_s[cols].to_dict(), index=range(n_minus_k + 1))
    # alt column is incremented using index
    obs_df[alt_type] = obs_df[alt_type] + obs_df.index
    return obs_df
