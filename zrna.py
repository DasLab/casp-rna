import numpy as np

def calc_z_scores(a, lower_is_better=False):
    '''
    Calculate the z-score for each model in a.

    Args:
        a (array): array of scores
        lower_is_better (bool): whether lower is better
    Returns:
        array of z-scores
    '''
    if lower_is_better: a=-a
    ans = (a - np.nanmean(a)) / np.nanstd(a)
    # The number of models per puzzle
    print("# of models: {}".format(a.size))

    # What is the mean TM-score
    print("The TM-score mean is {}".format(np.nanmean(a)))
    # What is std for TM-scores
    print("The TM-score std is {}".format(np.std(a)))

    return ans

def remove_outliers(a, thresh=-2):
    """
    Remove outliers from a by setting them to nan.

    Args:
        a (array): array of scores
        thresh (float): threshold for outliers
    Returns:
        array of scores with outliers set to nan
    """
    a[a < thresh] = np.nan
    return a

def calc_double_z_scores(raw, thresh=-2, lower_is_better=False):
    """
    Calculate the double pass z-score as described in the paper for each model in raw.

    Args:
        raw (array): array of scores
        thresh (float): threshold for outliers
        lower_is_better (bool): whether lower is better
    Returns:
        array of double pass z-scores
    """
    raw = np.copy(raw)
    first_pass = calc_z_scores(raw, lower_is_better=lower_is_better)
    filtered = remove_outliers(first_pass,thresh=thresh)

    raw[np.isnan(filtered)] = np.nan

    second_pass = calc_z_scores(raw, lower_is_better=lower_is_better)
    second_pass[second_pass < thresh] = thresh
    second_pass[np.isnan(second_pass)] = thresh

    return second_pass

def calc_zrna(df, metric_weights):  
    '''
    Calculate the zrna score for each target in df.

    Args:
        df (DataFrame): contains the target, gr_code, and metric columns
        metric_weights (dict): maps metric name to weight (max or min)
    Returns:
        DataFrame of targets and their zrna score
    '''  
    return sum([df["z_" + metric] * metric_weights[metric] for metric in metric_weights])

def best_score_across_all_confs(df, metrics_is_lower_better):    
    '''
    Calculate the best score across all conformations for each target.

    Args:
        df (DataFrame): contains the target, gr_code, and metric columns
        metrics_is_lower_better (dict): maps metric name to whether lower is better
    '''
    return df.groupby(["target", "gr_code"]).agg(metrics_is_lower_better).reset_index()

def get_group_score(df, agg="sum>0", score="z_rna"):
    '''
    Credit: Adapated from Rachael's EM pipeline.

    Aggregate the scores of each group in gr_code

    Args:
        df (DataFrame): contain gr_code and score (argument specified) column
        agg (str): how to aggregate scores options: 'sum>0','mean','sum' (default 'sum>0')
        score (str): what column to aggregate
    Returns:
        DataFrame of groups and their summed score
    '''
    vals = df.groupby("gr_code")[score]
    if agg == "sum>0":
        return vals.apply(lambda col: col[col > 0].sum()).reset_index()
    elif agg == "mean":
        return vals.mean().reset_index()
    elif agg == "sum":
        return vals.sum().reset_index()