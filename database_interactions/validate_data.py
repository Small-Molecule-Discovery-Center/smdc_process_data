### helper functions to validate that a column to be inserted the db is correct and fix if not
# import pandas as pd


def validate_column(valid_values=[], test_values=[]):
    """ A function to check if the values submitted are actually in the list of valid options.
    Args:
        valid_values (list): values that are allowed
        test_values (list): values to be tested
    Returns:
        mismatch_values (dict): invalid values (key) and their alternatives (val)
    """
    
    mismatch_values={}
    for v in test_values:
        v=str(v)
        if v in valid_values:
            pass
        else:
            words=v.split(' ')
            poss=[]
            for word in words:
                poss.extend([x for x in valid_values if word.lower() in x.lower()])
            mismatch_values[v]=list(set(poss))
    if len(mismatch_values)>0:
        for value, matches in mismatch_values.items():
            print(f"{value} not in list, closest matches are:")
            for i, match in enumerate(matches):
                print(i, match)
    return mismatch_values


def fix_column(df, col='owner', mismatch_values={}, match_ids=[]):
    """ A function to check if the values submitted are actually in the list of valid options.
    Args:
        df (Pandas DataFrame): df to be fixed
        col (str): column of df to be fixed
        mismatch_values (dict): dict from validate_column
        match_ids (list): a list of indices for each mismatch's replacement.
    Returns:
        newdf (Pandas DataFrame): df with fixed column
    """
    newdf=df.copy()
    for i, (value, matches) in enumerate(mismatch_values.items()):
        if value == 'nan':
            continue
        idx=match_ids[i]
        newdf.loc[newdf[col]==value, col]=matches[idx]
    return newdf

