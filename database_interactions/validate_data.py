### helper functions to validate that a column to be inserted the db is correct and fix if not
import pandas as pd
pd.options.mode.copy_on_write = True
pd.set_option('future.no_silent_downcasting', True)
import os
import numpy as np
from mysql_cnx import create_mysql_cnx
from datetime import date

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
            if value=='nan':
                continue
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

def check_notebook(df, cnx=None):
    if cnx is None:
        cnx=create_mysql_cnx(host='avidddb.ucsf.edu', db='QB3_HiTS_AVIDD')
    
    query="""
    select IN_SDESC from INV_NOTEBOOK
    """
    notebooks=pd.read_sql(query, cnx)
    checks=df[df.notebook_name.notna()].notebook_name.unique()
    missing=set(checks)-set(notebooks.IN_SDESC)
    assert len(missing)==0, f"missing notebook(s)\n{df[df.notebook_name.isin(missing)][['owner','notebook_name']].drop_duplicates()}"
    print("All notebooks are in the system.")

def check_owners(df, cnx=None):
    if cnx is None:
        cnx=create_mysql_cnx(host='avidddb.ucsf.edu', db='QB3_HiTS_AVIDD')
    
    query="""
    select DM_FNAME, DM_LNAME from DIRECTORY_MEMBER
    """
    owners=pd.read_sql(query, cnx)
    owners=owners.DM_FNAME+' '+owners.DM_LNAME
    mismatches=validate_column(valid_values=owners.tolist(), test_values=df.owner.unique())
    if len(mismatches)>0:
        return mismatches
    else:
        print("Owners OK")

def check_salts(df, cnx=None):
    if cnx is None:
        cnx=create_mysql_cnx(host='avidddb.ucsf.edu', db='QB3_HiTS_AVIDD')
    
    query="""
    select distinct(ICL_SALT) from INV_COMP_LOT
    """
    salts=pd.read_sql(query, cnx)
    
    salts=salts[salts.ICL_SALT.notna()].ICL_SALT.tolist()
    mismatches=validate_column(valid_values=salts, test_values=df[df.salt.notna()].salt.unique())
    if len(mismatches)>0:
        return mismatches
    else:
        print("Salts OK")

def check_solvates(df, cnx=None):
    if cnx is None:
        cnx=create_mysql_cnx(host='avidddb.ucsf.edu', db='QB3_HiTS_AVIDD')
    
    query="""
    select distinct(ICL_SOLVATE) from INV_COMP_LOT
    """
    solvs=pd.read_sql(query, cnx)
    solvs=solvs[solvs.ICL_SOLVATE.notna()].ICL_SOLVATE.tolist()
    mismatches=validate_column(valid_values=solvs, test_values=df[df.solvate.notna()].solvate.unique())
    if len(mismatches)>0:
        return mismatches
    else:
        print("Solvates OK")

def check_vendors(df, cnx=None):
    if cnx is None:
        cnx=create_mysql_cnx(host='avidddb.ucsf.edu', db='QB3_HiTS_AVIDD')
    
    query="""
    select distinct(IV_NAME) from INV_VENDOR
    """
    vendors=pd.read_sql(query, cnx)
    vendors=vendors[vendors.IV_NAME.notna()].IV_NAME.tolist()
    mismatches=validate_column(valid_values=vendors, test_values=df[df.vendor_name.notna()].vendor_name.unique())
    if len(mismatches)>0:
        return mismatches
    else:
        print("Vendors OK")

def check_lot_types(df, cnx=None):
    if cnx is None:
        cnx=create_mysql_cnx(host='avidddb.ucsf.edu', db='QB3_HiTS_AVIDD')
    
    query="""
    select distinct(ILT_LDESC) from INV_LOT_TYPE
    """
    lots=pd.read_sql(query, cnx)
    lots=lots[lots.ILT_LDESC.notna()].ILT_LDESC.tolist()
    mismatches=validate_column(valid_values=lots, test_values=df[df.lot_type.notna()].lot_type.unique())
    if len(mismatches)>0:
        return mismatches
    else:
        print("Lot types OK")

def check_avidd_input(df):
    assert 'po_number' not in df.columns, "Check lot notes vs PO number for Andrii"
    assert len(df[(df.lot_type=='synthesized internally')&(df.notebook_name.isna())])==0, "Missing notebook name"
    assert len(df[(df.lot_type=='synthesized internally')&(df.notebook_page.isna())])==0, "Missing notebook page"
    assert len(df[df.vendor_name=='(NULL)'])==0, '(Null) vendor name'
    assert len(df[df.vendor_name.isna()])==0, "Missing Vendor Name"
    try: 
        df.PARENT_AVIDD=df.PARENT_AVIDD.astype(float)
    except Exception:
        "Parent AVIDD is not numeric"
    assert len(df[(~df.salt.isna())&(df.salt_num.isna())])==0, "Missing salt num"
    assert len(df[(~df.solvate.isna())&(df.solvate_num.isna())])==0, "Missing solvate num"
    assert len(df[df.lot_alias.isna()])==0, "Missing lot alias"
    try:
        df.notebook_page=df.notebook_page.astype(float)
    except Exception:
        "Notebook page is not numeric"
    assert len(df[df.date_prepared.isna()])==0, "Missing date prepared"
    assert len(df[df.alias1.isna() & df.alias1_value.notna()])==0, "Missing alias1 type"
    assert len(df[df.alias2.isna() & df.alias2_value.notna()])==0, "Missing alias2 type"
    assert len(df[df.alias3.isna() & df.alias3_value.notna()])==0, "Missing alias3 type (Target)"
    check_notebook(df)
    print("Input file looks good!")
    return df

def concat_reg_files(regdir="./reg_files"):
    regfiles=[]
    for root, dir, files in os.walk(regdir):
        for file in files:
            # if not file.endswith('csv'):
            regfiles.append(os.path.join(regdir, file))
    regfiles.sort()

    dfs=[]
    for file in regfiles:
        if "._" in file:
            continue
        if "~$" in file:
            continue
        if ".DS_Store" in file:
            continue
        
        if file.endswith('xlsx'):
            df=pd.read_excel(file, sheet_name='Template', header=1)
        elif file.endswith('csv'):
            df=pd.read_csv(file, header=1)
            if 'date_prepared' not in df.columns:
                df=pd.read_csv(file)
        df=df.rename(columns={'provider_name':'vendor_name', 'provider_id':'lot_alias'})
        dfs.append(df)
    df=pd.concat(dfs)
    df=df.dropna(how='all', axis='index')
    df=df.reset_index(drop=True)
    return df

def clean_reg_files(df, today=None):
    df=df[df.smiles.notna()]
    if today is None:
        today=date.today().strftime("%Y%m%d")
    df.loc[:,'date_prepared']=pd.to_datetime(df.date_prepared.astype(str), format="mixed")
    df.notebook_name=df.notebook_name.astype(str)
    df.loc[:,'notebook_name']=[x.replace('.0','') for x in df.notebook_name]
    df.loc[:,'notebook_name']=df.notebook_name.replace('nan', np.nan)
    # fix aliases and parent_avidd IDs
    aliases=['lot_alias','alias1_value','alias2_value','alias3_value']
    df[aliases]=df[aliases].astype(str)
    for col in aliases:
        df.loc[df[col].str.startswith('RLA'), col]=df.loc[df[col].str.startswith('RLA'), col].str.replace('RLA ','RLA').str.replace('RLA-','RLA').str.replace('RLA','RLA-')
        if '_value' in col:
            basecol=col.replace('_value','')
            df[basecol]=df[basecol].astype(str)
            df.loc[df[col].str.startswith('RLA-'), basecol]='RLA_ID'
            for i, row in df.loc[df[col].str.startswith('RLA-')].iterrows():
                if row[basecol].isnumeric():
                    row[basecol]='RLA'+row[basecol].astype(int).astype(str)
            df[basecol]=df[basecol].replace('nan',None)
    df[aliases]=df[aliases].replace('nan',None)
    df.loc[df.lot_type!="synthesized internally", 'analytical']=None
    df.loc[:,'PARENT_AVIDD']=df.PARENT_AVIDD.astype(str).str.replace('AVI-','').str.replace('NO','nan').str.replace('AVI','').astype(float)

    # fix well locations
    if 'rack_row' in df.columns:
        df.loc[:,'rack_col']=df['rack_col'].fillna(0)
        df.loc[:,'rack_well']=df.loc[:,'rack_row']+df.rack_col.astype(int).astype(str).str.zfill(2)
        df.loc[:,'vial_id']=df.vial_id.fillna(0).astype(int).replace('0',np.nan)
        df.loc[:,'rack_col']=df['rack_col'].replace(0,np.nan)
        df.loc[:,'vial_id']=df.vial_id.replace('out of stock', np.nan)
    df=df[df.project.notna()]
    # truncate notebook names
    df.notebook_name=[str(x)[0:20] for x in df.notebook_name]
    df.notebook_name=df.notebook_name.replace('nan', np.nan)
    return df


