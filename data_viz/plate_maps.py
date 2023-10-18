# import mysql.connector
# from sqlalchemy import create_engine
import pandas as pd
import os
import numpy as np
import hashlib
import time
import humanize
import matplotlib.pyplot as plt
import seaborn as sns 
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
sns.set_theme(rc={'figure.dpi': 150, 'figure.figsize': (5, 3.75)})

import sys
# sys.path.append('/mnt/c/Users/akpau/Repos/data_proc')
sys.path.append('/home/apaulson/repos/data-proc/')

from mysql_cnx import create_mysql_cnx


###############################################################################  
def get_hash(path, size):
    """ A helper function to generate the md5 checksum of a file.
    Args:
        path (str): path to file
        size (int): chunk size of file to read at a time
    Returns:
        hex digest of md5 sum.
    """
    m = hashlib.md5()
    with open(path, 'rb') as f:
        b = f.read(size)
        while len(b) > 0:
            m.update(b)
            b = f.read(size)
    return m.hexdigest()    


###############################################################################
def draw_plate_heatmap(barcode='OTPARP1234-Eco1', df=None, wells='Well', values='SMDC_ID', labels='SMDC_ID', center=None, vmin=None, vmax=None, width=25, outdir=None):
    """
    Draw a heatmap of the plate using a platemap excel style file as your df
    Args:
        df (Pandas DataFrame): data to pivot into heatmap; optional if barcode is provided.
        row_col (str): name of column containing row values
        col_col (str): name of column containing column values
        values (str): name of column containing heatmap values (floats)
        labels (str): name of column containing heatmap label values (strings)
        barcode (str): name of plate or barcode to look up
    Returns:
        None
    Side effects:
        Draws a plate heatmap with values as labels
    """
    if df is None:
        df=make_platemaps_from_barcode(barcode=barcode, outdir=None)
        
    df['Row']=df[wells].str[0]
    df['Col']=df[wells].str[1:]
    df.Col=df.Col.astype(int)
    cbar=False
    pl  =  df.pivot_table(columns='Col', index='Row', values=values, dropna=False)
    pl  =  pl.astype(float)
    df[labels]=df[labels].astype(str).replace('nan','')
    plabels=df.pivot_table(columns='Col', index='Row', values=labels, aggfunc=lambda x: '\n'.join(x))
    for col in plabels.columns:
        plabels[col]=plabels[col].str.replace('.0', '', regex=False)
    if labels=='SMDC_ID':
        plabels=pl.copy()
        plabels=plabels.fillna(-99).astype(int).astype(str).replace('-99', '', regex=True)
        cbar=False
        
    # fig, ax = plt.subplots(1, figsize=(width,width/5))
    if cbar:
        fig, (ax, cax) = plt.subplots(1,2, figsize=(width,width/3),  gridspec_kw={"width_ratios":[1, .01]})
    else:
        fig, ax = plt.subplots(1, figsize=(width,width/3))
        cax=None
    sns.heatmap(pl, annot=plabels, fmt='s', linewidth=1, linecolor='white', ax=ax, cbar_ax = cax, cbar=cbar, center=center, vmin=vmin, vmax=vmax, 
                annot_kws={"size": 7}
               )
    ax.set_title(f"{barcode} {values}");
    plt.tight_layout()
    if outdir is not None:
        plt.savefig(os.path.join(outdir, barcode+'.png'), bbox_inches='tight')
        plt.close()


###############################################################################
def make_platemaps_from_barcode(barcode='OTPARP1234-Eco1', outdir='./', cnx=None):
    """
    Create an excel file with each sheet showing a plate map with vendor alias, 
    smdc_id, and lot for downstream reference.
    Args:
        barcode (str): Unique plate barcode to pull from database
        outdir (path): Place to save excel file. If None, no file will be saved.
    Returns:
        df: Long form table of all relevant plate information.
    """
    # pull data from hits
    close=False
    if cnx is None:
        close=True
        cnx = create_mysql_cnx()
    
    query=f"""
    select IW_COORDINATES, ICL_ALIAS as Vendor_Alias, ICL_IC_ID as SMDC_ID, IA_VALUE as Common_Name, IA_IAT_ID as Alias_Type,
        ICL_LOT_NUM as Lot_Num, IW_CONC_UM as Conc_uM, IW_FINAL_CONC_NM as Conc_nM, IC_SMILES as SMILES
    from INV_WELL left join INV_COMP_LOT on IW_ICL_ID = ICL_ID left join INV_ALIAS on IA_IC_ID=ICL_IC_ID left join INV_COMPOUND on IC_COMPOUND_ID=ICL_IC_ID
    where IW_IP_ID = 
        (select IP_ID from INV_PLATE where IP_BARCODE = '{barcode}');
    """
    # with cnx.connect() as conn:
    df=pd.read_sql(query, cnx)
    if close:
        cnx.dispose()
    
    # filter for common name or keep nan
    wells=[]
    for well in df.IW_COORDINATES.unique():
        tmp=df[df.IW_COORDINATES==well]
        tmp=tmp[tmp.Alias_Type=='Common_Name'] # this will be one row
        if len(tmp)==0:
            tmp=pd.DataFrame(df[df.IW_COORDINATES==well].iloc[0]).T # this will be one row
            tmp['Common_Name']=np.nan
        wells.append(tmp)
    df=pd.concat(wells)
    try:
        df=df.drop(columns='Alias_Type')
    except:
        pass
        
    # format data
    df['Row']=df.IW_COORDINATES.str.slice(start=0, stop=1)
    df['Col']=df.IW_COORDINATES.str.slice(start=1).astype(int)
    df=df.rename(columns={'IW_COORDINATES':'Well'})
    df['Plate']=barcode
    cols = ['Vendor_Alias', 'SMDC_ID', 'Common_Name', 'Lot_Num']
    
    # return or save
    if outdir is None:
        return df
    
    # save as plate layouts
    else:
        with pd.ExcelWriter(os.path.join(outdir, f'{barcode}_plate_maps.xlsx')) as writer:
            for col in cols:
                if col in ['Vendor_Alias','Common_Name']:
                    df[col]=df[col].fillna('None')
                    aggfunc=lambda x: ' '.join(x)
                else:
                    df[col]=df[col].fillna(-99).astype(int)
                    aggfunc='mean'
                tmpdf=df.pivot_table(index='Row', columns='Col', values=col, aggfunc=aggfunc)
                tmpdf=tmpdf.replace('None', np.nan).replace(-99, np.nan)
                tmpdf.index.name=col
                tmpdf.to_excel(writer, sheet_name=col)
            df=df.replace('None', np.nan).replace(-99, np.nan)
            df.to_excel(writer, sheet_name='Long_Form')
        return df


###############################################################################
def lookup_lot_by_user(ids=['285699','130958'] ,users=[("Amanda","Paulson")]):
    """ Find the most recent (highest) lotnum that the user has created for their compound.
    Args:
        ids (list of strings): The list of SMDC IDs to look up. Integers, but in string form.
        users (list of string tuples): The list of users in the form ("FirstName","LastName")
    Returns:
        fullresults (DataFrame): the SMDC IDs and Lots as two columns in a df.
    """
    # pull data from hits
    cnx = create_mysql_cnx()
    res=[]
    for fname, lname in users:
        print('looking up users..')
        query=f"""select DM_ID from DIRECTORY_MEMBER where DM_FNAME='{fname}' and DM_LNAME='{lname}';"""
        # with cnx.connect() as conn:
        user=pd.read_sql(query, cnx)
        user=user.DM_ID.iloc[0]
        
        print('looking up compounds..')
        format_strings = ','.join(['%s'] * len(ids))

        query=f"""
        select IC_COMPOUND_ID as SMDC_ID, max(ICL_LOT_NUM) as Lot from INV_COMPOUND 
        inner join INV_COMP_LOT on IC_COMPOUND_ID=ICL_IC_ID
        where IC_COMPOUND_ID in (%s) and ICL_INS_BY={user}
        group by IC_COMPOUND_ID;
        """ % format_strings
        fullresults=pd.read_sql(query, cnx, params=tuple(ids))
        res.append(fullresults)
    cnx.dispose()
    
    # create df
    print('creating df..')
    fullresults=pd.concat(res).reset_index(drop=True)
    return fullresults


###############################################################################
def lookup_lot_by_lot_num(ids=['285699','130958'] ,lot_nums=['1','1']):
    """ Find the most recent (highest) lotnum that the user has created for their compound.
    Args:
        ids (list of strings): The list of SMDC IDs to look up. Integers, but in string form.
        lot_nums (list of integers): The list of lots to look up. Integers, but in string form.
    Returns:
        fullresults (DataFrame): the SMDC IDs and Lots as two columns in a df.
    """
    # pull data from hits
    cnx = create_mysql_cnx()
    res=[]
    for fname, lname in users:
        print('looking up users..')
        query=f"""select DM_ID from DIRECTORY_MEMBER where DM_FNAME='{fname}' and DM_LNAME='{lname}';"""
        # with cnx.connect() as conn:
        user=pd.read_sql(query, cnx)
        user=user.DM_ID.iloc[0]
        
        print('looking up compounds..')
        format_strings = ','.join(['%s'] * len(ids))

        query=f"""
        select IC_COMPOUND_ID as SMDC_ID, max(ICL_LOT_NUM) as Lot from INV_COMPOUND 
        inner join INV_COMP_LOT on IC_COMPOUND_ID=ICL_IC_ID
        where IC_COMPOUND_ID in (%s) and ICL_INS_BY={user}
        group by IC_COMPOUND_ID;
        """ % format_strings
        fullresults=pd.read_sql(query, cnx, params=tuple(ids))
        res.append(fullresults)
    cnx.dispose()
    
    # create df
    print('creating df..')
    fullresults=pd.concat(res).reset_index(drop=True)
    return fullresults


###############################################################################
# def lookup_plates_by_ico_id(id=135):
#     cnx = create_mysql_cnx()
#     query=f"""
#     select IP_ID from INV_PLATE where IP_ICO_ID in 


###############################################################################
def read_platemap_block(plate_map=None, block=0, value_name='SMDC_ID'):
    """ Read a block formatted 384w plate map from excel.
    Args:
        plate_map (DataFrame): the dataframe containing the block formatted data
        block (int): the numerical ID of the block to read, 0 at the top
    Returns:
        df (DataFrame): the new dataframe containing the long-format block data.
        df columns are row, col and value_name
    """
    cols = ['row']
    cols.extend(range(1,25))
    df = plate_map.iloc[17*block:17*block+16, 0:25]
    df.columns=cols
    df = df.melt(id_vars = ['row'], value_vars=range(1,25), ignore_index=True, 
                 value_name=value_name, var_name='col')
    df.col=df.col.astype(str).str.rjust(2, '0')
    return df

###############################################################################
def make_dr_platemaps(dp, row, users, basedir):
    """ Construct dr platemap from info provided in excel file.
    TODO: refactor this so the text can be printed in top function... or add logging...
    """
    
    cnx = create_mysql_cnx()
    print(row.Grid_num)
    
    # destination plate info
    #smdc ids
    dids = read_platemap_block(plate_map=dp[row.Grid_num[0:31]], block=0, value_name='SMDC_ID')
    dids=dids.replace('DMSO', np.nan)
    if sum(dids.SMDC_ID.astype(str).str.contains('-'))>0:
        dids[['SMDC_ID','Lot']]=dids.SMDC_ID.astype(str).str.split('-', expand=True)
        dids.SMDC_ID=dids.SMDC_ID.replace('nan', np.nan).replace(u'\xa0', np.nan)
        dids.SMDC_ID=dids.SMDC_ID.astype(float)
        dids.Lot=dids.Lot.astype(float)
        
    dids['Destination well']=dids.row+dids.col
    dids['Destination plate'] = row.Grid_num
    if row.Raw_assay_ID != 'multiple':
        dids['Assay ID'] = int(row.Raw_assay_ID)
    else:
        aids=read_platemap_block(plate_map=dp[row.Grid_num[0:31]], block=2, value_name='Assay ID')
        aids['Destination well']=aids.row+aids.col
        dids=dids.merge(aids[['Destination well','Assay ID']], how='left')
    if ('Lot_info' in row.index):
        if (row.Lot_info == 'use map'):
            aids=read_platemap_block(plate_map=dp[row.Grid_num[0:31]], block=2, value_name='Lot')
            aids['Destination well']=aids.row+aids.col
            dids=dids.merge(aids[['Destination well','Lot']], how='left')

    dums= read_platemap_block(plate_map=dp[row.Grid_num[0:31]], block=1, value_name='conc')
    dums['Destination well']=dums.row+dums.col
    dids=dids.merge(dums[['Destination well','conc']], how='left')

    # source plate info
    if isinstance(row.Source_plate, str):
        if 'map' in row.Source_plate:
            sps=read_platemap_block(plate_map=dp[row.Grid_num[0:31]], block=2, value_name='Source plate')
            sps['Destination well']=sps.row+sps.col
            dids=dids.merge(sps[['Destination well','Source plate']])
            source_plates=sps[~sps['Source plate'].isna()]['Source plate'].unique().tolist()
        elif 'lot' in row.Source_plate:
            aids=read_platemap_block(plate_map=dp[row.Grid_num[0:31]], block=2, value_name='Lot')
            aids['Destination well']=aids.row+aids.col
            dids=dids.merge(aids[['Destination well','Lot']], how='left')
            source_plates=[]
        else:
            source_plates = row.Source_plate.split(';')
            source_plates = [x.strip() for x in source_plates]
            source_plates = [x for x in source_plates if x!='']
        # else:
        #     source_plates = 
    if len(source_plates)>0:
        dfs=[]
        for plate in source_plates:
            df=make_platemaps_from_barcode(barcode=plate, outdir=None, cnx=cnx)
            df=df[~df.SMDC_ID.isna()]
            dfs.append(df)
        df=pd.concat(dfs)
        df.columns=['Source well', 'alias', 'SMDC_ID', 'Common Name','Lot', 'Conc_uM', 'Conc_nM', 'SMILES','Row', 'Col', 'Source plate']
        dids = dids.merge(df[['SMDC_ID', 'Lot', 'Source plate', 'Source well']], how='left')
    else:
        dids['Source plate']='NA'
        dids['Source well']='NA'
    
    if 'Lot' not in dids.columns:
            logging.warning("no Lot from plate map")
            dids['Lot']=np.nan
    if len(dids[(~dids.SMDC_ID.isna())&(dids.Lot.isna())])>0:
        ids=dids[~dids.SMDC_ID.isna()].SMDC_ID.astype(int).astype(str).unique().tolist()
        lots=lookup_lot_by_user(ids=ids, users=users)
        dids.loc[(~dids.SMDC_ID.isna())&(dids.Lot.isna()),'Lot']=dids.loc[(~dids.SMDC_ID.isna())&(dids.Lot.isna()),'SMDC_ID'].map(dict(zip(lots.SMDC_ID,lots.Lot)))
        dids['Source plate']='NA'
        dids['Source well']='NA'

    dids=dids.fillna('NA')


    # check for zero concentrations
    if 0 in dids.conc.unique():
        minconc=input(f"There are zero's in the concentration column! {dids.conc.sort_values(ascending=False).round(4).unique().tolist()}.\nWhat would you like to replace 0 with? (Enter NA to ignore 0's): ")
        try:
            minconc=float(minconc)
            dids.loc[dids.conc==0, 'conc']=minconc
        except:
            dids=dids[dids.conc!=0]
        print(f"Min conc of plate is now {dids.conc.min()}")
    dids=dids.rename(columns={'conc':'[compound] uM'})
    dids=dids[['SMDC_ID', 'Lot', 'Source plate', 'Source well', 'Destination plate', 
               'Destination well', '[compound] uM', 'Assay ID']]
    cnx.dispose()
    
    # integrity check
    if len(dids) not in [96, 384, 1536]:
        logging.warning(f"Warning: your plate has {len(dids)} wells. Was a compound pulled from two plates?")
    
    # fill empty with NA
    dids=dids.fillna('NA')
    return dids


###############################################################################
def make_dr_platemaps_and_folders(infile='', users=[("Amanda","Paulson")], basedir='./', dose_response=True):
    """
    Create an excel file with dose response plate information for DR assay upload 
    to HiTS. Each sheet in the input file should be labeled with the destination 
    plate barcode and have SMDC IDs and Concs. There should be a sheet labeled 'Info'
    with additional data including source plate barcodes, assay IDs and dates.
    Args:
        barcode (str): Unique plate barcode to pull from database      
        basedir (path): Place to move data to and create folders in.
    Returns:
        df: Long form table of all relevant plate information.
    """
    
    # read plate layout file
    if 'xlsx' in infile:
        dp=pd.read_excel(infile, engine='openpyxl', sheet_name=None)
    else:
        dp={}
        dp['Info']=pd.read_csv(infile)
    
    info=dp['Info']
    info=info.dropna(how='all', axis='columns')
    info=info.groupby([x for x in info.columns if 'Grid_num' not in x]).agg(lambda x: ','.join(x)).reset_index()
    
    for assay, date in zip(info.Raw_assay_ID, info.Date_of_screen):
        df=info[(info.Raw_assay_ID==assay)&(info.Date_of_screen==date)].reset_index(drop=True)
        for i, row in df.iterrows():
            row.Date_of_screen=str(int(row.Date_of_screen))
            # create folders
            outdir=os.path.join(basedir, str(int(row.Raw_assay_ID))+'_'+row.Assay_name, f'toload_{i}', row.Date_of_screen)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                print(f'create folder {outdir}')
            # make and save platemaps if dose response    
            if dose_response:
                dids=make_dr_platemaps(dp, row, users, basedir)
                dids.to_excel(os.path.join(outdir, f'{row.Grid_num}_DR_plate.xlsx'), index=False)
            
            # create or modify grid_num.txt
            gridnum=pd.DataFrame(row.Grid_num.split(','), columns=['grid_num'])
            if os.path.exists(os.path.join(outdir, 'grid_num.txt')):
                gn=pd.read_csv(os.path.join(outdir, 'grid_num.txt'))
                gn=pd.concat([gn,gridnum])
                gn=gn.drop_duplicates(keep='first')
            else:
                gn=gridnum
            gn.to_csv(os.path.join(outdir, 'grid_num.txt'), index=False)
    
            # move data files to folders
            print('Copying data...')
            starttime=time.time()
            row.Raw_data_file = row.Raw_data_file.replace('\\','/').replace('Z:','').replace('Y:','').replace('/SMDC','SMDC').replace('/shared','').replace('sharedSMDC','SMDC').replace('//montreal.ucsf.edu','').replace('smb:','')
            # origfile=os.path.join('/mnt/z', row.Raw_data_file)
            origfile=os.path.join('/mnt/mac/Volumes/Shared/', row.Raw_data_file)
            # print(origfile)
            fname=row.Raw_data_file.rsplit('/',1)[-1]
            # print(fname)
            # outfile=os.path.join(outdir.replace('/Shared','').replace('/shared',''), fname)
            outfile=os.path.join(outdir, fname)
            # print(outfile)
            # if os.path.exists(outfile):
            #     # cmd="test $(md5sum %s %s | awk '{print $1}' | uniq | wc -l) == 1" % (outfile, origfile)
            #     # res=os.system(cmd)
            #     # if res==1:
            #     # if get_hash(outfile, 65536) != get_hash(origfile, 65536):
            #     cmd=f"cp '{origfile}' '{outdir.replace('/Shared','')}/'"
            #     os.system(cmd)
            # else:
            # cmd=f"cp '{origfile}' '{outdir.replace('/Shared','').replace('/shared','')}/'"
            cmd=f"cp '{origfile}' '{outdir}/'"
            os.system(cmd)
            print(f'Done copying {humanize. precisedelta(starttime - time.time(), minimum_unit="seconds", format="%d")}')
    


# ###############################################################################
# inv_well_cols=['IW_SOURCE_ID','IW_IP_ID','IW_BARCODE','IW_WELL_INDEX','IW_COORDINATES','IW_ICL_ID','IW_VOL_UL','IW_CONC_UM','IW_FINAL_CONC_NM','IW_DILUTION_FACTOR','IW_FREEZETHAW','IW_INS_BY','IW_INS_DATE','IW_CONTAM','IW_NOTES']
# def make_well_upload_list(pm):
    
    


###############################################################################
def clean_up_folders(infile=None, basedir='./', delete_orig=False):
    """ After data is uploaded into HiTS, delete files from all_data (optional)
    and move each dated folder from toload to loaded.
    Args:
        infile (path): path to the plate map info file
        basedir (path): path to the directory for which to move data
        delete_orig (bool): whether to delete the original data from its location
    """
    # read data
    if 'xlsx' in infile:
        dp=pd.read_excel(infile, engine='openpyxl', sheet_name='Info')
    else:
        dp=pd.read_csv(infile)
    dp=dp[["Raw_assay_ID","Assay_name","Date_of_screen","Raw_data_file"]].drop_duplicates()
    # delete original file from all_data or similar folder
    if delete_orig:
        for rdf in dp.Raw_data_file:
            rdf = rdf.replace('\\','/').replace('Z:','').replace('Y:','').replace('/SMDC','SMDC').replace('/shared','').replace('sharedSMDC','SMDC').replace('//montreal.ucsf.edu','').replace('smb:','')
            rdf=os.path.join('/mnt/mac/Volumes/Shared/',rdf)
            print(rdf)
            if os.path.exists(rdf):
                cmd=f"rm {rdf}"
                os.system(cmd)
    
    # move each date from toload to loaded folders
    
    for assay, date in zip(dp.Raw_assay_ID, dp.Date_of_screen):
        df=dp[(dp.Raw_assay_ID==assay)&(dp.Date_of_screen==date)].reset_index(drop=True)
        for i, row in df.iterrows():
            outdir=os.path.join(basedir, str(row.Raw_assay_ID)+'_'+row.Assay_name)
            toload=os.path.join(outdir, f'toload_{i}')
            loaded=os.path.join(outdir, 'loaded')
            if not os.path.exists(loaded):
                os.mkdir(loaded)
            print(loaded)
            for datadir in os.listdir(toload):
                oldpath=os.path.join(toload, datadir)
                oldpath=oldpath.replace('(','\(').replace(')','\)').replace(' ','\ ')
                loaded = loaded.replace('(','\(').replace(')','\)').replace(' ','\ ')
                cmd=f"mv {oldpath} {loaded}/"
                os.system(cmd)
                cmd=f"rm -r {oldpath}"
                os.system(cmd)
        