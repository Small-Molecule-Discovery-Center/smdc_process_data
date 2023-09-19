# functions to process user input and output Echo and Bruker SPR plate map files
# Created 20230918 by @paulsonak amanda.paulson@ucsf.edu


###############################################################################
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pandas.api.types import CategoricalDtype
import string
import sys
import os
# sys.path.append('./')
from datetime import date


###############################################################################
# known, immutable values for 384 well Bruker SPR
rows=list(string.ascii_uppercase)[0:16]*24
cols=list(np.arange(1,25))
cols=[str(x).zfill(2) for x in cols]*16
cols.sort()
all_wells=[row+str(col).zfill(2) for row,col in zip(rows, cols)]

channel_row_dict={1: 'A', 2: 'C', 3: 'E', 4: 'G', 5: 'I', 6: 'K', 7: 'M', 8: 'O'}
row_channel_dict={'A': 1,'B': 1,'C': 2,'D': 2,'E': 3,'F': 3,'G': 4,'H': 4,'I': 5,'J': 5,'K': 6,'L': 6,'M': 7,'N': 7,'O': 8,'P': 8}

solvent_corr_wells=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12',
                    'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12',
                    'E01','E02','E03','E04','E05','E06','E07','E08','E09','E10','E11','E12',
                    'G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12',
                    'I01','I02','I03','I04','I05','I06','I07','I08','I09','I10','I11','I12',
                    'K01','K02','K03','K04','K05','K06','K07','K08','K09','K10','K11','K12',
                    'M01','M02','M03','M04','M05','M06','M07','M08','M09','M10','M11','M12',
                    'O01','O02','O03','O04','O05','O06','O07','O08','O09','O10','O11','O12']

blank_wells=['A01','A05','A09','A13','A17','A21','B01','B05','B09','B13','B17','B21',
             'C01','C05','C09','C13','C17','C21','D01','D05','D09','D13','D17','D21',
             'E01','E05','E09','E13','E17','E21','F01','F05','F09','F13','F17','F21',
             'G01','G05','G09','G13','G17','G21','H01','H05','H09','H13','H17','H21',
             'I01','I05','I09','I13','I17','I21','J01','J05','J09','J13','J17','J21',
             'K01','K05','K09','K13','K17','K21','L01','L05','L09','L13','L17','L21',
             'M01','M05','M09','M13','M17','M21','N01','N05','N09','N13','N17','N21',
             'O01','O05','O09','O13','O17','O21','P01','P05','P09','P13','P17','P21']

dmso_cmpd_wells=list(set(all_wells)-set(solvent_corr_wells)-set(blank_wells))
dmso_cmpd_wells.sort()
plate2_cmpd_wells=list(set(all_wells)-set(blank_wells))
plate2_cmpd_wells.sort()

solv_corr_samples=['Warmup']*2+['Solv Corr']*10
solv_corr_vols=[1000]*2+[500,750,1000,1250,1500]*2

half_row_idxs=[0,12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,216,228,240,252,264,276,288,
          300,312,324,336,348,360,372,384,396,408,420,432,444,456,468,480,492,504,516,528,540,552,564,
          576,588,600,612,624,636,648,660,672,684,696,708,720,732,744,756]

today = date.today().strftime("%Y%m%d")


###############################################################################
def create_spr_files(input='test.txt'):
    print(f"Processing {input}...\n")

    # get dmso info
    df=pd.read_excel(input)
    name=df.iloc[2,0]
    name=name.lower().replace(' ','_')
    dmso=df.iloc[3,0]
    if dmso.lower()=='yes': dmso=True
    else: dmso=False
    dmso_pct=df.iloc[4,0]
    backfill_p1=df.iloc[5,0]
    backfill_p2=df.iloc[6,0]
    name, dmso, dmso_pct, backfill_p1, backfill_p2
    
    # read in data and format columns
    df=pd.read_excel(input, header=10)
    df.columns=[x.replace(' ','_').lower().replace('(','').replace(')','').replace('sample','sample_id').replace('sample_id_name','sample_name').replace('2:800_','').replace('plate_type_ldv_or_pp','source_plate') for x in df.columns]
    df.sample_name.unique()
    
    # check solvent correction channels and wells
    solv_sample=[1,5,9,13,17,21,25,29]
    solv_sample_vals=df[(df.plate==1)&(df.sample_id.isin(solv_sample))].sample_name.unique().tolist()
    solv_sample_vals=set([str(x).lower() for x in solv_sample_vals])
    assert (dmso) and ((solv_sample_vals=={'solv corr'}) or (solv_sample_vals=={'solv corr','nan'})), "In plate 1, samples 1,5,9,13,17,21,25,29 must be reserved for solvent correction and be labeled 'Solv Corr' if you are using DMSO."               
    
    # if solv_samples have nans, add back 'Solv Corr'
    df.loc[(df.plate==1)&(df.sample_id.isin(solv_sample)), 'sample_name']='Solv Corr'
    
    # give empty wells a sample name
    df.loc[df.sample_name.isna(), 'sample_name']='empty_'+df.loc[df.sample_name.isna(), 'sample_id'].astype(str)
    
    # get compound dilution info
    cmpd_cols=['conc', 'transfer_volume', 'source_well','real_conc']
    dfs=[]
    for i, row in df.iterrows():
        if str(row.sample_name)=='nan':
            continue
        if str(row.sample_name).lower()!='solv corr':
            tmp=pd.concat([pd.DataFrame(row).T]*9)
            stock_conc=row.stock_conc_um
            max_conc=row.highest_assay_conc_um
            dil=row.dilution_factor
            samples=[row.sample_name]*9
            concs=[]
            nl=[]
            stock_wells=[]
            real_concs=[]
            for j in np.arange(0,9):
                # calculate dilution
                concs.append(max_conc*(1/dil**j))
                
                # calculate nl in units of 2.5nL - low conc volume
                # nl.append(2.5*round(concs[-1]*100/stock_conc*1000/2.5*400))
                nl.append(2.5*np.ceil(concs[-1]*100/stock_conc*1000/2.5*400))
                
                # add well locations and update diluted nL
                if nl[-1]>=1000:
                    stock_wells.append(row.high_conc_location)
                    # nl[-1]=2.5*round(concs[-1]*100/stock_conc*1000/2.5)
                    nl[-1]=2.5*np.ceil(concs[-1]*100/stock_conc*1000/2.5)
                    real_concs.append(stock_conc*(nl[-1]/1000)/100)
                else:
                    stock_wells.append(row.diluted_stock_location)
                    real_concs.append(stock_conc/400*(nl[-1]/1000)/100)
            # create data frame for one sample
            for col, val_list in zip(cmpd_cols, [concs, nl, stock_wells, real_concs]):
                tmp[col]=val_list
        else:
            continue
        dfs.append(tmp)
    
    # concatenate all samples
    cmpd_info=pd.concat(dfs)
    cmpd_info['Backfill']=1000-cmpd_info.transfer_volume
    assert len(cmpd_info[cmpd_info.transfer_volume>1000])==0, "There are dilutions with more than 1000nL planned ejection"
    
    # preserve compound order from original file
    compcat=CategoricalDtype(categories=cmpd_info.sample_name.unique(), ordered=True)
    cmpd_info.sample_name=cmpd_info.sample_name.astype(compcat)
    
    # sort by ascending compound concentration
    cmpd_info=cmpd_info.sort_values(['plate', 'sample_name','real_conc'])
    compcat
    
    # create destination plate
    dest = _make_empty_plate()
    # make empty copy
    empty=_make_empty_plate()
    # make plate2 copy
    dest2=_make_empty_plate()
    dest2['plate']=2
    dest['plate']=1
    
    # get channels used
    channels_used=cmpd_info[~cmpd_info.sample_name.str.contains('empty')].channel.unique()
    channels_used_rows=[channel_row_dict[x] for x in channels_used]
    # get rows used
    rows_used=cmpd_info[~cmpd_info.sample_name.str.contains('empty')].dest_row.unique()
    
    # get plate-specific subsets
    plate_corr_wells=[x for x in solvent_corr_wells if x[0] in channels_used_rows]
    plate_blank_wells=[x for x in blank_wells if x[0] in rows_used]

    # fill in solvent correction and blanks only for channels/rows used
    if dmso:
        # make solvent correction wells per channel
        dest.sample_name=dest.destination_well.map(dict(zip(plate_corr_wells, solv_corr_samples*len(channels_used))))
        # make solvent correction volumes
        dest.transfer_volume=dest.destination_well.map(dict(zip(plate_corr_wells, solv_corr_vols*len(channels_used))))
    
    # solvent correction plate 1 only, then add plate 2 and continue
    dest=pd.concat([dest,dest2])
    # map blank locations and volumes where there is no dmso info
    dest.loc[dest.sample_name.isna(), 'sample_name']=dest.loc[dest.sample_name.isna(), 'destination_well'].map(dict(zip(plate_blank_wells, ['Blank']*len(plate_blank_wells))))
    dest.loc[dest.sample_name=='Blank', 'transfer_volume']=1000
    dest.shape
    
    # get ctrl wells only
    dest_ctrl=dest[~dest.sample_name.isna()].copy()
    
    # add source well info for blanks and solvent wells
    dest_ctrl.source_well=dest_ctrl.row+str(backfill_p1).zfill(2)
    dest_ctrl.loc[dest_ctrl.plate==2, 'source_well']=dest_ctrl.loc[dest_ctrl.plate==2, 'row']+str(backfill_p2).zfill(2)
    dest_ctrl.source_plate='PP'
    
    # sample well list should be in order due to the way I constructed it
    if dmso:
        sample_wells=dmso_cmpd_wells+plate2_cmpd_wells
    else:
        sample_wells=plate2_cmpd_wells*2
    
    # add list in order down cmpd_info list which should also be in order
    cmpd_info['destination_well']=sample_wells[0:len(cmpd_info)]
    cmpd_info['row']=cmpd_info.destination_well.str[0]
    cmpd_info['col']=cmpd_info.destination_well.str[1:]
    
    # # cut blanks for empty rows
    # dest_ctrl=dest_ctrl[dest_ctrl.row.isin(cmpd_info.row.unique())]
    
    # make final destination df
    dest=pd.concat([dest_ctrl,cmpd_info[['plate','row','col','sample_name', 'source_plate', 'source_well', 'destination_well','transfer_volume','real_conc','Backfill']]]).reset_index(drop=True)
    dest=dest.sort_values(['plate','destination_well']).reset_index(drop=True)
    
    # add back completely empty wells
    for plate in dest.plate.unique():
        tmp=empty[~empty.destination_well.isin(dest[dest.plate==plate].destination_well)].copy()
        tmp['plate']=plate
        dest=pd.concat([dest, tmp])
    
    # create well idx and sort dest df
    dest=dest.sort_values(['plate','destination_well'])
    dest=dest.reset_index(drop=True).reset_index(names='well_idx')
    dest.loc[dest.sample_name.astype(str).str.contains('empty'), 'sample_name']=np.nan
    
    # scan through each half-row and delete blanks if not used
    for i in half_row_idxs:
        tmp=dest[dest.well_idx.between(i,i+12, inclusive='left')]
        if tmp.sample_name.unique().tolist()==['Blank',np.nan]:
            print(f"removing {tmp.destination_well.iloc[0]} blanks")
            for col in ['sample_name', 'source_plate', 'source_well','transfer_volume']:
                dest.loc[dest.well_idx.between(i,i+12, inclusive='left'), col]=np.nan
    
    # rename plate types
    if dmso:
        dest.source_plate=dest.source_plate.replace('PP','384PP_DMSO2').replace('LDV','384LDV_DMSO')
    else:
        dest.source_plate=dest.source_plate.replace('PP','384PP_AQ_CP').replace('LDV','384LDV_AQ_B2')
    
    # cut plate 2 if no real samples
    if dest[dest.plate==2].sample_name.unique().tolist()==[np.nan]:
        dest=dest[dest.plate==1]
    
    # reformat backfill info
    backfill=dest.melt(id_vars=['sample_name', 'source_plate', 'source_well', 'destination_well','row', 'col', 'plate',],
                       value_name='new_transfer_volume', var_name='new_sample_name')
    backfill=backfill[backfill.new_sample_name=='Backfill']
    backfill=backfill[~backfill.new_transfer_volume.isna()]
    backfill=backfill[backfill.new_transfer_volume!=0]
    
    # change backfill plate and well location
    if dmso:
        backfill.source_plate='384PP_DMSO2'
    else:
        backfill.source_plate='384PP_AQ_CP'
    backfill.source_well=backfill.row+str(backfill_p1).zfill(2)
    backfill.loc[backfill.plate==2, 'source_well']=backfill.loc[backfill.plate==2, 'row']+str(backfill_p2).zfill(2)
    backfill=backfill[['new_sample_name','source_plate', 'source_well', 'destination_well', 'new_transfer_volume', 'row',
           'col', 'plate',]]
    backfill.columns=[x.replace('new_','') for x in backfill.columns]
    
    # draw
    sns.set_context('paper')
    for plate in dest.plate.unique():
        df=dest[dest.plate==plate].copy()
        # conc
        df.loc[df.sample_name.isin(['Warmup', 'Solv Corr', 'Blank']), 'real_conc']=1/1000000
        df.real_conc=df.real_conc.fillna(110)
        df.loc[df.sample_name.astype(str).str.contains('empty'), 'real_conc']=110
        df['log_conc']=np.log10(df.real_conc/1000000)
        df=df.rename(columns={'log_conc':''})
        # sample names
        df.sample_name=df.sample_name.fillna('')
        df.loc[df.sample_name.str.contains('empty'), 'sample_name']=''
        barcode=f'{today}_{name}_plate_'+str(plate)
        _draw_plate_heatmap(df=df, wells='destination_well', labels='sample_name',
                              values='', barcode=barcode, width=25, outdir="./")
        print(f"\nPlate map image saved as {barcode}.png")
    
    # add backfill info to echo plate df
    dest=pd.concat([dest, backfill])

    # remove nans
    dest_echo=dest[~dest.sample_name.isna()].copy()
    
    # remove zpadding
    dest_echo.destination_well=dest_echo.row+dest_echo.col.astype(int).astype(str)
    dest_echo['source_row']=dest_echo.source_well.str[0]
    dest_echo['source_col']=dest_echo.source_well.str[1:]
    dest_echo.source_well=dest_echo.source_row+dest_echo.source_col.astype(int).astype(str)
    
    # rename columns
    dest_echo=dest_echo[['sample_name', 'source_plate', 'source_well', 'destination_well',
           'transfer_volume', 'plate']]
    dest_echo.columns=['Sample Name','Source Plate Type','Source Well','Destination Well','Transfer Volume', 'plate']

    # save echo files for each plate
    for plate_num in dest_echo.plate.unique():
        plate=dest_echo[dest_echo.plate==plate_num]
        plate=plate.drop(columns='plate')
        plate[plate['Source Plate Type'].str.contains('PP')].to_csv(f"{today}_{name.lower().replace(' ','_')}_PP_echo_plate_{plate_num}.csv", index=False, encoding='utf-8-sig')
        print(f"\nPP echo plate {plate_num} saved to {today}_{name}_PP_echo_plate_{plate_num}.csv")
        if len(plate[plate['Source Plate Type'].str.contains('LDV')])>0:
            plate[plate['Source Plate Type'].str.contains('LDV')].to_csv(f"{today}_{name.lower().replace(' ','_')}_LDV_echo_plate_{plate_num}.csv", index=False, encoding='utf-8-sig')
            print(f"LDV echo plate {plate_num} saved to {today}_{name}_LDV_echo_plate_{plate_num}.csv")
    
    # create spr file
    spr=dest[dest.sample_name!='Backfill'].copy()
    spr['mw']=spr.sample_name.map(dict(zip(cmpd_info.sample_name, cmpd_info.mw_da)))
    spr['smiles']=spr.sample_name.map(dict(zip(cmpd_info.sample_name, cmpd_info.smiles_optional)))
    spr.loc[spr.sample_name.isin(['Warmup','Blank']), 'real_conc']=0
    spr.loc[~spr.real_conc.isna(), 'unit']=u"\u03bcM" # mu
    
    # dmso solvent corr percents
    if dmso:
        dmso_pct_range=list(np.arange(dmso_pct-0.5, dmso_pct+0.75, 0.25))
        spr.loc[spr.sample_name=='Solv Corr', 'real_conc']=dmso_pct_range*(int(len(spr[spr.sample_name=='Solv Corr'])/5))
        spr.loc[spr.sample_name=='Solv Corr', 'unit']="%"
    spr=spr[['sample_name','real_conc','unit','mw','smiles','plate']]
    
    # save bruker files
    for plate_num in spr.plate.unique():
        spr[spr.plate==plate_num].drop(columns=['plate']).to_csv(f"{today}_{name}_bruker_plate_{plate_num}.csv", index=False, header=False, encoding='utf-8-sig')
        print(f"\nBruker plate {plate_num} saved to {today}_{name}_bruker_plate_{plate_num}.csv")
    print('\n')


###############################################################################
def _make_empty_plate():
    '''Make empty echo destination plates
    Returns:
        A dataframe with echo columns and 384 rows with destination wells filled in
    '''
    empty=pd.DataFrame(columns=['sample_name','source_plate', 'source_well', 'destination_well', 'transfer_volume'])
    empty['row']=rows
    empty['col']=cols
    empty['destination_well']=all_wells
    empty=empty.sort_values('destination_well')
    return empty
    

###############################################################################
def _draw_plate_heatmap(barcode='OTPARP1234-Eco1', df=None, wells='Well', values='SMDC_ID', labels='SMDC_ID', width=25, outdir=None):
    """
    Draw a heatmap of the plate using a platemap excel style input as your df
    Args:
        df (Pandas DataFrame): data to pivot into heatmap; optional if barcode is provided.
        wells (str): name of column containing well IDs in the format [A-P]\d\d or [A-P]\d
        values (str): name of column containing heatmap values (floats)
        labels (str): name of column containing heatmap label values (strings)
        barcode (str): name of plate or barcode to look up
    Returns:
        None
    Side effects:
        Draws a plate heatmap with values as labels
    """
        
    df['Row']=df[wells].str[0]
    df['Col']=df[wells].str[1:]
    df.Col=df.Col.astype(int)

    pl  =  df.pivot_table(columns='Col', index='Row', values=values, dropna=False)
    pl  =  pl.astype(float)
    df[labels]=df[labels].astype(str).replace('nan','')
    plabels=df.pivot_table(columns='Col', index='Row', values=labels, aggfunc=lambda x: '\n'.join(x))
    for col in plabels.columns:
        plabels[col]=plabels[col].str.replace('.0', '', regex=False)

    fig, ax = plt.subplots(1, figsize=(width,width/3))
    sns.heatmap(pl, annot=plabels, fmt='s', linewidth=1, linecolor='white', ax=ax, cbar=False, annot_kws={"size": 7})
    ax.set_title(f"{barcode}");
    plt.tight_layout()
    if outdir is not None:
        plt.savefig(os.path.join(outdir, barcode+'.png'), bbox_inches='tight')
        plt.close()

