# mass spec tethering analysis code
# 11/26/2025 AKP
# For SMDC
import os
import numpy as np
import pandas as pd
import sys

try:
    sys.path.append('/Users/apaulson/repos/data-proc/')
    import plate_maps as pm
except ImportError:
    pass

import logging
logger = logging.getLogger()
import matplotlib.pyplot as plt
import seaborn as sns
import shutil

def extract_rpt_data(infile, samples_only=False, outdir="./"):
    """Extracts data from .rpt file and saves individual sample data to ./samples/ directory.
    
    Args:
        infile (str): Path to the .rpt file.
        HA_ID (int): High-throughput assay ID.
        outdir (str): Output directory to save intermediate files.
    Returns:
        pd.DataFrame: DataFrame containing extracted sample data."""
        
    file=infile
    print(file)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    cmd=f"grep -n '^Sample\t[0-9]' '{file}' > {os.path.join(outdir, 'samples.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^Well\t' '{file}' > {os.path.join(outdir, 'wells.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^FileName\t' '{file}' > {os.path.join(outdir, 'fnames.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^BPM\t' '{file}' > {os.path.join(outdir, 'bpms.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^BPI\t' '{file}' > {os.path.join(outdir, 'bpis.txt')}"
    os.system(cmd)

    cmd=f"wc -l < '{file}' > {os.path.join(outdir, 'linenum.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^;Mass' '{file}' > {os.path.join(outdir, 'massstart.txt')}"
    os.system(cmd)

    search='}'
    cmd=f"grep -n {search} '{file}' > {os.path.join(outdir, 'massend.txt')}"
    os.system(cmd)

    samples = pd.read_csv(os.path.join(outdir, 'samples.txt'), sep=':', header=None)
    wells = pd.read_csv(os.path.join(outdir, 'wells.txt'), sep=':', header=None)
    fnames = pd.read_csv(os.path.join(outdir, 'fnames.txt'), sep=':', header=None)
    bpms = pd.read_csv(os.path.join(outdir, 'bpms.txt'), sep=':', header=None)
    bpis = pd.read_csv(os.path.join(outdir, 'bpis.txt'), sep=':', header=None)

    samples.columns=['linenum','Sample']
    wells.columns=['linenum','Plate','Well']
    fnames.columns=['linenum','Fname']
    bpms.columns=['linenum','BPM']
    bpis.columns=['linenum','BPI']

    samples.Sample=samples.Sample.str.replace('Sample\t','').astype(int)
    wells=wells[wells.Well.notna()]
    wells.Plate=wells.Plate.str.replace('Well\t','').astype(int)

    fnames.Fname=fnames.Fname.str.replace('FileName\t','')
    bpms.BPM=bpms.BPM.str.replace('BPM\t','').astype(float)
    bpis.BPI=bpis.BPI.str.replace('BPI\t','').astype(float)

    wells['Sample']=np.nan
    fnames['Sample']=np.nan
    bpms['Sample']=np.nan
    bpis['Sample']=np.nan

    starts=pd.read_csv(os.path.join(outdir, 'massstart.txt'), sep=':', header=None)
    starts.columns=['linenum','headers']

    ends=pd.read_csv(os.path.join(outdir, 'massend.txt'), sep=':', header=None)
    ends.columns=['linenum','value']

    starts['Sample']=np.nan
    ends['Sample']=np.nan

    nl=pd.read_csv(os.path.join(outdir, 'linenum.txt'), header=None)
    linenum=nl.iloc[0][0]+1

    nsamples=len(samples)

    for i in range(nsamples):
        if i==nsamples-1:
            linerange=range(samples.linenum.loc[i], linenum)
        else:
            linerange=range(samples.linenum.loc[i], samples.linenum.loc[i+1]+1)
        sample=samples.Sample.loc[i]
        wells.loc[wells.linenum.isin(linerange), 'Sample']=sample
        fnames.loc[fnames.linenum.isin(linerange), 'Sample']=sample
        bpis.loc[bpis.linenum.isin(linerange), 'Sample']=sample
        bpms.loc[bpms.linenum.isin(linerange), 'Sample']=sample
        starts.loc[starts.linenum.isin(linerange), 'Sample']=sample
        ends.loc[ends.linenum.isin(linerange), 'Sample']=sample

    endrows=[]
    for i in range(len(starts)):
        if i==len(starts)-1:
            linerange=range(starts.linenum.loc[i], linenum)
        else:
            linerange=range(starts.linenum.loc[i], starts.linenum.loc[i+1]+1)
        tmp=pd.DataFrame(ends[ends.linenum.isin(linerange)].iloc[0]).T
        endrows.append(tmp)
    endrows=pd.concat(endrows)

    df=samples.merge(wells, how='left', on='Sample', suffixes=['_sample','_well'])
    df=df.merge(fnames, how='left', on='Sample')
    df=df.merge(bpis, how='left', on='Sample', suffixes=['_fname','_bpi'])
    df=df.merge(bpms, how='left', on='Sample')
    df=df.merge(starts, how='left', on='Sample', suffixes=['_bpm','_start'])
    df=df.merge(endrows, how='left', on='Sample')
    df=df[['Sample', 'Plate', 'Well', 'Fname', 'BPI', 'BPM', 'linenum_start', 'linenum']]
    df.columns=['Sample', 'Plate', 'Well', 'Fname', 'intensitymax', 'massatmax', 'linenum_start', 'linenum']
    df.linenum_start=df.linenum_start.replace(np.nan, -99).astype(int)

    if not os.path.exists(os.path.join(outdir, 'samples')):
        os.makedirs(os.path.join(outdir, 'samples'))
    
    file=infile
    
    with open(file, "rb") as f:
        content=f.readlines()

    for i, row in df.iterrows():
        if row.linenum_start==-99:
            continue
        samplefile=os.path.join(outdir, 'samples', f'Sample_{str(row.Sample).zfill(3)}_data.txt')
        with open(samplefile, "wb") as f:
            f.writelines(content[row.linenum_start-1:row.linenum-1])
    
    if samples_only:
        os.remove(os.path.join(outdir, 'samples.txt'))
        os.remove(os.path.join(outdir, 'wells.txt'))
        os.remove(os.path.join(outdir, 'fnames.txt'))
        os.remove(os.path.join(outdir, 'bpms.txt'))
        os.remove(os.path.join(outdir, 'bpis.txt'))
        os.remove(os.path.join(outdir, 'linenum.txt'))
        os.remove(os.path.join(outdir, 'massstart.txt'))
        os.remove(os.path.join(outdir, 'massend.txt'))
        df[['Sample', 'Plate', 'Well', 'Fname', 'intensitymax', 'massatmax',]].to_csv(os.path.join(outdir, 'samples', 'extracted_samples.csv'), index=False)

        samples_dir = os.path.join(outdir, 'samples')
        zip_base = os.path.join(outdir, 'samples')
        archive_path = shutil.make_archive(zip_base, 'zip', root_dir=samples_dir)
        print(f"Extracted and saved {len(df)} samples to the Samples folder; zipped folder is at {archive_path}.")
        return df[['Sample', 'Plate', 'Well', 'Fname', 'intensitymax', 'massatmax',]]
    return df


def calculate_adduct_and_cap_ranges(HA_ID, barcode, pmirm=False, cmirm=True, rmirm=True, scan_range=5, masses=[65036, 65215]):
    """Calculates adduct and cap mass ranges for tethering analysis.
    
    Args:
        HA_ID (int): High-throughput assay ID.
        barcode (str): Plate barcode.
        pmirm (bool): Whether protein mass is in reacted mass form.
        cmirm (bool): Whether cap mass is in reacted mass form.
        rmirm (bool): Whether reductant mass is in reacted mass form.
        scan_range (float): Scan range for mass binning.
        masses (list): List of protein masses to analyze.
    Returns:
        pd.DataFrame: DataFrame containing adduct and cap mass ranges."""

    # lookup info in db for assay and plate
    query1=f"""SELECT
    (SELECT HAPA_VALUE
    FROM HTS_ASSAY_PARAMETER
    INNER JOIN HTS_ASSAY ON HA_ID = HAPA_HA_ID
    WHERE HA_ID = {HA_ID} AND HAPA_HAPAT_ID = 57 AND HAPA_STATUS_ID = 1) AS 'site_count',
    (SELECT HAPA_VALUE
    FROM HTS_ASSAY_PARAMETER
    INNER JOIN HTS_ASSAY ON HA_ID = HAPA_HA_ID
    WHERE HA_ID = {HA_ID} AND HAPA_HAPAT_ID = 53 AND HAPA_STATUS_ID = 1) AS 'protein_mass',
    (SELECT HAPA_VALUE
    FROM HTS_ASSAY_PARAMETER
    INNER JOIN HTS_ASSAY ON HA_ID = HAPA_HA_ID
    WHERE HA_ID = {HA_ID}  AND HAPA_HAPAT_ID = 54 AND HAPA_STATUS_ID = 1) AS 'reductant_mass',
    (SELECT HA_ID
    FROM HTS_ASSAY
    WHERE HA_ID = {HA_ID}) AS HA_ID"""
    query2=f"""
    SELECT IP_ID, IW_ID, ICL_IC_ID, ICL_ID, GROUP_CONCAT(ICAF_FEATURE SEPARATOR ";") AS adduct_mass
    FROM INV_WELL
    LEFT OUTER JOIN INV_COMP_LOT ON IW_ICL_ID = ICL_ID
    LEFT OUTER JOIN INV_PLATE ON IP_ID = IW_IP_ID
    LEFT OUTER JOIN INV_COMP_ANNOT ON ICL_IC_ID = ICA_IC_ID AND ICA_ICAT_ID =  8
    LEFT OUTER JOIN INV_COMP_ANNOT_FEATURE ON ICAF_ICA_ID = ICA_ID
    WHERE IP_BARCODE = '{barcode}'
    GROUP BY IW_ID"""
    query3=f"""
    SELECT IP_ID, IW_ID, IW_COORDINATES as Well, ICL_IC_ID, ICL_ID, GROUP_CONCAT(ICAF_FEATURE SEPARATOR ";") AS cap_mass
    FROM INV_WELL
    LEFT OUTER JOIN INV_COMP_LOT ON IW_ICL_ID = ICL_ID
    LEFT OUTER JOIN INV_PLATE ON IP_ID = IW_IP_ID
    LEFT OUTER JOIN INV_COMP_ANNOT ON ICL_IC_ID = ICA_IC_ID AND ICA_ICAT_ID =  13
    LEFT OUTER JOIN INV_COMP_ANNOT_FEATURE ON ICAF_ICA_ID = ICA_ID
    WHERE IP_BARCODE = '{barcode}'
    GROUP BY IW_ID"""

    cnx = pm.create_mysql_cnx()
    nums=pd.read_sql(query1, cnx)
    adds=pd.read_sql(query2,cnx)
    caps=pd.read_sql(query3,cnx)
    cnx.dispose()

    addcaps1=adds.merge(caps)
    addcaps2=addcaps1.copy()
    hydrogen = 1.007825

    if masses is None:
        masses = [nums.protein_mass.astype(float).iloc[0]]
        addcaps_list=[addcaps1]
    elif len(masses)==1:
        addcaps_list=[addcaps1]
    else:
        addcaps_list=[addcaps1, addcaps2]
    
    dfs=[]
    for addcaps, mass in zip(addcaps_list, masses):

        protein_mass=mass
        site_count = nums.site_count.astype(float).iloc[0]
        reductant_mass = nums.reductant_mass.astype(float).iloc[0]

        fix_hydrogen=False
        if pmirm:
            protein_mass=protein_mass+hydrogen # add hydrogen to get free protein mass
            fix_hydrogen=True
        if not cmirm:
            caps.cap_mass=caps.cap_mass-hydrogen # remove hydrogen to get reacted mass
        if not rmirm:
            reductant_mass=reductant_mass-hydrogen # remove hydrogen to get reacted mass

        # free protein
        free_protein_min=protein_mass-scan_range
        free_protein_max=protein_mass+scan_range

        if fix_hydrogen:
            # remove hydrogen again
            protein_mass=protein_mass-hydrogen

        addcaps['protein_mass']=protein_mass
        addcaps['free_protein_min']=free_protein_min
        addcaps['free_protein_max']=free_protein_max
        
        # throw error if adduct_mass not defined
        logger.info(f"Adduct mass is not defined for {sum(addcaps.adduct_mass.isna())} rows.")
        addcaps=addcaps[~addcaps.adduct_mass.isna()]
        addcaps=addcaps[~addcaps.adduct_mass.astype(str).str.contains(';')]
        # addcaps.loc[addcaps.ICL_IC_ID==1128111, 'adduct_mass']=297.78

        # protein-adduct
        addcaps['protein_1adduct_min']=protein_mass+addcaps.adduct_mass.astype(float)-scan_range
        addcaps['protein_1adduct_max']=protein_mass+addcaps.adduct_mass.astype(float)+scan_range
        #protein-reductant
        addcaps['protein_1reductant_min']=protein_mass+reductant_mass-scan_range
        addcaps['protein_1reductant_max']=protein_mass+reductant_mass+scan_range
        #protein-cap
        addcaps['protein_1cap_min']=protein_mass+addcaps.cap_mass.astype(float)-scan_range
        addcaps['protein_1cap_max']=protein_mass+addcaps.cap_mass.astype(float)+scan_range

        # two binding sites
        if site_count >1:
            logger.info("Calculating 2 binding site adducts")
            # TODO: add checks for adduct mass existence
            # protein + 2 adduct
            addcaps['protein_2adduct_min']=protein_mass + 2*addcaps.adduct_mass.astype(float) - scan_range
            addcaps['protein_2adduct_max']=protein_mass + 2*addcaps.adduct_mass.astype(float) + scan_range
            # protein + adduct + cap
            addcaps['protein_1adduct_1cap_min'] = protein_mass + addcaps.adduct_mass.astype(float) + addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_1adduct_1cap_max'] = protein_mass + addcaps.adduct_mass.astype(float) + addcaps.cap_mass.astype(float) + scan_range
            # protein + adduct + reductant
            addcaps['protein_1reductant_1adduct_min'] = protein_mass + reductant_mass + addcaps.adduct_mass.astype(float) - scan_range
            addcaps['protein_1reductant_1adduct_max'] = protein_mass + reductant_mass + addcaps.adduct_mass.astype(float) + scan_range
            # protein + 2 reductant
            addcaps['protein_2reductant_min'] = protein_mass + 2*reductant_mass - scan_range
            addcaps['protein_2reductant_max'] = protein_mass + 2*reductant_mass + scan_range
            # protein + reductant + cap
            addcaps['protein_1reductant_1cap_min'] = protein_mass + reductant_mass + addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_1reductant_1cap_max'] = protein_mass + reductant_mass + addcaps.cap_mass.astype(float) + scan_range
            # protein + 2 cap
            addcaps['protein_2cap_min'] = protein_mass + 2*addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_2cap_max'] = protein_mass + 2*addcaps.cap_mass.astype(float) + scan_range
            
        if site_count>2:
            logger.info("Calculating 3 binding site adducts")    
            # protein + 3 adduct
            addcaps['protein_3adduct_min']=protein_mass + 3*addcaps.adduct_mass.astype(float) - scan_range
            addcaps['protein_3adduct_max']=protein_mass + 3*addcaps.adduct_mass.astype(float) + scan_range
            
            # protein + adduct + cap
            addcaps['protein_2adduct_1cap_min'] = protein_mass + 2*addcaps.adduct_mass.astype(float) + addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_2adduct_1cap_max'] = protein_mass + 2*addcaps.adduct_mass.astype(float) + addcaps.cap_mass.astype(float) + scan_range
            
            addcaps['protein_1adduct_2cap_min'] = protein_mass + addcaps.adduct_mass.astype(float) + 2*addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_1adduct_2cap_max'] = protein_mass + addcaps.adduct_mass.astype(float) + 2*addcaps.cap_mass.astype(float) + scan_range
            
            # protein + adduct + reductant
            addcaps['protein_2reductant_1adduct_min'] = protein_mass + 2*reductant_mass + addcaps.adduct_mass.astype(float) - scan_range
            addcaps['protein_2reductant_1adduct_max'] = protein_mass + 2*reductant_mass + addcaps.adduct_mass.astype(float) + scan_range
            
            addcaps['protein_1reductant_2adduct_min'] = protein_mass + reductant_mass + 2*addcaps.adduct_mass.astype(float) - scan_range
            addcaps['protein_1reductant_2adduct_max'] = protein_mass + reductant_mass + 2*addcaps.adduct_mass.astype(float) + scan_range
            
            # protein + adduct + cap + reductant
            addcaps['protein_1adduct_1cap_1reductant_min'] = protein_mass + reductant_mass + addcaps.cap_mass.astype(float) + addcaps.adduct_mass.astype(float) - scan_range
            addcaps['protein_1adduct_1cap_1reductant_max'] = protein_mass + reductant_mass + addcaps.cap_mass.astype(float) + addcaps.adduct_mass.astype(float) + scan_range
            
            # protein + 3 reductant
            addcaps['protein_3reductant_min'] = protein_mass + 3*reductant_mass - scan_range
            addcaps['protein_3reductant_max'] = protein_mass + 3*reductant_mass + scan_range
            
            # protein + reductant + cap
            addcaps['protein_2reductant_1cap_min'] = protein_mass + 2*reductant_mass + addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_2reductant_1cap_max'] = protein_mass + 2*reductant_mass + addcaps.cap_mass.astype(float) + scan_range
            
            addcaps['protein_1reductant_2cap_min'] = protein_mass + reductant_mass + 2*addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_1reductant_2cap_max'] = protein_mass + reductant_mass + 2*addcaps.cap_mass.astype(float) + scan_range
            
            # protein + 3 cap
            addcaps['protein_3cap_min'] = protein_mass + 3*addcaps.cap_mass.astype(float) - scan_range
            addcaps['protein_3cap_max'] = protein_mass + 3*addcaps.cap_mass.astype(float) + scan_range
        
        dfs.append(addcaps)
    addcaps=pd.concat(dfs)
    
    return addcaps, masses, site_count


def analyze_rpt_data(df, addcaps, masses, site_count, scan_range, assay_class="HTS", outdir="./"): 
    """Analyzes extracted .rpt data for tethering analysis.
    
    Args:
        df (pd.DataFrame): DataFrame containing extracted sample data.
        addcaps (pd.DataFrame): DataFrame containing adduct and cap mass ranges.
        masses (list): List of protein masses to analyze.
        site_count (int): Number of binding sites on the protein.
        outdir (str): Output directory containing sample data files.
    Returns:
        pd.DataFrame: DataFrame with analyzed results."""   
    
    df[['string']]=''
    df[['protein','protein_adduct','protein_secondary','protein_2adduct','protein_3adduct','well_sum','well_cnt']]=0

    df[['row','col']]=df.Well.str.split(',', expand=True)
    df.Well=df.row+df.col.str.zfill(2)
    addcaps['Plate']=df.Plate.iloc[0]

    df=df.merge(addcaps, how='outer')
    df=df.dropna(subset=['Sample'])
    logger.info(f"Analyzing data for {len(df)} samples.")
    
    if len(masses)==2:
        df2=df[df.protein_mass==masses[1]]
        df1=df[df.protein_mass==masses[0]]
    
        df1=df1.reset_index(drop=True)
        df2=df2.reset_index(drop=True)
    else:
        df1=df
        df1=df1.reset_index(drop=True)
        df2=None

    rows=[]
    for i, row in df1.iterrows():
        if row.linenum_start==-99:
            continue
        # read raw data
        file=f'{os.path.join(outdir, "samples")}/Sample_{str(row.Sample).zfill(3)}_data.txt'
        try:
            dat=pd.read_csv(file, sep='\t')
        except Exception as e:
            logger.warning(f"Error reading {file}; {e}")
            continue
        dat.columns=['mass', 'intensity']
        
        # create string annot of dat
        dat['string']=pd.cut(dat.intensity, [-1., 10.5, 30., 50., 70., 100., 200.], labels=['','.','-','i','l','!'])
        row['string']=dat['string'].astype(str).str.cat(sep='')
        
        # add protein intensities, remove from dat
        # protein with no adduct
        tmp=dat[(dat.mass.between(row.protein_1cap_min, row.protein_1cap_max))|(dat.mass.between(row.free_protein_min, row.free_protein_max))|(dat.mass.between(row.protein_1reductant_min, row.protein_1reductant_max))]
        row.protein=tmp.intensity.sum()
        dat=dat[~dat.index.isin(tmp.index)]
        row2=None
        # now, count and remove the second protein mass data from dat
        if df2 is not None:    
            row2=df2[df2.Sample==row.Sample].iloc[0]
            # protein with no adduct
            tmp=dat[(dat.mass.between(row2.protein_1cap_min, row2.protein_1cap_max))|(dat.mass.between(row2.free_protein_min, row2.free_protein_max))|(dat.mass.between(row2.protein_1reductant_min, row2.protein_1reductant_max))]
            row.protein=row.protein+tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
        
        # add (Single) tethered intensities, remove from dat
        tmp=dat[dat.mass.between(row.protein_1adduct_min, row.protein_1adduct_max)]
        row.protein_adduct=tmp.intensity.sum()
        dat=dat[~dat.index.isin(tmp.index)]
        
        if df2 is not None:    
            # add (Single) tethered intensities, remove from dat
            tmp=dat[dat.mass.between(row2.protein_1adduct_min, row2.protein_1adduct_max)]
            row.protein_adduct=row.protein_adduct+tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
        
        # add secondary labeling intensities, remove from dat
        if site_count>1:
            # 2caps, 2reds, 1cap1red (not labeled)
            tmp=dat[(dat.mass.between(row.protein_2cap_min, row.protein_2cap_max))|(dat.mass.between(row.protein_2reductant_min, row.protein_2reductant_max))|(dat.mass.between(row.protein_1reductant_1cap_min, row.protein_1reductant_1cap_max))]
            row.protein=row.protein+tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:    
                # 2caps, 2reds, 1cap1red (not labeled)
                tmp=dat[(dat.mass.between(row2.protein_2cap_min, row2.protein_2cap_max))|(dat.mass.between(row2.protein_2reductant_min, row2.protein_2reductant_max))|(dat.mass.between(row2.protein_1reductant_1cap_min, row2.protein_1reductant_1cap_max))]
                row.protein=row.protein+tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]
            
            # adduct-cap, adduct-red (single labeled at either spot)
            tmp=dat[(dat.mass.between(row.protein_1adduct_1cap_min,row.protein_1adduct_1cap_max))|(dat.mass.between(row.protein_1reductant_1adduct_min,row.protein_1reductant_1adduct_max))]
            row.protein_adduct=row.protein_adduct+tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:    
                # adduct-cap, adduct-red (single labeled at either spot)
                tmp=dat[(dat.mass.between(row2.protein_1adduct_1cap_min,row2.protein_1adduct_1cap_max))|(dat.mass.between(row2.protein_1reductant_1adduct_min,row2.protein_1reductant_1adduct_max))]
                row.protein_adduct=row.protein_adduct+tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]

            # adduct-adduct (double labeled at max 2 spots)
            tmp=dat[(dat.mass.between(row.protein_2adduct_min,row.protein_2adduct_max))]
            row.protein_2adduct=tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:    
                # adduct-adduct (double labeled at max 2 spots)
                tmp=dat[(dat.mass.between(row2.protein_2adduct_min,row2.protein_2adduct_max))]
                row.protein_2adduct=tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]
        
        if site_count>2:    
            # 3 caps, 3 reds, 2 cap 1 red, 2 red 1 cap (no label)
            tmp=dat[(dat.mass.between(row.protein_3cap_min,row.protein_3cap_max))|(dat.mass.between(row.protein_3reductant_min, row.protein_3reductant_max))|
                    (dat.mass.between(row.protein_2reductant_1cap_min, row.protein_2reductant_1cap_max))|(dat.mass.between(row.protein_1reductant_2cap_min, row.protein_1reductant_2cap_max))]
            row.protein=row.protein+tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:    
                # 3 caps, 3 reds, 2 cap 1 red, 2 red 1 cap (no label)
                tmp=dat[(dat.mass.between(row2.protein_3cap_min,row2.protein_3cap_max))|(dat.mass.between(row2.protein_3reductant_min, row2.protein_3reductant_max))|
                        (dat.mass.between(row2.protein_2reductant_1cap_min, row2.protein_2reductant_1cap_max))|(dat.mass.between(row2.protein_1reductant_2cap_min, row2.protein_1reductant_2cap_max))]
                row.protein=row.protein+tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]
            
            # 1adduct-2cap, 1adduct-2red, 1adduct-1cap-1red (single)
            tmp=dat[(dat.mass.between(row.protein_1adduct_2cap_min, row.protein_1adduct_2cap_max))|(dat.mass.between(row.protein_2reductant_1adduct_min, row.protein_2reductant_1adduct_max))|
                    (dat.mass.between(row.protein_1adduct_1cap_1reductant_min,row.protein_1adduct_1cap_1reductant_max))]
            row.protein_adduct=row.protein_adduct+tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:
                # 1adduct-2cap, 1adduct-2red, 1adduct-1cap-1red (single)
                tmp=dat[(dat.mass.between(row2.protein_1adduct_2cap_min, row2.protein_1adduct_2cap_max))|(dat.mass.between(row2.protein_2reductant_1adduct_min, row2.protein_2reductant_1adduct_max))|
                        (dat.mass.between(row2.protein_1adduct_1cap_1reductant_min,row2.protein_1adduct_1cap_1reductant_max))]
                row.protein_adduct=row.protein_adduct+tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]
            
            # 2adduct-1cap, 2adduct-1red (double)
            tmp=dat[(dat.mass.between(row.protein_2adduct_1cap_min, row.protein_2adduct_1cap_max))|(dat.mass.between(row.protein_1reductant_2adduct_min, row.protein_1reductant_2adduct_max))]
            row.protein_2adduct=tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:
                # 2adduct-1cap, 2adduct-1red (double)
                tmp=dat[(dat.mass.between(row2.protein_2adduct_1cap_min, row2.protein_2adduct_1cap_max))|(dat.mass.between(row2.protein_1reductant_2adduct_min, row2.protein_1reductant_2adduct_max))]
                row.protein_2adduct=tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]

            # 3 adduct (triple)
            tmp=dat[(dat.mass.between(row.protein_3adduct_min, row.protein_3adduct_max))]
            row.protein_3adduct=tmp.intensity.sum()
            dat=dat[~dat.index.isin(tmp.index)]
            
            if df2 is not None:
                # 3 adduct (triple)
                tmp=dat[(dat.mass.between(row2.protein_3adduct_min, row2.protein_3adduct_max))]
                row.protein_3adduct=tmp.intensity.sum()
                dat=dat[~dat.index.isin(tmp.index)]
        
        # add max intensities if elsewhere, remove from dat (labeling at max intensity spot ignoring where it should be labeled)
        # massatmax should be the same for either protein mass, it just probably is already removed bc it's the other peak
        tmp=dat[dat.mass.between(row.massatmax-scan_range, row.massatmax+scan_range)]
        row.protein_secondary=tmp.intensity.sum()
        dat=dat[~dat.index.isin(tmp.index)]

        # add residual signal (background or other labeling)
        row.well_sum=dat.intensity.sum()
        row.well_cnt=len(dat)
        rows.append(pd.DataFrame(row).T)
        
    df=pd.concat(rows)

    # rename columns
    df=df.rename(columns={
        'string':'intensityprint',
        'protein':'intensityprotein',            # sum of protein unlabeled peak (s) - unlabled single, unlabeled cap cap, unlabeled red red, unlabeled cap red
        'protein_adduct':'intensitytethered',    # sum of protein single labeled peak (s) - single label only, single + reductant, single + cap
        'protein_secondary':'intensitysecondary',# sum of protein at maximum intensity peak regardless of where it's expected
        'protein_2adduct':'intensitydoubletethered', # sum of protein at double labeled peaks
        'protein_3adduct':'intensitytripletethered', # sum of protein at triple labeled peaks
    })

    well_sum=df.well_sum.sum()                   # sum of all residual signal across all wells
    well_cnt=df.well_cnt.sum()                   # sum of number of measured points contributing to residual signal across all wells

    # calculate total (expected) intensity of all labeled or unlabeled protein
    df['intensitytotal']=df.intensityprotein+df.intensitytethered+df.intensitydoubletethered+df.intensitytripletethered

    # calc percentages

    # real labeling at expected masses (ignoring any secondary peak) (labeled protein / all protein)
    df['percentlabeled']=(df.intensitytethered+df.intensitydoubletethered+df.intensitytripletethered)/df.intensitytotal*100

    # only filter if AssayClass is 'HighThroughputAssay' not dose response
    if assay_class == 'HTS':
        df['note']=''
        df.loc[df.percentlabeled>=100,'note']='Percent labeled > 100%'
        df.loc[df.percentlabeled<0.1,'note']='Percent labeled < 0.1%'

    # secondary mass that doesn't represent expected mass peaks (secondary peak / all protein when you include the secondary peak)
    df['percentsecondary']=df.intensitysecondary/(df.intensitytotal+df.intensitysecondary)*100

    # single labeling at expected mass
    df['percentsinglelabeled']=df.intensitytethered/(df.intensitytotal)*100

    # percent double labeled at expected masses
    if site_count>1:
        df['percentdoublelabeled']=df.intensitydoubletethered/(df.intensitytotal)*100
        df['percenttriplelabeled']=df.intensitytripletethered/(df.intensitytotal)*100

    # ignore secondary peak if its % labeling is less than real labeling
    df.loc[df.percentsecondary<df.percentlabeled, 'percentsecondary']=0
    df.loc[df.percentsecondary<df.percentlabeled, 'spectrashift']=0

    # otherwise note spectra shift (delta in mass of max intensity peak vs protein alone peak
    df.loc[df.percentsecondary>=df.percentlabeled, 'spectrashift']=df.loc[df.percentsecondary>=df.percentlabeled, 'massatmax'].astype(float)-df.loc[df.percentsecondary>=df.percentlabeled, 'protein_mass'].astype(float)

    # what percentage of peaks were binned meaningfully?
    # ratio of expected mass intensity to well total sum
    # intensity of all peaks noted / all peaks noted plus extra signal
    df['signalsignificance']=100*(df.intensitytotal+df.intensitysecondary)/(df.well_sum+df.intensitytotal+df.intensitysecondary)
    # TODO: remove properties before insertion into db? "protein", "protein_adduct", "protein_secondary", "protein_2adduct", "well_sum", "mass", "intensity"
    return df


def prepare_for_plotting(df, addcaps, barcode):
    """Prepares analyzed data for plotting by merging with platemap information.
    
    Args:
        df (pd.DataFrame): DataFrame containing analyzed results.
        addcaps (pd.DataFrame): DataFrame containing adduct and cap mass ranges.
        barcode (str): Plate barcode.
    Returns:
        pd.DataFrame: Merged DataFrame ready for plotting."""
    
    mergecols=['IP_ID', 'IW_ID', 'ICL_IC_ID', 'ICL_ID','protein_mass']
    dropcols=addcaps.columns.difference(mergecols)

    df1=df.copy()
    df1=df1.drop(columns=dropcols)
    dfs=[df1]

    masses=addcaps.protein_mass.unique().tolist()
    if len(masses)==2:
        df2=df1.copy()
        dfs.append(df2)
    for tmp_df, mass in zip(dfs, masses):
        tmp_df.protein_mass=mass
    
    dfm=pd.concat(dfs)
    dfm=dfm.merge(addcaps, how='outer')

    pl=pm.make_platemaps_from_barcode(barcode=barcode, outdir=None)
    if pl.Conc_uM.unique()[0] is None:
        pl.Conc_uM=pl.Conc_nM/1000
    pl=pl.rename(columns={'Plate':'Barcode'})
    dfm=dfm.merge(pl, how='outer', on='Well')
    dfm['percentlabeled']=dfm['percentlabeled'].astype(float)
    dfm['hmlab']=dfm.percentlabeled.round(1)
    dfm=dfm.dropna(subset='Sample')
    return dfm


colors=['gray','lightpink','palevioletred','mediumvioletred','darkmagenta']

def plot_bpi(dfm, smdc, scan_range, outdir="./trace_images", save=True):
    """Plots BPI traces with adduct and cap mass annotations.

    Args:
        dfm (pd.DataFrame): DataFrame containing merged data for plotting.
        smdc (int): Compound ID for plotting.
        outdir (str): Directory to save the plots.
    Returns:
        list: List of labels used in the plots."""
    
    df1=dfm[dfm.protein_mass==dfm.protein_mass.unique()[0]]
    dfs=[df1]

    if len(dfm.protein_mass.unique())==2:
        df2=dfm[dfm.protein_mass==dfm.protein_mass.unique()[1]]
        dfs.append(df2)
    
    smdc=int(smdc)
    xmin=0
    xmax=dfm.protein_mass.max()+3

    for well in df1[df1.ICL_IC_ID==smdc].Well:
        barcode=df1[df1.Well==well].Barcode.iloc[0]
        sample=df1[df1.Well==well].Sample.iloc[0]
        conc=df1[df1.Well==well].Conc_uM.iloc[0]
        conc_unit='uM'
        if conc is None:
            conc=df1[df1.Well==well].Conc_nM.iloc[0]
            conc_unit='nM'
        pct=df1[df1.Well==well].hmlab.iloc[0]
    
        dat=pd.read_csv(os.path.join(outdir, 'samples', f'Sample_{str(sample).zfill(3)}_data.txt'), sep='\t')
        dat.columns=['Mass','BPI']
        dat=dat.astype(float)
        dat.head(2)
        
        colors=['gray','lightpink','palevioletred','mediumvioletred','darkmagenta']
        labels=['protein only','single','double','triple']

        fig1, ax1 = plt.subplots(1, figsize=(15,3))
        sns.lineplot(data=dat, x='Mass', y='BPI', color='navy', lw=1)

        fig2, ax2 = plt.subplots(1, figsize=(15,3))
        sns.lineplot(data=dat, x='Mass', y='BPI', color='navy', lw=1)        

        for i, df in enumerate(dfs):
            barcode=df[df.Sample==sample].Barcode.iloc[0]
            adduct_mass=df[df.Sample==sample].adduct_mass.astype(float).iloc[0] 
            protein_mass=df[df.Sample==sample].protein_mass.astype(float).iloc[0]
            xmin=protein_mass-adduct_mass-5*scan_range
            xmax=protein_mass+3*adduct_mass+5*scan_range
            for mult in [0,1,2,3]:
                ax1.axvspan(protein_mass+mult*adduct_mass-scan_range, protein_mass+mult*adduct_mass+scan_range, color=colors[mult], label=labels[mult], zorder=0, alpha=0.5)

            ax2.axvspan(protein_mass-scan_range, protein_mass+scan_range, color='gray', label='protein_only', zorder=0, alpha=0.5)
            protnums=df.loc[df['Sample']==sample,[x for x in df.columns if 'min' in x]].T
            protnums=protnums.sort_values(by=protnums.columns[0])
            protnums=protnums[~protnums[protnums.columns[0]].isna()]
            cols=protnums.index.tolist()
            cols=[x.replace('_min','') for x in cols]
            cols.remove('free_protein')
            
            pal=sns.color_palette("viridis", len(cols))
            
            for j, lab in enumerate(cols):
                ax2.axvspan(df[df.Sample==sample][f'{lab}_min'].astype(float).iloc[0],
                        df[df.Sample==sample][f'{lab}_max'].astype(float).iloc[0],
                        label=lab, color=pal[j], zorder=-j, alpha=0.5)

        
        for ax in [ax1, ax2]:
            ax.set_xlim(xmin, xmax)
            ax.legend()
            handles, labels=ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())


        if save:
            os.makedirs(os.path.join(outdir, 'trace_images'), exist_ok=True)

        if not np.isnan(conc):
            ax1.set_title(f'{barcode} {well} {smdc} {conc} {conc_unit} labeling assessment ({pct}% labeled)');
            if save:
                fig1.savefig(os.path.join(outdir, 'trace_images', f'{barcode}_{well}_{smdc}_{conc}_{conc_unit}_mults.png'))
                plt.close(fig1)
        else:
            ax1.set_title(f'{barcode} {well} {smdc} labeling assessment ({pct}% labeled)');
            if save:
                fig1.savefig(os.path.join(outdir, 'trace_images', f'{barcode}_{well}_{smdc}_mults.png'))
                plt.close(fig1)
            

        if not np.isnan(conc):
            ax2.set_title(f'{barcode} {well} {smdc} {conc} {conc_unit} labeling assessment ({pct}% labeled)');
            if save:
                fig2.savefig(os.path.join(outdir, 'trace_images', f'{barcode}_{well}_{smdc}_{conc}_{conc_unit}_full.png'))
                plt.close(fig2)
        else:
            ax2.set_title(f'{barcode} {well} {smdc} labeling assessment ({pct}% labeled)');
            if save:
                fig2.savefig(os.path.join(outdir, 'trace_images', f'{barcode}_{well}_{smdc}_full.png'))
                plt.close(fig2)