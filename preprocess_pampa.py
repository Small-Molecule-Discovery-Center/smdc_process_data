# preprocess pampa explorer data for upload to SMDC HiTS
# created 20230317 by @paulsonak amanda.paulson@ucsf.edu

import pandas as pd
import numpy as np
import os

def preprocess_pampa(infile='', source_plate="SMDC123", donor_well_conc=50):
    destination_plate=infile.split('/')[-1].replace('_results.xlsx','').replace('_results.xls','').replace('_Results.xlsx','').replace('_Results.xls', '').replace('_Results_processed.csv','')
    print('Destination plate:', destination_plate)
    
    # path to what it will be saved as
    outfile=os.path.join("./", destination_plate+'_processed.csv')
    
    # use pandas to read file
    try:
      pampa=pd.read_excel(infile, header=1)
    except:
      pampa=pd.read_csv(infile, index_col=0)
    # get ride of empty rows
    pampa=pampa[~pampa.pI.isna()]
    # reset the index column for ease of next steps
    pampa=pampa.reset_index(drop=True)
    # fill the rows from Sample column that are empty with 'nolab'
    pampa.Sample=pampa.Sample.fillna('nolab')
    # gut check - how many samples were measured - including 'nolab'?
    print('Unique samples pre:', pampa.Sample.nunique())
    
    
    # loop through each row and fill it in
    
    # initialize sample variable as None
    sample=None
    # go through each row
    for i, row in pampa.iterrows():
        # if Sample column is NOT nolab, ie has a real label,
        # set the sample variable to be what Sample is for that row
        if row.Sample!='nolab':
            sample=row.Sample
        # otherwise, set the Sample column for that row to what the sample variable is
        # in this case this will replace 'nolab' with the most recent previous true sample label
        else:
            pampa.loc[i,'Sample']=sample
    
    # gutcheck - is the number of samples 1 less than before, since all the 'nolabs' are gone?
    print('Unique samples post:', pampa.Sample.nunique())
    
    # rename other columns
    pampa=pampa.rename(columns={
        'Sample':'SMDC_ID',
        'Pe Well':'Destination well',
    })

    # add additional plate map columns
    i=0
    if i==0:
      try:
        pampa[['SMDC_ID','Lot']]=pampa.SMDC_ID.astype(str).str.split('-', expand=True)
      except:
        pampa['Lot']=np.nan
      pampa.Lot=pampa.Lot.astype(float)
      i=1
    
    pampa['Source plate']=source_plate
    pampa['Source well']=np.nan
    pampa['Destination plate']=destination_plate
    pampa['[compound] uM']=donor_well_conc

    # rename controls to SMDC_IDs and add lots
    ctrl_dict={
        'Theophylline':254802,
        "Verapamil":131810,
        "Corticosterone":1076478,
        'DMSO':np.nan
    }
    ctrl_lot_dict={
        'Theophylline':2,
        "Verapamil":13,
        "Corticosterone":2,
    }
    
    pampa.loc[pampa.SMDC_ID=='Theophylline', '[compound] uM']=250
    
    for smdc in ctrl_lot_dict:
      pampa.loc[pampa.SMDC_ID==smdc, 'Lot']=ctrl_lot_dict[smdc]
    pampa=pampa.replace({'SMDC_ID':ctrl_dict})
    
    pampa.SMDC_ID=pampa.SMDC_ID.astype(float)
    
    # fix Pe and notes column
    pampa.loc[pampa['P(10-6cm/s)']=='equilibrated','BCS code']='HIGH_EQ'
    pampa['P(10-6cm/s)']=pampa['P(10-6cm/s)'].replace('undetected','').replace('equilibrated','')

    
    # save new data as csv
    pampa.to_csv(outfile)
    print(f"Data saved to {outfile}")