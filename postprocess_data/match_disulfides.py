# postprocess disulfide tethering data for upload to AVIDD CDD Vault
# created 20240226 by @paulsonak amanda.paulson@ucsf.edu

import pandas as pd
import datetime

matches=pd.read_csv('/content/smdc_process_data/data/SMDC_disulfides_AVIDD.csv', index_col=0)

def match_disulfides(infile=''):   
    print("Matching AVIDD ID's...")
    # path to what it will be saved as
    outfile=(f'./{str(datetime.date.today())}_disulfide_for_avidd.csv')
    
    # use pandas to read file
    data=pd.read_csv(infile, header=1)
    # filter empty rows
    data=data[~data.Compound.isna()]

    # update match plate list with data plate names for merging
    plate_dict={}
    for plate in matches.plate.unique():
        plate_dict[plate]=[x for x in data.Plate.unique() if x.startswith(plate)][0]
    matches.plate=matches.plate.map(plate_dict)

    # merge data with matched ID's by plate/well
    newdata=data.merge(matches[['Molecule-Batch ID','plate','well']], left_on=['Plate','Well'], right_on=['plate','well'])

    # save new data as csv
    newdata[['Plate','Well','Molecule-Batch ID', 'percentlabeled']].to_csv(outfile, index=False)

    print(f'Saved matched data to {outfile}')