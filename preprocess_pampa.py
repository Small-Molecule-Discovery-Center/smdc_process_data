# preprocess pampa explorer data for upload to SMDC HiTS
# created 20230317 by @paulsonak amanda.paulson@ucsf.edu

def preprocess_pampa(infile=''):
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
    
    # look at new data
    pampa.head(12)
    
    # save new data as csv
    pampa.to_csv(outfile)
    print(f"Data saved to {outfile}")