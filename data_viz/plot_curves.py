### fit and plot curves for IC50s ########################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import math
from scipy.optimize import curve_fit 
import warnings
###################################################################################################
colors=[(0.212395, 0.359683, 0.55171), (0.153364, 0.497, 0.557724), (0.122312, 0.633153, 0.530398), (0.288921, 0.758394, 0.428426), (0.626579, 0.854645, 0.223353),(0.275191, 0.194905, 0.496005)]
sns.set_palette(colors)
pal=sns.color_palette(colors)

def plot_pct_inhib_curves(compound, raw, curve, targs=['PARP1','PARP2'], save_fig=False, graphcol='intensity_inhibition', fitcrit='R^2', sharey=False, save_lab='HiTS', save_dir='./curve_fits'):
    ''' Plot and save IC50 curves on top of raw data for visualization across various targets.
    Raw data frame should at minimum have columns
    [Compound, Target, <response column graphcol>, logconc, Plate, Well_Type]. 'solubility' optional.
    Curve data frame should have columns
    [Compound, Target, 'bottom','top','LogIC50','hill','RMSE',<fitcrit col>,'IC50','relation','pIC50','hitc','coff','notes']
    '''
    curvecols=['bottom','top','LogIC50','hill','RMSE',fitcrit,'IC50','pIC50','hitc','coff']
    sns.set_context('poster')
    rows=int(np.ceil(len(targs)/4))
    cols=int(np.ceil(len(targs)*2/rows))
    realcurvecols=['Compound']
    for col in curvecols:
        if col in curve.columns:
            curve.loc[:,col]=curve[col].astype(float)
            realcurvecols.append(col)
        # except:
        #     # warnings.warn(f'{col} column does not exist or was unable to be cast as float')
        #     pass
    
    fig, ax = plt.subplots(rows,cols,figsize=(cols*6,rows*10),sharey=sharey,gridspec_kw={'width_ratios': [5,1]*int(np.ceil(cols/2))})
    ax=ax.ravel()
    i=0
    minmin=0
    maxmax=0
    for targ in targs:
        # get curve data
        ppc=curve[(curve.Compound==compound) & (curve.Target==targ)].copy()
        note=ppc[['notes']]
        ppc=ppc[realcurvecols]
        if len(ppc)>1:
            ppc.IC50=ppc.IC50.astype(float)
            ppc=ppc.groupby('Compound').mean()
        if len(ppc)>0:
            for col in ['top','bottom','hill','pIC50',fitcrit,'IC50',]:
                ppc[col]=ppc[col].astype(float)
            if 'relation' in ppc.columns:
                ppc.relation=ppc.relation.replace(np.nan,'=')
            else:
                ppc['relation']='='
            if 'coff' in ppc.columns:
                coff = ppc.iloc[0].coff
            else:
                coff=None
        else:
            coff=None
        note=list(set(note.notes.dropna()))
        ppc['notes']='\n'.join(note)
        
        # get raw data
        ppr=raw[(raw.Compound==compound) & (raw.Target==targ)].copy()
        if len(ppr)==0:
            continue
        if 'RealPlate' in ppr.columns:
            plateNegCtrl=raw[(raw.RealPlate.isin(ppr.RealPlate.unique()))&(raw.Well_Type=='negativeControl')&(raw.Target==targ)]
        else:
            plateNegCtrl=raw[(raw.Plate.isin(ppr.Plate.unique()))&(raw.Well_Type=='negativeControl')&(raw.Target==targ)]
        if coff is None:
            vals=plateNegCtrl[graphcol].values
            coff=3*np.median(np.absolute(vals - np.median(vals))) # 3*BMAD
        
            
        # get axes bounds       
        logconcmin=ppr.logconc.min()-0.5
        logconcmax=ppr.logconc.max()+0.5
        pctmin = min(ppr[graphcol].min(), plateNegCtrl[graphcol].min(), coff)
        pctmax = max(ppr[graphcol].max(), plateNegCtrl[graphcol].max(), coff)
        pctmin = pctmin-0.05*pctmax
        pctmax = pctmax+0.05*pctmax
        
        # create curve if it's a hit and curve fit succeeded
        if len(ppc)>0:
            if (ppc.hitc.iloc[0]==1) and ('succeeded' in ppc.notes.iloc[0]):
                # try:
                x = np.linspace(logconcmin, logconcmax, 1000)
                y = ppc.bottom.iloc[0] + ((ppc.top.iloc[0] - ppc.bottom.iloc[0])/(1+10**((-ppc.pIC50.iloc[0] - x)*float(ppc.hill.iloc[0]))))
                ax[i].plot(x,y,c=pal[-1])
                # except:
                #     pass
                pctmin=min(pctmin, ppc.bottom.iloc[0])
                pctmax=max(pctmax, ppc.top.iloc[0])
                pctmin = pctmin-0.05*pctmax
                pctmax = pctmax+0.05*pctmax
           
        # create scatterplot whether it's a hit or not
        if 'RealPlate' in ppr.columns:
            sns.scatterplot(data=ppr, x='logconc', y=graphcol, hue='RealPlate', ax=ax[i], legend=False)
        else:
            ax[i].scatter(x=ppr.logconc.values, y=ppr[graphcol].values, s=50)
        if not isinstance(compound,str):
            try:
                compound=int(compound)
            except:
                pass
        if 'solubility' in ppc.columns:
            ax[i].set_title(f"{compound} {targ} IC50{ppc.relation.iloc[0]}{round(ppc.IC50.iloc[0], 4)} uM\n{ppc.notes.iloc[0]}\nSoluble: {ppc['solubility'].iloc[0]}; \nHITC: {int(ppc.hitc.iloc[0])};{fitcrit}={round(ppc[fitcrit].iloc[0], 4)}")
        elif len(ppc)==0:
            ax[i].set_title(f"{compound} {targ}\n No Fits")
        else:
            ax[i].set_title(f"{compound} {targ} IC50{ppc.relation.iloc[0]}{round(ppc.IC50.iloc[0], 4)} uM\n{ppc.notes.iloc[0]}\nHITC: {int(ppc.hitc.iloc[0])}; {fitcrit}={round(ppc[fitcrit].iloc[0], 4)}")
        ax[i].set_xlabel('log(Concentration in M)')
        ax[i].set_xlim(logconcmin, logconcmax)
        ax[i].set_ylim(pctmin, pctmax)
        ax[i].hlines(coff, xmin=logconcmin, xmax=logconcmax, color=pal[-2], zorder=0)
        
        if 'RealPlate' in plateNegCtrl.columns:
            sns.histplot(y=graphcol, data=plateNegCtrl, bins=np.arange(pctmin,pctmax,((pctmax-pctmin)/22)), hue='RealPlate', ax=ax[i+1], legend=False)
        else:
            ax[i+1].hist(plateNegCtrl[graphcol], orientation='horizontal', range=(pctmin,pctmax), bins=22)
        ax[i+1].set_ylim(pctmin, pctmax)
        ax[i+1].set_xlabel('negCtrlDist')
        i+=2
        
        # if sharey is true, use maxmax
        if pctmin < minmin:
            minmin=pctmin
        if pctmax > maxmax:
            maxmax=pctmax
        
    # finish creating plot
    if sharey:
        ax[0].set_ylim(minmin, maxmax)
    ax[0].set_ylabel(graphcol)
    fig.patch.set_facecolor('white')
    fig.tight_layout()
    
    # save plot if save
    if save_fig:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, f'{compound}_{save_lab}_curvefit.png'))
        plt.close()

def plot_sc_data(compound, raw, curve, targs=['PARP1','PARP2'], save_fig=False, graphcol='intensity_inhibition', save_lab='HiTS', save_dir='./curve_fits'):
    ''' Plot and save IC50 curves on top of raw data for visualization across various targets.
    Raw data frame should at minimum have columns
    [Compound, Target, <response column graphcol>, logconc, Plate, Well_Type]. 'solubility' optional.
    Curve data frame should have columns
    [Compound, Target, 'bottom','top','LogIC50','hill','RMSE',<fitcrit col>,'IC50','pIC50','hitc','coff','notes']
    '''
    sns.set_context('poster')
    rows=int(np.ceil(len(targs)/4))
    cols=int(np.ceil(len(targs)*2/rows))
    fig, ax = plt.subplots(rows,cols,figsize=(cols*6,rows*10),sharey=True,gridspec_kw={'width_ratios': [5,1]*int(np.ceil(cols/2))})
    ax=ax.ravel()
    i=0

    # minmax=raw[raw.Compound==compound]
    # if graphcol!='fold_change':
    #     pctmin = minmax[graphcol].min()-5
    #     pctmax = minmax[graphcol].max()+5
    #     pctmin = min(pctmin, -105)
    #     pctmax = max(pctmax, 105)
    # else:
    #     pctmin = minmax[graphcol].min()-0.5
    #     pctmax = minmax[graphcol].max()+0.5

    for targ in targs:
        # get curve data
        ppc=curve[(curve.Compound==compound) & (curve.Target==targ)].copy()
        # ppc=ppc[['Compound', 'Target', 'bottom','top','LogIC50','hill','RMSE',<fitcrit col>,'IC50','pIC50','hitc','coff','notes'
        if len(ppc)>1:
            note=list(set(ppc.notes.tolist()))[0]
            ppc.IC50=ppc.IC50.astype(float)
            note=set(ppc.notes.tolist())
            ppc=ppc.groupby('Compound').mean()
            ppc['notes']=str(note)
        coff = ppc.coff.iloc[0]

        # get raw data
        ppr=raw[(raw.Compound==compound) & (raw.Target==targ)].copy()
        if len(ppr)==0:
            continue
        plateNegCtrl=raw[(raw.Plate.isin(ppr.Plate.unique()))&(raw.Well_Type=='negativeControl')&(raw.Target==targ)]
        
        # get axes bounds       
        logconcmin=ppr.logconc.min()-0.5
        logconcmax=ppr.logconc.max()+0.5
        
        pctmin = min(ppr[graphcol].min(), plateNegCtrl[graphcol].min())
        pctmax = max(ppr[graphcol].max(), plateNegCtrl[graphcol].max())
        pctmin = pctmin-0.05*pctmax
        pctmax = pctmax+0.05*pctmax
        
        # create curve if it's a hit and curve fit succeeded
        if (ppc.hitc.iloc[0]==1) and ('succeeded' in ppc.notes.iloc[0]):
            try:
                x = np.linspace(logconcmin, logconcmax, 1000)
                y = ppc.bottom.iloc[0] + ((ppc.top.iloc[0] - ppc.bottom.iloc[0])/(1+10**((-ppc.pIC50.iloc[0] - x)*float(ppc.hill.iloc[0]))))
                ax[i].plot(x,y,c='cornflowerblue')
            except:
                pass

        # create scatterplot whether it's a hit or not
        ax[i].scatter(x=ppr.logconc.values, y=ppr[graphcol].values, s=50)
        if 'solubility' in ppc.columns:
            ax[i].set_title(f"{int(compound)} {targ}\nIC50={round(ppc.IC50.iloc[0], 4)} uM HITC: {int(ppc.hitc.iloc[0])}; {ppc.notes.iloc[0]}\nSoluble: {ppc['solubility'].iloc[0]}; RMSE={round(ppc['RMSE'].iloc[0], 4)}")
        else:
            ax[i].set_title(f"{int(compound)} {targ}\nIC50={round(ppc.IC50.iloc[0], 4)} uM HITC: {int(ppc.hitc.iloc[0])}; {ppc.notes.iloc[0]}\nRMSE={round(ppc['RMSE'].iloc[0], 4)}")
        ax[i].set_xlabel('log(Concentration in M)')
        ax[i].set_xlim(logconcmin, logconcmax)

        ax[i].hlines(coff, xmin=logconcmin, xmax=logconcmax, color='orange', zorder=0)
        ax[i+1].hist(plateNegCtrl[graphcol], orientation='horizontal', range=(pctmin,pctmax), bins=22)
        ax[i+1].set_xlabel('negCtrlDist')
        i+=2

    # finish creating plot
    ax[0].set_ylim(pctmin, pctmax)
    ax[0].set_ylabel(graphcol)
    fig.patch.set_facecolor('white')
    fig.tight_layout()

    # save plot if save
    if save_fig:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, f'{compound}_{save_lab}_curvefit.png'))
        plt.close()
        
        
###################################################################################################
# from stack overflow
def ModelAndScatterPlot(xData, yData, color='blue', ax=None, label=None, init=[100, 0, -5.5, 1],
                       bounds=([80,-20,-np.inf,-np.inf], [120, 20, np.inf,np.inf])):
    if ax is None:
        f = plt.figure(figsize=(5,5))
        ax = f.add_subplot(111)

    # first the raw data as a scatter plot
    ax.plot(xData, yData,  'D', color=colors[1])

    fittedParameters, pcov = curve_fit(hillfunc, xData, yData, init, bounds=bounds)
    
    # create data for the fitted equation plot
    xModel = np.linspace(min(xData), max(xData), num=100)
    yModel = hillfunc(xModel, *fittedParameters)

    # now the model as a line plot
    ax.plot(xModel, yModel, color=colors[1], label=label)

    # ax.set_xlabel('X Data') # X axis data label
    # ax.set_ylabel('Y Data') # Y axis data label
    # axes.set_ylim(-20,105)
    return ax

###################################################################################################
    
    
def hillfunc(x, Bottom, Top, LogIC50, HillSlope): 
    return (Bottom + (Top-Bottom)/(1+10**((LogIC50-x)*HillSlope)))


def fit_scipy_curves(cmpds, raw, graphcol, targs=['PARP1','PARP2'], init=np.array([100, -20, -5.5, 1]), 
                          bounds=([80,-20,-np.inf,-np.inf], [120, 20, np.inf,np.inf]), bcomm=''):
    ''' fit IC50 curves to raw concentration-response data. Raw data frame should at minimum have columns
    [Compound, Target, <response column graphcol>, logconc]
    '''
    allparams = {}

    for cmpd in cmpds:
        for i, targ in enumerate(targs):
            ppr=raw[(raw.Compound==cmpd) & (raw.Target==targ)].copy()
            xData=ppr.logconc
            yData = ppr[graphcol].values

            # these are the same as the scipy defaults
            # initialParameters = np.array([init[top], init[bottom], np.log10(float(init[ic50])/1000000), init[hill]])
            initialParameters = init

            # do not print unnecessary warnings during curve_fit()
            warnings.filterwarnings("ignore")

            try:
                # curve fit the test data
                fittedParameters, pcov = curve_fit(hillfunc, xData, yData, initialParameters, bounds=bounds)

                modelPredictions = hillfunc(xData, *fittedParameters) 

                absError = modelPredictions - yData

                SE = np.square(absError) # squared errors
                MSE = np.mean(SE) # mean squared errors
                RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
                Rsquared = 1.0 - (np.var(absError) / np.var(yData))
                IC50 = 10**fittedParameters[2]*1000000
                pIC50= -np.log10(IC50/1000000)

                fittedParameters=list(fittedParameters)
                fittedParameters.extend([RMSE, Rsquared, IC50, pIC50])
                fittedParameters.append(f'Fit succeeded with SciPy{bcomm}')
                allparams[targ+"_"+str(cmpd)]=fittedParameters

            except:
                # print(f'{targ} {cmpd} fit failed')
                allparams[targ+"_"+str(cmpd)]=[np.nan]*8+[f'Fit FAILED with SciPy{bcomm}']

    df=pd.DataFrame.from_dict(allparams).T
    df.columns=['top','bottom','LogIC50','hill','RMSE','R^2','IC50','pIC50','notes']
    # df.columns=['SciPy_'+col for col in df.columns]
    df.hill=df.hill*-1
    df=df.reset_index()
    df=df.rename(columns={'index':'Compound'})
    df[['Target','Compound']]=df.Compound.str.split('_', expand=True)
    # df.Compound=df.Compound.astype(int)
    df['relation']=np.nan
    df['hitc']=[1 if 'succeeded' in x else 0 for x in df.notes]

    return df


###################################################################################################


def tcpl_mc3_to_raw(mc3, graphcol='percent_inhibition', fix_raw=False):
    """Prep TCPL output to be used in plotting functions.
    
    Args:
        mc3: tcpl MC3 df to be prepped
    
    Returns:
        raw: modified df
    """
    raw=mc3.copy()
    raw.columns=['Compound', 'chid', 'casn', 'chnm', 'dsstox_substance_id', 'code', 'aeid',
       'Target', 'm0id', 'm1id', 'm2id', 'm3id', 'logc', graphcol, 'cndx', 'Well_Type',
       'Plate', 'rowi', 'coli', 'repi', 'resp_unit', 'conc_unit']
    raw['Conc_uM']=10**raw.logc
    raw['logconc']=np.log10(raw.Conc_uM/1000000)
    raw.Well_Type=raw.Well_Type.str.replace('n','negativeControl')
    nconcdict={}
    for plate in raw.Plate.unique():
        nconcdict[plate]=raw[raw.Plate==plate].cndx.max()
    raw['nconc']=raw.Plate.map(nconcdict)
    nplatedict={}
    for compound in raw.Compound.unique():
        nplatedict[compound]=raw[raw.Compound==compound].Plate.nunique()
    raw['nplate']=raw.Compound.map(nplatedict)
    raw['RealPlate']=raw.Plate
    raw.Plate = raw.Target
    if fix_raw:
        for targ in raw.Target.unique():
            for compound in raw.Compound.unique():
                if compound in [207,210]:
                    continue
            # print(compound, raw[raw.Compound==compound].nconc.unique().tolist())
                if (raw[(raw.Compound==compound)&(raw.Target==targ)].nconc.unique().tolist()==[3,9]):
                    raw=raw.drop(raw[(raw.Compound==compound)&(raw.Target==targ)&(raw.nconc==3)].index)
    return raw
    
def tcpl_mc5_to_curve(mc5):
    """Prep TCPL output to be used in plotting functions.
    
    Args:
        mc5: tcpl MC5 df to be prepped
    
    Returns:
        curve: modified df
    """
    curve=mc5.copy()
    curve['AC50_uM']=10**curve.modl_ga
    curve['pIC50']=-np.log10(curve.AC50_uM/1000000)
    curve=curve[['m5id', 'flag', 'spid', 'chid', 'casn', 'aenm', 'm4id', 'bmad', 'resp_max', 'resp_min',
       'max_mean', 'max_mean_conc', 'max_med', 'max_med_conc', 'logc_max',
       'logc_min', 'nconc', 'npts', 'nrep', 'nmed_gtbl', 'hitc', 'modl', 'fitc', 'coff', 'actp', 'modl_er',
       'modl_tp', 'modl_ga', 'modl_gw', 'modl_la', 'modl_lw', 'modl_rmse',
       'modl_prob', 'modl_acc', 'modl_acb', 'modl_ac10', 'resp_unit',
       'conc_unit', 'AC50_uM', 'pIC50']]
    curve.columns=['m5id', 'flag', 'Compound', 'chid', 'casn', 'Target', 'm4id', 'bmad', 'resp_max', 'resp_min',
           'max_mean', 'max_mean_conc', 'max_med', 'max_med_conc', 'logc_max',
           'logc_min', 'nconc', 'npts', 'nrep', 'nmed_gtbl', 'hitc', 'modl', 'fitc', 'coff', 'actp', 'modl_er',
           'top', 'modl_ga', 'hill', 'modl_la', 'modl_lw', 'RMSE',
           'modl_prob', 'modl_acc', 'modl_acb', 'modl_ac10', 'resp_unit',
           'conc_unit', 'IC50', 'pIC50']
    curve['bottom']=0
    curve['notes']=['TCPL fit succeeded' if x==1 else np.nan for x in curve.hitc]
    return curve