#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A. Paulson
Created on Thu Oct 17 15:28:23 2024
LCMS data analysis functions
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def read_lcms_results(infile, kind="spectrum"):
    """Read LCMS results from Excel file for given kind of data (spectrum or chromatogram)."""

    df=pd.read_excel(infile, sheet_name=None)
    dfs=[]
    for sheet in df.keys():
        if sheet=="Summary":
            continue
        tmp=pd.read_excel(infile, sheet_name=sheet, header=1)
        tmp.columns=[x.replace("\n","_").replace("(","").replace(")","").replace(" ","_").replace("/","").replace("%","pct").replace(":","",).lower() for x in tmp.columns]
        if kind=="spectrum":
            tmp=tmp.drop(0)
        if kind=="chromatogram":
            tmp=tmp.dropna(how='all')
            tmp=_process_chromatogram_df(tmp)
        tmp['sample_id']=sheet
        # ensure 'sample_id' is the first column
        cols = tmp.columns.tolist()
        if 'sample_id' in cols:
            cols = ['sample_id'] + [c for c in cols if c != 'sample_id']
            tmp = tmp[cols]
        dfs.append(tmp)
    df=pd.concat(dfs, axis=0).reset_index(drop=True)
    return df


def _process_chromatogram_df(df2):
    """Process LCMS chromatogram dataframe to segment by UV channels and sum absorption."""
    
    seps = ['TUVC','254','220','280','320']
    positions = df2['time_peak_maximum_msminutes'].str.contains('|'.join(seps)).dropna().index.tolist()

    segments = []
    for i, p in enumerate(positions):
        start = p
        end = positions[i+1] if i+1 < len(positions) else len(df2)
        label = df2.loc[p]['time_peak_maximum_msminutes']
        try:
            label=[x for x in seps if x in label]
            label=label[0]
        except Exception:
            pass
        tmp=df2.loc[start+1:end-1].copy()#.reset_index(drop=True)
        # tmp=tmp.set_index('time_peak_maximum_msminutes')
        # tmp.columns=[f"{label}_{col}" for col in tmp.columns]
        # segments[label] = tmp
        tmp['wavelength'] = label
        segments.append(tmp)
    df2=pd.concat(segments, axis=0)
    # df2=segments[seps[0]]
    # for sep in seps[1:]:
    #     if sep in segments:
    #         df2 = df2.join(segments[sep], how='outer')
    # df2['sum_maximum_absorption_mau'] = df2[[col for col in df2.columns if 'maximum_absorption_mau' in col]].sum(axis=1)
    return df2


def plot_lcms_spectrum(info_df, data_df, sample_id="04036CB-INO2-P12_1", scan_range=5, save_dir=None):
    """Plot LCMS spectrum with expected mass and common adducts highlighted."""
    
    smdc_id=info_df[info_df.sample_id==sample_id].SMDC_ID.astype(int).iloc[0]
    compound_mass=info_df[info_df.sample_id==sample_id].MolWt.iloc[0]
    fig, ax = plt.subplots(1, figsize=(15,4))
    sns.lineplot(data=data_df[data_df.sample_id==sample_id], x="mass_peak_maximum_mz", y="maximum_intensity_cs", ax=ax)
    ax.axvline(compound_mass, color='red', linestyle='--', label=smdc_id)
    ax.axvspan(compound_mass-scan_range, compound_mass+scan_range, color='yellow', alpha=0.3, label=f"Expected mass ({int(compound_mass)}±{scan_range})")
    ax.axvspan(130-scan_range, 130+scan_range, color='gray', alpha=0.3, label="ACN (130±5)")
    ax.axvspan(157-scan_range, 157+scan_range, color='gray', alpha=0.3, label="DMSO dimer (157±5)")
    # oxygen adduct
    o_adduct=compound_mass+16
    # ax.axvspan(o_adduct-scan_range, o_adduct+scan_range, color='cyan', alpha=0.3, label="Oxygen adduct (+16)")
    h2o_adduct=compound_mass+18
    # ax.axvspan(h2o_adduct-scan_range, h2o_adduct+scan_range, color='orange', alpha=0.3, label="H2O/NH4+ adduct (+18)")
    na_adduct=compound_mass+23
    ax.axvspan(na_adduct-scan_range, na_adduct+scan_range, color='blue', alpha=0.3, label="Na+ adduct (+23)")
    k_adduct=compound_mass+39
    # ax.axvspan(k_adduct-scan_range, k_adduct+scan_range, color='purple', alpha=0.3, label="K+ adduct (+39)")
    dmso_adduct=compound_mass+79
    ax.axvspan(dmso_adduct-scan_range, dmso_adduct+scan_range, color='green', alpha=0.3, label="DMSO adduct (+79)")
    ax.set_title(f"LCMS Spectrum for {sample_id}: {smdc_id} (MolWt: {round(compound_mass, 3)})")
    ax.set_xlabel("m/z")
    ax.set_ylabel("Max Intensity (CS)")
    ax.set_xlim(100,600)
    ax.legend()
    if save_dir is not None:
        plt.savefig(f"{save_dir}/{sample_id}_lcms_spectrum.png", dpi=300)
        plt.close()


def plot_uv_absorption_nm(sample_id="04036CB-INO2-P12_1",data_df=None, amp_col='maximum_absorption_mau', mu_col='wavelength_peak_maximum_nm', fwhm_col='peak_resolution', save_dir=None):
    """Plot individual Gaussians and their sum for UV absorption wavelengths.
    
    Parameters:
    sample_id: str, identifier for the sample
    data_df: pd.DataFrame, dataframe containing amplitude, mu, and fwhm columns
    amp_col: str, column name for amplitudes
    mu_col: str, column name for means (mu)
    fwhm_col: str, column name for full width at half maximum (FWHM)
"""
    
    x = np.linspace(start=200, stop=400, num=2000)
    gaussians, mus=_process_gaussians(data_df=data_df, amp_col=amp_col, mu_col=mu_col, fwhm_col=fwhm_col, x=x)

    # Plot individual curves
    plt.figure(figsize=(10, 5))
    for i, g in enumerate(gaussians):
        plt.plot(x, g, '--', label=f'Wavelength {mus[i]:.1f} nm')

    # Envelope (max at each x)
    envelope = np.max(gaussians, axis=0)

    sum_gaussians = np.sum(gaussians, axis=0)
    plt.plot(x, sum_gaussians, 'g-', label='Sum of Gaussians')

    # plt.plot(x, envelope, 'r-', label='Envelope (Curtain Shape)')
    plt.xlabel('x')
    plt.ylabel('Absorption')
    plt.xlim(200, 400)
    plt.legend()
    plt.title(f'Absorption Spectrum {sample_id}')
    plt.tight_layout()
    if save_dir is not None:
        plt.savefig(f"{save_dir}/{sample_id}_uv_absorption_nm.png", dpi=300)
        plt.close()
    else:
        plt.show()


def _process_gaussians(data_df, amp_col='maximum_absorption_mau', mu_col='time_peak_maximum_msminutes', fwhm_col='peak_resolution', x=None):
    amplitudes = data_df[amp_col].astype(float).tolist()
    mus = data_df[mu_col].astype(float).tolist()
    fwhms = data_df[fwhm_col].astype(float).tolist()

    # Convert FWHM to standard deviation (sigma)
    sigmas = [fwhm / (2 * np.sqrt(2 * np.log(2))) for fwhm in fwhms]
    gaussians = []

    for amp, mu, sigma in zip(amplitudes, mus, sigmas):
        gaussians.append(amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)))

    return gaussians, mus


# def plot_uv_absorption_time(sample_id="04036CB-INO2-P12_1",data_df=None, amp_col='maximum_absorption_mau', mu_col='time_peak_maximum_msminutes', fwhm_col='peak_resolution', wavelength_col='wavelength', save_dir=None):
#     """Plot individual Gaussians and their sum for UV absorption wavelengths.
    
#     Parameters:
#     sample_id: str, identifier for the sample
#     data_df: pd.DataFrame, dataframe containing amplitude, mu, and fwhm columns
#     amp_col: str, column name for amplitudes
#     mu_col: str, column name for means (mu)
#     fwhm_col: str, column name for full width at half maximum (FWHM)
# """
#     gaussians=[]
#     x = np.linspace(start=0, stop=9.2, num=2000)
#     mus_all=[]
#     for wavelength in data_df[wavelength_col].unique():
#         subset_df = data_df[data_df[wavelength_col] == wavelength]
#         sw_gaussians, mus=_process_gaussians(data_df=subset_df, amp_col=amp_col, mu_col=mu_col, fwhm_col=fwhm_col, x=x)
#         # Envelope (max at each x)
#         envelope = np.max(sw_gaussians, axis=0)
#         # Sum
#         sum_gaussians = np.sum(sw_gaussians, axis=0)
#         # gaussians.append(envelope)
#         gaussians.extend(sw_gaussians)
#         mus_all.extend(mus)
#     # Plot individual curves
#     plt.figure(figsize=(10, 5))
#     for i, g in enumerate(gaussians):
#         plt.plot(g, '--', label=f'Wavelength {mus_all[i]} nm')


#     # sum_gaussians = np.sum(gaussians, axis=0)
#     # plt.plot(x, sum_gaussians, 'g-', label='Sum of Gaussians')
#     envelope = np.max(gaussians, axis=0)
#     plt.plot(x, envelope, 'r-', label='Envelope (Curtain Shape)')

#     # plt.plot(x, envelope, 'r-', label='Envelope (Curtain Shape)')
#     plt.xlabel('x')
#     plt.ylabel('Absorption')
#     plt.xlim(0, 9)
#     # plt.yscale('log')
#     plt.legend()
#     plt.title(f'Absorption Spectrum {sample_id}')
#     plt.tight_layout()
#     if save_dir is not None:
#         plt.savefig(f"{save_dir}/{sample_id}_uv_absorption_time.png", dpi=300)
#         plt.close()
#     else:
#         plt.show()

def plot_uv_vs_time(data_df, sample_id="04036CB-INO2-P12_1", percent=False,
                       amp_col='maximum_absorption_mau',
                       mu_col='time_peak_maximum_msminutes',
                       wavelength_col='wavelength', save_dir=None, ax=None):
    """Lollipop plot: points at (mu_col, amp_col) with vertical lines from y=0 to y=amp_col,
    colored by wavelength_col.
    """
    # Prepare dataframe
    df = data_df[[mu_col, amp_col, wavelength_col]].copy()
    
    wavelengths = sorted(df[wavelength_col].unique())
    palette = dict(zip(wavelengths, sns.color_palette(n_colors=len(wavelengths))))
    if percent:
        new_amp_col = 'percent_max_abs'
    else:
        new_amp_col = amp_col
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 3))
    dfs=[]
    for w in wavelengths:
        sub = df[df[wavelength_col] == w].copy()
        sub['percent_max_abs'] = (sub[amp_col] / sub[amp_col].max()) * 100
        color = palette[w]
        # vertical lines
        ax.vlines(x=sub[mu_col], ymin=0, ymax=sub[new_amp_col], color=color, alpha=0.7, linewidth=1.5, zorder=-1)
        # endpoints
        ax.scatter(sub[mu_col], sub[new_amp_col], color=color, edgecolor='black', s=50, label=str(w))
        if percent:
            ax.plot(sub[mu_col], sub[new_amp_col], color='gray', alpha=0.5, zorder=-2)

        dfs.append(sub)
    df = pd.concat(dfs, axis=0)
    df=df.sort_values(by=mu_col)
    ax.set_xlabel(mu_col)
    ax.set_ylabel(new_amp_col)
    ax.set_title(f'UV Absorption Lollipop Plot: {sample_id}')
    # set sensible limits with a small margin
    x_min, x_max = df[mu_col].min(), df[mu_col].max()
    x_pad = 0.1
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    y_max = df[new_amp_col].max()
    if percent:
        ax.set_ylim(0, y_max * 1.1)
    else:
        ax.set_yscale("log")
    ax.legend(title=wavelength_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    if save_dir is not None:
        plt.savefig(f"{save_dir}/{sample_id}_uv_absorption_lollipop.png", dpi=300)
        plt.close()
    else:
        plt.show()
    return ax

def plot_tic_vs_time(sample_id="04036CB-INO2-P12_1", data_df=None, ax=None, save_dir=None):
    """Plot Total Ion Chromatogram (TIC) over UV absorption data."""
    if ax is None:
        fig, ax1 = plt.subplots(figsize=(12, 4))
    else:
        ax1 = ax

    # Plot UV absorption
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('UV Absorption (mAU)', color='tab:blue')
    sns.lineplot(data=data_df, x='time_peak_maximum_msminutes', y='maximum_absorption_mau', ax=ax1, color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # Create a second y-axis for TIC
    ax2 = ax1.twinx()
    ax2.set_ylabel('Total Ion Chromatogram (TIC)', color='tab:red')
    sns.lineplot(data=data_df, x='time_peak_maximum_msminutes', y='total_ion_chromatogram_tic', ax=ax2, color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    plt.title(f'TIC and UV Absorption for {sample_id}')
    plt.tight_layout()

    if save_dir is not None:
        plt.savefig(f"{save_dir}/{sample_id}_tic_uv_absorption.png", dpi=300)
        plt.close()
    else:
        plt.show()



