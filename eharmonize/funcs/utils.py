import pandas as pd
import pkg_resources
import numpy as np
import json
import os
import sys 
from datetime import datetime

class logstr():
    
    def __init__(self):
        now = datetime.now()
        dt_string = now.strftime("%B %d, %Y %H:%M:%S")

        self.message = "\n\nThis is a text file log for the eharmonize program.\n"
        self.run_settings = {
            "user": os.environ["USER"],
            "date": dt_string,
            "command": " ".join(sys.argv),
            }

    def stdO_file(self, message):
        self.message += message
        print(message)

    def to_file(self, outfile):
        log2write = "Runtime Information\n"
        for k, v in self.run_settings.items():
            log2write += "* %s: %s\n" %(k.title(), v)

        log2write += self.message

        with open(outfile, 'w') as f:
            f.write(log2write)

    # untested
    def abort(self, errortype, message, outfile):
      self.message += message
      self.to_file(outfile)
      raise errortype(message)

def input_check(dfI, covars):
    for covar in covars:
        if covar not in dfI.columns:
            raise KeyError("%s is a required column. Please check columns to ensure it is included/spelled correctly." %covar)

    cohort_info = {}
    sex_values = dfI.Sex.fillna("NA").value_counts()
    for k, v in sex_values.items():
        if k == 0:
            cohort_info["N_Females"] = v
        elif k == 1:
            cohort_info["N_Males"] = v
        else:
            cohort_info["N_Sex_%s" % k] = v
    
    site_values = dfI.SITE.fillna("NA").value_counts()
    for k, v in site_values.items():
        cohort_info["N_SITE_%s" %k] = v

    if "Dx" in dfI.columns:
        control_values = dfI.Dx.fillna("NA").value_counts()
        for k,v in control_values.items():
            if k == "case":
                cohort_info["N_cases"] = v
            elif k == "control":
                cohort_info["N_controls"] = v
            else:
                cohort_info["N_Dx_%s" %k] = v
    else:
        cohort_info["N_controls"] = dfI.shape[0]

    return cohort_info

def load_version(reference, metric):
    with open(pkg_resources.resource_filename(__name__, '../data/reference_meta.json'), 'r') as f:
        pkg_meta = json.load(f)
    pkg_settings = pkg_meta[reference]

    dfR = pd.read_csv(pkg_resources.resource_filename(__name__, '../data/%s' %pkg_settings[metric]["filepath"]))

    return pkg_settings, dfR

def age_check(dfI, pkg_settings, metric="FA"):
    outMessage = "\nChecking age range of input data\n"
    age_min = pkg_settings[metric]["age_min"]
    too_young = dfI.Age < age_min
    if too_young.sum() > 0:
        outMessage += "\nFound %i subjects younger than %i years old. Excluding them.\n" %(too_young.sum(), age_min)

    age_max = pkg_settings[metric]["age_max"]
    too_old = dfI.Age > age_max
    if too_old.sum() > 0:
        outMessage += "\nFound %i subjects older than %i years old. Excluding them.\n" %(too_old.sum(), age_max)

    age_excludes = (too_young | too_old)
    return age_excludes, outMessage

def column_match(dfI, roi_columns, metric):
    rois = roi_columns.str.split("_%s" %metric).str[0]
    roi_dict = {}
    tapetum = False
    outMessage = "\nChecking ROI names of input data\n"

    for c in dfI.columns:
        if c.startswith("TAP"):
            tapetum = True
        if c in ["AverageMD", "AverageAD", "AverageRD"]:
            new_name = c.replace("Average", "AverageFA_")
        else:
            new_name = c.replace("-", ".").replace("/", ".")
        if new_name in roi_columns:
            roi_dict[c] = new_name
        elif new_name in rois:
            roi_dict[c] = new_name + "_%s" % metric

    if tapetum:
        # TODO: actually switch names away from ENIGMA to JHU
        outMessage += "\nDetected Tapetum columns, which are not included in the ENIGMA-DTI version\n"
        unc_cols = dfI.columns[dfI.columns.str.startswith("UNC")]
        unc_replace = unc_cols.str.replace("UNC", "IFO")
        if len(unc_cols) > 0:
            outMessage += "Renaming UNC columns with IFO to match ENIGMA reference\n"
            for c, r in zip(unc_cols, unc_replace):
                new_name = r.replace("-", ".")
                if new_name in roi_columns:
                    roi_dict[c] = new_name
                elif new_name in rois:
                    roi_dict[c] = new_name + "_%s" % metric
        tap_cols = dfI.columns[dfI.columns.str.startswith("TAP")]
        tap_replace = tap_cols.str.replace("TAP", "UNC")
        if len(tap_cols) > 0:
            outMessage += "Renaming TAP columns with UNC to match ENIGMA reference\n"
            for c, r in zip(tap_cols, tap_replace):
                new_name = r.replace("-", ".")
                if new_name in roi_columns:
                    roi_dict[c] = new_name
                elif new_name in rois:
                    roi_dict[c] = new_name + "_%s" % metric

    dfi = dfI.rename(columns=roi_dict)

    return dfi, roi_dict, outMessage

def qc_images(rois2use, dfR, dfI, dfO, qcdir):
    import matplotlib.pyplot as plt
    import matplotlib
    import seaborn as sns

    matplotlib.use('Agg')
    font = {'size'   : 16}
    #       'family' : 'normal',
    #       'weight' : 'bold',

    matplotlib.rc('font', **font)

    linestyles = ['solid', 'dotted', 'dashed', 'dashdot']  

    age_min = dfI.Age.min()
    age_max = dfI.Age.max()
    dfR_filt = dfR.loc[(dfR.Age < age_max) & (dfR.Age > age_min)]
    if ("Dx" in dfI.columns) and ("Dx" in dfO.columns):
        dfI.loc[dfI.Dx == "control", "Data"] = "Raw - Control"
        dfI.loc[dfI.Dx == "case", "Data"] = "Raw - Case"
        dfO.loc[dfO.Dx == "control", "Data"] = "Harmonized - Control"
        dfO.loc[dfO.Dx == "case", "Data"] = "Harmonized - Case"
    else:
        dfI.loc[:, "Data"] = "Raw"
        dfO.loc[:, "Data"] = "Harmonized"
    plot_data = pd.concat([dfI, dfO], ignore_index=True)
    plot_data.dropna(subset=["Age"], inplace=True)
    plot_data.loc[:, "Age_int"] = plot_data["Age"].astype(int)

    lutfile = pkg_resources.resource_filename(__name__, '../data/ENIGMA_LUT.txt')
    lut = pd.read_csv(lutfile, delimiter="\t", index_col="Abbreviation")

    for roi in rois2use:
        # ax0 = sns.lmplot(plot_data, x="Age", y=roi, seed=0,
        #         row="SITE", hue="Data", lowess=True,
        #         line_kws={"ls": linestyles[:len(plot_data.Data.unique())]}, 
        #         scatter=True, scatter_kws={"s": 5, "alpha": 0.5},
        #         height=5, aspect=2)
        # ax = ax0.axes[0,0]
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        sns.lineplot(data=dfR, x="Age", y=roi,
                alpha = 0.5, ax=ax, color="black")
        sns.scatterplot(plot_data, x="Age", y=roi, 
                hue="SITE",style="Data", alpha=0.5, ax=ax)
        sns.lineplot(data=plot_data, x="Age_int", y = roi,
                hue="SITE", style="Data", ax=ax)

        # To have lines and histograms 
        # fig, ax = plt.subplots(1, 2 figsize=(10, 5))    
        # ax[0].plot(dfR_filt.Age, dfR_filt[roi],
        #         'k.', alpha = 0.5)
        # ax[0].set_title="ROI vs Age (years)"
        # ax[1].hist(dfR_filt[roi], bins=50, alpha=0.5)
        # sns.histplot(plot_data, x=roi, bins=10, hue="SITE",
        #         ax=ax[1])

        ax.set_title("ROI vs Age (years)")
        ax.set_xlabel("Age (years)")
        # TODO: this is specific to ENIGMA DTI
        # will break for other modalities (e.g., FreeSurfer)
        roi_key = roi.split("_")[0]
        if roi_key in lut.index:
            ax.set_ylabel(lut.loc[roi_key, "FullName"])
        else:
            ax.set_ylabel(roi)
        ax.legend(ncol=2, bbox_to_anchor=(1.1, 1.05))
        plt.tight_layout()
        plt.savefig(os.path.join(qcdir, "%s.png" %roi))
    return

