#!/usr/bin/python3

from pathlib import Path
import os, sys, re
import glob, shutil, math
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


program_name = Path(__file__).name
program_version = "2.1d"
project_name = sys.argv[-1]

"""
HONZO'S UPDATE NOTES
v1  - completed read stats summary and mapping summary tabulation and plots
    - started development of variant calling summaries
v2  - implemented object classes for each sample variant calling to keep sample-specific data
    - will have to eventually use this object class for all sections instead of just variant calling
    - removed OOP in favour of nested dictionaries, which can be converted into DF easier
    - indicated version 2.0d for Docker - path variables are reset wrt docker container
    
"""

# Changeable parameters
int_maxgenes = 10 # top genes to plot on stacked barchart
plots = ["summary", "mapping", "variant"]
#plots = ["variant"]

# inputs
dir_project = Path("/fastq")
pth_rsvlog = dir_project / "RNASeqVariant.log"

# outputs
str_pdfout = "rsv_postprocess.pdf"
str_summarytable = "ReadSummary_all.txt"
str_mappingtable = "MapSummary_all.txt"
str_varianttable = "VariantSummary_all.txt"

plt.rcParams["font.family"] = "Montserrat"

global nglight, ngmedium, ngdark
nglight = ["#E7EBF4","#F9E8E9","#FFF9E7","#EFFAF7","#F4F4F6"] # generic fill color
ngmedium = ["#407BFF","#FF6066","#FFEFB3","#5ADAAB","#8D9CAA"] # used for abundance plots
ngdark = ["#30509D","#BE444F","#EAA017","#2E9F99","#606A74"] # generic edge color


print("\n[rv-p] This is Honzo's {0} for project: {1}".format(program_name, project_name))
if ("-h" or "--help") in sys.argv or len(sys.argv) < 2:
    print("USAGE: python {0} -p <project name>".format(program_name))

# ============================================[ FUNCTIONS ]===========================================

def double_reads(readcount, decider):
    if decider: return readcount * 2
    else: return readcount

def extract_sv(fullname):
    return re.search('.*(SV\d{5}[\-\_]\d{4}).*', fullname).group(1)

def get_colors(n):
    # using a 4-colour panel + grey, generate a colour gradient of len(lst_values)
    n -= 1 # because the last one is 'other', and will be assigned grey regardless
    colors = []
    # manually create a step-by-step gradient, then convert back to hex
    steps = math.ceil(n / 5) + 1
    for cindex in range(len(ngmedium)-1):
        r1, g1, b1 = tuple(int(ngmedium[cindex][1:][i:i+2], 16) for i in (0, 2, 4))
        r2, g2, b2 = tuple(int(ngmedium[cindex+1][1:][i:i+2], 16) for i in (0, 2, 4))
        rdelta, gdelta, bdelta = (r2-r1)/steps, (g2-g1)/steps, (b2-b1)/steps
        for _ in range(steps):
            colors.append((r1, g1, b1))
            r1 += rdelta
            g1 += gdelta
            b1 += bdelta
    #colors.append((r1, g1, b1))

    for i in range(len(colors)):
        colors[i] = tuple([int(j) for j in colors[i]])
        colors[i] = ['#%02x%02x%02x' % colors[i]]

    # [[color1], [color2], [color3], ...] > [color1, color2, color3]
    colors = [c[0] for c in colors]
    print("[rv-p] {0}-color gradient created for a list of size: {1}".format(len(colors), n))
    return colors

def set_aes(ax, title=None, ylabel=None, xlabel=None):
    ax.grid(which='major', color='#E2E2E2', linestyle='--', alpha=0.3)
    ax.set_facecolor('#F5F5F5')
    for side in ["top","right"]:
        ax.spines[side].set_visible(False)

    ax.set_title(title, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.tick_params(axis='x', labelrotation=90)
    ax.legend(framealpha=0.7)
    return ax


# ============================================[ WORKFLOW ]============================================

print("[rv-p] Starting workflow...")
os.chdir(dir_project)

if os.path.isdir("postProcess"):
    print("[Warning] postProcess folder already exists. Removing now.")
    shutil.rmtree("postProcess")
os.mkdir("postProcess")

# Creating stacked barplot for read summary
with PdfPages("postProcess/"+str_pdfout) as pdf:
    if "summary" in plots:
        # -------------------------------[ Read Stats Summary ]---------------------------------------

        print("[rv-p] Getting read statistics from readsummary files...")
        lst_readfiles = sorted(glob.glob("RSV_OUTPUT/*/SV*_S*_readsummary.txt"))
        lst_dfsummaries = [pd.read_table(f, index_col=0) for f in lst_readfiles]
        df_summary = pd.concat(lst_dfsummaries, axis=1)
        df_summary.to_csv("postProcess/{0}".format(str_summarytable), sep="\t")

        arr_todouble = np.frompyfunc(lambda idx: '(x2)' in idx, 1, 1)(df_summary.index)
        lst_indexnames = [x.replace("(x2)", "(SE)").strip() for x in list(df_summary.index)]
        df_processed = pd.DataFrame(index = lst_indexnames)

        # This is the dataframe we'll be actively working on
        for splname in df_summary:
            df_processed[splname] = list(map(double_reads, df_summary[splname], arr_todouble))

        fig, ax = plt.subplots(figsize=(11.0, 8.5))
        df_summary.columns = [c.split("_S")[0] for c in df_summary.columns]
        colors = get_colors(len(df_processed.index))
        
        for i, (rowname, row) in enumerate(df_processed.iterrows()):
            # label is the type of read, for example label = Raw FASTQ (SE)
            ax.bar(list(df_summary.columns), row, color=colors[i], width=0.9, label=rowname)

        set_aes(ax, title="Summary of Mapping Stats", ylabel="Number of Reads", xlabel="Sample Name")
        plt.tight_layout()
        pdf.savefig()
        plt.close()


    if "mapping" in plots:
        # ------------------------------[ Mapping Stat Summary ]--------------------------------------

        print("[rv-p] Getting read mapping summary from indexStat files...")
        lst_readfiles = sorted(glob.glob("RSV_OUTPUT/*/*_genome_indexStats.txt"))
        lst_dfsummaries = []
        
        # There are unmapped pairs and unmapped mates. In either case, we want to classify them as unmapped reads
        # We take all unmapped counts from the idxstats file and combine them into an index called "Unknown"
        for readfile in lst_readfiles:
            splname = extract_sv(readfile)
            df_readfile = pd.read_table(readfile, index_col=0, usecols=[0,2,3], names=["gene",splname,"unmapped"])
            df_readfile[splname]["*"] = np.sum(df_readfile["unmapped"])
            df_readfile.rename(index={"*":"Unknown"}, inplace=True)
            lst_dfsummaries.append(df_readfile[splname])

        df_summary = pd.concat(lst_dfsummaries, axis=1)
        df_summary.to_csv("postProcess/{0}".format(str_mappingtable), sep="\t")
        #print(df_summary)

        # sort genes from highest to lowest count
        df_summary["total_gene"] = df_summary.sum(axis=1)
        df_summary.sort_values(by="total_gene", ascending=0, inplace=True)
        df_summary.drop("total_gene", axis=1, inplace=True)

        if len(df_summary.index) > int_maxgenes:
            # Get non top-X genes, sum the columns into one line, add it to the df containing top-X genes   
            df_topgenes = df_summary.iloc[:int_maxgenes]
            df_botgenes = df_summary.iloc[int_maxgenes:]

            if "Unknown" in df_topgenes.index:
                print("[rv-p] Combining Unknown with Other labels")
                df_botgenes = df_botgenes.append(df_topgenes.loc["Unknown"])
                df_topgenes.drop("Unknown", inplace=True)

            arr_botgenes = df_botgenes.sum(axis=0).rename("Other")
            df_summary = df_topgenes.append(arr_botgenes)

        else:
            # Move Unknown to the back
            arr_unknown = df_summary.loc["Unknown"]
            df_summary.drop("Unknown", inplace=True)
            df_summary = df_summary.append(arr_unknown)

        # Generate plots
        fig, ax = plt.subplots(figsize=(11.0, 8.5))
        colors = get_colors(len(df_summary.index))
        names = df_summary.columns

        for i, (rowname, row) in enumerate(df_summary.iterrows()):
            if i == 0:
                ax.bar(names, row, width=0.9, label=rowname, color=colors[i])
                np_bottom = np.array(df_summary.iloc[i])
            elif rowname == "Unknown":
                ax.bar(names, row, width=0.9, label=rowname, color='#C8C8C8', bottom=np_bottom)
                np_bottom = np_bottom + np.array(df_summary.iloc[i])
            else:
                ax.bar(names, row, width=0.9, label=rowname, color=colors[i], bottom=np_bottom)
                np_bottom = np_bottom + np.array(df_summary.iloc[i])

        set_aes(ax, title="Summary of Top Genes", ylabel="Counts", xlabel="Sample Name")
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    if "variant" in plots:
        # ---------------------------------[ Variant Counts Summary ]-------------------------------------

        # we're going to read through the summary files instead of entire VCF records
        print("[rv-p] Getting variant statistics from VCF summaries")
        lst_readfiles = sorted(glob.glob("RSV_OUTPUT/*/*_genome_variantStats.txt"))

        samples = {}
        for f, readfile in enumerate(lst_readfiles):
            splname = extract_sv(readfile)
            #print("[rv-p] Getting variant counts for sample:", splname)

            with open(readfile, "rt") as streamreader:
                for line in streamreader:
                    if line.strip() == "# ST, Substitution types:":
                        streamreader.readline()
                        variants = {}
                        for i in range(12):
                            #print(streamreader.readline().strip())
                            fields = streamreader.readline().strip().split()
                            variants[fields[2]] = int(fields[3])

                        samples[splname] = variants
                        break

        # Summarize variant call counts into a dataframe
        df_variants = pd.DataFrame(samples)
        df_variants.to_csv("postProcess/{0}".format(str_varianttable), sep="\t")
        df_variants = df_variants.loc[~(df_variants==0).all(axis=1)]

        # Plot variants (counts)
        fig, ax = plt.subplots(figsize=(11.0, 8.5))
        colors = get_colors(12)
        names = df_variants.columns

        for i, (rowname, row) in enumerate(df_variants.iterrows()):
            if i == 0:
                ax.bar(names, row, width=0.9, label=rowname, color=colors[i])
                np_bottom = np.array(df_variants.iloc[i])
            else:
                ax.bar(names, row, width=0.9, label=rowname, color=colors[i], bottom=np_bottom)
                np_bottom = np_bottom + np.array(df_variants.iloc[i])

        set_aes(ax, title="Variant Calls", ylabel="Counts", xlabel="Sample Name")
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # Get variant distributions as a plot
        df_variantpcts = df_variants.div(df_variants.sum(axis=0), axis=1)
        df_variantpcts.fillna(0, inplace=True)
        df_variantpcts.dropna(axis=0, how='any', inplace=True)
        #df_variantpcts = df_variantpcts.loc[~(df_variantpcts==0).all(axis=1)]
        # print(df_variantpcts)

        # if a sample has no variants, we want to make the whole bar plot grey instead of having no values
        spls_novariants = [n for n, p in samples.items() if sum(p.values()) == 0]
        if len(spls_novariants) > 0:
            print("[rv-p] No variants detected in:", spls_novariants)
            arr_novariants = [1.0 if n in spls_novariants else 0.0 for n in df_variants.columns]
            df_variantpcts.loc["None"] = arr_novariants
        
        # Plot variant distributions
        fig, ax = plt.subplots(figsize=(11.0, 8.5))
        names = df_variants.columns

        for i, (rowname, row) in enumerate(df_variantpcts.iterrows()):
            if i == 0:
                ax.bar(names, row, width=0.9, label=rowname, color=colors[i])
                np_bottom = np.array(df_variantpcts.iloc[i])
            elif rowname == "None":
                ax.bar(names, row, width=0.9, label=rowname, color='#C8C8C8', bottom=np_bottom)
                np_bottom = np_bottom + np.array(df_variantpcts.iloc[i])
            else:
                ax.bar(names, row, width=0.9, label=rowname, color=colors[i], bottom=np_bottom)
                np_bottom = np_bottom + np.array(df_variantpcts.iloc[i])

        set_aes(ax, title="Variant Call Distribution", ylabel="Percentage (%)", xlabel="Sample Name")
        plt.tight_layout()
        pdf.savefig()
        plt.close()

print("[rv-p] {0} completed.".format(program_name))

