import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from gwaslab.info.g_Log import Log

#STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
#1	TssA	Active TSS	Red	255,0,0
#2	TssAFlnk	Flanking Active TSS	Orange Red	255,69,0
#3	TxFlnk	Transcr. at gene 5' and 3'	LimeGreen	50,205,50
#4	Tx	Strong transcription	Green	0,128,0
#5	TxWk	Weak transcription	DarkGreen	0,100,0
#6	EnhG	Genic enhancers	GreenYellow	194,225,5
#7	Enh	Enhancers	Yellow	255,255,0
#8	ZNF/Rpts	ZNF genes & repeats	Medium Aquamarine	102,205,170
#9	Het	Heterochromatin	PaleTurquoise	138,145,208
#10	TssBiv	Bivalent/Poised TSS	IndianRed	205,92,92
#11	BivFlnk	Flanking Bivalent TSS/Enh	DarkSalmon	233,150,122
#12	EnhBiv	Bivalent Enhancer	DarkKhaki	189,183,107
#13	ReprPC	Repressed PolyComb	Silver	128,128,128
#14	ReprPCWk	Weak Repressed PolyComb	Gainsboro	192,192,192
#15	Quies	Quiescent/Low	White	255,255,255

color_dict={
    "E1": np.array([255,0,0]),
    "E2": np.array([255,69,0]),
    "E3": np.array([50,205,50]),
    "E4": np.array([0,128,0]),
    "E5": np.array([0,100,0]),
    "E6": np.array([194,225,5]),
    "E7": np.array([255,255,0]),
    "E8": np.array([102,205,170]),
    "E9": np.array([138,145,208]),
    "E10":np.array([205,92,92]),
    "E11":np.array([233,150,122]),
    "E12":np.array([189,183,107]),
    "E13":np.array([128,128,128]),
    "E14":np.array([192,192,192]),
    "E15":np.array([255,255,255])
}

color_dict_i={
    1: np.array([255,0,0]),
    2: np.array([255,69,0]),
    3: np.array([50,205,50]),
    4: np.array([0,128,0]),
    5: np.array([0,100,0]),
    6: np.array([194,225,5]),
    7: np.array([255,255,0]),
    8: np.array([102,205,170]),
    9: np.array([138,145,208]),
    10:np.array([205,92,92]),
    11:np.array([233,150,122]),
    12:np.array([189,183,107]),
    13:np.array([128,128,128]),
    14:np.array([192,192,192]),
    15:np.array([255,255,255])
}


def _plot_chromatin_state(region_chromatin_files, 
                          region_chromatin_labels,
                          region, 
                          fig, 
                          ax,
                          xlim_i,
                          fontsize = 12,
                          font_family = "Arial",
                          log=Log(),
                          verbose=True):
    '''
    files : a list of numbers
    '''    
    target_chr = region[0]
    target_start = region[1]
    target_end = region[2]

    offset_i = xlim_i[0] - region[1]

    ax.set_ylim([-0.05,0.1*len(region_chromatin_files)-0.05])
    ax.set_xlim([offset_i+target_start,offset_i+target_end])

    px_for_01 = ax.transData.transform([0,0])[1] - ax.transData.transform([0,0.1])[1]

    point=fig.dpi/72
    points_for_01 = px_for_01*72 / fig.dpi

    # each tissue
    for i,file in enumerate(region_chromatin_files):
        log.write(" -Loading : {}".format(file), verbose=verbose)
        enh = pd.read_csv(file,sep="\t",header=None)
        enh.columns=["ID","START","END","STATE"]
        enh["CHR"] =  enh["ID"].str.extract(r"chr([0-9]+)").astype("float").astype("Int64")
        enh["STATE_i"] =  enh["STATE"].str.extract(r"([0-9]+)_*").astype("float").astype("Int64")
        enh_in_region = (enh["CHR"] == target_chr) & ((enh["END"] > target_start) & (enh["START"]<target_end))
        df =enh.loc[enh_in_region,["STATE_i","START","END"]].sort_values("STATE_i",ascending=False)
        log.write("  -Number of records in specified region: {}".format(len(df)), verbose=verbose)
        # each block
        for index, row in df.iterrows():
            color=color_dict_i[row["STATE_i"]]
            ax.plot([offset_i + row["START"] ,offset_i + row["END"]],
                    [i*0.1,i*0.1],
                    c=color/255,linewidth=points_for_01,solid_capstyle="butt",rasterized=True)
    
    ## add stripe label
    if len(region_chromatin_labels) == len(region_chromatin_files):
        ax.set_yticks([i*0.1 for i in range(len(region_chromatin_labels))], region_chromatin_labels, fontsize=fontsize, family=font_family)
    else:
        ax.set_yticks(ticks=[])
    
    #ax.set_xticks(ticks=[])
    ax.invert_yaxis()
    return fig