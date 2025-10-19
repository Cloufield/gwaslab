import seaborn as sns
import matplotlib.pyplot as plt
from gwaslab.io.io_process_args import _update_args
from gwaslab.io.io_process_args import _update_arg
from gwaslab.g_Log import Log

def _plot_associations(associations, 
                       values="Beta",
                        fig_args=None,
                        log=Log(),
                        verbose=True,
                        sort="P-value",
                        fontsize=12,
                        font_family="Arial",
                        cmap=None,
                        ylabel="Traits",
                        xlabel="rsID - Gene Name",
                        **args               
                       ):
    
    log.write("Start to create heatmap for associations...", verbose=verbose)
    
    cmap = _update_arg(cmap, "RdBu")

    log.write(f" -Value to plot: {values}", verbose=verbose)
    log.write(f" -Total number of associations :{len(associations)} ", verbose=verbose)
    log.write(f" -Sorting associations by :{sort} ", verbose=verbose)

    associations = associations.sort_values(by=[sort]).drop_duplicates(subset=["rsID","GWASCATALOG_TRAIT"])

    log.write(f" -Keeping unique rsID-traits pairs :{len(associations)} ", verbose=verbose)

    associations = associations.dropna(subset=[values])
    log.write(f" -Keeping associations without NA in {values} :{len(associations)} ", verbose=verbose)

    log.write(f" -Total number of unique variants for plotting :{associations['rsID'].nunique()} ", verbose=verbose)
    log.write(f" -Total number of unique traits for plotting:{associations['GWASCATALOG_TRAIT'].nunique()} ", verbose=verbose)

    matrix_beta = associations.pivot_table(index=['rsID','gene.geneName'], 
                                    columns='GWASCATALOG_TRAIT', 
                                    values=[values])

    matrix_beta = matrix_beta.astype("float64")
    
    matrix_beta.columns = matrix_beta.columns.droplevel(0)

    height = max(2, len(matrix_beta)//2)
    width = max(2, len(matrix_beta.columns)//2)

    fig_args = _update_args(fig_args, dict(figsize=(width,height)))
    fig,ax = plt.subplots(**fig_args)
    
    sns.heatmap(matrix_beta.T, annot=True, fmt=".2f",cmap="RdBu",ax=ax)
    

    ax.set_xlabel(xlabel,fontsize=fontsize, fontfamily=font_family)
    ax.set_ylabel(ylabel,fontsize=fontsize, fontfamily=font_family)


    log.write("Finished creating heatmap for associations.", verbose=verbose)
    return fig, ax