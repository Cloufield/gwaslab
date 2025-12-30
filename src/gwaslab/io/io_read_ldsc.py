from typing import List, Optional, Dict, Any
import re
import pandas as pd
import numpy as np

def read_ldsc(filelist: Optional[List[str]] = None, mode: str = "h2") -> pd.DataFrame:
    """
    Read LDSC output files and parse heritability or genetic correlation results.
    
    Parameters
    ----------
    filelist : list of str, optional
        List of paths to LDSC output files
    mode : str, default "h2"
        Mode to parse: "h2" for heritability, "rg" for genetic correlation
        
    Returns
    -------
    pd.DataFrame
        Parsed LDSC results
    """
    if filelist is None:
        filelist = []
    #h2 mode
#####################################################################    
    is_liab = False
    if mode=="h2":
        summary = pd.DataFrame(columns = ['Filename', 'h2_obs', 'h2_se','Lambda_gc','Mean_chi2','Intercept','Intercept_se',"Ratio","Ratio_se"])
        
        for index, ldscfile in enumerate(filelist):
            print("Loading file "+str(index+1)+" :" + ldscfile +" ...")
            
            row={}
            
            with open(ldscfile,"r") as file:
                row["Filename"]=ldscfile.split("/")[-1]
                line=""
                while not re.compile('^Total Observed scale h2').match(line):
                    line = file.readline()
                    if "h2_liab" in line:
                        is_liab = True
                    if not line: break
                        
                try:
                    ## first line h2 se
                    objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(line)
                    row["h2_obs"]=objects[1]
                    row["h2_se"]=objects[2]

                    ##next line lambda gc

                    objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(file.readline())
                    row["Lambda_gc"] = objects[1]
                    ##next line Mean_chi2

                    objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(file.readline())
                    row["Mean_chi2"]=objects[1]
                    ##next line Intercept

                    objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(file.readline())
                    row["Intercept"]=objects[1]
                    row["Intercept_se"]=objects[2]
                    ##next line Ratio
                    
                    lastline=file.readline()
                    if re.compile('NA').findall(lastline):
                        row["Ratio"]="NA"
                        row["Ratio_se"]="NA"
                    elif re.compile('<').findall(lastline):
                        row["Ratio"]="Ratio < 0"
                        row["Ratio_se"]="NA"
                    else:
                        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+').findall(lastline)
                        row["Ratio"]=objects[1]
                        row["Ratio_se"]=objects[2]
                except:
                        row["h2_obs"]="NA"
                        row["h2_se"]="NA"
                        row["Lambda_gc"] = "NA"
                        row["Mean_chi2"]="NA"
                        row["Intercept"]="NA"
                        row["Intercept_se"]="NA"
                        row["Ratio"]="NA"
                        row["Ratio_se"]="NA"

            #summary = summary.append(row,ignore_index=True)
            row = pd.DataFrame([row], columns = summary.columns)
            summary = pd.concat([summary, row], ignore_index=True)
            if is_liab == True:
                summary = summary.rename(columns={
                    "h2_obs":"h2_liab",
                    "h2_se":"h2_liab_se"
                })
###############################################################################
    if mode=="rg":
        summary = pd.DataFrame(columns = ['p1',
                                          'p2',
                                          'rg',
                                          'se' ,
                                          'z','p',
                                          'h2_obs','h2_obs_se',
                                          'h2_int','h2_int_se',
                                          'gcov_int','gcov_int_se']
                               )

        for index, ldscfile in enumerate(filelist):
            print("Loading file "+str(index+1)+" :" + ldscfile +" ...")
            
            with open(ldscfile,"r") as file:
                line=""
                while not re.compile('^Summary of Genetic Correlation Results').match(line):
                    line = file.readline()
                    if not line: break
                line = file.readline() # header        
                
                if "h2_liab" in line:
                        is_liab = True

                line = file.readline() #line1
                
                ## first line h2 se
                while line.strip():
                    row = line.split()
                    row_series = pd.DataFrame([row], columns = summary.columns)
                    #summary = summary.append(row_series, ignore_index=True)
                    summary = pd.concat([summary, row_series], ignore_index=True)
                    line = file.readline()
        summary = summary.loc[summary["rg"]!="NA",:].copy() 
        summary[['rg','se' ,'z','p','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se']]  = summary[['rg','se' ,'z','p','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se']].astype("float32")  
    if is_liab == True:
            summary = summary.rename(columns={
                "h2_obs":"h2_liab",
                "h2_se":"h2_liab_se"
            })          
    return summary   


def read_popcorn(filelist: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Read Popcorn output files and parse cross-population genetic correlation results.
    
    Parameters
    ----------
    filelist : list of str, optional
        List of paths to Popcorn output files
        
    Returns
    -------
    pd.DataFrame
        Parsed Popcorn results
    """
    if filelist is None:
        filelist = []
    #h2 mode
#####################################################################    
    summary = pd.DataFrame(columns = ["Filename", 'sfile1', 'sfile2', 'mode', 'pg', 'pg_se','pg_z','pg_p', 'h1^2', 'h1^2_se','h1^2_z','h1^2_p', 'h2^2', 'h2^2_se','h2^2_z','h2^2_p'])

    for index, ldscfile in enumerate(filelist):
        print("Loading file "+str(index+1)+" :" + ldscfile +" ...")

        row={}
        try:
            with open(ldscfile,"r") as file:
                row["Filename"]=ldscfile.split("/")[-1]
                line=""
                while not re.compile('^Invoking command').match(line):
                    line = file.readline()
                    if not line: break


                        ## first line h2 se
                objects = re.compile('--sfile1 ([^\s]+) --sfile2 ([^\s]+)[ /n]').findall(line)
                row["sfile1"]=objects[0][0]
                row["sfile2"]=objects[0][1]

                #while not re.compile(r'^Jackknife iter:').match(line):
                #    line = file.readline()
                #    print(line)
                #    if not line: break
                while not re.compile(r'P \(Z\)').findall(line.strip()):
                    line = file.readline()
                    if not line: break

                #objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(file.readline())
                objects = file.readline().split()
                row["h1^2"] = objects[1]
                row["h1^2_se"] = objects[2]
                row["h1^2_z"] = objects[3]
                row["h1^2_p"] = objects[4]

                objects = file.readline().split()
                row["h2^2"] = objects[1]
                row["h2^2_se"] = objects[2]
                row["h2^2_z"] = objects[3]
                row["h2^2_p"] = objects[4]

                objects = file.readline().split()
                row["mode"] = objects[0]
                row["pg"] = objects[1]
                row["pg_se"] = objects[2]
                row["pg_z"] = objects[3]
                row["pg_p"] = objects[4]
        except:
            continue

        #summary = summary.append(row,ignore_index=True)
        row = pd.DataFrame([row], columns = summary.columns)
        summary = pd.concat([summary, row], ignore_index=True)
    return summary

def read_greml(filelist: List[str] = []) -> pd.DataFrame:
    summary = pd.DataFrame(columns = ["Filename", 'Sum of V(G)/Vp', 'SE','Pval', 'n'])
    for index, hsqfile in enumerate(filelist):
        print("Loading file "+str(index+1)+" :" + hsqfile +" ...")
        row={}
        try:
            with open(hsqfile,"r") as file:
                row["Filename"]=hsqfile.split("/")[-1]
                line=""
                while not re.compile('Sum of').match(line):
                    line = file.readline()

                objects = line.strip().split("\t")
                row["Sum of V(G)/Vp"]=objects[1]
                row["SE"]=objects[2]

                while not re.compile(r'Pval').findall(line.strip()):
                    line = file.readline()
                    if not line: break
                #objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(file.readline())
                
                objects = line.split()
                row["Pval"] = objects[1]
                print(row["Pval"])
                
                while not re.compile(r'^n').findall(line.strip()):
                    line = file.readline()
                    if not line: break
                #objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(file.readline())
                objects = line.split()
                row["n"] = objects[1]            
                print(row["n"])
        except:
            continue
        row = pd.DataFrame([row], columns = summary.columns)
        summary = pd.concat([summary, row], ignore_index=True)
    return summary

def parse_ldsc_summary(ldsc_summary: str) -> pd.DataFrame:
    is_liab = False
    if "Liability" in ldsc_summary:
        is_liab = True
    lines = ldsc_summary.split("\n")
    
    columns = ['h2_obs', 'h2_se','Lambda_gc','Mean_chi2','Intercept','Intercept_se',"Ratio","Ratio_se","Catagories"]
    
    summary = pd.DataFrame(columns = columns)
    
    row={}
    
    try:
        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[0])
        row["h2_obs"]=objects[1]
        row["h2_se"]=objects[2]

        ##check categories
        if len(lines) == 6:
            objects = re.compile('  -Categories:(.+)').findall(lines[1])
            row["Catagories"] = objects[0].strip()
            lines.pop(1)
        else:
            row["Catagories"] = "NA"

        ##next line lambda gc
        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[1])
        row["Lambda_gc"] = objects[1]
        ##next line Mean_chi2

        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[2])
        row["Mean_chi2"]=objects[1]
        ##next line Intercept

        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[3])
        row["Intercept"]=objects[1]
        row["Intercept_se"]=objects[2]
        ##next line Ratio

        if re.compile('NA').findall(lines[4]):
            row["Ratio"]="NA"
            row["Ratio_se"]="NA"
        elif re.compile('<').findall(lines[4]):
            row["Ratio"]="Ratio < 0"
            row["Ratio_se"]="NA"
        else:
            objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+').findall(lines[4])
            row["Ratio"]=objects[1]
            row["Ratio_se"]=objects[2]
    except:
        row["h2_obs"]="NA"
        row["h2_se"]="NA"
        row["Lambda_gc"] = "NA"
        row["Mean_chi2"]="NA"
        row["Intercept"]="NA"
        row["Intercept_se"]="NA"
        row["Ratio"]="NA"
        row["Ratio_se"]="NA"
        row["Catagories"] = "NA"

    #summary = summary.append(row,ignore_index=True)
    row = pd.DataFrame([row], columns = summary.columns)
    if is_liab == True:
        row = row.rename(columns={
            "h2_obs":"h2_liab",
            "h2_se":"h2_liab_se"
        })
    return row

def parse_partitioned_ldsc_summary(ldsc_summary: str) -> pd.DataFrame:
    is_liab = False
    if "Liability" in ldsc_summary:
        is_liab = True
    summary = pd.DataFrame(columns = ['h2_obs', 'h2_se','Lambda_gc','Mean_chi2','Intercept','Intercept_se',"Ratio","Ratio_se"])
    lines = ldsc_summary.split("\n")
    row={}
    try:
        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[0])
        row["h2_obs"]=objects[1]
        row["h2_se"]=objects[2]

        ##next line lambda gc
        row["Categories"] = lines[1].replace("-Categories: ", "")

        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[2])
        row["Lambda_gc"] = objects[1]
        ##next line Mean_chi2

        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[3])
        row["Mean_chi2"]=objects[1]
        ##next line Intercept

        objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+|NA').findall(lines[4])
        row["Intercept"]=objects[1]
        row["Intercept_se"]=objects[2]
        ##next line Ratio

        if re.compile('NA').findall(lines[5]):
            row["Ratio"]="NA"
            row["Ratio_se"]="NA"
        elif re.compile('<').findall(lines[5]):
            row["Ratio"]="Ratio < 0"
            row["Ratio_se"]="NA"
        else:
            objects = re.compile('[a-zA-Z\s\d]+:|[-0-9.]+[e]?[-0-9.]+').findall(lines[5])
            row["Ratio"]=objects[1]
            row["Ratio_se"]=objects[2]
    except:
        row["h2_obs"]="NA"
        row["h2_se"]="NA"
        row["Lambda_gc"] = "NA"
        row["Mean_chi2"]="NA"
        row["Intercept"]="NA"
        row["Intercept_se"]="NA"
        row["Ratio"]="NA"
        row["Ratio_se"]="NA"

    #summary = summary.append(row,ignore_index=True)
    row = pd.DataFrame([row], columns = summary.columns)
    if is_liab == True:
        row = row.rename(columns={
            "h2_obs":"h2_liab",
            "h2_se":"h2_liab_se"
        })
    return row