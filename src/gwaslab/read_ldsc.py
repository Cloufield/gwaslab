import re
import pandas as pd
import numpy as np

def read_ldsc(filelist=[],mode="h2"):
#h2 mode
#####################################################################    
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
                    if not line: break
                        
                
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
            #summary = summary.append(row,ignore_index=True)
            row = pd.DataFrame([row], columns = summary.columns)
            summary = pd.concat([summary, row], ignore_index=True)
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
    return summary   


def read_popcorn(filelist=[]):
#h2 mode
#####################################################################    
    summary = pd.DataFrame(columns = ["Filename", 'sfile1', 'sfile2', 'mode', 'pg', 'pg_se','pg_z','pg_p', 'h1^2', 'h1^2_se','h1^2_z','h1^2_p', 'h2^2', 'h2^2_se','h2^2_z','h2^2_p'])
    
    for index, ldscfile in enumerate(filelist):
        print("Loading file "+str(index+1)+" :" + ldscfile +" ...")
        
        row={}
        
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
            
            while not re.compile('^Jackknife iter:').match(line):
                line = file.readline()
                if not line: break
            ##next line lambda gc
            line = file.readline()
            
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

        #summary = summary.append(row,ignore_index=True)
        row = pd.DataFrame([row], columns = summary.columns)
        summary = pd.concat([summary, row], ignore_index=True)
    return summary   