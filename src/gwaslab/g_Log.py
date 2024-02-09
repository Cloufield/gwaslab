import time
class Log():
    def __init__(self):
        self.log_text=str(time.strftime('%Y/%m/%d %H:%M:%S'))+ " " + "Sumstats Object created."+ "\n"
    
    def write(self,*message,end="\n",show_time=True, verbose=True):
        if show_time is True:
            if verbose: print(str(time.strftime('%Y/%m/%d %H:%M:%S')),*message,end=end)
            self.log_text = self.log_text + str(time.strftime('%Y/%m/%d %H:%M:%S')) + " " + " ".join(map(str,message)) + end
        else:
            if verbose: print(*message,end=end)
            self.log_text = self.log_text + " ".join(map(str,message)) + end
    
    def warning(self,*message,end="\n",show_time=True, verbose=True):
        self.write(" #WARNING! {}".format(" ".join(map(str,message))), 
                   end=end, 
                   show_time=show_time,
                   verbose=verbose)

    def show(self):
        print(self.log_text)
    def save(self,path,verbose=True):
        with open(path,"w") as f:
            if verbose: print(str(time.strftime('%Y/%m/%d %H:%M:%S')) + " " + " -Save log file to : ", path)
            f.write(self.log_text)


    def log(self,*message,end="\n",show_time=True, verbose=True):
        if show_time is True:
            if verbose: print(str(time.strftime('%Y/%m/%d %H:%M:%S')),*message,end=end)
            self.log_text = self.log_text + str(time.strftime('%Y/%m/%d %H:%M:%S')) + " " + " ".join(map(str,message)) + end
        else:
            if verbose: print(*message,end=end)
            self.log_text = self.log_text + " ".join(map(str,message)) + end