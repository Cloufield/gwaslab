import time
class Log():
    def __init__(self):
        self.log_text=str(time.ctime(time.time()))+ " " + "Sumstats Object created."+ "\n"
    def write(self,*message,end="\n",show_time=True):
        if show_time is True:
            print(str(time.ctime(time.time())),*message,end=end)
            self.log_text = self.log_text + str(time.ctime(time.time())) + " " + " ".join(map(str,message)) + end
        else:
            print(*message,end=end)
            self.log_text = self.log_text + " ".join(map(str,message)) + end
    def show(self):
        print(self.log_text)
    def save(self,path):
        with open(path,"w") as f:
            print(str(time.ctime(time.time())) + " " + " -Save log file to : ", path)
            f.write(self.log_text)