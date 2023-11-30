from gwaslab.g_Sumstats import Sumstats

class SumstatsT( ):
    def __init__(self, sumstatsObject):
        
        print(isinstance(sumstatsObject, Sumstats))

        self.snp_info_cols = []
        self.stats_cols =[]
        self.other_cols=[]

        for i in sumstatsObject.data.columns:
            if i in ["SNPID","rsID","CHR","POS","EA","NEA"]:
                self.snp_info_cols.append(i)
            elif i in ["BETA","SE","P","N","EAF","Z","T","OR","OR95L","OR95U"]:
                self.stats_cols.append(i)
            else:
                self.other_cols.append(i)

        self.snp_info = sumstatsObject.data.loc[:,self.snp_info_cols]

    def match(self, sumstatsObject):
        pass

            
        
