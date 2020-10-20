#This generates a substitution matrix and returns a score for
#a given  csv matrix (I took the title off)

import numpy as np

class sub_matrix:
    def __init__(self,mat):
        self.prot_array=np.asarray(mat[0,])
        self.sub_scores=mat[1:,]
    
    def check_score(self,prot,prot2):
        pos1 = int(np.where(self.prot_array==prot)[0])
        pos2 = int(np.where(self.prot_array==prot2)[0])
        return int(self.sub_scores[pos1,pos2])