import numpy as np

class global_alignment:

    def __init__(self,sMat,cseq,rseq,gap_open=0,semiglobal=False):
        self.gap_open = gap_open
        out = np.zeros([len(cseq)+1,len(rseq)+1]) #values
        dirct = np.empty([len(cseq)+1,len(rseq)+1], dtype=object) #direction
        out[0,0] = 0
        for idx in range(1,out.shape[1]):
            out[0,idx] = out[0,idx-1]+gap_open
        for idx2 in range(1,out.shape[0]):
            out[idx2,0] = out[idx2-1,0]+gap_open
        for col in range(1,out.shape[1]):
            for row in range(1,out.shape[0]):
                n0 = out[row,col-1]+gap_open
                n1 = out[row-1,col]+gap_open
                diag = out[row-1,col-1]
                n2 = sMat.check_score(cseq[row-1],rseq[col-1])+diag
                arr = ['h','v','d']
                arr2 = [n0,n1,n2]
                mx = np.amax(arr2)
                out[row,col] = mx
                dirct[row,col]=''.join(np.take(arr, np.where(arr2==mx))[0])
        
        self.out = out.astype(int)
        print(out)
        if semiglobal==False:
            dirct[dirct==None] = "*"
            col1 = int(np.amax(out[out.shape[0]-1,:]))
            col2 = int(np.amax(out[:,out.shape[1]-1]))
            mx = max(col1,col2)
            arr1 = np.empty((0,0),dtype=str)
            arr2 = np.empty((0,0),dtype=str)
            #it = [np.where(out==mx)[0][-1],np.where(out==mx)[1][-1]]
            it = [out.shape[0]-1,out.shape[1]-1]
            for c in range(len(cseq) - it[0]):
                arr1 = np.append(arr1,cseq[len(cseq)-(c+1)])
                arr2 = np.append(arr2,"-")
                dirct[len(cseq)-c,len(rseq)] = "|"
            for r in range(len(rseq) - it[1]):
                arr1 = np.append(arr1,"-")
                arr2 = np.append(arr2,rseq[len(rseq)-(r+1)])
                dirct[len(cseq)-r,len(rseq)] = "-"
            while it[0] > 0 and it[1] > 0:
                if "d" in dirct[it[0],it[1]]:
                    arr1 = np.append(arr1, cseq[it[0]-1])
                    arr2 = np.append(arr2, rseq[it[1]-1])
                    dirct[it[0],it[1]] = '\\'
                    it[0] = it[0]-1
                    it[1] = it[1]-1
                elif "h" in dirct[it[0],it[1]]:
                    arr1 = np.append(arr1,"-")
                    arr2 = np.append(arr2, rseq[it[1]-1])
                    dirct[it[0],it[1]] = "-"
                    it[1] = it[1]-1
                else:
                    arr1 = np.append(arr1, cseq[it[0]-1])
                    arr2 = np.append(arr2, "-")
                    dirct[it[0],it[1]] = "|"
                    it[0] = it[0]-1
            for l in range(len(cseq) - len(arr1)):
                arr1 = np.append(arr1,cseq[l])
                arr2 = np.append(arr2,"-")
                dirct[l+1,0] = "|"
            for l in range(len(rseq) - len(arr2)):
                arr1 = np.append(arr1,"-")
                arr2 = np.append(arr2,rseq[l])
                dirct[l+1,0] = "-"

            self.arr1 = np.flip(arr1)
            self.arr2 = np.flip(arr2)
            self.dirct = dirct
            self.mx = mx

        else: 
            dirct[dirct==None] = "*"
            col1 = int(np.amax(out[out.shape[0]-1,:]))
            col2 = int(np.amax(out[:,out.shape[1]-1]))
            mx = max(col1,col2)
            arr1 = np.empty((0,0),dtype=str)
            arr2 = np.empty((0,0),dtype=str)
            it = [np.where(out==mx)[0][-1],np.where(out==mx)[1][-1]]
            #it = [out.shape[0]-1,out.shape[1]-1]
            for c in range(len(cseq) - it[0]):
                arr1 = np.append(arr1,cseq[len(cseq)-(c+1)])
                arr2 = np.append(arr2,"-")
                dirct[len(cseq)-c,len(rseq)] = "|"
            for r in range(len(rseq) - it[1]):
                arr1 = np.append(arr1,"-")
                arr2 = np.append(arr2,rseq[len(rseq)-(r+1)])
                dirct[len(cseq)-r,len(rseq)] = "-"
            while it[0] > 0 and it[1] > 0:
                if "d" in dirct[it[0],it[1]]:
                    arr1 = np.append(arr1, cseq[it[0]-1])
                    arr2 = np.append(arr2, rseq[it[1]-1])
                    dirct[it[0],it[1]] = '\\'
                    it[0] = it[0]-1
                    it[1] = it[1]-1
                elif "h" in dirct[it[0],it[1]]:
                    arr1 = np.append(arr1,"-")
                    arr2 = np.append(arr2, rseq[it[1]-1])
                    dirct[it[0],it[1]] = "-"
                    it[1] = it[1]-1
                else:
                    arr1 = np.append(arr1, cseq[it[0]-1])
                    arr2 = np.append(arr2, "-")
                    dirct[it[0],it[1]] = "|"
                    it[0] = it[0]-1
            for l in range(len(cseq) - len(arr1)):
                arr1 = np.append(arr1,cseq[l])
                arr2 = np.append(arr2,"-")
                dirct[l+1,0] = "|"
            for l in range(len(rseq) - len(arr2)):
                arr1 = np.append(arr1,"-")
                arr2 = np.append(arr2,rseq[l])
                dirct[l+1,0] = "-"

            self.arr1 = np.flip(arr1)
            self.arr2 = np.flip(arr2)
            self.dirct = dirct
            self.mx = mx