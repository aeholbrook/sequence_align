import numpy as np

class local_alignment:

    def __init__(self,sMat,cseq,rseq, gap_open):
        self.sMat=sMat #the sub matrix we'll be using
        self.cseq=cseq #the column sequence
        self.rseq=rseq #the row sequence
        self.gap_open=gap_open #the gap score. I called it gap_open in case I wanted to implement affine, which I didn't
        out = np.zeros([len(cseq)+1,len(rseq)+1]) #values
        dirct = np.empty([len(cseq)+1,len(rseq)+1], dtype=object) #direction
        out[0,0] = 0
        arr = ['h','v','d'] #I kno there are easier ways of doing this but i wanted to get better at using numpy
        for idx in range(1,out.shape[1]):
            out[0,idx] = out[0,idx-1]
        for idx2 in range(1,out.shape[0]):
            out[idx2,0] = out[idx2-1,0]
        for col in range(1,out.shape[1]):
            for row in range(1,out.shape[0]):
                n0 = out[row,col-1]+gap_open
                n1 = out[row-1,col]+gap_open
                diag = out[row-1,col-1]
                n2 = sMat.check_score(cseq[row-1],rseq[col-1])+diag
                arr2 = [n0,n1,n2]
                mx = np.amax(arr2)
                if mx > 0: out[row,col] = mx 
                else: mx = 0
                dirct[row,col]=''.join(np.take(arr, np.where(arr2==mx))[0])
            dirct[dirct==None] = "*"
        mx = np.amax(out)
        arr1 = np.empty((0,0),dtype=str)
        arr2 = np.empty((0,0),dtype=str)
        it = [np.where(out==mx)[0][0],np.where(out==mx)[1][0]]
        for c in range(len(cseq) - it[0]):
            arr1 = np.append(arr1,cseq[len(cseq)-(c+1)])
            arr2 = np.append(arr2,"-")
        for r in range(len(rseq) - it[1]):
            arr1 = np.append(arr1,"-")
            arr2 = np.append(arr2,rseq[len(rseq)-(r+1)])
        while it[0] > 0 and it[1] > 0 and out[it[0],it[1]] > 0:
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
                print(arr2)
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
        self.out = out
        self.dirct = dirct
        self.mx = mx