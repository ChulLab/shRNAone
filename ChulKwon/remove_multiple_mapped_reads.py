# remove_multiple_mapped_reads.py
# kwonschul

import sys

def remove_multiple_mapped_reads(sInSamFile, sOutSamFile):

    with open(sOutSamFile, 'w') as f:
        # Find FLAG:256 (multiple mapped) which has a FLAG:16 (reverse strand) mate 
        # --> Keep them.
        # Other FLAG:256 will be removed
        
        # First for loop to find the list
        
        lFlag16 = []
        lFlag256 = []
        for line in open(sInSamFile):            
            if line.startswith('@'):
                pass
            else:                
                field = line.split('\t')
                sQname = field[0]
                sFlag = field[1]
                
                if sFlag == '16':
                    lFlag16.append(sQname)
                elif sFlag == '256':
                    lFlag256.append(sQname)                
        lIntersect16and256 = list(set(lFlag16) & set(lFlag256))
        dIntersect16and256 = {k:[] for k in lIntersect16and256}
        
        # Second for loop to write a new SAM file
        
        for line in open(sInSamFile):            
            if line.startswith('@'):
                f.write(line)
            else:
                field = line.split('\t')
                sQname = field[0]
                sFlag = field[1]
                
                if sFlag == '0':
                    f.write(line)
                elif sFlag == '256':
                    if sQname in dIntersect16and256: # This 256 should be kept.
                        f.write(line)
                    else:
                        pass
                else:
                    pass # FLAG: 4, FLAG: 272, FLAG: 16 will be removed   
                        

if __name__ == '__main__':
    sInSamFile = sys.argv[1]
    sOutSamFile = sys.argv[2]
    remove_multiple_mapped_reads(sInSamFile, sOutSamFile)
