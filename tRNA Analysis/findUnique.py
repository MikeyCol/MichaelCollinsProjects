import sys
class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
       
        yield header,sequence



class tRNAseq:
    '''
    Represents the tRNA sequence of a specific tRNA 
    instantiation: 
    thisSeq = tRNAseq ('header', 'sequence')

    The functions in this class assume that they will be used only in the order that they are used in the main function of this file
    That is, findPowerSet --> findUnique --> findEssential

    the index name keeps track of the location of each tRNAseq in the shared lists. The way it is modified assumes you follow the ordering.
    '''
    essentialList = [] # list of essential sets of tRNA
    uniqueList = [] # list of unique but not essential sets of tRNA
    powerSetList = [] # list of all powersets generated from tRNA

    def __init__ (self, seq, header):
        self.header = header.replace(' ', '') # string reprentation of the header with the spaces removed
        self.seq = seq # string representation of the sequence
        self.index = -2 # used to find where the specific sequece is within the shared lists (essentialList, uniqueList, powerSetList)
       

    def findPowerSet(self):
        '''
        The first step in the algorithm, this function finds the powerset of each sequential substring in the sequence.
        ie:
        seq = 'abcd'
        returns ['a' 'ab' 'abc' 'abcd' 'b' 'bc' 'bcd' 'c' 'cd' 'd']

        it also appends the generated set to powerSetList and sets self.index to its location in powerSetList
        '''
        self.seq = self.seq.replace('-','')
        self.seq = self.seq.replace('_','')
        i = 0
        j = 0
        
        tRNASet = set()  # set to be returned
        for i in range(0,len(self.seq)): # iterates through all characters in sequence
            for j in range(i, len(self.seq)+1): # for each character it iterates through each character in front of it and adds that substring to the list       
                tRNASet.add(self.seq[i:j])
        
        
        self.index = len(tRNAseq.powerSetList)
        tRNAseq.powerSetList.append(tRNASet)

        return tRNASet    

    def findUnique(self):
        '''
        Second step in the algorithm, this function islotes the substrings in each power set in powerSetList that are only found in one power set.
        Returns a set containing the unique sub sequences in the power set

        call this only after findPowerSet has been run for each desired tRNAseq
        '''
        # creates a list of sets without the current tRNAset
        othertRNAsets = tRNAseq.powerSetList.copy()
        othertRNAsets.pop(self.index)
        # computes the union of those sets
        otherSeqs = set()
        for tRNAsequence in othertRNAsets:
            otherSeqs = otherSeqs.union(tRNAsequence)
         
        # computes the difference of the tRNAseq's power set and the union of all other tRNAseq's powersets.
        # this creates a list containing the unique subsequences
        unique = set()
        unique = unique.union(tRNAseq.powerSetList[self.index].difference(otherSeqs))
        self.index = len(tRNAseq.uniqueList) # index updated
        tRNAseq.uniqueList.append(unique) # appends set to uniqueList
        return unique

    def findEssential(self):
        '''
        Third and final step in the algorithm, this function eliminates the redundant substrings in uniqueList
        for instance if 'GTCA' is unique, so will every string following it, 'GTCAA' 'GTCAAT' .. 'GTCA....'
        This function removes these non-essential subsequences

        This function works by sorting the set of unique sub sequences using pyhtons default string sorting algorithm so that all nonessential strings are 
        placed just after the smallest version of itself.

        It would sort the previous strings into: ['GTCA' 'GTCAA' 'GTCAAT' ... 'GTCA....' ]
        it then iterates through the list and saves the first string as the 'essentialSeq' 
        this is used to compare to the next item in the list to see if it has the same essential component
        ie: 'GTCA'
        if it has a different beginning sequence it registers it as an essential subsequence and replaces 'essentialSeq' with the new subsequecne

        This chatches all strings in front of the essential seq but misses those behind it.
        To catch these, the strings are sorted by thier inverse and the same iteration is run again over the newly sorted list

        it then returns a set that consists of the essential subsequences 
        '''
        sortedUnique = sorted(tRNAseq.uniqueList[self.index]) # Lexigraphically sorts the set in uniqueList cooresponding to this tRNAseq
        essentialSeq = prevSeq = '---' # sets initial values so that it is allways replaced
        essential = set() # set to return
        for thisSeq in sortedUnique: # iterates through all sorted subsequences
            if thisSeq[0:len(essentialSeq)] != essentialSeq: # checks beggining for the essential subsequence and resets essentialSeq and appends to essentail if its different
                essentialSeq = thisSeq 
                essential.add(essentialSeq)
            
            prevSeq = thisSeq



        # This portion works exactly the same as the previous loop but checks the ending of the string instead of the beginning
        sortedEssential = sorted(essential, key=sortByInvert) # sorts the list by the inverted versions of each string
        
        essentialSeq = prevSeq = '---'
        essential = set()
        for thisSeq in sortedEssential:
            if thisSeq[len(thisSeq)-len(essentialSeq):] != essentialSeq:
                essentialSeq = thisSeq
                essential.add(essentialSeq)
            
            prevSeq = thisSeq

        self.index = len(tRNAseq.essentialList)
        tRNAseq.essentialList.append(essential)
        return essential

    def printSeq(self):
        '''
        prints the results of findPowerSet, findUnique and findEssential in a readable format
        '''
        print(self.header)
        print(self.seq)
        toPrint = [] # list of strings to print out 
        for sequence in tRNAseq.essentialList[self.index]: # iterates through the set of essential subsequences for the current tRNAseq
            position = 0 # keeps track of the loctation in the original sequence
            while(position < len(self.seq)): # iterates through each position in the sequence
                position = self.seq.find(sequence, position) # finds the subsequence in the original sequence
                if position < 0: # breaks when position is not found
                    break
                toPrint.append(((position)*'.')+self.seq[position:position+len(sequence)]) # appends the subsequence to toPrint with '.' to align the subsequecne with its location in the original sequence
                position += 1
                
        # sorts and prints the list of essential subsequences
        toPrint.sort(key=sortByDot) 
        for ordered in toPrint:
            print(ordered)


def sortByDot(str): # sorts by number of '.' characters in string
    return str.count('.')
def sortByHeader(tRNAseq): # sorts by the header of a tRNAseq object
    return tRNAseq.header
def sortByInvert(str): # sorts by inverted versions of strings
    return str[::-1]

def main(inCL=None):
    '''
    takes input file front stdin, assuemes input file is in the Fata format
    runs each function on the list of tRNAseqs and prints results to stdout
    '''
    reader = FastAreader()
    tRNAs = [] # this will contain all 22 tRNAseq objects
    headseq = list(reader.readFasta())
    for head, seq in headseq:

        tRNAs.append(tRNAseq(seq, head))
        tRNAs[-1].findPowerSet()

    for tRNA in tRNAs:
        tRNA.findUnique()
        tRNA.findEssential()
            
       
    tRNAs.sort(key=sortByHeader)
    for tRNA in tRNAs:
        tRNA.printSeq()
        
if __name__ == "__main__":
    main()  



