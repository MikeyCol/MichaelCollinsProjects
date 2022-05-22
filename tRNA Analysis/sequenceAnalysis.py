class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        self.protein = protein
        self.aaCompDict = self.aaComposition()

    def aaCount (self):
        aaFinal = 0;
        for aa in self.protein:
            if aa in self.aa2mw:
                aaFinal += 1
        return aaFinal


    def pI (self):
        '''
        returns the pI of the amino acid as a float to the second decimal
        '''
        zeroPH = 0 # at termination of loop, this is the pH that creates a charge in the amino acid closest to zero 
        for pH100 in range(0, 1400+1): # iterates over all pH values and finds the closest to zero
            pH = pH100 / 100
            if abs(self._charge_(pH)) < abs(self._charge_(zeroPH)):
                # print(self._charge_(pH))
                # print(self._charge_(zeroPH))
                zeroPH = pH
        return zeroPH


    def aaComposition (self) :
        '''
        used by _init_ 
        returns a dicitonary containig the number of times every amino acid appears in the protein
        Ex:
            'A' is the key to the number of Alaniens
            'G' is the ket the the nubmer of Glycenes
            etc.
        '''
        aaCompDict = {
        'A': 0,  'G': 0,  'M': 0, 'S': 0, 'C': 0,
        'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
        'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
        'W': 0,  'F': 0, 'L': 0, 'R': 0, 'Y': 0
        }
        
        for aa in self.protein: # iterates trough the protein's amino acid sequence and increments each 
            if aa in aaCompDict:
                aaCompDict[aa] += 1  
        return aaCompDict
    

    def _charge_ (self, pH):
        '''
        Finds the net charge of the protein at the given pH and retruns it as a float
        '''
        charge = (((10**self.aaNterm)/(10**self.aaNterm + 10**pH)) - ((10**pH)/(10**self.aaCterm + 10**pH))) # represents the charge of the protein at given pH
        for aa in self.aaCompDict: # iterates through the amino acid Composition dictionary 
            if aa in self.aa2chargePos: # The charges of the amino acids are totaled 
                charge += self.aaCompDict[aa]*((10**self.aa2chargePos[aa])/(10**self.aa2chargePos[aa] + 10**pH))
            elif aa in self.aa2chargeNeg:
                charge -= self.aaCompDict[aa]*((10**pH)/(10**self.aa2chargeNeg[aa] + 10**pH))

        return charge
    

    def molarExtinction (self):
        '''
        returns the Molar Extinction coeficeint  for the protein using the formula found by Gill and von Hippel 
        '''
        return self.aaCompDict['Y']*self.aa2abs280['Y'] + self.aaCompDict['W']*self.aa2abs280['W'] + self.aaCompDict['C']*self.aa2abs280['C']



    def massExtinction (self):
        '''
        Finds the mass extinction coefieceint for the protein by dividing its molar extinction coeficient by its molecular wieght
        '''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0



    def molecularWeight (self):
        '''
        Finds the molecular weight of the protein by summing the weight of each amino acid and subtracting by mass lost to water
        '''
        aaMW = self.mwH2O
        for aa in self.aaCompDict:
            aaMW += self.aaCompDict[aa]*(self.aa2mw[aa] - self.mwH2O)
                
        return aaMW



class NucParams:
    
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        
        inString = inString.replace(" ","").upper()
        self.codonSeq = [inString[i:i+3] for i in range(0, len(inString), 3)]  

        
    def addSequence (self, inSeq):
        '''
        Takes a codon sequence as a string and appends each codon to the codonSeq list as individual 3 character strings
        '''
        inSeq = inSeq.replace(" ","").upper() # removes whitespace and lowercase
        seqToAdd = [inSeq[i:i+3] for i in range(0, len(inSeq), 3)] # breaks up input into 3 character strings
        for i in range(0, len(seqToAdd)):
            self.codonSeq.append(seqToAdd[i])
        


    def aaComposition(self):
        '''
        Retruns a dictionary containing the count of each amino acid in the protein
        '-' represents codons that do not code for proteins (START / STOP codons)
        '''
        self.aaComp = {
        'A': 0,  'G': 0,  'M': 0, 'S': 0, 'C': 0,
        'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
        'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
        'W': 0,  'F': 0, 'L': 0, 'R': 0, 'Y': 0, 
        '-': 0
        }
        for codon in self.codonSeq:
            codon = codon.replace('T', 'U') # changes DNA to RNA
            if self.rnaCodonTable[codon] == '-': 
                self.aaComp[self.rnaCodonTable[codon]] += 1 # increments noncoding count
            elif codon in self.rnaCodonTable:
                self.aaComp[self.rnaCodonTable[codon]] += 1 # incremtns count of respective aa
          
        return self.aaComp
    
    def nucComposition(self):
        '''
        Returns a dictionary containing the count of each nucleic acid found in the codon sequence
        Takes DNA and RNA, counts number of 'N' as well 
        '''
        self.nucComp = {'A': 0, 'C': 0, 'G': 0,'U': 0, 'T': 0, 'N': 0} 
        for codon in self.codonSeq:  # iterates through codonSeq and counts each nucleic acid
            for nuc in codon:
                if nuc in self.nucComp:
                    self.nucComp[nuc] += 1
        return self.nucComp
      
    def codonComposition(self):
        '''
        much like aaComposition and nucComposition, codonComposition returns a dictionary containing the count of each codon
        '''
        self.codonComp =  {}
        for codon in self.codonSeq:
            codon = codon.replace('T', 'U') # changes all input to RNA
            if codon in self.rnaCodonTable:  # checks to see if the codon has been found, increments counter if already found and creates a new element if it has not been found yet
                if codon in self.codonComp:
                    self.codonComp[codon] += 1
                else:
                    self.codonComp[codon] = 1

        return self.codonComp

    def nucCount(self):
        '''
        returns the number of total number of nucleic acids in the sequence
        '''
        return sum(self.nucComposition().values())





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
                print(line)
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










########################################################################

'''
Psuedo-Code for OrfFinder

Generate reverse complement of Start/Stop codons
Iterate through the Fasta file storing every open reading frame in the form of orfs[rf1,...,rf6][orfs]

orfs are lists of positions where the first element of the list is a start codon if it is on the top strand, 
and a stop codon if its on the bottom strand, and vice versa. 
For example, [1, 6] means there is one start/stop codon at position 1 and one at position 6

Handle boundry cases where there is no start/stop codon to finish an orf by 
ending the orf at the boundry of the gene. 

Sort all reading frames first by length of the orf(high->low), then by initial position(low->high)
print each orf in the list, skipping any that dont meet the conditions of minGene length and longestGene
    

'''
########################################################################
def func(x):
    return (-x[0], x[2]) # storts orfs by length, then by start position, used by OrfFinder

def OrfFinder(CMDline):
    
    start = CMDline.args.start
    stop = CMDline.args.stop
    minGene = CMDline.args.minGene
    longestGene = CMDline.args.longestGene
    RCstart = [] 
    RCstop = []
    
    for i in range(0, len(start)):  # generates reverse complement of start and stop codons
        toAdd = ''
        for j in range(0, 3):
            if start[i][j] == 'A':
                toAdd += 'T'
                
            elif start[i][j] == 'T':
                toAdd += 'A'
                
            elif start[i][j] == 'G':
                toAdd += 'C'
                
            elif start[i][j] == 'C':
                toAdd += 'G'   
        #print(toAdd)
        RCstart.append(toAdd[::-1]) 

    for i in range(0, len(stop)):
        toAdd = ''
        for j in range(0, 3):
            if stop[i][j] == 'A':
                toAdd += 'T'
                
            elif stop[i][j] == 'T':
                toAdd += 'A'
                
            elif stop[i][j] == 'G':
                toAdd += 'C'
                
            elif stop[i][j] == 'C':
                toAdd += 'G'    
        #print(toAdd)   
        RCstop.append(toAdd[::-1]) 
    
    #print(start)
    #print(RCstart)
    #print(stop)
    #print(RCstop)
    reader = FastAreader()
    for head, seq in reader.readFasta():
        
        # 2d Lists that will hold each gene in the sequence as a 2 element list containing the position of the start and stop codon for the gene.
        #rf1-3 represent the three rfs on the top strand and 4-6 represent the three on the bottom strand
        # orderedOrfs will hold all of the rfs in another 2d list, with each list being of the form [length of rf, reading frame, start, stop]
        rf1 = [[]]
        rf2 = [[]] 
        rf3 = [[]]
        rf4 = []
        rf5 = []
        rf6 = []
        orderedOrfs = []
        # iterates through all headers and sequences, main loop of program
        '''
        The loop runs through the enitre sequence, looking at a 3character substring that shifts one charater each iteration
        The substring is checked to see if it is a start or stop codon
        if it is, the position is saved
        at the end of the loop rf1-6 will hold all orfs in thier respective rf. 
        The lists are of the form [start codon 1, start codon 2, ..., start codon n, stop codon], start codons that overlap are put into the same list
        for the reverse complement rfs the lists are of the form [stop, start 1, ..., start n]
        '''

       
        
        for i in range(0,len(seq)-1):
            if seq[i:i+3] in start: # places start codon in the apropreate reading frame
                rf = i%3 # determines reading frame
                if rf == 0:
                    rf1[-1].append(i+1)
                    
                elif rf == 1:
                    rf2[-1].append(i+1)
                    
                elif rf == 2:
                    rf3[-1].append(i+1)

            elif seq[i:i+3] in stop: # places stop codon in the apropreate reading frame
                rf = i%3 # determines reading frame
                if rf == 0: 
                    if len(rf1) > 0: 
                        if len(rf1[-1]) == 1: # removes redundant stop codons
                            rf1.pop(-1)
                        rf1[-1].append(i+3)
                        rf1.append([])
                    else:
                        rf1[-1].append(i+3)
                        rf1.append([])
                    
                elif rf == 1:
                    if len(rf2) > 0: 
                        if len(rf2[-1]) == 1: # removes redundant stop codons
                            rf2.pop(-1)
                        rf2[-1].append(i+3)
                        rf2.append([])
                    else:
                        rf2[-1].append(i+3)
                        rf2.append([])
                        
                elif rf == 2:
                    if len(rf3) > 0: 
                        if len(rf3[-1]) == 1: # removes redundant stop codons
                            rf3.pop(-1)
                        rf3[-1].append(i+3)
                        rf3.append([])
                    else:
                        rf3[-1].append(i+3)
                        rf3.append([])

            # Cases for reverse complement
            elif seq[i:i+3] in RCstart:
                rf = (len(seq) - i) % 3
                if rf == 0:
                    if len(rf4) > 0 and isinstance(rf4[-1], list):
                        rf4[-1].append(i+3)
                    else: # handles case where there is no stop codon before the start codon
                        rf4.append([1, i+3])
                        

                elif rf == 1:
                    if len(rf5) > 0 and isinstance(rf5[-1], list):
                        rf5[-1].append(i+3)
                    else: # handles case where there is no stop codon before the start codon
                        rf5.append([1, i+3])
                       

                elif rf == 2:
                    if len(rf6) > 0 and isinstance(rf6[-1], list):
                        rf6[-1].append(i+3)
                    else: # handles case where there is no stop codon before the start codon
                        rf6.append([1, i+3])
                        

            elif seq[i:i+3] in RCstop:
                rf = (len(seq) - i) % 3
                if rf == 0:
                    if len(rf4) > 0 and len(rf4[-1]) < 2: # removes redundant stop codons
                        rf4.pop(-1)
                    rf4.append([i+1])
                   
                    
                elif rf == 1:
                    if len(rf5) > 0 and len(rf5[-1]) < 2: # removes redundant stop codons
                        rf5.pop(-1)
                    rf5.append([i+1])
                   
                    
                elif rf == 2:
                    if len(rf6) > 0 and len(rf6[-1]) < 2: # removes redundant stop codons
                        rf6.pop(-1)
                    rf6.append([i+1])
                    
        


        rfs = [rf1,rf2,rf3,rf4,rf5,rf6] # list of reading frames to be iterated through
        
        # handles boundry cases for the top strand and puts rfs of desired length into orderrfs
        for i in range(0, 3):
            # case where seq has no start/stop codons
            if len(rfs[i][0]) == 0:
                rfs[i][0] = [1, len(seq)]
            # case where there is a start codon without a stop after it
            elif len(rfs[i][-1] )> 0:
                rfs[i][-1].append(len(seq))
            # case where there is a stop codon without a start before it
            elif len(rfs[i][0]) == 1:
                rfs[i][0] = [1, rfs[i][0][0]]
            
            for j in range(0, len(rfs[i])): # iterates through all orfs in the current rf
                '''
                if longestGene: # Takes the end of the lists to ensure longest gene in orf is chosen
                    if len(rfs[i][j]) > 1 and rfs[i][j][-1]-rfs[i][j][0] + 1 > minGene: # checks if gene is greater than desired length and Handles stop codons in between rfs
                        orderedOrfs.append([rfs[i][j][-1] - rfs[i][j][0] + 1, i+1, rfs[i][0], rfs[i][-1]])
            
                else:
                '''
                for k in range(0, len(rfs[i][j])-1): # iterates through each start and stop codon in the current orf
                    if len(rfs[i][j]) > 1 and rfs[i][j][-1]-rfs[i][j][k] + 1 > minGene: # checks if gene is greater than desired length and Handles stop codons in between rfs
                        orderedOrfs.append([rfs[i][j][-1] - rfs[i][j][k] + 1, i+1, rfs[i][j][k], rfs[i][j][-1]]) # adds gene to orderedOrfs
        
        # handles boundry cases for the reverse complement strand and puts rfs of desired length into orderrfs
        for i in range(3, 6):
            
            # case where seq has no start/stop codons
            if len(rfs[i]) == 0:
                rfs[i].append([1, len(seq)])
            # case where there is a stop codon with no start after it
            elif len(rfs[i][-1]) == 1:
                rfs[i][-1].append(len(seq))
            
            for j in range(0, len(rfs[i])): # iterates through all sets of start and stop codons in the current rf
                '''
                if longestGene: # Takes the end of the lists to ensure longest gene in orf is chosen
                        if len(rfs[i][j]) > 1 and rfs[i][j][-1]-rfs[i][j][0] + 1 > minGene: # checks if gene is greater than desired length and Handles stop codons in between rfs
                            orderedOrfs.append([rfs[i][j][-1] - rfs[i][j][0] + 1, i+1, rfs[i][0], rfs[i][-1]])
                else:
                '''
                for k in range(0, len(rfs[i][j])): # iterates through each start and stop codon in the current set
                    if len(rfs[i][j]) > 1 and abs(rfs[i][j][0]-rfs[i][j][k]) + 1 > minGene: # checks if gene is greater than desired length and Handles stop codons in between rfs
                        orderedOrfs.append([abs(rfs[i][j][0]-rfs[i][j][k]) + 1, -(i-2), rfs[i][j][0], rfs[i][j][k]]) # adds gene to orderedOrfs 
                        
        orderedOrfs.sort(key=func) # sorts the orfs first by length, then by starting position
        
        # writes finished rfs to stdout
        
        print(head+'\n')
        for rf in orderedOrfs:
            if rf[1] < 0:
                print('{:=} {:5d}..{:5d} {:5d}\n'.format(rf[1],rf[2],rf[3],rf[0]))
            else:
                print('+{:=} {:5d}..{:5d} {:5d}\n'.format(rf[1],rf[2],rf[3],rf[0]))
        
        
   