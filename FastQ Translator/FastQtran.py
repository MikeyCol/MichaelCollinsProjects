#!/usr/bin/env python3
# Name: Michael Collins 1603172
# Date: 6/7/21
# Partners: None
import sys
class FastQreader:

    ''' 
    Define objects to read FastQ files.
    
    instantiation: 
    thisReader = FastQreader ('testTiny.file')
    usage:
    for head, seq, qulity in thisReader.readFasta():
        print (head,seq,quality)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFastq (self):
        ''' Reads an entire FastQ record and returns the header squence and quality lines separately'''
        header = ''
        sequence = ''
        quality = ''

        with self.doOpen() as fileH:
            # These will hold the header, sequence and quality lines respectively
            header = '' 
            sequence = ''
            quality = ''
            
            # skip to first fastq header
            line = fileH.readline()
            while not line.startswith('@') :
                line = fileH.readline()
                
            header = line[1:].rstrip()
            sequence = ''.join(next(fileH).rstrip().split()).upper()
            sequence = sequence.replace('*', 'N')
            sequence = sequence.replace('.', 'N')

            for line in fileH:
                if line.startswith ('@'):
                    yield header,sequence,quality
                    header = line[1:].rstrip()
                    sequence = ''.join(next(fileH).rstrip().split()).upper()
                    sequence = sequence.replace('*', 'N')
                    sequence = sequence.replace('.', 'N')
                    quality = ''

                elif line.startswith('+'):
                    quality = ''.join(next(fileH).strip().split())
                    
                    
        yield header,sequence,quality

def FastQtran(CMDline):
    '''
    Determines the type of translation defined by the parameters, PHRED33 to PHRED64 for example
    Uses populateTransDict and translate to execute that translation
    '''
    inStyle = CMDline.args.inputStyle
    outStyle = CMDline.args.outputStyle

    
    reader = FastQreader()
    # These if statements determeine what type of input the input file is in, and then print out the translation 
    if len(inStyle) == 5: # narrows inStyle to P33in or P64in
        populateTransDict(int(inStyle[1:3]),int(outStyle[1:3])) # poplulates transDict with PHRED to PHRED with the choseen offsets
        for head,seq,quality in reader.readFastq():
            print('@'+head)
            print(seq)
            print('+'+head)
            print(translate(quality))
        
    elif inStyle == 'PBin': # handels case when input is set to PBin
        try:
            B = CMDline.args.bShift
        except ValueError:
            print('offset B not set\nUse -B = offset')
            quit()
        populateTransDict(B,int(outStyle[1:3])) # populates transDict with PHREDB to PHRED33 or 64 
        for head,seq,quality in reader.readFastq():
            print('@'+head)
            print(seq)
            print('+'+head)
            print(translate(quality))


    else : # Handles solexa to PHRED translation
        populateTransDict(int(outStyle[1:3]),-1)
        for head,seq,quality in reader.readFastq():
            print('@'+head)
            print(seq)
            print('+'+head)
            print(translate(quality))



import math
transDict = {} # Used in translate, contains key value pairs corresponding to the quality score translation needed
def populateTransDict(offset1, offset2):
    '''
    Takes two offest arguments(integers) and populates transDict with the appropreate translation, be it PHRED33 to PHRED64, solexa to PHRED33, etc
    '''
    if offset2 < 0:  # handles solexa to PHRED translations
        for i in range(-5,63):
            transDict[i+64] = int(round(10*math.log10(10**(i/10)+1))) + offset1

    # Two separate conidtionals are here to make sure extrainous dictionary elements are not created 
    elif offset2 < offset1: # handles PHRED to PHRED
        
        for i in range(offset2,127):
            value = i + offset2 - offset1
            if value > 126:
                value = 126
            elif value < 33:
                value = 33
            transDict[i] = value

    else: # also handles PHRED to PHRED
        for i in range(offset1,127):
            value = i + offset2 - offset1
            if value > 126:
                value = 126
            elif value < 33:
                value = 33
            transDict[i] = value
    
    return

def translate(quality):   
    '''
    Runs the quality line through the popluated translation dictionary, 
    replaceing each character with its corresponding value in the dictionary

    returns a string containing the translated quality line

    '''
    out = []
    for char in quality:
        newChar = transDict.get(ord(char))
        out.append(chr(newChar))

    return ''.join(out)
