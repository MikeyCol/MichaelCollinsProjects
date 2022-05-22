#!/usr/bin/env python3 
# Name: Michael Collins 1603172
# Date: 6/7/21
# Partners: None
import FastQtran as FQT
class FastQin() :
    '''
    Command Line program that runs FastQtran with specific outputs and inputs

    attributes:
    -in --inputStyle
        This defines the type of input by type of quality score
        P33in represents PHRED33 (33 offset), this is the defualt option
        P664in represents PHRED64 (64 offset)
        PBin repersetns a PHRED quality encoding with an offset of B which is determined by the parameter -B --bShift
        P64SOlin repersetns the solexa encoding (64 offset)

    -out --outStyle
        defines the type of output by quality score
        P33out gives an output in PHRED33 quality score encoding, this is the defualt option
        P64out gives an output in PHRED64 quality score encoding

    -B --bShift
        defeines the offset of B if PBin is used 
    
    EX:
    ./FastQin -in=P64SOlin -out=P33out <input.file> output.file

    '''
    
    def __init__(self, inOpts=None) :
        '''
        Command line program that handles input/output for FastQtran.py
        
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - ', 
                                             epilog = 'Program epilog - ', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input> output'
                                             )
        self.parser.add_argument('-in', '--inputStyle', choices = ['P33in', 'P64in','PBin','P64SOlin'], default='P33in', help = '')
        self.parser.add_argument('-out', '--outputStyle', choices = ['P33out', 'P64out'], default = 'P33out', help = '')
        self.parser.add_argument('-B', '--bShift', type = int, action = 'store', default = 0, help = '')
       
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)



def main(inCL=None):
    if inCL is None:
        myCommandLine = FastQin()
    else :
        myCommandLine = FastQin(inCL)

    FQT.FastQtran(myCommandLine)
    

    
if __name__ == "__main__":
    main() 