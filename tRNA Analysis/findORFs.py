import sequenceAnalysis as seqA
from operator import itemgetter



class CommandLine() :
    '''
    Command Line program

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of ob-1ect instantiated from CommandLine.
    For example, if myCommandLine is an ob-1ect of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - Finds orfs in a Fasta file', 
                                             epilog = 'Program epilog - ', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an rf')
        self.parser.add_argument('-mG', '--minGene', type=int,  default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)



def main(inCL=None):
   
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)

    seqA.OrfFinder(myCommandLine)

    
if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN