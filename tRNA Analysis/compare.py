import fileinput
import sys
'''
The Text class houses the compare function which compares two text files for equality.

The program prints the first inconsistent line the text files and the line number, 
ex: 
File1:The dog walked down the road. 

File2:The dog ran down the street. 

At line: 6
'''
class Text:

    def __init__ (self, fname1='', fname2='' ):
        '''contructor: saves attribute fname '''
        self.fname1 = fname1
        self.fname2 = fname2

    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def compare (self):
        with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2:
            i = 0;
            for l1, l2 in zip(f1, f2):
                i += 1
                if l1 != l2:
                    return 'File 1:\n' + l1 + '\n' + 'File 2:\n' + l2 + '\nAt line: {0}'.format(i)
        return 'No inconsistencies found'


def main(inCL=None):
   texts = Text()
   print(texts.compare())
    

if __name__ == "__main__":
    main()  
