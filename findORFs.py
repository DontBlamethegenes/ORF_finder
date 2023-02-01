#!/usr/bin/env python3
# Name: Jose Alberto Aldapa (jaldapa)
# Group Members: “None”
'''
Design:
-read one sequence at a time
-instantiate a list of frames(-3,-2,-1,2,3), each with 3 lists of start stop and length of gene
-while loop from 0 to 2(n) to read frames one through 3
    -instantiate a list of startCodons
    -make list of codons starting from position (n) called seq
    -for codon in seq
        -codonNumber = 0
        -if  ['ATG' | 'TTG' | 'GTG']
            -push n + (codonNumber times 3) + index of codon to first list of reading frame 1
        -elif  ['TAA' | 'TGA' | 'TAG']
            -push location of last character of stop codon
            -calculateDistance() and return to list of stops in frame (n+1)
            -erase values from startCodons
Once list of list is filled out, a print method is used to make sence of data and organize it
into a list of tuples, each tuple being an open reading frame.
'''
import sequenceAnalysis

########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('inFile', action='store', help='input file name')
        self.parser.add_argument('outFile', action='store', help='output file name')
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), default=100,
                                 action='store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-t', '--stop', action='append', default=['TAG', 'TGA', 'TAA'], nargs='?',
                                 help='stop Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


########################################################################
# Main
# Here is the main program
#
#
########################################################################



def main(inCL=None):
    '''
    Find genes from a fasta file input and output to a designated txt file
    using command line class and the myOrf class
    '''
    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    print(myCommandLine.args)  # prints the arguments in command line
    inFile = myCommandLine.args.inFile  # has the input file name
    outFile = myCommandLine.args.outFile  # has the output file name
    longestGene = myCommandLine.args.longestGene  # is  True if only the longest Gene
    start = myCommandLine.args.start  # is a list of start codons
    stop = myCommandLine.args.stop  # is a list of stop codons
    minGene = myCommandLine.args.minGene# is the minimum Gene length to include

    myReader = sequenceAnalysis.FastAreader(inFile)
    myOrf = sequenceAnalysis.OrfFinder(start, stop, longestGene, minGene, outFile)
    open(outFile, 'w').close()
    for head, seq in myReader.readFasta():
        myOrf.findGenes(seq, head)





#######

if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN

