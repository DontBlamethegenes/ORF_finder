class OrfFinder:
    '''Read a sequence and find open reading frames of all frames and print them out to a text file. Use list of lists
    of lists as initial container then a list of tuples as the final container.'''

    def __init__(self, start, stop, longest, geneLength, outFile):
        '''instantiates lists of starts, stops and distances sorted by list of frames'''
        self.frameList = [[[], [], []], [[], [], []], [[], [], []], [[], [], []], [[], [], []], [[], [], []]]
        self.startList = start
        self.stopList = stop
        self.longestGene = longest
        self.minGeneLength = geneLength
        self.outPutFile = outFile
    def findGenes(self, seq, head):
        '''Find valid start and stop codons of an open reading frame and calculate the length of Orf.'''
        n = 0  # keeps track of what frame is being read
        startCodons = self.startList  # List of start codons. default=['ATG'] edit: should work with ATG' | 'TTG' | 'GTG'
        stopCodons = self.stopList    # List of stop codons. default=['TAG', 'TGA', 'TAA']
        seqHead = head
        reverseSeq = seq[::-1]
        revCompSeq = ""
        for base in reverseSeq:
            if base == "A":
                revCompSeq += "T"
            elif base == "T":
                revCompSeq += "A"
            elif base == "G":
                revCompSeq += "C"
            elif base == "C":
                revCompSeq += "G"

        while n <= 5:  # loops through frames
            if n <= 2:  # makes a list of codons in the sequence
                z = n
                listCodons = [seq[(i + z):(i + z) + 3] for i in
                              range(0, len(seq), 3)]
            elif n >= 3:  # makes a list of codons in the reverse complement sequence
                z = n-3
                listCodons = [revCompSeq[(i + z):(i + z) + 3] for i in
                              range(0, len(seq), 3)]
            if n >= 3:
                listStart = [len(seq)-(z)]
            else:
                listStart = [(n+1)]
            codonNumber = 0  # keeps track of codon being read
            for codon in listCodons:
                if codon in startCodons:

                        if n >= 3:
                            if len(seq) - (z + codonNumber * 3) != len(seq)-(z):# fixes issue that happens when seq starts with ATG
                                listStart.append(len(seq) - (z + codonNumber * 3))
                        else:
                            if (z + codonNumber * 3 + 1) != (n + 1):  # fixes issue that happens when seq starts with ATG
                                listStart.append(z + codonNumber * 3 + 1)  # adds location of start codon to list

                elif codon in stopCodons and len(listStart) != 0:
                    if n >= 3:  # Handles negative frames
                        if self.longestGene:  # Only saves the longest Orf per stop codon.
                                longestStartIndex = listStart[0]
                                stopIndex = len(seq) - (z + codonNumber * 3 + 2)
                                lengthOfOrf = 1 + longestStartIndex - stopIndex

                                self.frameList[n][1].append(longestStartIndex)
                                self.frameList[n][0].append(stopIndex)
                                self.frameList[n][2].append(lengthOfOrf)
                        else:
                            for start in listStart:
                                    stopIndex = len(seq) - (z + codonNumber * 3 + 2)

                                    self.frameList[n][1].append(start)
                                    self.frameList[n][0].append(stopIndex)
                                    self.frameList[n][2].append(1 + start - stopIndex)
                    else:
                        if self.longestGene:
                            if listStart[0] == (n+1) and len(listStart) == 1:
                                self.frameList[n][0].append(1)
                                self.frameList[n][1].append(z + codonNumber * 3 + 3)
                                self.frameList[n][2].append((z + codonNumber * 3 + 3))
                            else:
                                self.frameList[n][0].append(listStart[0])
                                self.frameList[n][1].append(z + codonNumber * 3 + 3)
                                self.frameList[n][2].append((z + codonNumber * 3 + 3 - listStart[0] + 1))
                        else:
                            for start in listStart:

                                self.frameList[n][0].append(start)
                                self.frameList[n][1].append(z + codonNumber*3 + 3)
                                self.frameList[n][2].append((z + codonNumber*3 + 3)-start + 1)
                    listStart = []
                codonNumber += 1
            if len(listStart) != 0:
                for start in listStart:
                    if n >= 3:
                        if start == (len(seq)-(z)):
                            self.frameList[n][1].append(len(seq))
                            self.frameList[n][2].append(len(seq))
                        else:
                            self.frameList[n][1].append(start)
                            self.frameList[n][2].append(start)
                        self.frameList[n][0].append(1)
                    else:
                        if start == (n+1):
                            self.frameList[n][0].append(1)
                            self.frameList[n][2].append(len(seq))
                        else:
                            self.frameList[n][0].append(start)
                            self.frameList[n][2].append(len(seq) - start + 1)
                        self.frameList[n][1].append(len(seq))


            n += 1
        self.printList(seqHead)
        self.frameList = [[[], [], []], [[], [], []], [[], [], []], [[], [], []], [[], [], []], [[], [], []]]

    def printList(self, seqHead):
        '''Print start/stop index, length and frame of each orf to a designated txt file'''
        listOfOrfs = []
        x = 0
        while x <= 5:  # itterates through 6 frames
            mappedValues1 = zip(self.frameList[x][0], self.frameList[x][1], self.frameList[x][2])  # makes a list of tuple
            for values in mappedValues1: # adds corresponding frame to each tuple
                if x >= 3:
                    values = values + (-(x-2),)
                elif x <= 2:
                    values = values + ((x+1),)
                listOfOrfs.append(values)
            x += 1
        finalList = sorted(listOfOrfs, key=lambda x: (-x[2], x[3]), )  # sorts list by length and frame
        file1 = open(self.outPutFile, "a")  # opens file in append mode
        file1.write("{}{}".format(seqHead, "\n"))
        for tup in finalList:
            if tup[2] >= self.minGeneLength: #sets the minimum Orf length
                if tup[3] > 0:
                    file1.write('{:<1s}{:<6d}{:<4d}{:<1s}{:>9d}{:>9d}{}'.format("+", tup[3], tup[0], "..", tup[1], tup[2], "\n"))
                else:
                    file1.write('{:<7d}{:<4d}{:<1s}{:>9d}{:>9d}{}'.format(tup[3], tup[0], "..", tup[1], tup[2], "\n"))
        file1.close()