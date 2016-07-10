#!/usr/bin/env python3

#By Michael Weinstein, 2016... still testing this version.  Please credit if used.

class TwoBitDNA(object):
    
    def __init__(self):
        import sys
        self.byteOrder = sys.byteorder
        self.dataStructureTemplate = [["Q", 64, 0],   #lists available data structures and their bitlengths, ORDER LARGE TO SMALL
                                      ["I", 32, 0],
                                      ["H", 16, 0],
                                      ["B", 8, 0]]
        self.dataStructureBitMap = {}
        for structure in self.dataStructureTemplate:  #building a map of structures to bitlengths based on the above
            self.dataStructureBitMap[structure[0]] = structure[1]
        self.encodingKey = {"A":"00","C":"01","G":"10","T":"11"}  #mapping base to bit value
        self.decodingKey = {}  
        for base in list(self.encodingKey.keys()):  #reverse map of the above
            self.decodingKey[self.encodingKey[base]] = base
        self.cachedEncodingData = {}
        
    def validSeq(self, dna):   #Makes sure that the DNA seq submitted is upper case and only valid letters.  4-bit coding of degenerate sequence may be available in the future, but not planned.
        for letter in dna:
            if letter not in "ATGC":
                return False
        return True
    
    def encode(self, dna):  #turns DNA sequence into a 2 bit encoded value returned as an integer
        dna = dna.upper()  #sequence is handled as upper case
        if not self.validSeq(dna):
                raise RuntimeError("DNA sequence for 2 bit conversion must only consist of ATGC")
        packBits = self.calculatePackBits(dna)
        bitString = ""
        for letter in dna:
            bitString = self.encodingKey[letter] + bitString   #builds the binary as a text string of 1 and 0 with the first bases in the least significant digit
        return int(bitString, 2)  #return the value as an integer... this WILL lose any A values on the end, so we will need to know the intended sequence length and restore those values
    
    def calculatePackBits(self, dna):  #figures out how many pack bits will be needed to get the total bitlength up to a multiple of 8
        seqBits = self.calculateSeqBits(dna)  #this will work either for int or seq and give the appropriate result, no need to check type
        if seqBits % 8 == 0:  #keeps us from getting back an 8 when we do not need packing bits
            return 0
        else:
            return 8 - (seqBits % 8)
                    
    def calculateSeqBits(self, dna):  #calculates how many bits will be needed to encode the DNA itself
        try:
            dna = int(dna)
            return dna * 2
        except ValueError:
            pass
        return len(dna) * 2
    
    def calculateBitLength(self, dna): #determines how many bits will be required based on sequence length (counted from a passed sequence or an integer) and returns a tuple of how many sequence bytes are needed and how many extra packing bytes are needed to get it to a multiple of 8 for storage
        try:
            dna = int(dna)
            dna = "T" * dna
        except ValueError:
            pass
        if len(dna) in self.cachedEncodingData:
            return (self.cachedEncodingData[len(dna)]["seqBits"], self.cachedEncodingData[len(dna)]["packBits"])
        return(self.calculateSeqBits(dna), self.calculatePackBits(dna))
    
    def getPackingPattern(self, dna):  #this will generate the string for packing the data into a struct library data type
        try:
            dna = int(dna)
            dna = "T" * dna
        except ValueError:
            pass
        if len(dna) in self.cachedEncodingData:  #if we already did the work, use the cached pattern
            return self.cachedEncodingData[len(dna)]["packPattern"]  #use work we already did
        seqBits, packBits = self.calculateBitLength(dna)
        totalBitLength = seqBits + packBits
        remainingBitLength = totalBitLength
        dataStructures = self.dataStructureTemplate.copy()  #getting a copy of the datastructures template.  Using a copy to change the counts at the last index
        for i in range(0, len(dataStructures)):  #note that by this point, we should be an even multiple of 8 bits due to packing bits added and remainders should not be a problem here
            if not remainingBitLength:
                break   #if we have accounted for all bits needed, stop trying to account for them
            structureBits = dataStructures[i][1]
            if remainingBitLength >= structureBits:
                dataStructures[i][2] = remainingBitLength // structureBits  #the only time we should ever be using a structure type more than once is for the largest structure for very long DNA sequences.  The length should be a multiple of 8 and our structures are 8*8, 8*4, 8*2, and 8*1
                remainingBitLength -= dataStructures[i][2] * structureBits
        packPattern = ""
        for structure in self.dataStructureTemplate:  #build up the actual packing pattern string for a struct object
            packPattern += structure[0] * structure[2]
        self.cachedEncodingData[len(dna)] = {"packPattern":packPattern, "packBits":packBits, "seqBits":seqBits}  #cache work we just did
        return packPattern
                
    def pack(self, dna):  #the business end of the process... use bit shifting to get bits of interest and subtracting to remove segments already added
        dna = dna.upper()
        if not self.validSeq(dna):
                raise RuntimeError("DNA sequence for 2 bit conversion must only consist of ATGC")
        values = []  #this will ultimately return a list of int values that should correspond to the packing pattern
        packingPattern = self.getPackingPattern(dna) #if we already cached it, this will be much faster, otherwise we cache it now
        seqBits, packBits = self.calculateBitLength(dna)
        encodedDna = self.encode(dna)
        encodedDna = encodedDna << packBits  #add in the packing bits at the least significant end so they can be shifted out later
        remainingEncoded = encodedDna  #tracking the bits we have left to account for
        for i in range(0, len(packingPattern)):
            bitsNotOfInterest = 0  #counting bits at the right (least significant) end of the pattern that will not be added to the pack this iteration
            for structure in packingPattern[i+1:]:  #get the packing pattern structures we will NOT be handling this iteration
                bitsNotOfInterest += self.dataStructureBitMap[structure]  #and add their bit lengths to the bits not of interest
            value = remainingEncoded >> bitsNotOfInterest   #the value to add for this iteration is the bits of interest with the ones not of interest shifted off
            remainingEncoded -= value << bitsNotOfInterest  #now we chop off the bits we just accounted for by shifting the value bits back into position and subtracting them out (this will zero out the most significant bits that were added)
            values.append(value)  #and we add the new value to our packing list
            # print(bin(encodedDna)[2:] + "|")  #these may be commented out, but they just allow for tracking of the pack building on STDOUT
            # for number in values:
            #     print(bin(number)[2:] + "|", end = "")
            # print(bin(remainingEncoded)[2:])
        return values
            
    def unpack(self, values, seqLength):  #reverse of the above, returns the integer representation of the DNA.  Since A bases (00's) at the most significant end may have been lost, we need to know the intended sequence length in order to get them back
        packingPattern = self.getPackingPattern(seqLength)
        seqBits, packBits = self.calculateBitLength(seqLength)
        if not len(values) == len(packingPattern):  #this should never happen, and if it does it should raise the alarm bells
            raise RuntimeError("Number of values and packing pattern must be the same length")
        encodedDna = 0  #this just initializes an integer for us to add the packed parts to
        for i in range(0, len(packingPattern)):
            bitsNotAdded = 0  #tracking the trailing bits on the lest significant end that will be covered by other packed values on a later iteration
            for structure in packingPattern[i+1:]:
                bitsNotAdded += self.dataStructureBitMap[structure]  #adding the bit length for each packed value we will not be handling with this iteration so we can bit shift appropriately
            encodedDna += values[i] << bitsNotAdded  #bit shifting the value being handled on this current iteration into the appropriate position and adding it
        encodedDna = encodedDna >> packBits  #shifting off any packing bits that may have been added
        return encodedDna
    
    def decode(self, encodedDna, seqLength):  #can take either the unpacked single int or the packed list
        if type(encodedDna) == list:  #if the passed value was a list of packed values directly from a the struct, unpack them
            encodedDna = self.unpack(encodedDna, seqLength)
        seqBits, packBits = self.calculateBitLength(seqLength)
        binaryDna = bin(encodedDna)[2:]  #gets the binary value as a string without the first two characters that will be a prefix of no significance
        missingLeadingZeros = seqBits - len(binaryDna)  #figure out how many (if any) A bases (00's) were lost from the most significant end in order to replace them
        binaryDna = "0" * missingLeadingZeros + binaryDna  #actually replace them
        humanReadableSequence = ""
        for i in range(len(binaryDna) - 2, -1, -2):  #take binary values from the string as pairs starting from the high index numbered end
            humanReadableSequence += self.decodingKey[binaryDna[i:i+2]]  #map the binary (as a string) back to the letter and add it
        return humanReadableSequence
    
    def countMismatches(self, seq1, seq2, shortestSeqLength, endClip = 0): #can either take the unpacked sequence (preferred) or a list or tuple containing packed seq and seq length in that order
        if type(seq1) == list or type(seq1) == tuple:
            seq1 = self.unpack(seq1[0],seq1[1])
        if type(seq2) == list or type(seq2) == tuple:
            seq2 = self.unpack(seq2[0],seq2[1])  #things get tricky now because we may have lost any leading zeros (A's) from the encoding
        variations = seq1 ^ seq2  #doing an xor on the two sequences will generate one or two 1 values where they differ and a pair of zeros where they match
        #peekvar = bin(variations)
        checkBits = int("11", 2)  #just a binary 3f, but by shifting it over two positions, we can do an AND operation that will return at least a single 1 value (evaluates to True) if there was a mismatch between the sequences at that doublet
        basesToCheck = shortestSeqLength - endClip
        mismatches = 0
        if not variations: #handles the rare perfect match quickly
            return 0
        for i in range(0, basesToCheck):  #will check either the entire length of the shortest of the two sequences or that length minus the endClip
            #peekCheck = bin(checkBits)
            #peekRes = bin(variations & checkBits)
            if variations & checkBits:  #this does the AND operation on the moving check bits and variations to test each bit doublet for mismatching
                mismatches += 1
            checkBits = checkBits << 2
        return mismatches
    
    def withinMismatchTolerance(self, seq1, seq2, shortestSeqLength, mismatchTolerance, endClip = 0): #does the exact same thing as mismatch counting, except this will stop if a mismatch threshold is reached and return False, otherwise it will return True.  Optimized operation for filtering.
        if type(seq1) == list or type(seq1) == tuple:
            seq1 = self.unpack(seq1[0],seq1[1])
        if type(seq2) == list or type(seq2) == tuple:
            seq2 = self.unpack(seq2[0],seq2[1])  
        variations = seq1 ^ seq2
        #peekvar = bin(variations)
        checkBits = int("11", 2)  #just a binary 3
        basesToCheck = shortestSeqLength - endClip
        mismatches = 0
        if not variations: #handles the rare perfect match quickly
            return True
        for i in range(0, basesToCheck):
            #peekCheck = bin(checkBits)
            #peekRes = bin(variations & checkBits)
            if variations & checkBits:
                mismatches += 1
                if mismatches > mismatchTolerance:
                    return False
            checkBits = checkBits << 2
        return True
    
#begin testcode for debugging
seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
seq2 = "ATGTTCCTACGTACGTACCT"
twobit = TwoBitDNA()
encoded = twobit.encode(seq)
packingPattern = twobit.getPackingPattern(seq)
pack = twobit.pack(seq)
unpacked = twobit.unpack(pack, len(seq))
decoded = twobit.decode(unpacked, len(seq))
encoded2 = twobit.encode(seq2)
mismatches = twobit.countMismatches(encoded, encoded2, len(seq2))
withinTolerance = twobit.withinMismatchTolerance(encoded, encoded2, len(seq), 3, 3)


print("something")

