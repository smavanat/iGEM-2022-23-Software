import json
import string

#File for storing output
outputFile = open("miRNA.txt", "w")

#Dictionary for storing user input
miRNADict = {}

#dictionary defining complimentary base pairing (no wobble pairing present)
baseDict = {"A":"U", "C":"G"}

#List of bases
baseList = ["A", "U", "C", "G"]

#list of and regions
andRegions = ("AAAAAA", "UUUUUU", "CCCCCC", "GGGGGG")

#computes complimentary strand 
def compStrand(rnaStrand):
    compliment = ""
    for i in rnaStrand:
        for key, value in baseDict.items():
            if i == key:
                compliment += value
            elif i == value:
                compliment += key
    return compliment

#reverses string
def reverse(rnaStrand):
  return rnaStrand[::-1]


#computes and gate
def andGateGen(sequenceDict, startIndex):#The way we do it recursively is by asking for a start index as an input, then replace the indexes of the dictionary with index, index+1
    endIndex = len(sequenceDict)
    if(startIndex == endIndex -1):
        print("Done")
    else:
        #splitting each miRNA
        miRNA1, miRNA2 = compStrand(list(sequenceDict.values())[startIndex] [:len(list(sequenceDict.values())[startIndex+1])//3]), compStrand(list(sequenceDict.values())[1] [:len(list(sequenceDict.values())[1])//3])
        print("length of spliced Strand1: ", len(miRNA1))
        print("length of spliced Strand2: ", len(miRNA2))
        outputFile.write(f'\nlength of spliced Strand1: {len(miRNA1)}')
        outputFile.write(f'\nlength of spliced Strand2: {len(miRNA2)}')

        #index is incrimented later to update the andRegion such that its adjacent base doesn't match itself
        index = 0

        #this sets the default value of the andRegion to AAAAAA
        andRegion = andRegions[index] 

        #checking to see if the andRegion borders an identical base and changing it to correct this by incrimenting the index
        #this only needs to be done twice as in the "worst case scenario", ALL of the below code will run, making the anndRegion different to BOTH of its neighboring bases
        if (miRNA1[-1] == andRegion[1]):
            index += 1
            andRegion = andRegions[index] 
            if (miRNA2[0] == andRegion[1]):
                index += 1
                andRegion = andRegions[index]
        print(miRNA1 + andRegion + miRNA2)
        outputFile.write("\n" + miRNA1 + andRegion + miRNA2)
        andGateGen(miRNADict, startIndex+ 1)
        return((miRNA1 + andRegion + miRNA2))


#Ensures correct inputs
def checkCorrectCharacters(characterList, userInput:string):
    userInput = userInput.replace(" ","")
    for element in userInput:
        if(characterList.count(element)==0):
            break
        else: 
            return True
    return False


#Handles user input
def start():
    print("Input the number of miRNA strands")
    numMiRNA = int(input())
    i = 0
    while(i < numMiRNA):
        print("input the miRNA and its sequence")
        miRNA = input()#User input must be in the same form as this: name,base sequence. Base sequence needs to be in all caps
        index = miRNA.split(",")
        if(checkCorrectCharacters(baseList, index[1])):
            miRNADict.update({index[0]: index[1]})
            i+=1
        else:
            print("Incorrect input")
    json.dump(miRNADict, outputFile)
    andGateGen(miRNADict, 0)
start()
