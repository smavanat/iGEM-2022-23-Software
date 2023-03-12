import json
import string
import textwrap

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

def splitString(miRNAInput):
    storageArray = []
    miRNAInput = miRNAInput.replace(" ","")
    n = len(miRNAInput)//3
    part1 = miRNAInput[:n]
    part2 = miRNAInput[2*n:]
    storageArray.append(compStrand(part1))
    storageArray.append(compStrand(part2))
    return storageArray

#computes and gate
def andGateGen(sequenceDict):#The way we do it recursively is by asking for a start index as an input, then replace the indexes of the dictionary with index, index+1
    #endIndex = len(sequenceDict)
    storageDict = sequenceDict.values()
    storageArray = []
    miRNAPartArray = []

    for element in storageDict:
        storageArray = splitString(element)
        miRNAPartArray.append(storageArray[0])
        miRNAPartArray.append(storageArray[1])

    del miRNAPartArray[0]
    del miRNAPartArray[len(miRNAPartArray)-1]

    for element in miRNAPartArray: 
        if(miRNAPartArray.index(element) % 2==0):

            #index is incrimented later to update the andRegion such that its adjacent base doesn't match itself
            index = 0

            #this sets the default value of the andRegion to AAAAAA
            andRegion = andRegions[index]

            miRNA1 = miRNAPartArray[miRNAPartArray.index(element)]
            miRNA2 = miRNAPartArray[miRNAPartArray.index(element) + 1]
            #checking to see if the andRegion borders an identical base and changing it to correct this by incrimenting the index
            #this only needs to be done twice as in the "worst case scenario", ALL of the below code will run, making the anndRegion different to BOTH of its neighboring bases
            if (miRNA1[-1] == andRegion[1]):
                index += 1
            andRegion = andRegions[index] 
            if (miRNA2[0] == andRegion[1]):
                index += 1
                andRegion = andRegions[index]
            print(miRNA1 + andRegion + miRNA2)
        else:
            continue
    
    print("Done")
#        return((miRNA1 + andRegion + miRNA2))


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
    andGateGen(miRNADict)
start()