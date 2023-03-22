import string
import math as m
from nupack import *
import matplotlib.pyplot as plt

#defining the toehold switch complex analysis model
thsModel = Model(material = "rna", celsius = 37)

#Filse for storing output
outputFile = open("triggerComplex.txt", "w")
outputFile2 = open("toeholdSwitch.txt", "w")

#Dictionary for storing user input
miRNADict = {}

#dictionary defining complimentary base pairing (no wobble pairing present)
baseDict = {"A":"U", "C":"G"}

#List of bases
baseList = ["A", "U", "C", "G"]

#list of and regions
andRegions = ("AAAAAA", "UUUUUU", "CCCCCC", "GGGGGG")

#complex structure in DU+
complexStructure = []

#strands to be put together in the file
strands = []

#list of domain names for NUPACK
domains = list((string.ascii_lowercase))

#list of complex.seq
seq = []

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

#splits each miRNAInput into thirds and adds those to the storageArray
def splitString(miRNAInput):
    storageArray = []
    miRNAInput = miRNAInput.replace(" ","")
    n = round(len(miRNAInput)//3)
    part1 = miRNAInput[:n]
    part2 = miRNAInput[2*n+1:]
    part3 = miRNAInput[n:2*n+1]#modified
    storageArray.append(compStrand(part1))
    storageArray.append(compStrand(part2))
    storageArray.append(part3)#Modified
    return storageArray

#computes and gate
def andGateGen(sequenceDict):
    storageDict = sequenceDict.values()
    storageArray = []
    miRNAPartArray = []
    global unpairedMiddleArray
    unpairedMiddleArray = []#modified
 
    for element in storageDict:
        storageArray = splitString(element)
        miRNAPartArray.append(storageArray[0])
        miRNAPartArray.append(storageArray[1])
        unpairedMiddleArray.append(storageArray[2])#modified

    #makes sure code runs on the end of 1 strand and the beginning of the other
    del miRNAPartArray[0]
    del miRNAPartArray[len(miRNAPartArray)-1]
    middleIndex = 0
    for element in miRNAPartArray:

        if(miRNAPartArray.index(element) % 2==0):
            #index is incrimented later to update the andRegion such that its adjacent base doesn't match itself
            index = 0

            #this sets the default value of the andRegion to AAAAAA
            andRegion = andRegions[index]

            global miRNA1, miRNA2

            miRNA1 = reverse(miRNAPartArray[miRNAPartArray.index(element)])

            miRNA2 = miRNAPartArray[miRNAPartArray.index(element) + 1]            

            #defining a repeating unit for the complex
            duplex1 = "D" + str(len(miRNA1)) + "+"
            complexAnd = "U" + str(len(andRegion))
            duplex2 = "D" + str(len(miRNA2)) + "+"
            unpairedMiddle = "U" + str(len(unpairedMiddleArray[middleIndex]))#modified

            repeatUnit = [duplex1, complexAnd, duplex2, unpairedMiddle]       
            
            
            andGateDict = {}
            andGateArray = []
            for gate in andRegions:
                if (reverse(miRNA2)[0] and miRNA1[-1] == gate[1]):
                    energy = structure_energy(strands = [gate, compStrand(gate)], 
                                             structure = "D6+", model = thsModel)
                    andGateArray.append(energy)
                    andGateDict[energy] = gate
            andGateArray = sorted(andGateArray)
            andRegion = andGateDict[andGateArray[0]]
            
            andGate = miRNA1 + andRegion + reverse(miRNA2)

            strands.append(andGate)
            middleIndex+=1#modified

            #adding repeat unit to an empty array
            complexStructure.extend(repeatUnit)
        else:
            continue
    print("Done")

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
    andGateGen(miRNADict)
start()

unpairedStartLength = len(list(miRNADict.values())[0]) - round((len(list(miRNADict.values())[0])//3))

#this was the annoying thing I was talking about ages ago where due to divisibility by 3 resulting in a number that may round up or down depending on the remainder. fixed by brute force
if (unpairedStartLength % 3 == 2):
    unpairedStartLength = unpairedStartLength + 1 
elif (unpairedStartLength % 3 == 1):
    unpairedStartLength = unpairedStartLength - 1
else:
    unpairedStartLength = unpairedStartLength

#converting unpairedStart and End to DU+ notation and adding to the beginning and end of complexStructure respectively
unpairedStart = "U" + str(unpairedStartLength)
unpairedEnd = "U" + str(len(list(miRNADict.values())[-1]) - len(miRNA2))
complexStructure[-1] = unpairedEnd
complexStructure.insert(0, unpairedStart)

#adds the base sequences of each miRNA and AND gate to strands to be used later
for value in miRNADict.values():
    strands.append(value)

domainCodes = []

#this writes the text to be put in the txt file for assigning domains to the DU+ structure, eg
#domain a = AUCGACGACUGCCGUA
#domain b = UUCGUAGCUAGCUGA etc...
for strand, domain in zip(strands, domains):
    domainCode = "domain " + domain + " = " + strand
    domainCodes.append(domainCode)

#this takes letters of the alphabet in order until there is 1 for each miRNA and AND gate
domainsUsed = domains[0:len(domainCodes)]

#the domainsUsed list is split into 2 halves, one with letters to be assigned to miRNAs and the other to AND gates
mid = int((len(domainsUsed) + 1)/2)
miRNAs = list(domainsUsed[mid-1:])
andGates = list(domainsUsed[0:mid-1])

seq.extend(miRNAs)

integer = 0

#reordering the miRNAs and AND gates to be in a pattern of: miRNA, AND gate, miRNA, AND gate, miRNA etc...
while (integer < len(andGates)):
    for domain in seq:
        if((seq.index(domain)%2 ==0) and (seq.index(domain) + 1 < len(seq))):
            seq.insert(seq.index(domain) + 1, andGates[integer])
            integer += 1
        else:
            continue

#writes all the calculated info to a text file
outputFile.write(f'# miRNA complex structure in DU+ notation \nstructure complex = {"".join(complexStructure)}\n')
outputFile.write(f'\n# base sequences for each RNA strand')
for code in domainCodes:
    outputFile.write(f"\n" + code)
outputFile.write(f'\n \n# mapping RNA strands to the complex')
outputFile.write(f'\ncomplex.seq = ')
for sequence in seq:
    outputFile.write(sequence + " ")

#toehold switch stuff
#assigning start and end of toeholdRegion
toeholdStart = list(miRNADict.values())[0][0:unpairedStartLength]
toeholdEnd = list(miRNADict.values())[-1][len(miRNA2):]

toeholdAndGates = []

#every strand but the first and last one which are already accounted for in toeholdStart and End
middleStrands = strands[0:-1]
del middleStrands[len(middleStrands)//2]

#checks the rest of the strands to see if they are and gates
for loopStrand in middleStrands:
    count = 1
    for char in range(0, len(loopStrand)):
        if (loopStrand[char - 1] == loopStrand[char]):
            count+=1
            if (count == 6):
                toeholdAndGates.append(count*loopStrand[char])
            else:
                continue
        else:
            count = 1

del unpairedMiddleArray[0]
del unpairedMiddleArray[-1]
 
integer = 0

#reordering sequences to form the toehold region
while (integer < len(unpairedMiddleArray)):
    for domain in toeholdAndGates:
        if((toeholdAndGates.index(domain)%2 ==0) and (toeholdAndGates.index(domain) + 1 < len(toeholdAndGates))):
            toeholdAndGates.insert(toeholdAndGates.index(domain) + 1, unpairedMiddleArray[integer])
            integer += 1
        else:
            continue

toeholdRegion = []

#forming the toehold region
toeholdRegion.insert(0, toeholdStart)
toeholdRegion.extend(toeholdAndGates)
toeholdRegion.append(toeholdEnd)
toeholdRegion = "".join(str(i) for i in toeholdRegion)

#this is the unpaired region of the toehold switch where the miRNA complex binds
toehold = round((len(list(miRNADict.values())[0])//2))

#the paired section of the toehold switch that unzips until the start codom
thsDuplex = len(toeholdRegion) - toehold

#the entire switch structure in DU+ including the hairpin loop
thsStructure = ("U"+ str(toehold) + "D"+ str(thsDuplex) + "(U3D5(U15)U3)U21")

#sequence of the toehold switch
toeholdSwitchSequence = (compStrand(toeholdRegion) + "GGAUUUGCAAAAAAAAGAGGAGAGUAAAAUG" + reverse(toeholdRegion)[0:thsDuplex] + "AACCUGGCGGCAGCGCAAAAG")

#writing the text file for ths
outputFile2.write(f'# toehold switch structure in DU+ notation \n')
outputFile2.write(f'structure toeholdSwitch = ' + thsStructure + '\n')
outputFile2.write(f'\n# base sequences for each section of the toehold switch \n')
outputFile2.write(f'domain a = ' + compStrand(toeholdRegion) + '\n')
outputFile2.write(f'domain b = GGAUUUGCAAAAAAAAGAGGAGA N5 AUG \n')
outputFile2.write(f'domain c = ' + reverse(toeholdRegion)[0:thsDuplex] + '\n')
outputFile2.write(f'domain d = AACCUGGCGGCAGCGCAAAAG \n')
outputFile2.write(f'# toeholdSwitch sequence = ' + toeholdSwitchSequence + '\n')
outputFile2.write(f'\n# mapping RNA domains to the toehold switch \n')
outputFile2.write(f'toeholdSwitch.seq = a b c d')

#creating separate Strand object for toehold switch
thsObject = Strand(toeholdSwitchSequence, domains[len(domainsUsed)])

#attribute dict is a dictionary pairing each strand to a corresponding letter, these do not need to be the same letters used in the text files as these are for analysing the MFE of complexes
attributeDict = {}

#strand dict is the dictionary pairing each strand object to a corresponding concentration in M to then be passed onto a tube
strandDict = {}

#note to self:
# strand_n = Strand("ACGUAGCUAGCUAGCAUC", name="n")
# strand_n = Strand(Sequence, name="")


#passing attribute dict to the Strand class' object constructor and outputting that to a strand dict
for strand, letter in zip(strands, domainsUsed):
    attributeDict[strand] = letter 
    name = letter 
    strandDict[Strand(strand, name = name)] = 5e-6

t1 = Tube(strands = strandDict, complexes = SetSpec(max_size = len(strands)), name = "tube 1")

tube_result = tube_analysis(tubes = [t1], compute = ["pairs", "mfe"], model = thsModel)

sampled_structures = sample(strands = [strands[1], strands[0], strands[2]], num_sample = 5, model = thsModel)