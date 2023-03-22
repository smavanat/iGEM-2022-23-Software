import string 
import math as m
#from nupack import *
#import matplotlib.pyplot as plt

#thsModel = Model(material='rna', celsius=37)

#File for storing output
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
    part3 = miRNAInput[n:2*n+1]
    storageArray.append(compStrand(part1))
    storageArray.append(compStrand(part2))
    storageArray.append(part3)
    return storageArray

#computes and gate
def andGateGen(sequenceDict):
    storageDict = sequenceDict.values()
    storageArray = []
    miRNAPartArray = []
    global unpairedMiddleArray
    unpairedMiddleArray = []
    
    #Takes each of the thirds and adds the first and last to the miRNA part array (used for designing the and gates) and the middle is used to determine the linker sequence
    for element in storageDict:
        storageArray = splitString(element)
        miRNAPartArray.append(storageArray[0])
        miRNAPartArray.append(storageArray[1])
        unpairedMiddleArray.append(storageArray[2])

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

            #checking to see if the andRegion borders an identical base and changing it to correct this by incrimenting the index
            #this only needs to be done twice as in the "worst case scenario", ALL of the below code will run, making the anndRegion different to BOTH of its neighboring bases
            #if (miRNA1[-1] == andRegion[1]):
            #    index += 1
            #    andRegion = andRegions[index]
            #if (reverse(miRNA2)[0] == andRegion[1]):
            #    index += 1
            #    andRegion = andRegions[index]
            #if (reverse(miRNA2)[0] and miRNA1[-1] == andRegion[1]):
            #    index += 1
            #    andRegion = andRegions[index]            

            #defining a repeating unit for the complex using DU+ notation.
            duplex1 = "D" + str(len(miRNA1)) + "+"
            complexAnd = "U" + str(len(andRegion))
            duplex2 = "D" + str(len(miRNA2)) + "+"
            unpairedMiddle = "U" + str(len(unpairedMiddleArray[middleIndex]))

            repeatUnit = [duplex1, complexAnd, duplex2, unpairedMiddle]  
            
            andGateDict = {}
            andGateArray = []
            for gate in andRegions:
                if (reverse(miRNA2)[0] and miRNA1[-1] == gate[1]):
                    #energy = structure_energy(strands = [storageDict[miRNAIndex], f'{compStrand(miRNA1)}{gate}{compStrand(miRNA2)}', storageDict[miRNAIndex]], 
                    #                         structure = repeatUnit1, model = thsModel)
                    #energy = structure_energy(strands = [gate, compStrand(gate)], 
                    #                         structure = "D6+", model = thsModel)
                    andGateArray.append(gate)
                    #andGateDict[energy] = gate
            print(andGateArray)
            #miRNAIndex +=1
            andGateArray = sorted(andGateArray)
            andRegion = andGateArray[0]

            andGate = miRNA1 + andRegion + reverse(miRNA2)

            strands.append(andGate)
            middleIndex+=1

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

#Calculating unpaired start of linker sequence
unpairedStartLength = len(list(miRNADict.values())[0]) - round((len(list(miRNADict.values())[0])//3))

if (unpairedStartLength % 3 == 2):
    unpairedStartLength = unpairedStartLength + 1 
elif (unpairedStartLength % 3 == 1):
    unpairedStartLength = unpairedStartLength - 1
else:
    unpairedStartLength = unpairedStartLength

#Transposing the start codon into DU+ notation
unpairedStart = "U" + str(unpairedStartLength)
unpairedEnd = "U" + str(len(list(miRNADict.values())[-1]) - len(miRNA2))
complexStructure[-1] = unpairedEnd
complexStructure.insert(0, unpairedStart)

for value in miRNADict.values():
    strands.append(value)

domainCodes = []

for strand, domain in zip(strands, domains):
    domainCode = "domain " + domain + " = " + strand
    domainCodes.append(domainCode)

domainsUsed = domains[0:len(domainCodes)]

mid = int((len(domainsUsed) + 1)/2)

#Putting all the domains in one array in order
miRNAs = list(domainsUsed[mid-1:])
andGates = list(domainsUsed[0:mid-1])

seq.extend(miRNAs)

integer = 0

#Placing elements from andGates between each term in miRNAs
while (integer < len(andGates)):
    for domain in seq:
        if((seq.index(domain)%2 ==0) and (seq.index(domain) + 1 < len(seq))):
            seq.insert(seq.index(domain) + 1, andGates[integer])
            integer += 1
        else:
            continue

#Outputting all the information to a text file
outputFile.write(f'# miRNA complex structure in DU+ notation \nstructure complex = {"".join(complexStructure)}\n')
outputFile.write(f'\n# base sequences for each RNA strand')
for code in domainCodes:
    outputFile.write(f"\n" + code)
outputFile.write(f'\n \n# mapping RNA strands to the complex')
outputFile.write(f'\ncomplex.seq = ')
for sequence in seq:
    outputFile.write(sequence + " ")

#toehold switch stuff
toeholdStart = list(miRNADict.values())[0][0:unpairedStartLength]
toeholdEnd = list(miRNADict.values())[-1][len(miRNA2):]
toeholdAndGates = []

middleStrands = strands[0:-1]
del middleStrands[len(middleStrands)//2]

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

while (integer < len(unpairedMiddleArray)):
    for domain in toeholdAndGates:
        if((toeholdAndGates.index(domain)%2 ==0) and (toeholdAndGates.index(domain) + 1 < len(toeholdAndGates))):
            toeholdAndGates.insert(toeholdAndGates.index(domain) + 1, unpairedMiddleArray[integer])
            integer += 1
        else:
            continue

toeholdRegion = []

toeholdRegion.insert(0, toeholdStart)
toeholdRegion.extend(toeholdAndGates)
toeholdRegion.append(toeholdEnd)
toeholdRegion = "".join(str(i) for i in toeholdRegion)

toehold = round((len(list(miRNADict.values())[0])//2))
thsDuplex = len(toeholdRegion) - toehold

thsStructure = ("U"+ str(toehold) + "D"+ str(thsDuplex) + "(U3D5(U15)U3)U21")

toeholdSwitchSequence = (compStrand(toeholdRegion) + "GGAUUUGCAAAAAAAAGAGGAGAGUAAAAUG" + reverse(toeholdRegion)[0:thsDuplex] + "AACCUGGCGGCAGCGCAAAAG")

#Outputting stuff to the text files.
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

strands.append(toeholdSwitchSequence)
domainsUsed.append(domains[len(domainsUsed)])

attributeDict = {}
strandDict = {}

#for strand, letter in zip(strands, domainsUsed):
#    attributeDict[strand] = letter 
#    name = letter 
#    strandDict[Strand(strand, name = name)] = 5e-6
    
#t1 = Tube(strands = strandDict, complexes = SetSpec(max_size = len(strandDict)), name="t1")
#tube_result = tube_analysis(tubes=[t1], compute=['mfe'], model=thsModel)