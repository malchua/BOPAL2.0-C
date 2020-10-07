from Event import Event
import numpy as np
import globals
import copy

printToConsole = False

#############################################
### Multiple Sequence Alignment Functions ###
#############################################

######################################################
# findOrthologsByMultiSequenceAlignment
# Parameters: Two orthologous operons, along with the neighbour events we are trying to align with
# Description: Finds the neighbour event that includes operon1. We then use the neighbour operon in that event and perform multiple sequence alignment.
######################################################
def findOrthologsByMultiSequenceAlignment(operon1, operon2, event, neighborEvents):
    if globals.printToConsole:
        print('Computing multiple sequence alignment matrix for: {%s, %s}...' % (', '.join(map(str, operon1.sequence)), ', '.join(map(str, operon2.sequence))))

    score = -1
    for neighbourEvent in neighborEvents:
        if operon1.fragmentIndex == neighbourEvent.fragmentDetails1.fragmentIndex and neighbourEvent.genome1Name != neighbourEvent.genome2Name:
            # print "Neighbour Events:"
            # print neighbourEvent.toString()
            # print "Genome1: " + neighbourEvent.genome1Name
            # print "Fragment1 Id: " + str(neighbourEvent.fragmentDetails1.fragmentIndex)
            # print "Sequence: " + '[%s]' % ', '.join(map(str, neighbourEvent.fragmentDetails1.sequence))
            # print "Genome2: " + neighbourEvent.genome2Name
            # print "neighbour fragment Id: " + str(neighbourEvent.fragmentDetails2.fragmentIndex)
            # print "Sequence: " + '[%s]' % ', '.join(map(str, neighbourEvent.fragmentDetails2.sequence))
            # print "Operon info:"
            # print "operon1 Id: " + str(operon1.fragmentIndex)
            # print "Sequence: " + '[%s]' % ', '.join(map(str, operon1.sequence))
            # print "operon2 Id: " + str(operon2.fragmentIndex)
            # print "Sequence: " + '[%s]' % ', '.join(map(str, operon2.sequence))
            # print " "
            score, event = performMultiSequenceAlignment(operon1.sequence, operon2.sequence, neighbourEvent.fragmentDetails2.sequence, event)
            

    # print "Traversed neighbour events"

    #initialize the matrix to store the global alignment scores
    #alignmentMatrix = [[[ 0.0 for x in range(0, len(strain2.genomeFragments))] for y in range(0, len(strain1.genomeFragments))] for z in range(0, len(strain3.genomeFragments))]
    #eventMatrix = [[None for x in range(0, len(strain2.genomeFragments))] for y in range(0, len(strain1.genomeFragments))]

    return score, event

######################################################
# performMultiSequenceAlignment
# Parameters: Three orthologous operons which we are aligning
# Description: Performs multiple sequence alignment on the three operons
######################################################
def performMultiSequenceAlignment(operon1, operon2, operon3, event):

    # Initialize the score and direction matrices for the alignment
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1, len(operon3)+1))
    dirMatrix = np.zeros((len(operon1)+1, len(operon2)+1, len(operon3)+1), dtype=int)

    # Intialize edges of the matrices
    for a in range(1, len(operon1)+1):
        scoreMatrix[a][0][0] = scoreMatrix[a-1][0][0] + calculateSumOfPairsScore(operon1[a-1], '-', '-')
        dirMatrix[a][0][0] = 7
    for a in range(1, len(operon2)+1):
        scoreMatrix[0][a][0] = scoreMatrix[0][a-1][0] + calculateSumOfPairsScore('-', operon2[a-1], '-')
        dirMatrix[0][a][0] = 6
    for a in range(1, len(operon3)+1):
        scoreMatrix[0][0][a] = scoreMatrix[0][0][a-1] + calculateSumOfPairsScore('-', '-', operon3[a-1])
        dirMatrix[0][0][a] = 5

    # Initialize three of the faces of the matrices
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], '-')
            scoreMatrix[a][b][0], dirMatrix[a][b][0] = findMax(scoreMatrix, a, b, 0, operon1, operon2, operon3)
    for a in range(1, len(operon1)+1):
        for c in range(1, len(operon3)+1):
            spScore = calculateSumOfPairsScore(operon1[a-1], '-', operon3[c-1])
            scoreMatrix[a][0][c], dirMatrix[a][0][c] = findMax(scoreMatrix, a, 0, c, operon1, operon2, operon3)
    for b in range(1, len(operon2)+1):
        for c in range(1, len(operon3)+1):
            spScore = calculateSumOfPairsScore('-', operon2[b-1], operon3[c-1])
            scoreMatrix[0][b][c], dirMatrix[0][b][c] = findMax(scoreMatrix, 0, b, c, operon1, operon2, operon3)

    # Perform the multiple sequence alignment
    print ""
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            for c in range(1, len(operon3)+1):
                # spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], operon3[c-1])
                # print str(spScore) + " ",
                scoreMatrix[a][b][c], dirMatrix[a][b][c] = findMax(scoreMatrix, a, b, c, operon1, operon2, operon3)
            # print ""
        # print "\n"

    # Compute the number of events that occured between the operons
    # event = alignmentTraceback(scoreMatrix, dirMatrix, operon1, operon2, operon3, event)

    np.set_printoptions(precision=3, linewidth=np.inf)
    if printToConsole:
        print scoreMatrix
        print dirMatrix
        print ""

    print operon1
    print operon2
    print operon3
    print ""

    # Compute the number of events that occured between the operons
    event = traceback(scoreMatrix, dirMatrix, operon1, operon2, operon3, event)

    return scoreMatrix[len(operon1)][len(operon2)][len(operon3)], event

######################################################
# calculateSumOfPairsScore
# Parameters: Three genes were are using to calculate a sum of pairs score
# Description: Calculates the sum of pairs score using only three genes
######################################################
def calculateSumOfPairsScore(gene1, gene2, gene3):
    sopScore = 0.0
    numPairings = 3.0

    sopScore += compareGenes(gene1, gene2)
    sopScore += compareGenes(gene1, gene3)
    sopScore += compareGenes(gene2, gene3)

    return sopScore/numPairings

######################################################
# compareGenes
# Parameters: Two genes we are comparing
# Description: Produces a score between two genes based on their similarity
######################################################
def compareGenes(gene1, gene2):
    score = 0.0

    if gene1 != '-' and gene2 != '-':
        if gene1.split('_')[0].strip() == gene2.split('_')[0].strip():
            if gene1.strip() == gene2.strip():
                score = globals.match
            else:
                score = globals.codonCost
        else:
            score = globals.substitutionCost
    elif gene1 == '-' and gene2 == '-':
        # Score for a gap alignment
        score = 0
    else:
        score = globals.deletionCost

    return score

######################################################
# findMax
# Parameters: The score matrix along with indices for the current position in the matrix. The operons being compared are also included.
# Description: Finds the max score of all scores calculated from each possible direction of the score matrix. Also returns the direction of the max score.
######################################################
def findMax(scoreMatrix, a, b, c, operon1, operon2, operon3):
    maxScore = -999.0

    # Match in gene1,gene2,gene3
    if a != 0 and b != 0 and c != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], operon3[c-1])
        if maxScore < scoreMatrix[a-1][b-1][c-1] + spScore:
            maxScore = scoreMatrix[a-1][b-1][c-1] + spScore
            direction = 1
    # Match in gene1,gene2. Gap in gene 3
    if a != 0 and b != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], '-')
        if maxScore < scoreMatrix[a-1][b-1][c] + spScore:
            maxScore = scoreMatrix[a-1][b-1][c] + spScore
            direction = 2
    # Match in gene1,gene3. Gap in gene 2
    if a != 0 and c != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], '-', operon3[c-1])
        if maxScore < scoreMatrix[a-1][b][c-1] + spScore:
            maxScore = scoreMatrix[a-1][b][c-1] + spScore
            direction = 3
    # Match in gene2,gene3. Gap in gene 1
    if b != 0 and c != 0:
        spScore = calculateSumOfPairsScore('-', operon2[b-1], operon3[c-1])
        if maxScore < scoreMatrix[a][b-1][c-1] + spScore:
            maxScore = scoreMatrix[a][b-1][c-1] + spScore
            direction = 4
    # Gap in gene1, gene2
    if c != 0:
        spScore = calculateSumOfPairsScore('-', '-', operon3[c-1])
        if maxScore < scoreMatrix[a][b][c-1] + spScore:
            maxScore = scoreMatrix[a][b][c-1] + spScore
            direction = 5
    # Gap in gene1, gene3
    if b != 0:
        spScore = calculateSumOfPairsScore('-', operon2[b-1], '-')
        if maxScore < scoreMatrix[a][b-1][c] + spScore:
            maxScore = scoreMatrix[a][b-1][c] + spScore
            direction = 6
    # Gap in gene2, gene3
    if a != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], '-', '-')
        if maxScore < scoreMatrix[a-1][b][c] + spScore:
            maxScore = scoreMatrix[a-1][b][c] + spScore
            direction = 7

    if maxScore == -999.0:
        print "ERROR: No max score found. Direction was not chosen. Exiting..."
        sys.exit(0)

    return maxScore, direction

######################################################
# traceback
# Parameters: The score matrix and direction matrix, along with the operons being compared.
# Description: Produces an optimal alignment of the three operons.
######################################################
def traceback(matrix, dirMatrix, operon1, operon2, operon3, event):

    i = len(operon1)
    j = len(operon2)
    k = len(operon3)

    match = 0
    codonMismatch = 0
    mismatch = 0
    substitution = 0

    # Tracks index and genes for codon mismatches in both strains
    codonMismatchIndexesStrain1 = []
    codonMismatchIndexesStrain2 = []
    codonMismatchGenesStrain1 = []
    codonMismatchGenesStrain2 = []

    # Tracks substitution indexes and genes for both strains
    substitutionIndexesStrain1 = []
    substitutionIndexesStrain2 = []
    substitutionGenesStrain1 = []
    substitutionGenesStrain2 = []

    # Tracks the genes in a gap and the index of those genes in both strains note: details stored in an array of arrays
    operon1Gaps = []
    operon1Gap = []
    operon1GapIndexes = [] #This is used to determine the position of the genes with respect to the genome
    operon1GapIndex = []
    operon1ConsecutiveGap = False #Tracks consecutive gaps

    operon2Gaps = []
    operon2Gap = []
    operon2GapIndexes = []
    operon2GapIndex = []
    operon2ConsecutiveGap = False #Tracks consecutive gaps

    # track the alignment
    alignmentSequence1 = []
    alignmentSequence2 = []
    alignmentSequence3 = []

    # Tracks where the extra genes are from
    gap1Indexes = [] #This is used to determine where to insert the genes into the alignment
    gap2Indexes = []

    # Tracks which events are known because of the alignment with the neighbour
    operon1GapEvents = []
    operon1GapEvent = []
    operon2GapEvents = []
    operon2GapEvent = []
    
    selfDuplication = '' #Used only for self global alignment
    # selfPosition = event.fragmentDetails1.startPositionInGenome

    while i > 0 or j > 0 and k > 0:
        # Case 1: Perfect match
        if i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2) and operon1[i-1] == operon2[j-1]:
            # Self global alignment
            # if event.fragmentDetails1.isNegativeOrientation == False:
            #     selfDuplication = operon2[j-1] + ' ' + str((i-1) + selfPosition) + ', ' + selfDuplication
            # else:
            #     selfDuplication = selfDuplication + operon2[j-1] + ' ' + str(len(operon1) - (i-1) + selfPosition) + ', '

            match += 1
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])

            if dirMatrix[i][j][k] == 1:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            
        # Case 2: Codon mismatch
        elif i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            # Self global alignment
            # if event.fragmentDetails1.isNegativeOrientation == False:
            #     selfDuplication = '!' + operon2[j-1] + ' ' + str(-1) + ', ' + selfDuplication
            # else:
            #     selfDuplication = selfDuplication + '!' + operon2[j-1] + ' ' + str(-1) + ', '

            # Increment the Id counter to ensure Id id unique
            globals.codonMismatchId += 1
            
            codonMismatch += 1
            alignmentSequence1.insert(0, operon1[i-1] + '-#' + str(globals.codonMismatchId) + '#')
            alignmentSequence2.insert(0, operon2[j-1] + '-#' + str(globals.codonMismatchId) + '#')

            codonMismatchIndexesStrain1.append(i-1)
            codonMismatchGenesStrain1.append(operon1[i-1] + '-#' + str(globals.codonMismatchId) + '#')

            codonMismatchIndexesStrain2.append(j-1)
            codonMismatchGenesStrain2.append(operon2[j-1] + '-#' + str(globals.codonMismatchId) + '#')
            # alignmentSequence1.insert(0, operon1[i-1])
            # alignmentSequence2.insert(0, operon2[j-1])

            if dirMatrix[i][j][k] == 1:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            
        # Case 3: Substitution
        elif i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2):
            # Self global alignment
            # if event.fragmentDetails1.isNegativeOrientation == False:
            #     selfDuplication = '!' + operon2[j-1] + ' ' + str(-1) + ', ' + selfDuplication
            # else:
            #     selfDuplication = selfDuplication + '!' + operon2[j-1] + ' ' + str(-1) + ', '
                
            # Increment the Id counter to ensure the ID is unique
            globals.substitutionId += 1
            
            substitution += 1
            alignmentSequence1.insert(0, operon1[i-1] + '-@' + str(globals.substitutionId) + '@')
            alignmentSequence2.insert(0, operon2[j-1] + '-@' + str(globals.substitutionId) + '@')

            substitutionIndexesStrain1.append(i-1)
            substitutionGenesStrain1.append(operon1[i-1] + '-@' + str(globals.substitutionId) + '@')

            substitutionIndexesStrain2.append(j-1)
            substitutionGenesStrain2.append(operon2[j-1] + '-@' + str(globals.substitutionId) + '@')
            # alignmentSequence1.insert(0, operon1[i-1])
            # alignmentSequence2.insert(0, operon2[j-1])

            if dirMatrix[i][j][k] == 1:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            
        # Case 4: Mismatch- Gap in operon 2
        elif i > 0 and (dirMatrix[i][j][k] == 3 or dirMatrix[i][j][k] == 7):
            index = i-1
            mismatch += 1

            if dirMatrix[i][j][k] == 3:
                alignmentSequence3.insert(0, operon3[k-1])
                # Adding to operon1 to pre-empt the switch done later in code
                eventTitle = "loss" 
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")
                eventTitle = "dup"

            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon2ConsecutiveGap:
                operon2Gap.insert(0, operon1[index])
                operon2GapIndex.insert(0, index)
                operon2GapEvent.insert(0, eventTitle)

                operon2ConsecutiveGap = True
            else:
                if len(operon2Gap) > 0:
                    operon2Gaps.insert(0, operon2Gap)
                    operon2GapIndexes.insert(0, operon2GapIndex)
                    operon2GapEvents.insert(0, operon2GapEvent)
                operon2Gap = []
                operon2GapIndex = []
                operon2GapEvent = []
                operon2Gap.insert(0, operon1[index])
                operon2GapIndex.insert(0, index)
                operon2GapEvent.insert(0, eventTitle)
                gap2Indexes.insert(0, len(alignmentSequence2))
                operon2ConsecutiveGap = True
            # Only for debugging
            # alignmentSequence1.insert(0, operon1[i-1])
            # alignmentSequence2.insert(0, "-")

            i -= 1
            operon1ConsecutiveGap = False

        # Case 5: Mismatch - Gap in operon 1
        elif j > 0 and (dirMatrix[i][j][k] == 4 or dirMatrix[i][j][k] == 6):
            #Self global alignment
            # if event.fragmentDetails1.isNegativeOrientation == False:
            #     selfDuplication = '!' + operon2[j-1] + ' ' + str(-1) + ', ' + selfDuplication
            # else:
            #     selfDuplication = selfDuplication + operon2[j-1] + ' ' + str(-1) + ', '
                
            index = j - 1
            mismatch += 1

            if dirMatrix[i][j][k] == 4:
                alignmentSequence3.insert(0, operon3[k-1])
                # Adding to operon2 to pre-empt the switch done later in code
                eventTitle = "loss"
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")
                eventTitle = "dup"

            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon1ConsecutiveGap:
                operon1Gap.insert(0, operon2[index])
                operon1GapIndex.insert(0, index)
                operon1GapEvent.insert(0, eventTitle)

                operon1ConsecutiveGap = True
            else:
                if len(operon1Gap) > 0:
                    operon1Gaps.insert(0, operon1Gap)
                    operon1GapIndexes.insert(0, operon1GapIndex)
                    operon1GapEvents.insert(0, operon1GapEvent)
                operon1Gap = []
                operon1GapIndex = []
                operon1GapEvent = []
                operon1Gap.insert(0, operon2[index])
                operon1GapIndex.insert(0, index)
                operon1GapEvent.insert(0, eventTitle)
                gap1Indexes.insert(0, len(alignmentSequence1))
                operon1ConsecutiveGap = True
            # Only for debugging
            # alignmentSequence2.insert(0, operon2[j-1])
            # alignmentSequence1.insert(0, "-")

            j -= 1
            operon2ConsecutiveGap = False

        # Case 6: Mismatch - Gap in operon 1 and operon 2 (need to move operon3 index)
        else:
            if dirMatrix[i][j][k] == 5:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
                # alignmentSequence1.insert(0, "-")
                # alignmentSequence2.insert(0, "-")
            else:
                print "Error: Case 6 encountered but direction was not 5. Direction was " + str(dirMatrix[i][j][k])

    # print [gene.center(7) for gene in alignmentSequence1]
    # print [gene.center(7) for gene in alignmentSequence2]
    # print [gene.center(7) for gene in alignmentSequence3]

    event.selfDuplication = selfDuplication[0:(len(selfDuplication) - 2)] + ';' #Remove the last comma and space and add a semicolon 
    #Empty any remaining gaps
    if len(operon1Gap) > 0:
        operon1Gaps.insert(0, operon1Gap)
        operon1GapIndexes.insert(0, operon1GapIndex)
        operon1GapEvents.insert(0, operon1GapEvent)
        operon1Gap = []
        operon1GapIndex = []
        operon1GapEvent = []

    if len(operon2Gap) > 0:
        operon2Gaps.insert(0, operon2Gap)
        operon2GapIndexes.insert(0, operon2GapIndex)
        operon2GapEvents.insert(0, operon2GapEvent)
        operon2Gap = []
        operon2GapIndex = []
        operon2GapEvent = []

    #The indexes values need to be flipped b/c right now they're oriented from right to left
    if len(gap1Indexes) > 0:
        for x in range(0, len(gap1Indexes)):
            gap1Indexes[x] = len(alignmentSequence1) - gap1Indexes[x]
    if len(gap2Indexes) > 0:
        for x in range(0, len(gap2Indexes)):
            gap2Indexes[x] = len(alignmentSequence2) - gap2Indexes[x]

    #Need to swap the gap lists since the gaps refer to extra genes
    temp = operon1Gaps
    operon1Gaps = operon2Gaps
    operon2Gaps = temp

    temp = operon1GapIndexes
    operon1GapIndexes = operon2GapIndexes
    operon2GapIndexes = temp

    temp = operon1GapEvents
    operon1GapEvents = operon2GapEvents
    operon2GapEvents = temp

    temp = gap1Indexes
    gap1Indexes = gap2Indexes
    gap2Indexes = temp

    #Match, mismatch details
    event.setNumMatches(match)
    event.setNumMismatches(mismatch)

    #Codon Mismatch details
    event.setNumCodonMismatches(codonMismatch)
    event.setCodonMismatchGenesStrain1(codonMismatchGenesStrain1)
    event.setCodonMismatchGenesStrain2(codonMismatchGenesStrain2)
    event.setCodonMismatchIndexesStrain1(codonMismatchIndexesStrain1)
    event.setCodonMismatchIndexesStrain2(codonMismatchIndexesStrain2)

    #Substitution details
    event.setNumSubstitutions(substitution)
    event.setSubstitutionIndexesStrain1(substitutionIndexesStrain1)
    event.setSubstitutionIndexesStrain2(substitutionIndexesStrain2)
    event.setSubstitutionGenesStrain1(substitutionGenesStrain1)
    event.setSubstitutionGenesStrain2(substitutionGenesStrain2)

    #Alignment details
    event.setOperon1Alignment(alignmentSequence1)
    event.setOperon2Alignment(alignmentSequence2)

    #Gap details
    event.setOperon1Gaps(operon1Gaps)
    event.setOperon2Gaps(operon2Gaps)
    event.setOperon1GapPositions(operon1GapIndexes)
    event.setOperon2GapPositions(operon2GapIndexes)
    event.setOperon1GapIndexes(gap1Indexes)
    event.setOperon2GapIndexes(gap2Indexes)
    event.setOperon1GapEvents(operon1GapEvents)
    event.setOperon2GapEvents(operon2GapEvents)

    #Used for debugging
    #print('These are the operons being compared: %s, %s' %(operon1, operon2))
    #print('This is the resulting alignment: %s, %s' %(alignmentSequence1, alignmentSequence2))
    #print('These are the extra genes for operon 1: %s' %(operon1Gaps))
    #print('These are the indexes for extra genes in operon 1: %s' %(gap1Indexes))
    #print('These are the extra genes for operon 2: %s' %(operon2Gaps))
    #print('These are the indexes for extra genes in operon 2: %s' %(gap2Indexes))

    return event

def splitByEvents(gap, positions, events):
    newGaps = []
    newGap = []
    newPositions = []
    newPosition = []
    newEvents = []

    prevEvent = None

    for i in range(0, len(gap)):
        currEvent = events[i]

        if currEvent == prevEvent:
            newGap.append(gap[i])
            newPosition.append(positions[i])
        else:
            if len(newGap) > 0:
                newGaps.append(newGap)
                newPositions.append(newPosition)
                newEvents.append(prevEvent)
            newGap = []
            newPosition = []
            newGap.append(gap[i])
            newPosition.append(positions[i])

        prevEvent = currEvent

    if len(newGap) > 0:
        newGaps.append(newGap)
        newPositions.append(newPosition)
        newEvents.append(prevEvent)
        newGap = []
        newPosition = []

    # print "Splitting gap"
    # print gap
    # print positions
    # print events
    # print newGaps
    # print newPositions
    # print newEvents

    return newGaps, newPositions, newEvents


def main():
    global printToConsole
    printToConsole = False

    globals.initialize() #Initialize the globals file

    op1 = ["16S", "Mer_GGC", "Lys_GCA", "Ala_GCA", "23S", "5S", "16S", "Met_CAA"]
    op2 = ["16S", "Ile_AUC", "Arg_UCU", "Leu_AGU", "Ala_GCA", "23S", "5S"]
    op3 = ["16S", "Gln_CAA", "Arg_UCU", "Leu_AGU", "Arg_GGC", "Lys_GCA", "Ala_GCA", "23S", "5S", "Met_CAA"]

    event = Event(0)
    event.setScore(1.0)
    event.setDistance(abs(5 - 10))
    event.setFragmentDetails1(op1)
    event.setFragmentDetails2(op2)
    event.setGenome1Name("strain1")
    event.setGenome2Name("strain2")
    event.setTechnique('Multiple Sequence Alignment')

    score, event = performMultiSequenceAlignment(op1, op2, op3, event)

    print str(score)
    print event.operon1Gaps
    print event.operon1GapPositions
    print event.operon1GapIndexes
    print event.operon1GapEvents
    print " "
    print event.operon2Gaps
    print event.operon2GapPositions
    print event.operon2GapIndexes
    print event.operon2GapEvents

    splitByEvents(event.operon1Gaps[1], event.operon1GapPositions[1], event.operon1GapEvents[1])
    splitByEvents(event.operon2Gaps[0], event.operon2GapPositions[0], event.operon2GapEvents[0])

    op1 = ["16S", "Ile_AUC", "23S", "5S"]
    op2 = ["16S", "Ala_GCA", "23S", "5S"]
    op3 = ["16S", "Ala_GCA", "23S", "5S", "16S"]

    event = Event(0)
    event.setScore(1.0)
    event.setDistance(abs(5 - 10))
    event.setFragmentDetails1(op1)
    event.setFragmentDetails2(op2)
    event.setGenome1Name("strain1")
    event.setGenome2Name("strain2")
    event.setTechnique('Multiple Sequence Alignment')

    performMultiSequenceAlignment(op1, op2, op3, event)

    op1 = ["16S", "Ile_AUC", "23S", "16S"]
    op2 = ["16S", "Ala_GCA", "23S", "5S"]
    op3 = ["16S", "23S", "Ala_GCA", "5S", "16S"]

    event = Event(0)
    event.setScore(1.0)
    event.setDistance(abs(5 - 10))
    event.setFragmentDetails1(op1)
    event.setFragmentDetails2(op2)
    event.setGenome1Name("strain1")
    event.setGenome2Name("strain2")
    event.setTechnique('Multiple Sequence Alignment')

    performMultiSequenceAlignment(op1, op2, op3, event)

    op1 = ['16S', 'Ala_GCC', '23S', '5S']
    op2 = ['16S', 'Ala_GCA', '23S', '5S']
    op3 = ['16S', 'Ala_GCA', 'Ala_GCC', '23S', '5S']

    event = Event(0)
    event.setScore(1.0)
    event.setDistance(abs(5 - 10))
    event.setFragmentDetails1(op1)
    event.setFragmentDetails2(op2)
    event.setGenome1Name("strain1")
    event.setGenome2Name("strain2")
    event.setTechnique('Multiple Sequence Alignment')

    score, event = performMultiSequenceAlignment(op1, op2, op3, event)
    print str(score)
    print event.operon1Gaps
    print event.operon1GapPositions
    print event.operon1GapIndexes
    print event.operon1GapEvents
    print " "
    print event.operon2Gaps
    print event.operon2GapPositions
    print event.operon2GapIndexes
    print event.operon2GapEvents


if __name__ == "__main__":
    main()
