import os
import sys
import getopt
import time
from Bio import Phylo

printToConsole = False

### CONSTANTS ###
ORTHOALIGN_PATH =  "OrthoAlign/OrthoAlign/"; ##I recommend using an absolute path
ORTHOALIGN_EXEC = "Aligning"
DUPLOSS_PATH = "2-SPP/"  ##I recommend using an absolute path
DUPLOSS_EXEC = "duploss"

testSetDir = ""

### METHODS ###
def getFirstLineFromFile(filename):
    """
    Return the first line in the file.
    """
    try:
        f = open(filename, 'r')
    except IOError as e:
        sys.exit("Error opening file "+filename+': '+e)
    line = f.readline()
    f.close()
    return line

def remBracketsDistributeNegSign(genome):
    """
    Removes all square brackets (representing operons) and distributes the negative signs in front of operons to all genes inside it.
    Changes <o> and <t> for [o] and [t].
    Returns the new string representing the updated genome.
    """
    negOpSplit = genome.split(",-[")
    
    for i in range(1, len(negOpSplit)):  #skipping the first index --> genomes don't start with -[
        op = negOpSplit[i][:negOpSplit[i].index(']')]  #slice of the operon
        op = op.replace(",", ",-")
        op = '-' + op
        singletons = negOpSplit[i][negOpSplit[i].index(']') + 1:]
        #print "op" + str(i) + " = " + op
        negOpSplit[i] = op + singletons
        
    genome = ",".join(negOpSplit)
    genome = genome.replace("[", "")
    genome = genome.replace("]", "")
    genome = genome.replace("<o>", "[o]")
    return genome.replace("<t>", "[t]")
    
def getCostAndAncestorFromOutFile(outfile):
    """
    Opens either a duploss or orthoAlign output file, finds the cost and returns it.
    """
    try:
        f = open(outfile, 'r')
    except IOError as e:
        sys.exit("Error opening file "+outfile+': '+e)
        
    line = f.readline()
    while line:
        splitted = line.split("ost = ")
        if len(splitted) > 1:
            cost = float(splitted[1])
            break
        line = f.readline()
    
    line = f.readline()
    while line:
        if ">Ancestor" in line:  #the next line should be the ancestor
            ancestor = f.readline()
            break
        line = f.readline()
    
    f.close()
    return cost, ancestor

def traverseNewickTree(node, parentNode):
    leftSibling = None
    rightSibling = None

    #Check there's any descendants
    if len(node.clades) > 0:
        leftSibling = traverseNewickTree(node.clades[0], node)
        if len(node.clades) > 1:
            rightSibling = traverseNewickTree(node.clades[1], node)

    if node.name != None and len(node.name) > 0:
        print node.name
        genomeFile = testSetDir + '/' + node.name + '/sequence.txt'

        if(not os.path.exists(genomeFile)):
            sys.exit("Error: "+genomeFile+" does not exist.\n")
        genome = getFirstLineFromFile(genomeFile)

        #Removing all whitespaces
        genome = genome.replace(" ", "")
        #Removing square brackets and distributing negative signs inside operons, updating origin and terminus
        genome = remBracketsDistributeNegSign(genome)

        return genome

    #Case 1: Both siblings exist therefore we need to construct their ancestor
    if leftSibling != None and rightSibling != None:
        command = "java -classpath " + ORTHOALIGN_PATH + " " + ORTHOALIGN_EXEC + " -dt " + leftSibling + " " + rightSibling + " > " + orthoAlignOutFile
        os.system(command)
        
        orthoCost, orthoAncestor = getCostAndAncestorFromOutFile(orthoAlignOutFile)
        return orthoAncestor
    #Case 2: Only the left sibling exists so return it
    elif leftSibling != None:
        # print "Case 2"
        return leftSibling
    #Case 3: Only the right sibling exists so return it
    elif rightSibling != None:
        # print "Case 3"
        return rightSibling
    #Case 4: None of the siblings exist so return NULL
    else:
        # print "Case 4"
        return None


def runOrthoAlignOnTree(newickFileName, testSetDir):
    print "Running OrthoAlign on newickTree"
    global testSetDir = testSetDir
    newickTree = Phylo.read(newickFileName, 'newick')

    orthoAlignStartTime = time.time()
    result = traverseNewickTree(newickTree.clade, None)
    orthoAlignRunTime = time.time() - orthoAlignStartTime

    fileDirectory = testSetDir.split("/")
    runtimePath = "/".join(fileDirectory[:-1])
    with open(runtimePath + "/OrthoRuntimes.txt", "a+") as runtimeFile:
            runtimeFile.write("%f " % (orthoAlignRunTime))

##
## Main
##
if __name__ == '__main__':

    usage = "usage: Run_2-SPP_OrthoAlign -hfdo genome1 genome2 destinationFolderName\n"+\
            " -f genome1 and genome2 are filenames instead of genome strings.\n"\
            " -h print this message.\n"\
            " -d run DupLoss.\n"\
            " -o run OrthoAlign.\n"\

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhdo")
    except getopt.GetoptError:
        sys.exit(usage)

    input_files = False
    runDupLoss = False
    runOrthoAlign = False
    
    for (opt,val) in opts:
        if(opt == "-h"):
            sys.exit(usage)
        if(opt == "-f"):
            input_files = True
        if(opt == "-d"):
            runDupLoss = True
        if(opt == "-o"):
            runOrthoAlign = True
    
    if len(args) != 3:
        sys.exit(usage)

    if not runDupLoss and not runOrthoAlign:
        sys.exit(usage)
        
        
    if(input_files):
        file1 = args[0]
        file2 = args[1]
        if(not os.path.exists(file1)):
            sys.exit("Error: "+file1+" does not exist.\n")
        if(not os.path.exists(file2)):
            sys.exit("Error: "+file2+" does not exist.\n")
        genome1 = getFirstLineFromFile(file1)
        genome2 = getFirstLineFromFile(file2)
    else:
        genome1 = args[0]
        genome2 = args[1]
        
    #Output filenames with output prefix
    duplossOutFile = args[2] + "/duploss.out"
    orthoAlignOutFile = args[2] + "/orthoAlign.out"
        
    #Removing all whitespaces
    genome1 = genome1.replace(" ", "")
    genome2 = genome2.replace(" ", "")
    
    #Removing square brackets and distributing negative signs inside operons, updating origin and terminus
    genome1 = remBracketsDistributeNegSign(genome1)
    genome2 = remBracketsDistributeNegSign(genome2)
        
    #print "genome 1 = " + genome1 + "\ngenome 2 = " + genome2
    
    #Running OrthoAlign
    if runOrthoAlign:
        command = "java -classpath " + ORTHOALIGN_PATH + " " + ORTHOALIGN_EXEC + " -dt " + genome1 + " " + genome2 + " > " + orthoAlignOutFile
        orthoAlignStartTime = time.time()
        os.system(command)
        orthoAlignRunTime = time.time() - orthoAlignStartTime

        orthoCost, orthoAncestor = getCostAndAncestorFromOutFile(orthoAlignOutFile)
    
    #Running Duploss
    if runDupLoss:
        command = "python " + DUPLOSS_PATH + DUPLOSS_EXEC + " -eiqdt " + genome1 + " " + genome2 + " > " + duplossOutFile
        duplossStartTime = time.time()
        os.system(command)
        duplossRunTime = time.time() - duplossStartTime
        
        duplossCost, duplossAncestor = getCostAndAncestorFromOutFile(duplossOutFile)
    
    if printToConsole:
        if runDupLoss:
            print "Duploss cost = " + str(duplossCost)
            print "Duploss ancestor =    " + duplossAncestor
        if runOrthoAlign:
            print "OrthoAlign cost = " + str(orthoCost)
            print "OrthoAlign ancestor = " + orthoAncestor
        
    fileDirectory = args[2].split("/")
    runtimePath = "/".join(fileDirectory[:-1])
    if runOrthoAlign:
        with open(runtimePath + "/OrthoRuntimes.txt", "a+") as runtimeFile:
            runtimeFile.write("%f " % (orthoAlignRunTime))
        
    if runDupLoss:
        with open(runtimePath + "/DuplossRuntimes.txt", "a+") as runtimeFile:
            runtimeFile.write("%f " % (duplossRunTime))