import ast, sys
from ImageProcessor import PhyloTreeConstruction
from SQLiteRetriever import RecordRetrival
from SequenceCollector import collectTaxonomy

def Main(accessionList = None, fileHandle = '', prune='',
         JupyterNotebook= False,
         runGenomicContextFig = True,
         runIntronFig= True,
         runDomainFig= True,
         printLeafAccession= False,
         printLeafDesc= False):

    inputFile = 'InputSeqs.txt'
    JupyterNotebookInlineFig = JupyterNotebook


    V = RecordRetrival()

    if JupyterNotebookInlineFig == True:
        if accessionList != None:

            V.retrieveRecordsbyList(accessionList)

        elif accessionList == None and fileHandle == None:
            raise ValueError('Need to pass list of accession numbers to Main.py')

        elif fileHandle != '':
            V.retrieveFileRecords(fileHandle)

        pass

    elif JupyterNotebookInlineFig == False:
        V.retrieveFileRecords(inputFile)



    records = V.pullDBrecords(dbfile='Records.db')

    proteinAccession = [record[0][1] for record in records]
    proteinSeq = [record[0][2] for record in records]
    proteinDescription = [record[0][3] for record in records]
    genomicContext = [ast.literal_eval(record[0][12]) for record in records]
    parentDomains = [ast.literal_eval(record[0][13]) for record in records]
    introns = [record[0][14] for record in records]
    exonLength = [record[0][15] for record in records]
    taxonomy = [record[0][16] for record in records]
    commonNames = [record[0][17] for record in records]


    Phylo = PhyloTreeConstruction(proteinAccession,
                                  proteinSeq,
                                  proteinDescription,
                                  genomicContext,
                                  parentDomains,
                                  introns,
                                  exonLength,
                                  JupyterNotebookInlineFig,
                                  prune,
                                  printLeafAccession,
                                  printLeafDesc,
                                  commonNames)


    if runIntronFig == True:
        Phylo.renderingTreewithIntrons()

    if runGenomicContextFig == True:
        Phylo.renderingTreewithGenomicContext()

    if runDomainFig == True:
        Phylo.renderingTreewithDomains()


    return [proteinAccession,
            proteinDescription,
            genomicContext,
            parentDomains,
            introns,
            exonLength,
            commonNames]

if __name__ == '__main__':
    Main()