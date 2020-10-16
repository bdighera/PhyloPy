import re, itertools, os, subprocess, sys, math, dendropy, pprint, ast
from ImageProcessor import PhyloTreeConstruction
from SQLiteRetriever import RecordRetrival, SQLiteChecker

def Main():
    # Params
    # *********************************************************
    inputFile = 'InputSeqs.txt'
    runGenomicContextFig = False
    runIntronFig = True
    runDomainFig = True
    # **********************************************************
    V = RecordRetrival()

    V.retrieveRecords(inputFile)

    records = V.pullDBrecords(dbfile='Records.db')

    proteinAccession = [record[0][1] for record in records]
    proteinSeq = [record[0][2] for record in records]
    proteinDescription = [record[0][3] for record in records]
    genomicContext = [ast.literal_eval(record[0][12]) for record in records]
    parentDomains = [ast.literal_eval(record[0][13]) for record in records]
    introns = [record[0][14] for record in records]
    exonLength = [record[0][15] for record in records]

    Phylo = PhyloTreeConstruction(proteinAccession,
                                  proteinSeq,
                                  proteinDescription,
                                  genomicContext,
                                  parentDomains,
                                  introns,
                                  exonLength)


    if runIntronFig == True:
        Phylo.renderingTreewithIntrons()

    if runGenomicContextFig == True:
        Phylo.renderingTreewithGenomicContext()

    if runDomainFig == True:
        Phylo.renderingTreewithDomains()

if __name__ == '__main__':
    Main()