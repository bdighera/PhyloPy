from Bio import SeqIO
from time import sleep
from pprint import pprint
from ArguementHandler import ArguementParser
from SequenceCollector import SequeuceCollector, GenomicContext
from Sqlite import SQliteRecordInput, SQLiteChecker
import sys

'''
This is the main file for PhyloPy
Input: text file of fasta formatted protein records. Can also take text file of protein accession numbers
Output: Will run each protein accession indivudally through the entire pipeline and store in SQLite DB

Currently set to populate HSP70 table in records.db, if choosing to populate other/new table see documentation in sqlite.py
- table names must be changed manually
'''

def Main():
    # Initialize Argument parser class
    # DEPRECATED FUNCTION
    # However, function could be utilized for CLI interactions if needed in future
    Args = ArguementParser()

    # Use Argument.parse() method to retrieve input params from user
    # DEPRECATED
    # could be used again, retrieves input from the user via CLI
    proteinPath, cwd, rec, outfile, genomicPath, CDSPath, tree = Args.parse()

    # Input File name in CWD
    proteinPath = 'InputSeqs.txt'

    # sets the timer value between ping requests to the NCBI server
    timer = 6

    # Uses biopython class to parse the input file into various components including:
    # sequence, accession number, name
    # variable proteinSeqRecords consists of list of
    proteinSeqRecords = SeqIO.parse(proteinPath, 'fasta')

    # Iterates through each parsed record in the input file
    for proteinRecord in proteinSeqRecords:

        try:

            # Uses the record accession number to check if it is already in the SQL database, if not: continue.
            SQL = SQLiteChecker(proteinRecord.id)

            # If accession number is not located in the DB then it will run through the pipeline
            if SQL.checkRecords() != True:

                # Initialize Sequence Collector class with input of the working protein accession record
                Seq = SequeuceCollector(proteinRecord)

                # collects the protein ID for the current working protein
                # protein ID is needed in order to establish relationship between NCBI databases ie. Gene, Protein, etc.
                # runs timer after the method is executed in order to not ping NCBI server too fast - or they will SHUT YOU DOWN
                Seq.collectProteinIDs()
                sleep(timer)

                # Runs the method to collect the CDS and taxonomy
                # TO understand this method better (along with other methods of class Sequence Collector) visit SequenceCOllector.py in CWD
                CDS, Taxonomy = Seq.collectCDS()
                sleep(timer)

                # Collects the entire sequence of the gene and the gene ID
                Genomic, GeneID = Seq.collectGenomicSeq()
                sleep(timer)

                # Collects the intron phases and exon lengths - this method requires executible files from ./exec folder
                # Internal function -does not ping NCBI server
                IntronPhases, ExonLengths = Seq.collectIntrons(CDS, Genomic)

                # Collects the domains for the protein
                Domains = Seq.collectProteinDomains()
                sleep(timer)

                # Initialize Genomic Context class
                # Collects neighboring genes to parent protein, coding direction, and domains of those genes
                GC = GenomicContext(Genomic)

                # Collects record for  +/- 50k basepairs up/downstream of parent gene
                gcRecord = GC.fetchRecord()
                sleep(timer)

                # parses raw genomic context data collected from NCBI
                parsedGCrecord = GC.parseRecord(gcRecord)

                # Initialize the SQLite class with all of the data for the working protein
                SQL = SQliteRecordInput(Seq.proteinRecord,
                                        Seq.proteinID,
                                        CDS,
                                        Genomic,
                                        GeneID,
                                        parsedGCrecord,
                                        Domains,
                                        IntronPhases,
                                        ExonLengths,
                                        Taxonomy)

                # Uploads records to the DB
                # REMEMBER TO CHANGE THE NAME OF THE TABLE TO REFLECT WHERE YOU WANT THE RECORDS TOGO!
                SQL.uploadRecords()

            # If protein accession number is in the database then it will skip over and go to the next record
            else:
                print('Record ' + proteinRecord.id + ' is already in database. Proceeding to next record')

        # Errors in some of the sequences do occur on the NCBI server end. For this errors the DNAJC_Errors file is populated with those accession #s
        except IndexError:
            print('Index Error Occurred on: ' + str(proteinRecord.id))
            with open('DNAJC_Errors', 'a') as handle:
                handle.write(str(proteinRecord.id) + '\n')
                handle.close()

if __name__ == '__main__':
    Main()






