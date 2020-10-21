import sqlite3, sys
from Bio import SeqIO

class RecordRetrival():

    def __init__(self):
        pass

    #This method is used to retrieve desired records for visualization from input file
    #input file consists of file in cd with fasta formatting
    def retrieveFileRecords(self, filename):


        with open(filename, 'r') as handle:
            contents = handle.readlines()
            contents = [i.strip('\n') for i in contents]

            handle.close()

        self.inputAccessionList = contents

        return self.inputAccessionList

    def retrieveFastaRecords(self, filename):

        print('Ensure that input file is in fasta formatting')

        file = SeqIO.parse(filename, 'fasta')

        self.inputAccessionList = [item.name for item in file]

        return self.inputAccessionList

    #This method is used to take the input accession list and use it as primary key for sqlite db
    #pulls the db entry from the accession number of input file
    def pullDBrecords(self, dbfile):

        inputList = self.inputAccessionList

        #SQLiteChecker class takes the input accession list and the name of the dbfile
        #dbfile currently set to records.db - only current records file
        SQL = SQLiteChecker(proteinAccession= inputList, dbfile= dbfile)

        records = SQL.checkRecords()

        return records

    def retrieveRecordsbyList(self, accessionList):
        self.inputAccessionList = accessionList
        return self.inputAccessionList


class SQLiteChecker():
    def __init__(self, proteinAccession, dbfile):
        self.dbfile= dbfile
        self.connect = sqlite3.connect(dbfile)
        self.proteinAccessionList = proteinAccession

    # Checks to make sure that the record isn't already in the table before it runs it. Checks against the accession number of the protein fasta
    def checkRecords(self):
        #Will check all of the records in all of the tables

        dataList = []

        C = self.connect.cursor()

        for accession in self.proteinAccessionList:

            #SQL statement getting all of the table names ie. DNAJC, HSP70
            C.execute("SELECT name FROM sqlite_master WHERE type='table'")

            for tablename in C.fetchall():
                tablename = tablename[0]

                #Checking all of the tables for the accession number
                C.execute('SELECT * FROM {t} WHERE ProteinAccession= (?)'''.format(t=tablename), (accession,))

                data = C.fetchall()

                if data != []:
                    dataList.append(data)

        return dataList