import sqlite3, sys
from uuid import uuid4

'''
These classes control of the SQLite DB. Remember that there must be manual editing for the SQL table names
This is needed to make sure that the accessions are checked against/added to the correct table in the DB.
If you require a new table in the db then run this script and update the method makeNewTable of class makeNewTable to reflect the necessary table name
Remember to update that in the other two classes as well
'''

#Takes all of the records from the various collection procedures and formats them, and inputs them into the database
class SQliteRecordInput():

    def __init__(self, Protein, ProteinID, CDS, Genomic, GeneID, GC, Domains, IntronPhase, ExonLength, Taxonomy):
        self.conn = sqlite3.connect('Records.db')
        self.ProteinRecord = Protein
        self.ProteinSeq = str(Protein.seq)
        self.ProteinAccession = Protein.id
        self.ProteinDescription =Protein.description
        self.ProteinID = int(ProteinID)
        self.CDS = CDS
        self.CDSSeq = str(CDS.seq)
        self.CDSAccession = CDS.id
        self.CDSDescription = CDS.description
        self.Genomic = Genomic
        self.GenomicSeq = str(Genomic.seq)
        self.GenomicAccession = str(Genomic.id)
        self.GenomicDescription = Genomic.description
        self.GeneID = int(GeneID)
        self.GC = str(GC)
        self.Domains = str(Domains)
        self.Introns = str(IntronPhase)
        self.ExonLength = str(ExonLength)
        self.uuid = str(uuid4())
        self.taxonomy = str(Taxonomy)

    def uploadRecords(self):

        C = self.conn.cursor()

        # C.execute('''CREATE TABLE HSP70 (UUID TEXT PRIMARY KEY ,ProteinAccession TEXT UNIQUE, ProteinSequence TEXT, ProteinDescription TEXT,
        #                           ProteinID INTEGER, CDSAccession TEXT, CDSSeq TEXT, CDSDescription TEXT, GenomicAccession TEXT,
        #                         GenomicSeq TEXT, GenomicDescription TEXT, GeneID INTEGER, GenomicContext TEXT, Introns TEXT, ExonLength TEXT)''')

        try:
            C.execute('''INSERT INTO HSP70 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',(self.uuid, self.ProteinAccession, self.ProteinSeq, self.ProteinDescription, self.ProteinID,
                                                                                  self.CDSAccession, self.CDSSeq, self.CDSDescription, self.GenomicAccession,
                                                                                  self.GenomicSeq, self.GenomicDescription,
                                                                                  self.GeneID,
                                                                                  self.GC,
                                                                                    self.Domains,
                                                                                  self.Introns,
                                                                                  self.ExonLength,
                                                                                    self.taxonomy))
        except sqlite3.IntegrityError:
            pass
        self.conn.commit()
        self.conn.close()

#Checks to make sure that the database does not already contain a record before running the NCBI mining procedures
class SQLiteChecker():
    def __init__(self, proteinAccession):
        self.connect = sqlite3.connect('Records.db')
        self.proteinAccession = proteinAccession

    #Checks to make sure that the record isn't already in the table before it runs it. Checks against the accession number of the protein fasta
    def checkRecords(self):
        C = self.connect.cursor()

        C.execute('SELECT ProteinAccession FROM HSP70 WHERE ProteinAccession= (?)''',(self.proteinAccession,))

        data = C.fetchall()

        if not data:
            return False
        else:
            return True

#Makes a brand new table in the database so that different sets of proteins can be added
class makeNewTable():
    def __init__(self):
        self.connect = sqlite3.connect('Records.db')
    def makeNewTable(self):
        C = self.connect.cursor()

        '''
        Insert the desired table name where it says insert table name
        '''
        C.execute('''CREATE TABLE INSERTTABLENAME (UUID TEXT PRIMARY KEY ,ProteinAccession TEXT UNIQUE, ProteinSequence TEXT, ProteinDescription TEXT,
                                  ProteinID INTEGER, CDSAccession TEXT, CDSSeq TEXT, CDSDescription TEXT, GenomicAccession TEXT,
                                GenomicSeq TEXT, GenomicDescription TEXT, GeneID INTEGER, GenomicContext TEXT, ParentDomains TEXT,
                                Introns TEXT, ExonLength TEXT, Taxonomy TEXT)''')

        self.connect.commit()


if __name__ == '__main__':

    SQL = makeNewTable()
    SQL.makeNewTable()