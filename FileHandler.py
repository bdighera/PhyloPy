import os, re, sys
from Bio import SeqIO, Seq, SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class IntronFileHandler():

    def __init__(self):
        pass

    def getCDSTmpPath(self):
        return os.path.join('execs', 'tmp', 'CDS.tmp')

    def getGenomicTmpPath(self):
        return os.path.join('execs', 'tmp', 'Genomic.tmp')

    def getSpideyPath(self):
        return os.path.join(os.path.join('execs', 'spidey.linux.64'))

    def clearPreviousInput(self, path):
        with open(path, 'w') as fileHandler:
            fileHandler.write('')
            fileHandler.close()

    def writeNewInput(self, path, seq):
        with open(path, 'w') as fileHandler:
            fileHandler.write(seq)
            fileHandler.close()

    def spideyOutputParser(self, rawSpideyOutput):

        throwaway_lines = rawSpideyOutput[3:]
        number_of_introns = throwaway_lines[1].split()[3]
        introns = throwaway_lines[2:2 + int(number_of_introns)]
        intronPhases = []
        exonLengths = []

        if int(number_of_introns) > 1:

            for exon_bit in introns[:-1]:
                intron_location = exon_bit.split()[4].split()[0].split('-')
                intron_start = int(intron_location[0])
                intron_end = int(intron_location[1])
                intron_phase = (intron_end - intron_start) % 3

                intronPhases.append(intron_phase)
                # print 'Intron Phase(s): ' + str(intron_phase)
                exonLengths.append(exon_bit.split()[4].split())
                # print 'Exon Length(s): ' + str((exon_bit.split()[4].split()))


        else:
            return 'NONE', 'NONE'


        return intronPhases, exonLengths

class MSAfileHandler():

    def __init__(self):
        pass

    def getUnalignedTmpPath(self):
        return os.path.join('execs', 'tmp', 'unaligned.tmp')

    def getAlignedTmpPath(self):
        return os.path.join('execs', 'tmp', 'aligned.tmp')

    def clearPreviousInput(self, path):
        with open(path, 'w') as fileHandler:
            fileHandler.write('')
            fileHandler.close()

    def writeNewInputFile(self, path, seqRecords):
        SeqIO.write(seqRecords, path, 'fasta')

    def writeNewOutputFile(self):
        pass

    def getClustaloExecPath(self):
        return os.path.join('execs', "clustalo-1.2.0")
      
    def msa_FileCorrection(self):

        # the purpose of this function is to convert the MSA protein file into one that can include the description


        path = os.path.join('execs', 'tmp', 'aligned.tmp')
        msa_input_file = SeqIO.parse(path, 'fasta')

        genus_species_pattern = re.compile('\[\w+ \w+\]')
        accession_pattern = re.compile('[XNP_]+\d*.\d')

        sequence_list = []

        for item in msa_input_file:
            accession_number = item.description.replace(' ', '_').replace(':', '_').replace('[', '_').replace(']',
                                                                                                              '_').replace(
                'LOW_QUALITY_PROTEIN', '_')

            # accession_search = re.search(accession_pattern, str(item))
            # accession_number = accession_search.group()
            #
            # genus_species_search = re.search(genus_species_pattern, str(item))
            # genus_species = genus_species_search.group()
            #
            # new_msa_input = str('>' + accession_number + '_' + genus_species + '\n' + item.seq + '\n')
            new_msa_input = str('>' + accession_number + '\n' + item.seq + '\n')
            sequence_list.append(new_msa_input)

        with open(path, 'w') as clear_file:
            clear_file.write('')
            clear_file.close()
        for item in sequence_list:
            with open(path, 'a') as msa_file_rewrite:
                msa_file_rewrite.write(item)
                msa_file_rewrite.close()


class treeOBjFileHandler():

    def __init__(self):
        pass

    def getTreeOutputPath(self):
        return os.path.join('execs','tmp', 'unrootedTree.nwk')

    def getTreeInputPath(self):
        return os.path.join('execs', 'tmp', 'aligned.tmp')

    def getRootedTreePath(self):
        return os.path.join('execs','tmp', "rooted_tree.nwk")


    def getProteinAccession(self, descriptionRecords):
        # This method will use regular expression take the input of a list of protein descriptions and remove the accession number
        '''
        :input: list of protein descriptions ie. [NP_005338.1_HSPA5_Homo_sapians, etc.]
        :return:list of protein accessions ie. [NP_005338.1, etc. ]
        '''
        accession_pattern = re.compile('[XN][P][_ ]\d+.\d')

        recordList = descriptionRecords

        accessionList = []

        for record in recordList:
            accession_expression = re.search(accession_pattern, '%s' % (record))
            accession_expression_out = accession_expression.group()

            accessionList.append(accession_expression_out.replace(' ', '_'))

        return accessionList


class GCfileHandler():

    def __init__(self):
        pass

    def fetchGenes(self, gbEntrezResultObj):
        gene_name_list = []
        gene_start_list = []
        gene_end_list = []
        gene_id_list = []
        gene_direction_list = []
        protein_accession_dict = {}
        i = -1

        for gb_gene in gbEntrezResultObj['GBSeq_feature-table']:

            if gb_gene['GBFeature_key'] == 'CDS':

                for protein_info in gb_gene['GBFeature_quals']:

                    if protein_info['GBQualifier_name'] == 'gene':
                        name = protein_info['GBQualifier_value']

                    if protein_info['GBQualifier_name'] == 'protein_id':
                        protein_accession = protein_info['GBQualifier_value']

                        protein_accession_dict[name] = protein_accession

        for gb_gene in gbEntrezResultObj['GBSeq_feature-table']:

            if gb_gene['GBFeature_key'] == 'gene':

                for gb_info in gb_gene['GBFeature_quals']:

                    if gb_info['GBQualifier_name'] == 'gene':
                        name = str(gb_info['GBQualifier_value'])

                        try:
                            if protein_accession_dict[name]:

                                gene_name_list.append(name)
                                i += 1
                                start = int(gb_gene['GBFeature_intervals'][0]['GBInterval_from'])
                                gene_start_list.append(start)
                                end = int(gb_gene['GBFeature_intervals'][0]['GBInterval_to'])
                                gene_end_list.append(end)

                                if start < end:
                                    direction = "+"
                                    gene_direction_list.append(direction)
                                else:
                                    direction = "-"
                                    gene_direction_list.append(direction)

                                try:
                                    if 'GeneID:' in gb_info['GBQualifier_value']:
                                        gene_id = int(gb_info['GBQualifier_value'].replace('GeneID:', ''))
                                        gene_id_list.append(gene_id)
                                except:
                                    print("No GBQualifier_value")
                        except:
                            pass

        return gene_name_list, gene_start_list, gene_end_list, gene_direction_list, protein_accession_dict

#TODO: need to make this class work more efficiently
class ImageProcessingHandler():

    def __init__(self):
        pass

    def intron_fix(self, working_seq, uncalculated_intron):
        # aligned sequence counter is the number of the leaf - or the number of the protein being worked on
        # uncalculated intron is the intron output from spidey then divided by 3 and rounded either up or down (it is in context of aa)

        # bring in the msa file and get a list of all of the sequences

        MSA = MSAfileHandler()
        msa_aligned_file = MSA.getAlignedTmpPath()
        msa_seqs = [str(item.seq) for item in SeqIO.parse(msa_aligned_file, 'fasta')]
        msa_unparsed_accession = [str(item.id) for item in SeqIO.parse(msa_aligned_file, 'fasta')]

        msa_accession = self.accession_parser(msa_unparsed_accession)

        msa_dict = dict(zip(msa_accession, msa_seqs))

        if uncalculated_intron != None:
            i = 1
            n = 0
            aa_counter = 0
            gap_counter = 0
            old_introns = uncalculated_intron
            newPositions = []

            for old_intron in old_introns:

                while old_intron > i:

                    msa_sequence = str(msa_dict.get(working_seq))

                    char = msa_sequence[n]
                    n += 1
                    if char == '-':
                        old_intron += 1
                        gap_counter += 1

                    if char != '-':
                        aa_counter += 1
                        i += 1

                    continue

                newPositions.append(aa_counter + gap_counter)

            return newPositions, len(str(msa_dict.get(working_seq)))

        return None, len(str(msa_dict.get(working_seq)))

    def accession_parser(self, accession_description_list):

        accession_list = []
        for accession_description in accession_description_list:
            np_xp_pattern = re.compile('N[P]|X[P]')
            digits_pattern = re.compile('\d+.\d')

            np_xp_search_obj = re.search(np_xp_pattern, accession_description)
            digits_search_obj = re.search(digits_pattern, accession_description)

            np_xp_name = np_xp_search_obj.group()
            digits_name = digits_search_obj.group()
            final_accession = str(np_xp_name + '_' + digits_name)

            accession_list.append(final_accession)

        return accession_list

    def intronHomogenizer(self, exonLocation):
        pass



