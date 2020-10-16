from Bio import Entrez, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess, dendropy, sys, pprint, os
from randomcolor import RandomColor

from FileHandler import IntronFileHandler, MSAfileHandler, treeOBjFileHandler, GCfileHandler

class SequeuceCollector():

    def __init__(self, proteinRecord):
        self.proteinRecord = proteinRecord
        self.proteinID = ''

    def collectProteinIDs(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        proteinAccession = self.proteinRecord.id

        eSearch = Entrez.esearch(db='protein',term= proteinAccession, rettype='gb', api_key= '42c8b18ba1ceca33301e1fa5e582eed81509' )
        result = Entrez.read(eSearch, validate=False)

        proteinID = result['IdList'][0]

        self.proteinID = proteinID

        return proteinID

    def collectCDS(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        proteinID = self.proteinID

        generalEfetch= SeqIO.read(Entrez.efetch(db='protein', id=proteinID, rettype='gb', retmode='text', api_key='42c8b18ba1ceca33301e1fa5e582eed81509'), 'genbank')

        Taxonomy = generalEfetch.annotations['taxonomy']

        result = [result.qualifiers['coded_by'] for result in generalEfetch.features if result.type == 'CDS'][0][0]

        accession, codingRegion = result.split(':')
        startseq, endseq = codingRegion.split('..')

        recordEfetch = SeqIO.read(Entrez.efetch(db='nuccore', id=accession, seq_start=startseq, seq_top=endseq, rettype='fasta', api_key='42c8b18ba1ceca33301e1fa5e582eed81509'), 'fasta')

        return recordEfetch, Taxonomy

    def collectGenomicSeq(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        proteinID = self.proteinID

        elinkResult = Entrez.read(Entrez.elink(db='gene', dbfrom='protein', id=proteinID, api_key='42c8b18ba1ceca33301e1fa5e582eed81509'))
        geneID = elinkResult[0]['LinkSetDb'][0]['Link'][0]['Id']

        generalEfetch= Entrez.read(Entrez.efetch(db='gene', id=geneID, rettype='fasta', retmode='xml', api_key='42c8b18ba1ceca33301e1fa5e582eed81509'), validate=False)

        accession = generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_accession'] + '.' + generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_version']
        startseq = generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
        endseq = generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']

        genomeEfetch = SeqIO.read(Entrez.efetch(db='nuccore', id=accession, seq_start=int(startseq), seq_stop=int(endseq), rettype='fasta', api_key='42c8b18ba1ceca33301e1fa5e582eed81509'),'fasta')

        return genomeEfetch, geneID

    def collectIntrons(self, CDS, Genomic):

        CDSseq = str('>' + CDS.description + '\n' + CDS.seq)
        GenomicSeq = str('>' + Genomic.description + '\n' + Genomic.seq)

        I = IntronFileHandler()

        CDSfilePath = I.getCDSTmpPath()
        GenomicfilePath = I.getGenomicTmpPath()

        I.clearPreviousInput(CDSfilePath)
        I.clearPreviousInput(GenomicfilePath)

        I.writeNewInput(CDSfilePath, CDSseq)
        I.writeNewInput(GenomicfilePath, GenomicSeq)

        IntronCalculatorExec = I.getSpideyPath()

        intronProcessoutput = subprocess.Popen([IntronCalculatorExec,
                                           "-i", GenomicfilePath,
                                           "-m", CDSfilePath,
                                           "-p", "1"], stdout=subprocess.PIPE)

        return I.spideyOutputParser(intronProcessoutput.stdout.readlines())

    def collectProteinDomains(self):
        Entrez.email ='bdighera@csu.fullerton.edu'

        proteinAccession = self.proteinRecord.id

        e_fetch = SeqIO.parse(Entrez.efetch(db='protein', id=proteinAccession, retmax=1000, rettype='gb', retmode='fasta',api_key='42c8b18ba1ceca33301e1fa5e582eed81509'), 'gb')

        resultFeatures = [result.features for result in e_fetch][0]

        rawDomains = [result for result in resultFeatures if result.type == 'Region']

        domainLocation =[str(result.location).split('[')[1].split(']')[0] for result in rawDomains]
        domainName = [result.qualifiers['region_name'][0] for result in rawDomains]


        Domains = {domainName[i]:domainLocation[i] for i in range(len(domainLocation))}

        return Domains

class GenomicContext():
    def __init__(self, geneRecord):

        self.geneRecord = geneRecord
        self.timer = 0
        self.run_name = ''
        self.input_protein_accession_number = ['NP_005336.3']
        self.input_protein_sequence = []
        self.retrieved_protein_ids = []
        self.retrieved_full_cds = []
        self.parsed_cds_acession = []
        self.gene_accession_list = []
        self.gene_start_sequence_list = []
        self.gene_end_sequence_list = []
        self.gene = []
        self.intron_phase = []
        self.exon_lengths = []
        self.accession_dict_with_introns = {}
        self.number_protein_IDs = 0
        self.number_gene_IDs = 0
        self.number_CDS_IDs = 0
        self.msa_aligned_protein = []
        self.genomic_context = {}
        self.genomic_context_coords = ''
        self.genomic_context_motifs = {}
        self.trueGenestart = 0
        self.trueGeneend = 0
        self.trueGenedirection = ''
        self.genomic_context_gene_counter = 0
        self.protein_name = 'STRING'

    # Collects the genomic context of the parent protein
    def fetchRecord(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        geneRecord = self.geneRecord

        #Gets gene accession number
        accession = str(geneRecord.description.split(':')[0])
        #Gets the start sequence of the gene and moves 50k basepairs upstream
        startseq = int(geneRecord.description.split(' ')[0].split(':')[1].split('-')[0]) - 50000

        #Gets the end sequence of the gene and moves 50k basepairs downstream
        endseq = int(geneRecord.description.split(' ')[0].split(':')[1].split('-')[0]) + 50000

        # Begins the efetch call to acquire all the neighboring genes of the gene of interest (genomic context efetch call)

        GCfetch= Entrez.read(Entrez.efetch(db="nuccore",
                                                           id=accession,
                                                           seq_start=startseq,
                                                           seq_stop=endseq,
                                                           rettype='gb',
                                                           retmode='xml', validate=False,
                                                           api_key='42c8b18ba1ceca33301e1fa5e582eed81509'))

        return GCfetch

    # Goes through the Gene portion of the genomic context fetch, (Entrez Dictionary) and grabs the gene direction, name, and start/end point
    def parseRecord(self, GCrecord):

        #Working on rewriting code to be more pythonic, start is the next 3 lines:
        #rawCDSRecord = [record['GBFeature_quals'] for record in GCrecord[0]['GBSeq_feature-table'] if record['GBFeature_key'] == 'CDS'][0]
        #rawGeneRecord = [record for record in GCrecord[0]['GBSeq_feature-table'] if record['GBFeature_key'] == 'gene']
        #geneName = [record['GBQualifier_value'] for record in rawCDSRecord if record['GBQualifier_name'] == 'gene']

        gb_fetch = GCrecord[0]
        gene_name_list = []
        gene_start_list = []
        gene_end_list = []
        gene_id_list = []
        gene_direction_list = []
        protein_accession_dict = {}
        i = -1

        for gb_gene in gb_fetch['GBSeq_feature-table']:

            if gb_gene['GBFeature_key'] == 'CDS':

                for protein_info in gb_gene['GBFeature_quals']:

                    if protein_info['GBQualifier_name'] == 'gene':
                        name = protein_info['GBQualifier_value']

                    if protein_info['GBQualifier_name'] == 'protein_id':
                        protein_accession = protein_info['GBQualifier_value']

                        protein_accession_dict[name] = protein_accession

        for gb_gene in gb_fetch['GBSeq_feature-table']:

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

        #return gene_name_list, gene_start_list, gene_end_list, gene_direction_list, protein_accession_dict

        GC_List = []
        GC = GenomicContext(self.geneRecord)
        for i in range(len(gene_name_list)):

            completeDomains, domain = GC.fetchDomains(protein_accession_dict[gene_name_list[i]])

            GC_List.append({
                'gene_name':gene_name_list[i],
                'gene_start_seq':gene_start_list[i],
                'gene_end_seq':gene_end_list[i],
                'coding_direction':gene_direction_list[i],
                'protein_accession':protein_accession_dict[gene_name_list[i]],
                'domain': completeDomains[0].values()[0]

            })

        return GC_List

    # Domain Fetch for Genomic Context
    def fetchDomains(self, gc_accession):

        Complete_Domains = []
        domainNameList = []

        e_fetch = SeqIO.parse(
            Entrez.efetch(db='protein', id="%s" % gc_accession, retmax=1000, rettype='gb', retmode='fasta',
                          api_key='42c8b18ba1ceca33301e1fa5e582eed81509'), 'gb')

        for seq_record in e_fetch:
            domain_list = []
            accession_number = seq_record.id

            for i in range(len(seq_record.features)):
                if seq_record.features[i].type == 'Region':
                    domain_location = str(seq_record.features[i].location).split('[')[1].split(']')[0]
                    domain_name = str(seq_record.features[i].qualifiers['region_name'][0])
                    domainNameList.append(domain_name)

                    domain_list.append([domain_name, domain_location])

            Complete_Domains.append(dict([(accession_number, domain_list)]))

        rand_color = RandomColor()
        Domains = [domain for domain in set(domainNameList)]

        return Complete_Domains, Domains

    def domain_colors(self, completed_domains):

        parent_protein_domains = completed_domains

        # pprint.pprint(parent_protein_domains)

        domain_list = []

        # Goes in and pulls out all of the individual domains for all of the proteins in the genomic context. List is cleaned removing duplicate domains from GC proteins
        for child_protein_domains in parent_protein_domains:

            for each_domain_list in child_protein_domains:

                if each_domain_list['domain'] != []:

                    for i in range(len(each_domain_list['domain'])):
                        domain_list.append(each_domain_list['domain'][i])

        domain_list = list(dict.fromkeys(domain_list))

        rand_color = RandomColor()

        domains_dict_colors = {domain: rand_color.generate()[0] for domain in set(domain_list)}

        # Assigns the correct color to the domains. This makes it so all gc domains are assigned the same color based on their name
        for child_protein_domains in parent_protein_domains:

            for each_domain_list in child_protein_domains:

                if each_domain_list['domain'] != []:

                    color_list = []
                    for i in range(len(each_domain_list['domain'])):
                        domain = each_domain_list['domain'][i]
                        correct_color = domains_dict_colors.get(domain)
                        domain_color_pair = {domain: correct_color}

                        color_list.append(domain_color_pair)

                        each_domain_list['color'] = color_list

                else:

                    each_domain_list['color'] = []

        return parent_protein_domains

