from Bio.Align.Applications import ClustalOmegaCommandline
import itertools, os, subprocess, sys, math, dendropy, pprint, ast
from randomcolor import RandomColor
from FileHandler import MSAfileHandler, treeOBjFileHandler, ImageProcessingHandler
from IPython.display import display


from ete3 import Tree, SeqMotifFace, TreeStyle, TextFace

class PhyloTreeConstruction(object):

    def __init__(self, proteinAccession, proteinSeq, proteinDescription, GenomicContext, ParentDomains, Introns, ExonLenghts, JupyterNoteBookFigure, prune,
                 printLeafAccession, printLeafDesc, commonNames, GeneID):
        self.proteinAccessions = proteinAccession
        self.proteinSeqs = proteinSeq
        self.proteinDescs = proteinDescription
        self.GenomicContexts = GenomicContext
        self.parentDomains = ParentDomains
        self.Introns = Introns
        self.exonLengths = ExonLenghts
        self.JNfig = JupyterNoteBookFigure
        self.commonNames = commonNames
        self.collectMultipleSequencingAlignment()
        self.rootedTreeConstruction()
        self.prune = prune
        self.printLeafAccession = printLeafAccession
        self.printLeafDesc = printLeafDesc
        self.geneID = GeneID


    def rootedTreeConstruction(self):

        in_file = os.path.join('execs', 'tmp', 'aligned.tmp')
        out_file = os.path.join('execs', 'tmp', 'unrooted_tree.nwk')

        subprocess.call(["./execs/FastTree", "-out", out_file, in_file])

        if self.JNfig == False:
            print('\n' + subprocess.list2cmdline(["./execs/FastTree", in_file, ">", out_file]))

        rooted_tree = dendropy.Tree.get_from_path(out_file, schema='newick')
        rooted_tree.reroot_at_midpoint()
        rooted_tree.write_to_path(os.path.join('execs', 'tmp', "rooted_tree.nwk"), schema='newick')

        newick_rooted_file = open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'r')
        read_edit_newick = newick_rooted_file.read()

        # The tree generated sometimes has the [&R] region - if not stripped it will throw error. Try except handles if the [&R] is not generated
        try:

            stripped_tree = read_edit_newick.strip('\[&R\] ')
            with open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'w') as writeStrippedTree:
                writeStrippedTree.write('')
                writeStrippedTree.write(stripped_tree)

                with open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'w') as writeStrippedTree:
                    writeStrippedTree.write('')
                    writeStrippedTree.write(stripped_tree)

        except AttributeError:
            pass

        if self.JNfig == False:
            print('\n' + 'ROOTED TREE HAS BEEN CONSTRUCTED...')

    def collectMultipleSequencingAlignment(self):

        MSA = MSAfileHandler()

        in_file = MSA.getUnalignedTmpPath()
        out_file = MSA.getAlignedTmpPath()

        MSA.clearPreviousInput(in_file)
        MSA.clearPreviousInput(out_file)

        protein_description = self.proteinDescs
        protein_sequence = self.proteinSeqs
        common_names = self.commonNames

        with open(in_file, 'a') as msaprotein_writeFile:
            for i in range(len(protein_description)):

                protein = '\n'+ '>' + str(protein_description[i]) + str(common_names[i]).replace(' ', '_') + '\n' + str(protein_sequence[i])
                msaprotein_writeFile.write(protein)


        if self.JNfig == False:
            print('\n' + 'CREATING MULTIPLE SEQUENCING ALIGNMENT...')

        clustalomega_cline = ClustalOmegaCommandline(cmd=os.path.join('execs', "clustalo-1.2.0"),
                                                     infile=in_file,
                                                     outfile=out_file, verbose=True, auto=True, force=True)
        clustalomega_cline()
        MSA.msa_FileCorrection()

        if self.JNfig == False:
            print('\n' + 'MULTIPLE SEQUENCING ALIGNMENT HAS BEEN CREATED...')

    def constructTreeObj(self):

        tree = treeOBjFileHandler()

        MSAfilePath = tree.getTreeInputPath()
        unrootedTreePath = tree.getTreeOutputPath()

        subprocess.call(["./execs/FastTree", "-out", MSAfilePath, unrootedTreePath])

        rootedTree = dendropy.Tree.get_from_path(unrootedTreePath, schema='newick')

        rootedTree.reroot_at_midpoint()

    def assignDomainColors(self, Domains):

        #Iterates through the domains taken from the database, removes the duplicates, and assigns a random color which will be used in the final figure
        rawDomainNames = [domain.keys() for domain in Domains]
        domainNames = []
        for domain in rawDomainNames:
            for eachDomain in domain:
                domainNames.append(eachDomain)

        domainNames = list(dict.fromkeys(domainNames))

        randomColor = RandomColor()
        domainNameColors = {domain:randomColor.generate()[0] for domain in domainNames}

        return domainNameColors

    def renderingTreewithDomains(self):

        proteinAccessions = self.proteinAccessions
        parentDomains = self.parentDomains

        treeObj = treeOBjFileHandler()
        treeObj.getRootedTreePath()

        with open(treeObj.getRootedTreePath()) as nwkTreeFile:
            nwkTree = nwkTreeFile.read()
            dt = Tree(nwkTree)

        dts = TreeStyle()
        dts.title.add_face(TextFace('PhyloPy - Protein Ortholog Finding Tool by Bryan Dighera: Protein Domains', fsize= 16,), column= 0)
        dts.allow_face_overlap = True
        dts.show_leaf_name = True
        dts.show_branch_support = True

        leafNames = dt.get_leaf_names()

        accessionDomains = {proteinAccessions[i]: parentDomains[i] for i in range(len(leafNames))}

        domainColors = self.assignDomainColors(accessionDomains.values())

        domainMotifs = []

        # The leaf names contain the description so the accession must be stripped in order to index with protein accessions from db
        leafAccessionExtracted = treeObj.getProteinAccession([leaf for leaf in leafNames])

        for leaf in leafAccessionExtracted:

            domains = accessionDomains[leaf]
            domainLen = len(accessionDomains[leaf])

            domainName = [list(domains.keys())[i] for i in range(domainLen)]
            domainStart = [list(domains.values())[i].split(':')[0] for i in range(domainLen)]
            domainEnd = [list(domains.values())[i].split(':')[1] for i in range(domainLen)]

            domainColor = [domainColors[domain] for domain in domains]

            leafMotfis = [[int(domainStart[i].strip('<')), int(domainEnd[i].strip('>')), "<>", None, 12, "Black", domainColor[i], "arial|4|white|%s" % str(domainColor[i])] for i in range(domainLen)]

            domainMotifs.append(leafMotfis)

        IPH = ImageProcessingHandler()

        MSASeqLen = IPH.intron_fix(leafNames[0], None)[1]

        domainSeqFace = [SeqMotifFace(' ' * MSASeqLen, gapcolor="black", seq_format='line', scale_factor=1,
                                      motifs=domainMotifs[i]) for i in range(len(domainMotifs))]

        for i in range(len(domainSeqFace)):

            (dt & leafNames[i]).add_face(domainSeqFace[i], 0, "aligned")

        if self.JNfig == True:
            display(dt.render('%%inline', tree_style=dts))
        elif self.JNfig == False:
            dt.show(tree_style=dts)

    def renderingTreewithIntrons(self):

        proteinAccessions = self.proteinAccessions
        introns = self.Introns
        exonLengths = self.exonLengths

        treeObj = treeOBjFileHandler()
        treeObj.getRootedTreePath()

        with open(treeObj.getRootedTreePath()) as nwkTreeFile:
            nwkTree = nwkTreeFile.read()
            t = Tree(nwkTree)
            nwkTreeFile.close()

        ts = TreeStyle()
        ts.title.add_face(
            TextFace('PhyloPy - Protein Ortholog Finding Tool by Bryan Dighera: Intron Location and Phases',
                     fsize=16, ), column=0)
        ts.allow_face_overlap = True
        ts.show_leaf_name = True
        ts.show_branch_support = True

        leafNames = t.get_leaf_names()

        accessionIntrons = {proteinAccessions[i]: [introns[i], exonLengths[i]] for i in range(len(leafNames))}

        dummyIntronMotif = [[0, 400, "-", None, 12, "Black", "Black", None]]
        MSASeqlen = 0
        intronMotifs = []

        if self.printLeafDesc != False:
            for leaf in leafNames:
                print(leaf)


        #The leaf names contain the description so the accession must be stripped in order to index with protein accessions from db
        leafAccessionExtracted = treeObj.getProteinAccession([leaf for leaf in leafNames])

        for leaf in leafAccessionExtracted:  # Corrects introns, and builds intron motifs

            if self.printLeafAccession != False:
                print(leaf)

            #leaf = str(leaf.split('_')[0] + '_' + leaf.split('_')[1])
            intronPhases = accessionIntrons[leaf][0]
            exonLengths = accessionIntrons[leaf][1]


            if intronPhases and exonLengths != 'NONE':

                IPH = ImageProcessingHandler()

                intronPhases = ast.literal_eval(intronPhases)
                exonLengths = ast.literal_eval(exonLengths)

                exonLength = [math.floor(int(exonLengths[i][0].split('-')[1]) / 3) for i in range(len(exonLengths))]

                try:
                    exonLocation, MSASeqlen = IPH.intron_fix(leaf, exonLength)

                except IndexError:
                    print('Index Error Found at: %s. Please remove and rerun.' % leaf)

                recordMotifs = []

                for i in range(len(exonLocation)):

                    intronPhase = int(intronPhases[i]) % 3

                    if intronPhase == 0:
                        recordMotifs.append([exonLocation[i] - 5,
                                             exonLocation[i] + 5,
                                             "[]", None, 5, "Grey", "Grey", "arial|4|white|%s" % str(exonLocation[i])])

                    elif intronPhase == 1:
                        recordMotifs.append([exonLocation[i] - 5,
                                             exonLocation[i] + 5,
                                             "[]", None, 5, "Black", "Black",
                                             "arial|4|white|%s" % str(exonLocation[i])])

                    elif intronPhase == 2:
                        recordMotifs.append([exonLocation[i] - 5,
                                             exonLocation[i] + 5,
                                             "[]", None, 5, "Blue", "Blue", "arial|4|white|%s" % str(exonLocation[i])])

                    else:
                        recordMotifs.append(dummyIntronMotif)

                intronMotifs.append(recordMotifs)


            else:
                intronMotifs.append(dummyIntronMotif)  # Adds the intron mo #A

        intronSeqFace = [SeqMotifFace('-' * MSASeqlen, gapcolor="black", seq_format='line', scale_factor=1,
                                      motifs=intronMotifs[i]) for i in range(len(intronMotifs))]

        for i in range(len(intronSeqFace)):
            (t & leafNames[i]).add_face(intronSeqFace[i], 0, "aligned")

        intronPhases = 3

        ts.legend.add_face(
            SeqMotifFace("A" * 1, [[0, 80, "[]", None, 8, "Gray", 'Gray', "arial|4|white|%s" % str('Phase 0')]]),
            column=0)

        ts.legend.add_face(
            SeqMotifFace("A" * 1, [[0, 80, "[]", None, 8, "Black", 'Black', "arial|4|white|%s" % str('Phase 1')]]),
            column=0)

        ts.legend.add_face(
            SeqMotifFace("A" * 1, [[0, 80, "[]", None, 8, "Blue", 'Blue', "arial|4|white|%s" % str('Phase 2')]]),
            column=0)

        if self.JNfig == True:
            display(t.render('%%inline', tree_style=ts))
        elif self.JNfig == False:
            t.show(tree_style=ts)

    def renderingTreewithGenomicContext(self):

        #Set parent proteins and genomic context retrieved from SQLite into a variable
        proteinAccessions = self.proteinAccessions
        parentGC = self.GenomicContexts

        accessionGCdict = {proteinAccessions[i]:parentGC[i] for i in range(len(proteinAccessions))}

        #Largest length of the genomic context
        maxGClength = max([len(parentGC[i]) for i in range(len(parentGC))])

        #Strip all the domains from the entire datastructure so that all domains are stored in a single list
        GCdomains = itertools.chain(*[[parentGC[i][j]['domain'] for j in range(len(parentGC[i]))] for i in range(len(proteinAccessions))])
        GCdomains = list(itertools.chain(*GCdomains))
        GCdomains = list(itertools.chain(*GCdomains))[::2]


        #Assign each domain a color, as key value pair (dict), which will be assigned during motif construction
        rand_color = RandomColor()
        GCcolors = {GCdomains[i]:rand_color.generate()[0] for i in range(len(GCdomains))}


        treeObj = treeOBjFileHandler()
        treeObj.getRootedTreePath()

        with open(treeObj.getRootedTreePath()) as nwkTreeFile:
            nwkTree = nwkTreeFile.read()
            t = Tree(nwkTree)
            nwkTreeFile.close()

        ts = TreeStyle()
        ts.title.add_face(
            TextFace('PhyloPy - Protein Ortholog Finding Tool by Bryan Dighera: Genomic Context',
                     fsize=16, ), column=0)
        ts.allow_face_overlap = True
        ts.show_leaf_name = True
        ts.show_branch_support = True

        leafNames = t.get_leaf_names()

        GCMotifs = []

        # The leaf names contain the description so the accession must be stripped in order to index with protein accessions from db
        leafAccessionExtracted = treeObj.getProteinAccession([leaf for leaf in leafNames])

        geneIDs = {self.proteinAccessions[i]:self.geneID[i] for i in range(len(leafNames))}


        for j, leaf in enumerate(leafAccessionExtracted):

            try:
                # Function which takes the current leaf genomic context and aligns all working genes in same direction
                record = treeObj.fix_coding_direction(accessionGCdict[leaf], geneIDs[leaf])
                record = treeObj.fix_coding_alignment(record, geneIDs[leaf], maxGClength)

                coding_direction = [record[i]['coding_direction'] for i in range(len(record))]
                geneName = [record[i]['gene_name'] for i in range(len(record))]
                numberofDomains = [record[i]['domain'] for i in range(len(record))]
                flip = [record[i]['flip'] for i in range(len(record))]

                numberofGenes = len([record[i]['img_start'] for i in range(len(record))])
                start_gene_location = [record[i]['img_start'] for i in range(len(record))]
                end_gene_location = [record[i]['img_end'] for i in range(len(record))]

                #Code below is the old verion of gene collection and location
                #numberofGenes = len([record[i]['gene_start_seq'] for i in range(len(record))])
                #start_gene_location = [i for i in range(10, int(numberofGenes * 100 + 100), 100)]
                #end_gene_location = [i for i in range(90, int(numberofGenes * 100 + 200), 100)]

                recordMotifs = []

                for i in range(numberofGenes):

                    if i != None:

                        if coding_direction[i] == '-' and flip[i] == False:


                            genomic_context_motif = [start_gene_location[i], end_gene_location[i], "[]", 12, 12, "Black", "White", "arial|1|black|%s" % geneName[i]]
                            direction_motif = [start_gene_location[i], int(start_gene_location[i]) - 10, ">", 12, 12,
                                               "Black", "Black", None]

                            recordMotifs.append(genomic_context_motif)


                            start_domain_location = [i for i in range(start_gene_location[i], end_gene_location[i]-10, 5)]
                            end_domain_location = [i for i in range(start_gene_location[i]+10, end_gene_location[i], 5)]

                            #TODO: Make it so that the domain name is properly parsed into the color square, stripping @ | might not be the best option
                            domainMotif = [[start_domain_location[j], end_domain_location[j] - 5, "[]", 12, 12, GCcolors[k[0]],
                                            GCcolors[k[0]], "arial|1|black|%s" % k[0].split('|')[0]] for j,k in enumerate(numberofDomains[i])]

                            for motif in domainMotif:
                                recordMotifs.append(motif)

                            recordMotifs.append(direction_motif)

                        elif coding_direction[i] == '+' and flip[i] == True:

                            genomic_context_motif = [start_gene_location[i], end_gene_location[i], "[]", 12, 12, "Black", "White", "arial|1|black|%s" % geneName[i]]
                            direction_motif = [int(start_gene_location[i]), int(start_gene_location[i]) - 10, ">", 12, 12,
                                               "Black", "Black", None]

                            recordMotifs.append(genomic_context_motif)


                            start_domain_location = [i for i in range(start_gene_location[i], end_gene_location[i] - 10, 5)]
                            end_domain_location = [i for i in range(start_gene_location[i] + 10, end_gene_location[i], 5)]

                            #TODO: Make it so that the domain name is properly parsed into the color square, stripping @ | might not be the best option
                            domainMotif = [[start_domain_location[j], end_domain_location[j] - 5, "[]", 12, 12, GCcolors[k[0]],
                                            GCcolors[k[0]], "arial|1|black|%s" % k[0].split('|')[0]] for j, k in enumerate(numberofDomains[i])]

                            for motif in domainMotif:
                                recordMotifs.append(motif)

                            recordMotifs.append(direction_motif)


                        elif coding_direction[i] == '-' and flip[i] == True:

                            genomic_context_motif = [start_gene_location[i], end_gene_location[i], "[]", 12, 12, "Black", "White", "arial|1|black|%s" % geneName[i]]
                            direction_motif = [end_gene_location[i], int(end_gene_location[i]) + 10, ">", 12, 12,
                                               "Black", "Black", None]

                            recordMotifs.append(genomic_context_motif)


                            start_domain_location = [i for i in range(start_gene_location[i], end_gene_location[i]-10, 5)]
                            end_domain_location = [i for i in range(start_gene_location[i]+10, end_gene_location[i], 5)]

                            #TODO: Make it so that the domain name is properly parsed into the color square, stripping @ | might not be the best option
                            domainMotif = [[start_domain_location[j], end_domain_location[j] - 5, "[]", 12, 12, GCcolors[k[0]],
                                            GCcolors[k[0]], "arial|8|black|%s" % k[0].split('|')[0]] for j,k in enumerate(numberofDomains[i])]

                            for motif in domainMotif:
                                recordMotifs.append(motif)

                            recordMotifs.append(direction_motif)

                        elif coding_direction[i] == '+' and flip[i] == False:

                            genomic_context_motif = [start_gene_location[i], end_gene_location[i], "[]", 12, 12, "Black", "White", "arial|1|black|%s" % geneName[i]]
                            direction_motif = [end_gene_location[i], int(end_gene_location[i]) + 10, ">", 12, 12,
                                               "Black", "Black", None]

                            recordMotifs.append(genomic_context_motif)


                            start_domain_location = [i for i in range(start_gene_location[i], end_gene_location[i] - 10, 5)]
                            end_domain_location = [i for i in range(start_gene_location[i] + 10, end_gene_location[i], 5)]

                            #TODO: Make it so that the domain name is properly parsed into the color square, stripping @ | might not be the best option
                            domainMotif = [[start_domain_location[j], end_domain_location[j] - 5, "[]", 12, 12, GCcolors[k[0]],
                                            GCcolors[k[0]], "arial|1|black|%s" % k[0].split('|')[0]] for j, k in enumerate(numberofDomains[i])]

                            for motif in domainMotif:
                                recordMotifs.append(motif)

                            recordMotifs.append(direction_motif)
                    else:
                        dummyIntronMotif = [0, 0, "[]", None, 12, "White", "White", None]
                        recordMotifs.append(dummyIntronMotif)

                GCMotifs.append(recordMotifs)

            except IndexError:
                print('Genomic Context Index Error at Sequence: %s' % leaf)
        GCSeqFace = [SeqMotifFace(gapcolor='white', seq_format='line', scale_factor=1,
                                      motifs=GCMotifs[i]) for i in range(len(GCMotifs))]

        for i in range(len(GCSeqFace)):
            (t & leafNames[i]).add_face(GCSeqFace[i], 0, "aligned")


        if self.JNfig == True:
            ts.scale=120
            display(t.render('%%inline', tree_style=ts, dpi=800))
        elif self.JNfig == False:
            t.show(tree_style=ts)
