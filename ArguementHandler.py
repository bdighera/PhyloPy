import os, argparse, sys

class ArguementParser():

    def __init__(self):
        pass

    def parse(self):

            # Initializes the parsing function
            parser = argparse.ArgumentParser()

            # Adds the parsing arguments, sets default value, type, and help response

            parser.add_argument('-protein', default=bool(False), type=str, help='FASTA formatted protein list, .fa')
            parser.add_argument('-dir', default=os.getcwd(), type=str, help='Working Directory, Default current directory')
            parser.add_argument('-rec', default=bool(False), type=str,
                                help='Recursively runs for protein files in other directories - BETA')
            parser.add_argument('-outfile', default=os.getcwd(), type=str, help='Path to output file directory')
            parser.add_argument('-genomic', default=bool(False), type=str, help='FASTA formatted genomic list, .fa')
            parser.add_argument('-CDS', default=bool(False), type=str, help='FASTA formatted CDS list, .fa')
            parser.add_argument('-tree', default=bool(False), type=str, help='Newick Tree object, .nwk')


            args = vars(parser.parse_args())

            return args['protein'], args['dir'], args['rec'], args['outfile'], args['genomic'], args['CDS'], args['tree']

