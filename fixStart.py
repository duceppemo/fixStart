#!/usr/local/env python3

__author__ = 'duceppemo'
__version__ = '0.1'


from Bio import SeqIO
import os

class FixStart(object):

    def __init__(self, args):
        import multiprocessing

        # Define variables based on supplied arguments
        self.args = args
        self.assembly = args.assembly
        self.ref = args.ref
        self.out = args.output

        # Number of cpu
        self.cpus = int(multiprocessing.cpu_count())

        # Default output path is input assembly path
        if not self.out:
            self.out = os.path.dirname(self.assembly) + "/"

        self.position_list = list()  # contain tuples (contig, position, sense

        # run the script
        self.run()

    def run(self):
        self.check_dependencies()
        self.fix_start()

    def fix_start(self):
        if self.ref:
            self.find_dnaA(self.ref, self.assembly)
        else:
            self.create_default_ref()
            self.find_dnaA(self.ref, self.assembly)

        # parse blast xml file
        self.parse_blast_xml(self.out + 'blastn.xml')
        # reverse complement if needed

    def parse_blast_xml(self, xml):
        """
        test
        :param xml: Blast output file in xml format
        :return:
        """
        from Bio.Blast import NCBIXML

        with open(xml, 'r') as f:
            for record in NCBIXML.parse(f):
                if record.alignments:  # skip queries with no matches
                    print("QUERY: %s" % record.query[:60])
                    for align in record.alignments:
                        for hsp in align.hsps:
                            if hsp.expect < float('1e-10'):
                                print("MATCH: %s " % align.title[:60])
                                print(hsp.expect)

    def create_default_ref(self):
        """
        test
        :return: A fasta file in the output directory with the dnaA sequence of
                 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 (GCF_000006945)
        """

        default_ref_dnaA = ['>NC_003197.2_cds_NP_462738.1 dnaA Salmonella Typhimurium LT2',
                    'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACG',
                    'CCCGTTGCAGGCGGAACTGAGCGATAACACGCTGGCTTTGTATGCGCCAAACCGTTTTGTGCTCGATTGGGTAAGAGATA',
                    'AGTACCTCAACAATATCAATGGATTACTCAACACATTCTGCGGCGCGGATGCCCCACAACTGCGTTTTGAAGTGGGAACA',
                    'AAGCCCGTTACGCAAACGCTAAAAACGCCTGTGCATAACGTTGTCGCGCCTGCGCAGACAACAACGGCGCAGCCGCAGCG',
                    'CGTAGCGCCTGCGGCCCGTTCGGGCTGGGATAACGTACCAGCGCCAGCGGAGCCGACCTACCGCTCCAACGTCAATGTAA',
                    'AACATACATTTGATAACTTCGTCGAAGGTAAATCCAACCAACTGGCGCGCGCGGCGGCACGTCAGGTGGCGGATAATCCT',
                    'GGCGGCGCTTATAACCCGTTATTCCTCTATGGCGGCACCGGTCTGGGTAAAACTCACCTGCTGCACGCGGTGGGTAACGG',
                    'CATTATGGCGCGTAAACCCAACGCGAAAGTCGTGTATATGCACTCCGAGCGCTTTGTGCAGGACATGGTAAAAGCCCTGC',
                    'AAAATAACGCCATCGAAGAGTTTAAACGCTACTACCGCTCCGTTGACGCATTGCTGATCGACGATATTCAATTCTTCGCC',
                    'AATAAAGAACGATCCCAGGAAGAGTTTTTCCATACCTTTAACGCTCTGCTGGAAGGCAATCAGCAGATCATTTTGACGTC',
                    'GGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTCTCAAGTCCCGCTTTGGTTGGGGGCTGACAGTGGCGATCG',
                    'AACCGCCAGAACTGGAAACCCGCGTGGCGATCCTGATGAAAAAAGCGGACGAAAATGATATTCGTCTGCCAGGCGAAGTA',
                    'GCGTTCTTTATCGCCAAACGTCTACGCTCTAACGTGCGTGAACTGGAAGGCGCGTTGAACCGGGTGATTGCCAATGCCAA',
                    'CTTTACCGGCCGGGCGATTACTATCGACTTTGTGCGCGAAGCGCTGCGCGATTTACTGGCGTTGCAGGAAAAACTGGTCA',
                    'CCATCGACAATATTCAGAAGACGGTGGCCGAGTATTACAAAATCAAAATTGCGGATCTGCTTTCTAAGCGTCGATCTCGC',
                    'TCGGTAGCACGTCCGCGCCAGATGGCTATGGCGCTGGCAAAAGAGCTCACCAACCACAGTCTGCCGGAAATCGGCGATGC',
                    'GTTTGGCGGGCGCGACCACACTACCGTACTTCATGCCTGTCGTAAAATTGAGCAACTGCGTGAAGAAAGCCACGATATCA',
                    'AAGAAGATTTTTCGAATTTAATCAGAACATTGTCGTCGTGA']
        default_ref_dnaA = "\n".join(default_ref_dnaA)

        default_ref_repA = ['>NC_003277.2_cds_NP_490498.3 repA Salmonella Typhimurium LT2',
                            'GTGCAGGCAGAAGTGACTGACACCCGAACACTGTTCACTCATTACCGACAGGTCAAAAATCCGAATCCGGAATTCACGCC',
                            'GAGAGAAGGGAAAAAGACCCTGCCGTTCTGTCGTAAGCTGATGGCGAAAGCCGAAGGGTTCACGTCCCGTTTTGATTTTT',
                            'CCGTCCATGTGGCGTTCGTTCGTTCGCTGGGAAAGCGTCACCGGATGCCGCCTTTGCTGCGCCGTCGTGCCATCGATGCG',
                            'CTGCTTCAGGGGTTGTGCTTCCATTATGATCCACTGGCCAACCGTGTACAGAGATCCATCACCAATCTGGCTATAGAGTG',
                            'CGGTCTGGCCACTGAGTCAAAAAGTGGTAATCTGTCCATCACCCGCGCCACACGGGCGCTGAAGTTTATGGCAGAGCTGG',
                            'GGCTGATTACCTACCAGACCGAGTACGATCCGCAGATTGGCTGCAACATCCCGACCGATATCACCTTCACGCCGGCACTG',
                            'TTCTCCGCGCTCGACGTTTCTGACGTGGCTGTTATGGCCGCCCGCCGCAGTCGTGTCGAATGGGAAAATCAGCAGCGCAA',
                            'AAAGCAGAATCTGGAGCCGCTGGAGATGGATGAACTGATAGCGAAAGCCTGGCGATTTGTGCGCGAACGGTTCCGCAGCT',
                            'ACCAGTCTGAACGTAAACTGCACGGACTGAAACGCGCCCGGGCGCGGCGTGACGCGGACAGAAGCCGTAAGGACATCGAA',
                            'ACACTGGTGAAACAGCAGCTGACGCGGGAATACGCCAGCGGACGCTTTACCGGCGGGCTTGATGCCATGAAGCGGGAGCT',
                            'GCAGCGGCGCGTGAAAGAGCGCATGATGATGTCGCGGGGCAAAAACTACACGCGACTGACCATGGCCACCGTCCCCATAT',
                            'AA']
        default_ref_repA = "\n".join(default_ref_repA)

        out_ref_file = self.out + 'default_dnaA_repA.fasta'
        with open(out_ref_file, 'w') as out:
            out.write(default_ref_dnaA + "\n"
                      + default_ref_repA + "\n")
        self.ref = out_ref_file

    def find_dnaA(self, ref, query):
        """
        Check if dnaA sequence is present. Runs blastn.
        :param ref: fasta file with dnaA or repA sequences to look for
        :param query: assembly file in fasta format
        :return: a blast output file in xml format
        """
        import re

        self.make_blast_db(ref)
        self.run_blastn(ref, query)

        # Remove blastdb index files
        for f in os.listdir(self.out):
            pattern = os.path.basename(self.ref) + ".n*"
            if re.search(pattern, f):
                os.remove(os.path.join(self.out, f))
        # Parse blast ouput file

    def make_blast_db(self, ref):
        """
        Create a custom blast database from fasta file
        :param ref: Reference fasta file
        :return: The index files required for blast to run on the reference fasta file
        """
        import subprocess
        subprocess.run(['makeblastdb', '-in', ref, '-dbtype', 'nucl', '-parse_seqids', '-hash_index'])

    def run_blastn(self, ref_db, query):
        """
        Perform blastn using biopython
        :param ref: A fasta file for which "makeblastdb' has already been run
        :param query: The assembly, in fasta format
        :return: XML file
        """
        from Bio.Blast.Applications import NcbiblastnCommandline

        out_xml = self.out + 'blastn.xml'
        blastn_cline = NcbiblastnCommandline(db=ref_db, query=query, evalue='1e-10',
                                             outfmt=5, max_target_seqs=5,  out=out_xml,
                                             num_threads=self.cpus)
        (stdout, stderr) = blastn_cline()

        if stderr:
            raise Exception('There was a problem with the blast')

    def check_dependencies(self):
        pass


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Reorder assembly to have the dnaA gene first')
    parser.add_argument('-a', '--assembly', metavar='my_assembly.fasta',
                        required=True,
                        help='Assembly file in fasta format')
    parser.add_argument('-r', '--ref', metavar='reference_dnaA.fasta',
                        required=False,
                        help='Fasta file with dnaA sequence(s) of interest')
    parser.add_argument('-o', '--output', metavar='reordered_assembly.fasta',
                        required=False,
                        help='Output fasta file. Will default to "input_fixedstart.fasta" if omitted')

    # Get the arguments into an object
    arguments = parser.parse_args()

    FixStart(arguments)
