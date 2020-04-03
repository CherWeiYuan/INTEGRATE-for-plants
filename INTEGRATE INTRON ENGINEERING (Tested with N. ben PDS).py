from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from sys import exit


#########################
#User-defined parameters#
#########################

# Which exon to target? (To target first exon, target_exon = 1)
target_exon = 2

# Parameters
insertion_pos = 49 
target_seq_len = 32
pam_len = 2

# Cargo
# Remember to add stop codon; ensure no BsaI or BpiI sites; sequence provided here is eGFP)
# What frame is the last amino acid in exon in? If in frame 1, use EXON IN FRAME 1 sequence.
# Example: Last amino acid in exon is in frame 1, such as 'AAA'. Use EXON IN FRAME 1.
# Last amino acid in exon is in frame 2, such as 'A'. Use EXON IN FRAME 2.
# Last amino acid in exon is in frame 3, such as 'AA'. Use EXON IN FRAME 3.
# add hastag to undesired insert frame, leaving only one for python to use

'''EXON IN FRAME 1'''
cargo = Seq('TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAGgtaatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcagctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactcacggcatggacgagctgtacaagtaaggatccccgGAAATGACCTCTTCATGAAAATGTCGTAAACTTACTGAAGTTTAGACCATGAAGAGGCAATATCAATTTATGGGTGTGATAATTATCAATTTATGGGTGTAATTATCATTTTATGGTTGTATCAACA',
            IUPAC.unambiguous_dna)

'''EXON IN FRAME 2'''
cargo = Seq('TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAGgtaAGatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcagctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactcacggcatggacgagctgtacaagtaaggatccccgGAAATGACCTCTTCATGAAAATGTCGTAAACTTACTGAAGTTTAGACCATGAAGAGGCAATATCAATTTATGGGTGTGATAATTATCAATTTATGGGTGTAATTATCATTTTATGGTTGTATCAACA',
            IUPAC.unambiguous_dna)

'''EXON IN FRAME 3'''
cargo = Seq('TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAGgtaAatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcagctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactcacggcatggacgagctgtacaagtaaggatccccgGAAATGACCTCTTCATGAAAATGTCGTAAACTTACTGAAGTTTAGACCATGAAGAGGCAATATCAATTTATGGGTGTGATAATTATCAATTTATGGGTGTAATTATCATTTTATGGTTGTATCAACA',
            IUPAC.unambiguous_dna)

# Intron seq without the first 6 conserved bp, which is usually 'GTAAGT' of which only GT is crucial.
intron = Seq('TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAG')

# Reading Genbank / Fasta file
# For Genbank file: gb = SeqIO.read("Nicotiana benthamiana PDS (Niben101Scf01283g02002.1).gb", "gb")
# For Fasta file: gb = SeqIO.read("Nicotiana benthamiana PDS (Niben101Scf01283g02002.1).fasta", "fasta")
gb = SeqIO.read("Nicotiana benthamiana PDS (Niben101Scf01283g02002.1).fasta", "fasta")


#Manually adding exon locations (if genbank files does not have exons)

''' 
How to find positions? Go to Genbank browser for TAIR10 genome, hover cursor over the desired gene, click 'Genbank View'. 
Scroll down and you see something like this:
    
CDS             join(552..719,808..940,

Then,  write:
seqfeat = [SeqFeature(FeatureLocation(551, 719), type="exon"),
           SeqFeature(FeatureLocation(807, 940), type="exon")]

IMPT: Note that instead of 552, we write 551.
'''

manual_annotation = True #change to False if annotation is available in gb file. If False, ignore seqfeat variable (leave it as it is, doesn't affect the script)

#seqfeat for Arabidopsis thaliana PDS3
'''
seqfeat = [SeqFeature(FeatureLocation(551, 719), type="exon"),
           SeqFeature(FeatureLocation(807, 940), type="exon"),
           SeqFeature(FeatureLocation(1035, 1124), type="exon"),
           SeqFeature(FeatureLocation(1566, 1624), type="exon"),
           SeqFeature(FeatureLocation(1801, 1957), type="exon"),
           SeqFeature(FeatureLocation(2521, 2670), type="exon"),
           SeqFeature(FeatureLocation(2898, 3018), type="exon"),
           SeqFeature(FeatureLocation(3306, 3520), type="exon"),
           SeqFeature(FeatureLocation(3602, 3705), type="exon"),
           SeqFeature(FeatureLocation(3800, 3848), type="exon"),
           SeqFeature(FeatureLocation(4023, 4069), type="exon"),
           SeqFeature(FeatureLocation(4220, 4408), type="exon"),
           SeqFeature(FeatureLocation(4501, 4664), type="exon"),
           SeqFeature(FeatureLocation(4829, 4895), type="exon"),       
           ]
'''


#seqfeat for Nicotiana benthamiana PDS
seqfeat = [SeqFeature(FeatureLocation(877, 1099), type="exon"),
           SeqFeature(FeatureLocation(1187, 1320), type="exon"),
           SeqFeature(FeatureLocation(2283, 2372), type="exon"),
           SeqFeature(FeatureLocation(2944, 2998), type="exon"),
           SeqFeature(FeatureLocation(3239, 3372), type="exon"),
           SeqFeature(FeatureLocation(3662, 3811), type="exon"),
           SeqFeature(FeatureLocation(3974, 4094), type="exon"),
           SeqFeature(FeatureLocation(4258, 4472), type="exon"),
           SeqFeature(FeatureLocation(4564, 4667), type="exon"),
           SeqFeature(FeatureLocation(4936, 4984), type="exon"),
           SeqFeature(FeatureLocation(5853, 5899), type="exon"),
           SeqFeature(FeatureLocation(7032, 7220), type="exon"),
           SeqFeature(FeatureLocation(7325, 7488), type="exon"),       
           ]


###########
#Functions#
###########

# Extract exons from either manual_annotation (if True), or Genbank file (if manual_annotation = False)
def extract_exons(gb, seqfeat, manual_annotation): #returns a list of exons
    try: 
        exon_list = []
        if manual_annotation == True: #uses seqfeat to extract features
            for i in seqfeat:
                exon_list += [i.extract(gb.seq)]
            return exon_list
        elif manual_annotation == False: #uses annotation downloaded in gb file to extract features
            for feature in gb.features:
                if feature.type == "exon":
                    start = feature.location.start.position
                    end = feature.location.end.position #plus six accounts for first 6 nucleotide of intron, 
                    #such that the cargo is added after conserved 5' donor sites, which is usually 'GTAAGT'
                    exon_list += [gb[start:end].seq,]
            return exon_list
    except IndexError:
        print('Error: No such exon. Change number in target_exon variable.')
        exit()

# Global variable to indicate start of cargo sequence 
cargo_start = None

# Global variable to indicate end of intron sequence 
intron_end = None

# Generating a genomic sequence with inserted sequence
def genomic_with_insertion(gb, seqfeat, manual_annotation, target_exon, cargo):
    
    #finds position at the end of exon
    end_exon = gb.seq.find(extract_exons(gb, seqfeat, manual_annotation)[target_exon-1]) \
    + len(extract_exons(gb, seqfeat, manual_annotation)[target_exon-1]) 
    
    #position of insertion point
    global cargo_start
    cargo_start = end_exon + 6
    
    #insert cargo
    genomic_seq = gb.seq[0:cargo_start] + cargo + gb.seq[cargo_start:]
    
    #finding the end position of intron
    global intron_end
    intron_end = genomic_seq.find(intron) + len(intron)
    
    return genomic_seq

# Output here is a genomic sequence to be pasted in NetGene2 to find predicted donor and acceptor sites
# Another output is the desired intron-exon junctions, which can also be compared to NetGene2 output
def intron_pos(genomic_seq):
    if genomic_seq.count(intron) > 1:
        error = 'More than one intron sequences in genome sequence. Please check the sequences and run program again.'
        return error
    print('Desired donor splice sites, direct strand:')
    print('exon ^ intron')
    print(str(genomic_seq[cargo_start - 16: cargo_start - 6]) + '^' + str(genomic_seq[cargo_start - 6: cargo_start + 4]))
    print('')
    print('Desired acceptor splice sites, direct strand:')
    print('intron ^ exon')
    print(str(genomic_seq[intron_end - 10: intron_end]) + '^' + str(genomic_seq[intron_end: intron_end + 10]))

def find_PAM(genomic_seq, intron):
    start = cargo_start - insertion_pos - target_seq_len - pam_len 
    end = cargo_start - insertion_pos - target_seq_len  
    return genomic_seq[start:end]

#######################
#Running the programme#
#######################
genomic_seq = genomic_with_insertion(gb, seqfeat, manual_annotation, target_exon, cargo)

print('~~~~~~RESULTS~~~~~~~')
print(f'The PAM is {find_PAM(genomic_seq, intron)}')
print('The known optimal PAM sequences are GT, CC, CA and CT')
print('')

intron_pos(genomic_seq)
print('')

print('The genomic sequence with insert is')
print(str(genomic_seq))

print('')
print('Please conduct the following checks:')
print('1. Desired intron sequence is correctly spliced out.')
print('2. PAM is correct.')
print('3. The frame the last triplet in exon is in. Choose the correct insert sequence. This ensures the translated CDS is in-frame.')
print('4. No BsaI or BpiI sites.')









