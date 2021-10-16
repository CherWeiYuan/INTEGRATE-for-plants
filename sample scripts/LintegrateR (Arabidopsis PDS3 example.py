from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from sys import exit

#########################
#User-defined parameters#
#########################

# Parameters
insertion_pos = 49 
target_seq_len = 32
pam_len = 2
pams = ['GT','CC','CA','CT']

# Sequence of L-Cargo-R in three frames
cargo_frame1 = 'TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAGgtaatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcagctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactcacggcatggacgagctgtacaagtaaggatccccgGAAATGACCTCTTCATGAAAATGTCGTAAACTTACTGAAGTTTAGACCATGAAGAGGCAATATCAATTTATGGGTGTGATAATTATCAATTTATGGGTGTAATTATCATTTTATGGTTGTATCAACA'
cargo_frame2 = 'TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAGgtaAGatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcagctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactcacggcatggacgagctgtacaagtaaggatccccgGAAATGACCTCTTCATGAAAATGTCGTAAACTTACTGAAGTTTAGACCATGAAGAGGCAATATCAATTTATGGGTGTGATAATTATCAATTTATGGGTGTAATTATCATTTTATGGTTGTATCAACA'
cargo_frame3 = 'TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAGgtaAatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcagctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactcacggcatggacgagctgtacaagtaaggatccccgGAAATGACCTCTTCATGAAAATGTCGTAAACTTACTGAAGTTTAGACCATGAAGAGGCAATATCAATTTATGGGTGTGATAATTATCAATTTATGGGTGTAATTATCATTTTATGGTTGTATCAACA'

# Expected protein tag sequence
tag = Seq('MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK')

# Sequence of intron seq without the first 6 conserved bp, which is usually 'GTAAGT'. Should end in AG.
intron = 'TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAACTTAACCACAAAACAACCATATATTGATATCTCACAAAACAACCATAAGTTGATATTTTTGTGAATCGAGTATTTCAGCAAAACTACTGCAGTAAGgactctagaCTAATttatttcttgcAG'

# Uploading gene sequence
# For Genbank file: gb = SeqIO.read("Nicotiana benthamiana PDS (Niben101Scf01283g02002.1).gb", "gb")
# For Fasta file: gb = SeqIO.read("Nicotiana benthamiana PDS (Niben101Scf01283g02002.1).fasta", "fasta")
gb = SeqIO.read("Genomic PDS3 Arabidopsis thaliana chr4 8187190 to 8198320.gb", "gb")

# Manually adding exon locations
# manual_annotation = True if genbank files does not have exons labelled)
manual_annotation = True 

#Example: Arabidopsis PDS3
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
HELP

How to find positions? Go to Genbank browser for TAIR10 genome, 
hover cursor over the desired gene, click 'Genbank View'. 

Scroll down and you see:
    
CDS             join(552..719,808..940,

Then,  write:
seqfeat = [SeqFeature(FeatureLocation(551, 719), type="exon"),
           SeqFeature(FeatureLocation(807, 940), type="exon")]

IMPT: Note that instead of 552, we write 551.
'''

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

# Generating a genomic sequence to be pasted in NetGene2 to find predicted donor and acceptor sites
# and also the desired intron-exon junctions, which can also be compared to NetGene2 output
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

# Outputs a PAM sequence
def find_PAM(genomic_seq, intron):
    start = cargo_start - insertion_pos - target_seq_len - pam_len 
    end = cargo_start - insertion_pos - target_seq_len  
    return genomic_seq[start:end]

# Finds the 32-nt targeting sequence
def find_crRNA(genomic_seq):
    intron_pos_start = genomic_seq.find(intron) - 6
    crRNA_start = intron_pos_start - insertion_pos - target_seq_len
    crRNA_end = intron_pos_start - insertion_pos
    return genomic_seq[crRNA_start:crRNA_end]

# Finds CDS
def CDS(gb):
    x = ''
    for i in extract_exons(gb, seqfeat, manual_annotation):
        x += i
    return x

# Command line tool for users
def retrieve(exon_number, *arg):
    try:
        database = data[exon_number]
        for i in arg:
            if i == 'cargo': print(database['cargo'])
            if i == 'crrna': print(database ['crrna'])
            if i == 'pam': print(database['pam'])
            if i == 'seq': print(database['seq'])
    except KeyError:
        print("'Exon number ' + str(exon_number) is not a suitable target.")


#########################
#Executing the programme#
#########################

# Global variable to indicate start of cargo sequence 
cargo_start = None
# Global variable to indicate end of intron sequence 
intron_end = None

# Reading sequence of L-Cargo-R
cargo1 = Seq(cargo_frame1,IUPAC.unambiguous_dna)
cargo2 = Seq(cargo_frame2,IUPAC.unambiguous_dna)
cargo3 = Seq(cargo_frame3,IUPAC.unambiguous_dna)
cargo_lst = [(cargo1,1), (cargo2,2), (cargo3,3)]
intron = Seq(intron, IUPAC.unambiguous_dna)
data = {}

# Loop for exons, and nested loop for cargoes
for k in range(len(extract_exons(gb, seqfeat, manual_annotation))):
    target_exon = k + 1
    for cargo_type in cargo_lst:
        genomic_seq = genomic_with_insertion(gb, seqfeat, manual_annotation, 
                                             target_exon, cargo_type[0])
            
        #screen by PAM
        if find_PAM(genomic_seq, intron)in pams:
            pass
        else:
            continue
        
        #checking if tag, from end of intron to end of tag, is in frame
        
        #first, get exon seq
        exon_seq = extract_exons(gb, seqfeat, manual_annotation)[target_exon-1]

        #second, get seq after artificial intron
        intron_pos_end = genomic_seq.find(intron) + len(intron)
        length = len(genomic_seq[intron_pos_end: intron_pos_end + len(cargo_type[0])])
        while length%3 != 0:
            length += 1
        seq_after_artificial_intron = genomic_seq[intron_pos_end: intron_pos_end + length]
        
        #third, combine the two seqs
        exon_plus_tag = exon_seq + seq_after_artificial_intron
        
        #fourth, check if tag in frame. If yes, print results.
        if tag in exon_plus_tag.translate():
            print('--------------------Exon number: ' + str(target_exon) + '---------------------')
            
            #print('Cargo: Frame ' + str(cargo_type[1]))
            #print('')
            
            print(f'The 32-bp crRNA is {find_crRNA(genomic_seq)}')
            print('')
            
            print(f'The PAM is {find_PAM(genomic_seq, intron)}') 
            print('')
            intron_pos(genomic_seq)
            print('')
            
            #print('The genomic sequence with insert is')
            #print(str(genomic_seq))
            
            data[target_exon] = {'crrna': find_crRNA(genomic_seq),'pam': find_PAM(genomic_seq, intron), 'seq': genomic_seq}
        

        
print('')
print('***Please conduct the following checks:')
print('1. If desired intron sequence is predicted to be correctly spliced out by NetGene2.')
print('2. Cargo frame to be used.')
print('3. No BsaI or BpiI sites in your cargo.')
print('')
print('You may use retrieve(exon_number, *arg) to retrieve data,')
print("where arg is either 'crrna', 'pam','seq' or a combination.")
print("For example, use retrieve(10, 'seq') to output the genomic sequence with insertion,")
print("which is handy for checking the frame of cargo to use in Snapgene")
