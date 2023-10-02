from collections import defaultdict

DNA = "ataccgatt"
# TTACGA

def validate_dna(dna_seq):
    #covert to uppercase
    seqm =  dna_seq.upper()
    print(seqm)

    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")

    if valid == len(seqm):
        print("valid")
        return True
    else:
        print("invalid")
        return False
    

validate_dna(DNA)

def frequency (seq):
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] +=1
        else:
            dic[s] = 1
    return dic

def validate_dna(dna_seq):
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")

    if valid == len(seqm):
        print("Valid")
        return True
    else:
        print("Invalid")
        return False

def frequency(seq):
    dic = {} #dictionary

    for s in seq.upper():
        if s in dic:
            dic[s] += 1
        else:
            dic[s] = 1
    return dic

def count_percentage(freq_data):
    dna_length = sum(freq_data.values())
    print("DNA Length = ", str(dna_length))

    g, c = 0, 0

    if "G" in freq_data:
        percent_g = freq_data.get("G")/dna_length * 100

    if "C" in freq_data:
        percent_c = freq_data.get("C")/dna_length * 100

    return percent_g, percent_c

print("DNA = ", DNA)

validate_result = validate_dna(DNA)
if(validate_result == True):
    freq_data = frequency(DNA)
    print(freq_data)

    #count_percentage(freq_data)
    result_g, result_c = count_percentage(freq_data)
    print("G = ", round(result_g), "%")
    print("C = ", round(result_c), "%")
print("-----------------------------")

########## hw
# from collections import defaultdict

# DNA = "ataccgatt"

print("input DNA  :", DNA.upper()) 
input_dna = DNA.upper()

replacement_dict = {'T': 'A','A':'T','C':'G','G':'C'}
#Complement
complement = ''.join(replacement_dict.get(char, char) for char in input_dna)
print("complement :", complement)
# complement = "AATGCTAAT"
#mRNA
mRNA = complement.replace('T', 'U')
print("mRNA       :", mRNA)
print(" ")

list = []

# Split the string into 3-letter parts
for i in range(0, len(mRNA), 3):
    codon = mRNA[i:i+3]
    list.append(codon)

# aminoacid sequence
for i, codon in enumerate(list, start=1):
    print(f"Part {i}: {codon}")
    if codon == "UUU" or codon == "UUC":
        print("Aminoacid = Phe(F)")
    elif codon == "UUA" or codon == "UUG":
        print("Aminoacid = Leu(L)") 
    elif codon == "UCU" or codon == "UCC" or codon == "UCA" or codon == "UCG":
        print("Aminoacid = Ser (S)")
    elif codon == "UAU" or codon == "UAC":
        print("Aminoacid = Tyr (Y)")
    elif codon == "UAA" or codon == "UAG":
        print("STOP")
    elif codon == "UGU" or codon == "UGC":
        print("Aminoacid = Cys (C)")
    elif codon == "UGA":
        print("STOP")
    elif codon == "UGG":
        print("Aminoacid = Trp (W)")
    elif codon == "CUU" or codon == "CUC" or codon == "CUA" or codon == "CUG":
        print("Aminoacid = Leu (L)")
    elif codon == "CCU" or codon == "CCC" or codon == "CCA" or codon == "CCG":
        print("Aminoacid = Pro (P)")
    elif codon == "CGU" or codon == "CGC" or codon == "CGA" or codon == "CGG":
        print("Aminoacid = Arg (R)")
    elif codon == "AUU" or codon == "AUC" or codon == "AUA":
        print("Aminoacid = Ile (I)")
    elif codon == "AUG":
        print("Aminoacid = Met (M)")
    elif codon == "ACU" or codon == "ACC" or codon == "ACA" or codon == "ACG":
        print("Aminoacid = Thr (T)")
    elif codon == "AAU" or codon == "AAC":
        print("Aminoacid = Asn (N)")
    elif codon == "AAA" or codon == "AAG":
        print("Aminoacid = Lys (K)")
    elif codon == "AGU" or codon == "AGC":
        print("Aminoacid = Ser (S)")
    elif codon == "AGA" or codon == "AGG":
        print("Aminoacid = Arg (R)")
    elif codon == "GUU" or codon == "GUC" or codon == "GUA" or codon == "GUG":
        print("Aminoacid = Val (V)")
    elif codon == "GCU" or codon == "GCC" or codon == "GCA" or codon == "GCG":
        print("Aminoacid = Ala (A)")
    elif codon == "GAU" or codon == "GAC":
        print("Aminoacid = Asp (D)")
    elif codon == "GAA" or codon == "GAG":
        print("Aminoacid = Glu (E)")
    elif codon == "GGU" or codon == "GGC" or codon == "GGA" or codon == "GGG":
        print("Aminoacid = Gly (G)")
print(" ")

# Initialize a dictionary to store codon frequencies
codon_freq = defaultdict(int)

# Count the frequency of each codon
for codon in list:
    codon_freq[codon] += 1

# Print the codon frequencies
for codon, freq in codon_freq.items():
    print(f'{codon}: {freq}')
