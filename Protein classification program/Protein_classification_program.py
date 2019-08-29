"""
This is a program which will translates a gene sequence into a polypeptide. 
Then a codon bais will be predicted, the charge and structure of the protein will 
also be predicted.

By Thulani Tshabalala
   tshabalalaboyt@gmail.com
"""

File_Location= input("Input FASTA file location:  ") # will take in the file location from user

Codon_info = [ [ [ ['TTT',0], ['TTC',0] ], 'Phe'], [ [ ['TTA',0] ,['TTG',0], ['CTT',0], ['CTC',0], 
['CTA',0], ['CTG',0] ] ,'Lue'], [ [ ['ATT',0], ['ATC',0], ['ATA',0] ],  'Ile' ], [ [ ['ATG',0] ], 'Met' ], 
[ [ ['GTT',0], ['GTC',0], ['GTA',0], ['GTG',0] ], 'Val' ], [ [ ['TCT',0], ['TCC',0], ['TCA',0], ['TCG',0], 
['AGT',0], ['AGC',0] ], 'Ser' ], [ [ ['CTT',0], ['CCC',0], ['CCA',0], ['CCG',0] ], 'Pro' ], [ [ ['ACT',0], 
['ACC',0], ['ACA',0], ['ACG',0] ], 'Thr' ], [ [ ['GCT',0], ['GCC',0], ['GCA',0], ['GCG',0] ], 'Ala' ], 
[ [ ['TAT',0], ['TAC',0] ], 'Tyr' ], [ [ ['CAT',0], ['CAC',0] ], 'His' ], [ [ ['CAA',0], ['CAG',0] ], 'Gln' ],
[ [ ['AAT',0], ['ACC',0] ], 'Asn' ], [ [ ['AAA',0], ['AAG',0] ], 'Lys' ], [ [ ['GAT',0], ['GAC',0] ], 'Asp' ], 
[ [ ['GAA',0], ['GAG',0] ], 'Glu' ], [ [ ['TGT',0], ['TGC',0] ], 'Cys' ], [ [ ['TGG',0] ], 'Trp' ], [ [ ['CGT',0], 
['CGC',0], ['CGA',0], ['CGG',0], ['AGA',0], ['AGG',0] ], 'Arg' ], [ [ ['GGT',0], ['GGC',0], ['GGA',0], ['GGG',0] ], 'Gly' ],
[ [ ['TAA',0], ['TCA',0], ['TGA',0] ], 'STOP' ]  ]

AminoAcid_Symbols = [ ['Alanine','Ala','A'], ['Arginine','Arg','R'], ['Asparagine', 'Asn', 'N'], 
['Aspartic acid', 'Asp', 'D'], ['Cysteine', 'Cys','C'], ['Glutamic acid', 'Glu', 'E'], ['Glutamine', 'Gln', 'Q'], 
['Glycine', 'Gly','G'], ['Histidine','His',	'H'], ['Isoleucine', 'Ile',	'I'], ['Leucine', 'Leu', 'L'], 
['Lysine', 'Lys','K'], ['Methionine', 'Met', 'M'], ['Phenylalanine', 'Phe',	'F'], ['Proline', 'Pro', 'P'],
['Serine', 'Ser', 'S'], ['Threonine', 'Thr',	'T'], ['Tryptophan', 'Trp',	'W'], ['Valine', 'Val',	'V'] ]

def translation(gene, codon_info): # translations the gene sequence to an polypeptide
    stop = 0                                
    nucleotide = gene[0]                            
    codon_info_Return= codon_info
    polypeptide_Return = ''
    sequence = gene
    amino_Acid =''
    while('L' != nucleotide ):
        if('G'== nucleotide or 'T'== nucleotide or 'A'== nucleotide or 'C'== nucleotide):
            amino_Acid = amino_Acid + nucleotide   # Addes nucleotides together to form an amina acid
            stop= stop+1                           # stop to see if codon is formed
        else:
            print("Nucleotide sequence corrupted! \n >>> Program shutdown <<<")
            return ()
        if(stop==3):
            codon_Num = 0
            stop=0
            while(codon_Num<len(codon_info_Return)):
                boolStopper = 0
                list_aminoAcid = codon_info_Return[codon_Num] #list_ aminoAcid gives [List of codons, Amino Acid]
                Codon_Sequence = list_aminoAcid[0]     # codon_Sequence gives List of [codon, number]
                CodonType_Num=0
                while(CodonType_Num<len(Codon_Sequence)):
                    nucleotides_ForCodon= Codon_Sequence[CodonType_Num]  # nucleotides_ForCodon gives codon, number
                    if(amino_Acid==nucleotides_ForCodon[0]): # checks whether codons match with Amino Acid codon
                        nucleotides_ForCodon[1]= nucleotides_ForCodon[1]+1  # Increment the codon usage
                        polypeptide_Return = polypeptide_Return+ list_aminoAcid[1] # Adds amino acid to the polypeptide
                        boolStopper = 1
                        break
                    CodonType_Num= CodonType_Num+1
                if(boolStopper==1):
                    amino_Acid = ''
                    break
                codon_Num = codon_Num+1
        sequence = sequence[1:]             #Removes leading nucleotide in the sequence
        nucleotide = sequence[0]            
    return([polypeptide_Return,codon_info_Return])

def Codon_Bias (codon_Info, polypeptid): # Predict possible codon bais for an amino
    print('Bais is considered when codon is used more then 50 times')
    moveTho_codon_Info=0
    BiasOutput =''
    codon_Used = ''
    sequence = polypeptid
    aminoAcid_Codon = sequence[0:3]
    sequence = sequence[3:]
    while(moveTho_codon_Info<len(codon_Info)):
        aminoAcid_info = codon_Info[moveTho_codon_Info]
        AminoAcid = aminoAcid_info[1]
        if(aminoAcid_Codon==AminoAcid):
            if (sequence != ''):
                aminoAcid_Codon = sequence[0:3]
                sequence = sequence[3:]
            codons = aminoAcid_info[0]
            CodonType_Num = 0
            codonUsage = 0
            codonB = 0
            while(CodonType_Num<len(codons)):
                codonProperties= codons[CodonType_Num]  
                CodonUsage = codonProperties[1]
                if(CodonUsage>codonB):
                    codon_Used =''
                    codonB = CodonUsage
                    codon_Used = codonProperties[0] + ' used ' + str(CodonUsage)
                elif (CodonUsage>=codonB and codonB>50):
                    codon_Used = codon_Used + ', ' + codonProperties[0]
                CodonType_Num = CodonType_Num + 1
            if(codonB==0):
                AminoAcid = AminoAcid + ' -> No codon bais'
                BiasOutput = BiasOutput + AminoAcid + ' \n'
            else:
                BiasOutput = BiasOutput + AminoAcid + ' -> ' + codon_Used + ' \n'
            moveTho_codon_Info =  moveTho_codon_Info + 1 
            
        else:
            moveTho_codon_Info =  moveTho_codon_Info + 1 
    return(BiasOutput)

def Polypeptid_Charge (polypeptid):
    sequence = polypeptid
    aminoAcid = sequence[0:3]
    sequence = sequence[3:]
    Positive_Arg = 0 
    Positive_Lys = 0 
    Positive_His = 0 
    Negitive_Asn = 0
    Negitive_Gln = 0
    polypeptid_ChargePosition = ''
    OutputString =''
    while(aminoAcid != ''):
        if(aminoAcid=='Arg'):
            Positive_Arg = Positive_Arg + 1
            polypeptid_ChargePosition = polypeptid_ChargePosition + ' + '
        elif(aminoAcid == 'Lys'):
            Positive_Lys = Positive_Lys + 1
            polypeptid_ChargePosition = polypeptid_ChargePosition + ' + '
        elif(aminoAcid == 'His'):
            Positive_His = Positive_His + 1
            polypeptid_ChargePosition = polypeptid_ChargePosition + ' + '
        elif(aminoAcid=='Asn'):
            Negitive_Asn = Negitive_Asn + 1
            polypeptid_ChargePosition = polypeptid_ChargePosition + ' - '
        elif(aminoAcid =='Gln'):
            Negitive_Gln = Negitive_Gln +1
            polypeptid_ChargePosition = polypeptid_ChargePosition + ' - '
        else :
            count = 0
            while (count< len(AminoAcid_Symbols)):
                AminoAcid_info = AminoAcid_Symbols[count]
                if(AminoAcid_info[1]==aminoAcid):
                    AminoSym = AminoAcid_info[2]
                    polypeptid_ChargePosition = polypeptid_ChargePosition + AminoSym
                    break
                count = count + 1
        if(sequence ==""):
            break
        aminoAcid = sequence[0:3]
        sequence = sequence[3:]
    OutputString = OutputString + 'This polypeptide contains: Negative Amino Acids -> Asparagine ' + str(Negitive_Asn) +'\n' 
    OutputString = OutputString + '                                                -> Glutamine  ' + str(Negitive_Gln) + '\n'
    OutputString = OutputString + '                         : Positive Amino Acids -> Aspartic acid ' + str(Positive_Arg) + '\n'
    OutputString = OutputString + '                                                -> Lysine ' + str(Positive_Lys) + '\n'
    OutputString = OutputString +'                                                -> Histidine ' +  str(Positive_His) + '\n'
    OutputString = OutputString + 'Polypeptide and charge location \n'+ polypeptid_ChargePosition
    return (OutputString)

def predicted_structure (polypeptid):



try:
    with open(File_Location,'r') as  Gene_Sequence: #open the FASTA file and stores the gene sequence into GENE_Sequence
        gene_Information = Gene_Sequence.readlines() # Stores FASTA File heading
        nucleotide_List = Gene_Sequence.readlines() # stores file as list of nucleotides 
        nucleotide_sequence =''
        for index in nucleotide_List: # retrieve every nucleotide squence per line into single string of nucleotides
            sequence_Length = len(index)
            nucleotide_sequence = nucleotide_sequence + index[:sequence_Length-2]
        nucleotide_sequenceStop = nucleotide_sequence +'L'
finally:
    print("Error. The file is not found ! \n >>> Program shutdown <<<") # error for file not found








