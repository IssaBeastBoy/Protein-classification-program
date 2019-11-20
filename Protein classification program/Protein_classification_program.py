"""
This is a program which will translates a gene sequence into a polypeptide. 
Then a codon bais will be predicted, the charge and structure of the protein will 
also be predicted.

By Thulani Tshabalala
   tshabalalaboyt@gmail.com
"""
print ("This program functions for only prokaryotes and eukaryotes")
File_Location= input("Input FASTA file locationc containing the gene sequence:  ") # will take in the file location from user
Protein_source= input("For prokaryotes enter P and for eukaryotes enter E: ")

#Protein information
Alanine_Info= {'Name': 'Alanine', 'symbols':['Ala','A'], 'codon':['GCU', 'GCC', 'GCA', 'GCG']}
Arginine_Info= {'Name': 'Arginine', 'symbols':['Arg','R'], 'codon':['CGU', 'CGC', 'CGA', 'CGG']}
AsparticAcid_Info={'Name': 'Aspartic acid' , 'symbols':['Asp', 'D', 'codon':['GAU', 'GAC']}
Cysteine_Info={'Name': 'Cysteine' , 'symbols':['Cys','C'], 'codon':['UGU', 'UGC']}
GlutamicAcid_Info={'Name': 'Glutamic acid', 'symbols':['Glu', 'E'], 'codon':['GAA', 'GAG']}
Glutamine_Info={'Name': 'Glutamine', 'symbols':['Gln', 'Q'], 'codon':['CAA', 'CAG']}
Glycine_Info={'Name': 'Glycine', 'symbols':['Gly','G'], 'codon':['GGU', 'GGC', 'GGA', 'GGG']}
Histidine_Info={'Name': 'Histidine', 'symbols':['His', 'H'], 'codon':['CAU', 'CAC']}
Isoleucine_Info={'Name': 'Isoleucine', 'symbols':['Ile', 'I'], 'codon':['AUU', 'AUC', 'AUA']}
Leucine_Info={'Name': 'Leucine', 'symbols':['Leu', 'L'], 'codon':['CUU', 'CUC', 'CUA', 'CUG']}
Lysine_Info={'Name': 'Lysine', 'symbols':['Lys','K'], 'codon':['AAA', 'AAG']}
Methionine_Info={'Name': 'Methionine', 'symbols':['Met', 'M'], 'codon':['AUG']}
Phenylalanine_Info={'Name': 'Phenylalanine', 'symbols':['Phe', 'F'], 'codon':['UUU', 'UUC']}
Proline_Info={'Name': 'Proline', 'symbols':['Pro', 'P'], 'codon':['CCU', 'CCC', 'CCA', 'CCG']}
Serine_Info={'Name': 'Serine', 'symbols':['Ser', 'S'], 'codon':['UCU', 'UCC', 'UCA', 'UCG']}
Threonine_Info={'Name': 'Threonine', 'symbols':['Thr', 'T'], 'codon':['ACU', 'ACC', 'ACA', 'ACG']}
Tryptophan_Info={'Name': 'Tryptophan', 'symbols':['Trp', 'W'], 'codon':['UGG']}
Valine_Info={'Name': 'Valine', 'symbols':['Val', 'V'], 'codon':['GUU', 'GUC', 'GUA', 'GUG']}

Proteins_List =

def translation(gene, codon_info): # translations the gene sequence to an polypeptide
    stop = 0                                
    nucleotide = gene[0]                            
    codon_info_Return= codon_info
    polypeptide_Return = ''
    sequence = gene
    amino_Acid =''
    Liststru2 =''
    Liststru1 =''
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
            if(amino_Acid == 'TAA' or amino_Acid =='TCA' or amino_Acid == 'TGA'):
                print('Polypeptide sequence -> ' + polypeptide_Return + '\n')
                print ('Codon Bias for polypeptid -> \n' + Codon_Bias(codon_info_Return, polypeptide_Return))
                print('Polypeptide charged positions -> \n' + Polypeptid_Charge (polypeptide_Return)  )
                List_Ofstructures = predicted_structure(polypeptide_Return)
                Liststru1 = List_Ofstructures[0]
                Liststru2 = List_Ofstructures[1]
                print('\nPolypeptide predicted structure -> \n' + Liststru1 + '\n' + Liststru2)
                sequence = sequence[1:]
                return (translation(sequence, Codon_info))
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

def Codon_Bias (codon_Info, polypeptid): # Predict possible codon bais for an amino
    moveTho_codon_Info=0
    BiasOutput =''
    codon_Used = ''
    sequence = polypeptid
    aminoAcid_Codon = sequence[0:3]
    sequence = sequence[3:]
    end_Ofsequence = 0
    while(moveTho_codon_Info<len(codon_Info)):
        aminoAcid_info = codon_Info[moveTho_codon_Info]
        AminoAcid = aminoAcid_info[1]
        if(aminoAcid_Codon==AminoAcid):
            if (sequence != ''):
                aminoAcid_Codon = sequence[0:3]
                sequence = sequence[3:]
            else:
                end_Ofsequence = 1
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
            moveTho_codon_Info =  moveTho_codon_Info + 1 
            if(codonB==0):
                AminoAcid = AminoAcid + ' -> No codon bais'
                BiasOutput = BiasOutput + AminoAcid + ' \n'
                moveTho_codon_Info =0
                if(end_Ofsequence == 1):
                    return(BiasOutput)
            else:
                BiasOutput = BiasOutput + AminoAcid + ' -> ' + codon_Used + ' \n'
                moveTho_codon_Info = 0
                if(end_Ofsequence == 1):
                    return(BiasOutput)        
            
        else:
            moveTho_codon_Info =  moveTho_codon_Info + 1 

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
    OutputString = OutputString + '                                                -> Histidine ' +  str(Positive_His) + '\n'
    OutputString = OutputString + 'Polypeptide and charge location \n'+ polypeptid_ChargePosition
    return (OutputString)

def predicted_structure (polypeptid):
    print('\nSecondary predicted protein structures for polypeptid, base on their polarity.')
    sequence = polypeptid
    aminoAcid = ' '
    startOf_sequenceBeta = 1
    startOf_sequenceAlpha = 1
    close_Alpha = 1
    close_Beta = 1
    alphaSeq_increBool = 0
    beatSeq_increBool = 0
    endOf_sequence = 0
    Alpha_Helix = 'Predicted alpha helix structures: \n'
    Beta_Sheets = 'Predicted beta sheet strcutures: \n'
    while(aminoAcid != ''):
        aminoAcid = sequence[0:3]
        sequence = sequence[3:]
        if(aminoAcid == 'Ala' or aminoAcid =='Ile' or aminoAcid =='Leu' or aminoAcid =='Met' or aminoAcid == 'Phe' or aminoAcid =='Val' or aminoAcid =='Pro' or aminoAcid =='Gly' or aminoAcid =='Trp' or aminoAcid =='Tyr' or aminoAcid =='Met'):
            if(close_Beta == 0):
                Beta_Sheets = Beta_Sheets +'[ End of sequence: ' + str(endOf_sequence) + ' ] \n'
                close_Beta = 1
                startOf_sequenceBeta = endOf_sequence + 1
                beatSeq_increBool = 0
            if(alphaSeq_increBool == 0):
                Alpha_Helix = Alpha_Helix + '[ Start of sequence: ' + str(startOf_sequenceAlpha) + '] -> '
                alphaSeq_increBool = 1
                close_Alpha = 0
            Alpha_Helix = Alpha_Helix + aminoAcid
        if(aminoAcid =='Gln' or aminoAcid =='Asn' or aminoAcid =='His' or aminoAcid =='Ser' or aminoAcid =='Thr' or aminoAcid =='Tyr' or aminoAcid =='Cys'):
            if(close_Alpha == 0):
                Alpha_Helix = Alpha_Helix +'[ End of sequence: ' + str(endOf_sequence) + ' ] \n'
                close_Alpha = 1
                alphaSeq_increBool = 0
                startOf_sequenceAlpha = endOf_sequence + 1
            if(beatSeq_increBool == 0):
                Beta_Sheets = Beta_Sheets + '[ Start of sequence: ' + str(startOf_sequenceBeta) + '] -> ' 
                beatSeq_increBool = 1
                close_Beta = 0
            Beta_Sheets = Beta_Sheets + aminoAcid
        if ( aminoAcid==''):
            if(close_Alpha == 0):
                Alpha_Helix = Alpha_Helix + '[ End of sequence: ' + str(endOf_sequence) + '] \n'
            if(close_Beta == 0):
                Beta_Sheets = Beta_Sheets + '[ End of sequence: ' + str(endOf_sequence) + '] \n'
            return ([Beta_Sheets, Alpha_Helix])
        endOf_sequence = endOf_sequence + 1 

def transcription (protein_File): # transcriptes protein sequence to RNA sequence and save to 
    sequence = protein_File
    transcript_sequence = ''
    while (sequence !=  ''):
        nucleotide = sequence[:1]
        if(nucleotide == 'A' or nucleotide == 'a'):
            transcript_sequence = transcript_sequence + 'U'
        elif(nucleotide == 'T' or nucleotide == 't'):
             transcript_sequence = transcript_sequence + 'A'
        elif(nucleotide == 'G' or nucleotide == 'g'):
            transcript_sequence = transcript_sequence + 'C'
        elif(nucleotide == 'C' or nucleotide == 'c'):
            transcript_sequence = transcript_sequence + 'G'
        sequence = sequence[1:]
    transcript_File = open("transcripted RNA seqeunce.txt", "w+")
    return transcript_sequence
    
try:
    with open(File_Location,'r') as  Gene_Sequence: #open the FASTA file and stores the gene sequence into GENE_Sequence
        gene_Information = Gene_Sequence.readline() # Stores FASTA File heading
        nucleotide_sequence =''
        outFile = transcription(Gene_Sequence.read().replace('\n', ''))
        for line in  Gene_Sequence: # retrieve every nucleotide squence per line into single string of nucleotides
            line = line [:(len(line)-1)]
            nucleotide_sequence = nucleotide_sequence + line
        nucleotide_sequence = nucleotide_sequence +'L'
        translated = translation(nucleotide_sequence, Codon_info)
        
finally:
    print("Error. The file is not found ! \n >>> Program shutdown <<<") # error for file not found