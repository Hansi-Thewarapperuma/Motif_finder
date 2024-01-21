'''
Project 5 - Motif Finder
Input : FASTA file containing multiple sequences/ DNA sequence
Output : Sequence motifs

Author : Hansi Thewarapperuma
Date: 10/03/2022
'''


import SequenceMotifs_methods

# main method
if __name__ == '__main__':
    # obtain the user input fasta file
    fasta_file = input('Enter the FASTA file: ')

    # executing method 1 to split the fasta and insert the seq as values into the dictx
    dictx = SequenceMotifs_methods.SequenceMotifs.fasta_split(fasta_file)
    # print(dictx)

    # empty dict to extract only the DNA sequences from the above dict having all kind of seq types
    DNAseq_diction = {}
    print('\nSequence types in the FASTA file: \n')
    for key,value in dictx.items():
        seqType = SequenceMotifs_methods.SequenceMotifs.get_Seq_Type(SequenceMotifs_methods.SequenceMotifs, value)
        print(key,'-------------------------',seqType)

        # insert the seq into DNA dict that get the seq type as DNA
        if seqType == 'DNA':
            DNAseq_diction[key] = value
    print('\nDNA sequences extracted from the FASTA file: \n',DNAseq_diction,'\n')

    # executing method 2 to check for the availability of elements
    print('************ Availability of elements in the sequences ************')
    search_type = input("\nEnter the element type (AT content/ GC content/ TATA box/ CCAAT box) : ")
    for value in DNAseq_diction.values():
        # creating an object
        obj = SequenceMotifs_methods.SequenceMotifs(value)
        output1 = obj.search_fasta_file(DNAseq_diction, search_type)
    print('\nPresence of the ', search_type, ' : \n', output1)

    # executing method to find the availability of the user's input pattern in sequences
    print('\n************ Availability of entered pattern in the sequences ************\n')
    pattern = input("\nEnter the desired pattern: ")
    output2 = SequenceMotifs_methods.SequenceMotifs.search_fasta_file_for_pattern(DNAseq_diction,pattern)
    print('\nPresence of the entered pattern: \n',output2,'\n')

    # executing method to predict the availability of promoter region
    print('\n************ Availability of promoter region in the sequences ************\n')
    sequence = input("\nEnter the DNA sequence: ")
    out = SequenceMotifs_methods.SequenceMotifs.get_Seq_Type(SequenceMotifs_methods.SequenceMotifs,sequence)
    output3 = SequenceMotifs_methods.SequenceMotifs.predict_promoter_sequence(sequence)
    print(' \nThe sequence type of the entered sequence: \n',out)
    print('\nPresence of promoter sequence: \n',output3)


