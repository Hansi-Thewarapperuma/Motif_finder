'''
Project 5 - Motif Finder
class and methods

Author : Hansi Thewarapperuma
Date: 10/03/2022
'''


import re

class SequenceMotifs:

    # constructor method with input parameters
    def __init__(self, sequence):
        self.sequence = sequence


    # instance method to get the sequence type
    def get_Seq_Type(self, sequence):
        list1 = ['D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'V', 'W', 'Y']
        list2 = ['B','J','U','X','Z']
        # define an empty string
        type = ''

        if "U" in sequence:
            type = "mRNA"
            return type
        for base in sequence:
            if base in list1:
                type = "Protein"
                break
            elif base in list2:
                type = "INVALID sequence"
                break
            else:
                type = "DNA"
        return type

    # instance method to find whether TATA box is present
    def find_tata_box(self):
        tata_regex = re.compile('TATA[AT]A[AT]')
        if tata_regex.search(self.sequence.upper()):
            return True, tata_regex.search(self.sequence.upper()).span()
        else:
            return False
        # mo = tata_regex.finditer(self.sequence)
        # for item in mo:
            # return item.group()
            # return item.span()

    # instance method to find whether CCAAT bos is present
    def find_ccaat_box(self):
        ccaat_regex = re.compile('GGCCAATCT')
        if ccaat_regex.search(self.sequence.upper()):
            return True , ccaat_regex.search(self.sequence.upper()).span()
        else:
            return False

    # instance method to get AT content
    def calculate_at_content(self):
        at_count = 0
        for base in self.sequence.upper():
            if base == 'A' or base == 'T' :
                at_count += 1
        at_content = at_count/ len(self.sequence)
        return at_content
        # return float(at_count) / len(self.sequence)

    # instance method to get GC content
    def calculate_gc_content(self):
        gc_count = 0
        for base in self.sequence.upper():
            if base == 'G' or base == 'C':
                gc_count += 1
        gc_content = gc_count / len(self.sequence)
        return gc_content

    # fasta split method to seperate multiple sequences contained in a fasta file and present it in a dictionary
    @staticmethod
    def fasta_split(fasta_file):
        # read the fasta and store the sequences in a dictionary
        seq_dict = {}
        with open(fasta_file) as file:
            for line in file:
                # for the header line
                if line.startswith('>'):
                    # strip to remove unwanted characters
                    seq_name = line.strip()[1:]
                    seq_dict[seq_name] = ''
                else:
                    # concatenate the rest of the lines
                    seq_dict[seq_name] += line.strip()
        return seq_dict

    # method to find element types when a FASTA file containing multiple sequences are given as the input
    @staticmethod
    def search_fasta_file(seq_dict, search_type):
        # define empty dictionary
        results = {}
        for seq_name, seq in seq_dict.items():
            sequence = SequenceMotifs(seq)
            if search_type == 'TATA box':
                result = sequence.find_tata_box()
            elif search_type == 'CCAAT box':
                result = sequence.find_ccaat_box()
            elif search_type == 'AT content':
                result = sequence.calculate_at_content()
            elif search_type == 'GC content':
                result = sequence.calculate_gc_content()
            else:
                raise ValueError('This is an INVALID search type. Valid search types: TATA box / CCAAT box / AT content / GC content')
            results[seq_name] = result
        return results

    # method to find search for a user given sequence motif pattern in a file containing multiple FASTA sequences. The user must provide the sequence pattern as a regular expression, which should be taken as a parameter for the method
    @staticmethod
    def search_fasta_file_for_pattern(seq_dict, pattern):
        results = {}
        for seq_name, seq in seq_dict.items():
            if re.search(pattern.upper(), seq):
                results[seq_name] = True
            else:
                results[seq_name] = False
        return results

    # method to predict the presence of promoter sequence for user given DNA sequence
    @staticmethod
    def predict_promoter_sequence(sequence):
        tata_box_present = SequenceMotifs(sequence).find_tata_box()
        ccaat_box_present = SequenceMotifs(sequence).find_ccaat_box()
        gc_content = SequenceMotifs(sequence).calculate_gc_content()
        # GC islands (CpG islands) found at 5' region of gene has a GC content of approx. 50% or greater
        if tata_box_present or ccaat_box_present or gc_content >= 0.5:
            return True
        else:
            return False