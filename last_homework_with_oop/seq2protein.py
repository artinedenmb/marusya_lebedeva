import json

with open('rna_codon_table.json') as f:
    #######################################################################
    # TODO:                                                               #
    # Use json.load function to load contents of json file to a dict      #
    #######################################################################
    # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    CODON_TABLE = json.load(f)
    #codon_set =  set(CODON_TABLE.values())
    #codon_set.remove('Stop')


# Вспомогательный класс для вызова SequenceError
# raise SequenceError('Error text')
class SequenceError(Exception):
    def __init__(self, message: str) -> None:
        self.message = message
    
    def __str__(self) -> str:
        return f'{self.message}'

class Sequence(object):
    seq_type = None
    types = set(['DNA', 'RNA', 'Protein']) # all possible sequence types

    #######################################################################
    # TODO:                                                               #
    # Create sets: _prot_acids (a set of all amino acids using            #  
    # json codon_table), _dna_nucls (a set of all nucleotides in a DNA),  #
    # _rna_nucls (a set of all nucleotides in a RNA)                      #
    #######################################################################
    # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    _prot_acids =  set(CODON_TABLE.values())
    _prot_acids.remove('Stop')
    _dna_nucls = set('A', 'T', 'G', 'C')
    _rna_nucls = set('A', 'U', 'G', 'C')
    # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    def __init__(self, file_name: str) -> None:
        """
        Input:
        - file_name : FASTA file with sequence_name, sequence_type and 
        the sequence itself

        Output: None
        """
        #######################################################################
        # TODO:                                                               #
        # Using _parse and _check methods check that input sequence is        #
        # correct. Add an extra check that resulting sequence type corresponds#
        # to the given sequence_type. if everything is fine then create       #
        # sequence_attribute to store given sequence and change self.seq_type # 
        # to sequence_type, else raise Error                                  #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        self.file_name = file_name  
        self.seq_type, self.seq_itself   = self._parse(file_name)
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        
    def _parse(self, file_name: str) -> tuple[str, str]:
        """
        Input:
        - file_name : FASTA file with sequence_name, sequence_type and 
        the sequence itself

        Output: A tuple with (sequence, sequence_type)
        """
        #######################################################################
        # TODO:                                                               #
        # Open file_name and read its contents.                               #
        # Input file format:                                                  #
        # >sequence_name sequence_type                                        #
        # sequence                                                            #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        with open(self.file_name, "r") as file:
            lines = file.readlines()
            snt = lines.split()[1]

            for i, line in enumerate(lines)[1:]:
                line = line.strip() 
                lines[i] = line

            sn = ''.join(lines)

        return snt, sn
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def _check(self, string: str) -> bool:
        """
        Input:
        - string : sequence from FASTA file

        Output: A boolean value (True or False)
        """
        #######################################################################
        # TODO:                                                               #
        # Check that given type is in self.types                              #
        # Check that every element of given string is either in               #
        # self._prot_acids/self._dna_nucls/self._rna_nucls. If its true,      #
        # return True, else return False                                      #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        if self.seq_type in self.types:
            #choose needed set
            if self.seq_type == 'DNA':
                self.needed_set = self._dna_nucls
            if self.seq_type == 'RNA':
                self.needed_set = self._rna_nucls
            if self.seq_type == 'Protein':
                self.needed_set = self._prot_acids
            for i in len(self.seq_itself):
                if self.seq_itself[i] not in self.needed_set:
                    return False
            return True
        else:
            return False
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def hamming_distance(self, sequence: str) -> int:
        """
        Input:
        - sequence : another sequence of nucleotides

        Output: number of different letters 
        between sequence_attribute and given string sequence
        """
        #######################################################################
        # TODO:                                                               #
        # First, check that sequence_attribute and given string have the same #
        # length, if not raise Error.                                         #
        # If the length of strings is the same, loop over one of the strings  #
        # and count different letters.                                        #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        self.another_seq = sequence
        if len(self.seq_itself) != len(self.another_seq ):
            raise NameError('HiThere')
        else:
            def counter(p, n):
                if p != n:
                    return 1
                return 0
            
            symbs = list(map(counter, self.seq_itself, self.another_seq))
            return sum(symbs)
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def count_nucleotides(self) -> None:
        """
        Input: None

        Output: None
        """
        #######################################################################
        # TODO:                                                               #
        # Raise an error, since transrcibe method can work only with elements #
        # of DNA, RNA or Protein classes                                      #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        raise NameError('HiThere')
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def to_protein(self) -> None:
        """
        Input: None

        Output: None
        """
        #######################################################################
        # TODO:                                                               #
        # Raise an error, since transrcibe method can work only with elements #
        # of DNA, RNA or Protein classes                                      #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        raise NameError('HiThere')
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def transcribe(self) -> None:
        """
        Input: None

        Output: None
        """
        #######################################################################
        # TODO:                                                               #
        # Raise an error, since transrcibe method can work only with elements #
        # of DNA or RNA classes                                               #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        raise NameError('HiThere')
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

class DNA(Sequence):
    _type = 'DNA'
    def count_nucleotides(self) -> dict:
        """
        Input: None

        Output:
        - a dictionary with keys 'A', 'T', 'G', 'C' and their 
        corresponding amounts in sequence_attribute. 
        {'A': count_A, 'T': count_T, 'G': count_G, 'C': count_C}
        """
        #######################################################################
        # TODO:                                                               #
        # Counting 'A's, 'T's, 'G's, 'C's either by                           #
        # looping over sequence_attribute or using a standard string method.  #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        A, T, G, C = (self.seq_itself.count('A'), self.seq_itself.count('T'),
            self.seq_itself.count('G'), self.seq_itself.count('C'))
        d = {'A': A, 'T': T, 'G': G, 'C': C}
        return d
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def complement_dna(self) -> str:
        """
        Input: None

        Output: a NEW string, where all 'A's were changed to 'T's
        and vice versa, all 'C's changed to 'G's and vice versa
        """
        #######################################################################
        # TODO:                                                               #
        # Create a new empty string. Loop over sequence_attribute with        #
        # if-statements, while adding corresponding letters                   #
        # to the empty string.                                                #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        def apply(x):
            if x == 'A':
                x = 'T'
            elif x == 'T':
                x = 'A'
            elif x == 'C':
                x = 'G'
            else:
                x = 'C'
            return x
        
        listochek = list(map(apply, self.seq_itself))
        return ''.join(listochek)
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def transcribe(self) -> str:
        """
        Input: None

        Output: a NEW string, where all 'T's were changed to 'U's
        """
        #######################################################################
        # TODO:                                                               #
        # Create a new empty string. Loop over sequence_attribute with        #
        # if-statements, while adding corresponding letters                   #
        # to the empty string.                                                #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        return self.seq_itself.replace('T', 'U')
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def to_protein(self) -> str:
        """
        Input: None

        Output: a NEW string, where all nucleotides are replaced with codons
        """
        #######################################################################
        # TODO:                                                               #
        # First transcribe DNA sequence to RNA using transcribe() method.     #
        # Second find 'AUG' - start codon. If it is found use json codon table#
        # to transcribe every 3 nucleotides of sequence_attribute to a codon, #
        # if a stop-codon is met then stop the transcription, else transcribe #
        # untill the end of the sequence.                                     #
        # If start codon was not found - raise Error                          #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        if self.seq_itself.find('ATG') == -1:
            return False
        return True
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def restriction_slices(self) -> int:
        """
        Input: None

        Output: number of slices by EcoRI restrictase
        """
        #######################################################################
        # TODO:                                                               #
        # EcoRI restrictase slices 'GAATTC'/'CTTAAG' sequences, return        #
        # resulting number of slices                                          #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        return len(self.seq_itself.split('GAATTC', 'CTTAAG'))
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

class RNA(Sequence):
    _type = 'RNA'
    def count_nucleotides(self) -> dict:
        """
        Input: None

        Output:
        - a dictionary with keys 'A', 'U', 'G', 'C' and their 
        corresponding amounts in sequence_attribute. 
        {'A': count_A, 'U': count_U, 'G': count_G, 'C': count_C}
        """
        #######################################################################
        # TODO:                                                               #
        # Counting 'A's, 'U's, 'G's, 'C's either by                           #
        # looping over sequence_attribute or using a standard string method.  #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        A, U, G, C = (self.seq_itself.count('A'), self.seq_itself.count('U'),
            self.seq_itself.count('G'), self.seq_itself.count('C'))
        d = {'A': A, 'U': U, 'G': G, 'C': C}
        return d
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def transcribe(self) -> str:
        """
        Input: None

        Output: a NEW string, where all 'U's were changed to 'T's
        """
        #######################################################################
        # TODO:                                                               #
        # Create a new empty string. Loop over sequence_attribute with        #
        # if-statements, while adding corresponding letters                   #
        # to the empty string.                                                #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        return self.seq_itself.replace('U', 'T')
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def to_protein(self) -> str:
        """
        Input: None

        Output: a NEW string, where all nucleotides are replaced with codons
        """
        #######################################################################
        # TODO:                                                               #
        # First find 'AUG' - start codon. If it is found use json codon table #
        # to transcribe every 3 nucleotides of sequence_attribute to a codon, #
        # if a stop-codon is met then stop the transcription, else transcribe #
        # untill the end of the sequence.                                     #
        # If start codon was not found - raise Error                          #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        if self.seq_itself.find('AUG') == -1:
            return False
        
        self.protein_seq = ''
        
        start_place = self.seq_itself.index('AUG')
        for i in range(start_place, len(self.seq_itself), 3):
            cod = self.seq_itself[i:i+3]
            for key, item in CODON_TABLE.items():
                if key in ["UAA", "UAG", "UGA"]:
                    break
                if cod == key:
                    self.protein_seq += item

        return self.protein_seq
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

class Protein(Sequence):
    _type = 'Protein'
    #######################################################################
    # TODO:                                                               #
    # Create two sets of positive charge amino acids and negative charge  #
    # amino acids: _pos_acids and _neg_acids                              #
    #######################################################################
    # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    _pos_acids = set(['K', 'R', 'H'])
    _neg_acids = set(['D', 'E'])
    # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    def count_amino_acids(self) -> dict:
        """
        Input: None

        Output:
        - a dictionary with amino acids as keys and their corresponding
        amounts in sequence_attribute
        {'A': count_A, 'U': count_U, 'G': count_G, 'C': count_C ...}
        """
        #######################################################################
        # TODO:                                                               #
        # Using dictionary with codons as values, create a new dictionary     #
        # with amino acids and count them in sequence_attribute by            #
        # looping over sequenceattribute or using a standard string method.   #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        K, R, H, D, E = (self.seq_itself.count('K'), self.seq_itself.count('R'),
            self.seq_itself.count('H'), self.seq_itself.count('D'), self.seq_itself.count('E'))
        self.d_positive = {'K': K, 'R': R, 'H': H}
        self.d_negative = {'D': D, 'E': E}
        return self.d_positive, self.d_negative
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def to_protein(self) -> str:
        """
        Input: None

        Output:
        - sequence of amino acids
        {'A': count_A, 'U': count_U, 'G': count_G, 'C': count_C ...}
        """
        #######################################################################
        # TODO:                                                               #
        # Returning sequence_attribute is enough                              #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        return self.seq_itself
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    
    def charge(self) -> int:
        """
        Input: None

        Output:
        - resulting charge of sequence
        """
        #######################################################################
        # TODO:                                                               #
        # Loop over the amino acids and check if they are in _pos_acids or    #
        # _neg_acids. If in _pos_acids then charge increases by 1, if in      #
        # _neg_acids decreases by 1, else does not changes                    #
        #######################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        cnt_charger = sum(self.d_positive.values()) - sum(self.d_negative.values())
        return cnt_charger
        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****