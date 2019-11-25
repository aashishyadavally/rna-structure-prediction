"""
This script contains the :class: ``PredictStableStructure`` which computes
the maximum number of permissible letter pairs from the input sequence
of 4 letters (A, C, G, U), of which A-U, U-A, G-C, C-G are legal pairs.
Furthermore, it denotes the nucleotide bases which are predicted to be
paired with '{' and '}', and those predicted to not be paired with any
other by '.'. 

Author:
-------
Aashish Yadavally
"""

LEGAL_PAIR_MATRIX = [
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
]

PAIR_MATRIX_MAPPING = {
    ('A','U'): 0,
    ('U','A'): 1,
    ('G','C'): 2,
    ('C','G'): 3,
    ('G','U'): 4,
    ('U','G'): 5,
}


def read_structure(filename='sequence.txt'):
    """This function reads the sequence from 'sequence.txt' file
    for which to compute the maximum count of base pairs.

    Expected, default filename is 'sequence.txt'.

    Returns:
        sequence (str):
            RNA Sequence
    """
    with open(filename, 'r') as sequence_file:
        sequence = sequence_file.readline()
    return sequence


class PredictStableStructure:
    """Predicts the nucleotide bases in the input sequence which might
    form a pair with the legal nucleotide bases, representing the pairs
    with '{', '}'; further writing them to an output file ``output.txt``.

    Arguments:
        n_bases (list):
            List of nucleotide bases in the sequence.
        legal_pairs (dict):
            Dictionary of legal pairs with their corresponding
            gamma-scores.
        table (list):
            Matrix or table with the filled-in values according to
            the objective function.
        output (list):
            List of '{', '}' and '.', wherein, it denotes whether
            the nucleotide base is predicted to be paired or not.
    """
    def __init__(self, sequence, legal_pair_matrix, pair_matrix_mapping):
        """Initializes :class: ``PredictStructure``.

        Arguments:
            sequence (str):
                RNA Sequence
            legal_pair_matrix (list):
                Matrix mapping nucleotide pairs to nucleotide pairs.
            pair_matrix_mapping (dict):
                Mapping nucleotide-base pairs to their indices in
                legal_pair_matrix.
        """
        n_bases = list(sequence)
        n = len(n_bases)
        self.n_base_pairs = [tuple([n_bases[i], n_bases[i+1]]) if i != n-1 else tuple([n_bases[i], '']) for i in range(0, n)]
        self.legal_pair_matrix = legal_pair_matrix
        self.pair_matrix_mapping = pair_matrix_mapping
        self.base_pairs_braces = []
        initial_table = self.initialize_table()
        self.table = self.fill_table(initial_table)
        self.trace_non_crossings(0, len(self.n_base_pairs) - 1)
        self.output = self.decode_output()
        print('Printing DP-Table:')
        for row in self.table:
            print(row)
        print('Printing sequence:')
        print(sequence)
        print('Printing output:')
        print(self.output)


    def initialize_table(self):
        """Initializes the DP-Table, annotating all the entries with
        the value '-1', and further, assigning the value '0' to the
        entries part of the base case.

        Returns:
            table (list):
                Matrix or table of dimensions n*n, where 'n' is the
                length of the RNA sequence.
        """
        # Initialilzing DP-Table with '-1'
        table = [[-1 for i in range(len(self.n_base_pairs))] \
                    for j in range(len(self.n_base_pairs))]
        # Setting up base cases
        for i in range(len(self.n_base_pairs)):
            # Filling main diagonal with '0'
            table[i][i] = 0
            # Filling diagonal below main-diagonal with '0'
            if i > 0:
                table[i][i-1] = 0
        return table


    def is_legal(self, pair):
        """
        Arguments:
            pair (tuple):
                Nucleotide base pair.

        Returns:
            (bool):
                True or False, depending on whether (a, b) is a legal pair.
        """
        if pair in list(self.pair_matrix_mapping.keys()):
            return True
        else:
            return False


    def fill_table(self, table):
        """Fills the DP-table with values according to the following
        objective function:

        table[i][j] = max {
                            gamma-score + table[i+1][j-1]; if legal pair

                            table[i+1][j]; if base at 'i' doesn't have pair

                            table[i][j-1]; if base at 'j' doesn't have pair

                            max (i<k<j){
                                table[i][k] + table[k+1][j]                             
                                }   
                          }

        Arguments:
            table (list):
                Matrix or table with initialized values.

        Returns:
            table(list):
                Matrix of table with filled-in values according to the
                objective function.
        """
        n = len(self.n_base_pairs)
        diagonal_elements = [[i, i] for i in range(n)]

        for distance_from_diagonal in range(1, n):
            for i in range(n - distance_from_diagonal):
                row, col = diagonal_elements[i][0], diagonal_elements[i][1] + distance_from_diagonal
                if col - row > 0:
                    pairs = tuple([self.n_base_pairs[row], self.n_base_pairs[col]])
                    conditions = []

                    if self.is_legal(pairs[0]) and self.is_legal(pairs[1]):
                        pair_indices = [self.pair_matrix_mapping[pairs[0]], self.pair_matrix_mapping[pairs[1]]]
                        gamma_score = self.legal_pair_matrix[pair_indices[0]][pair_indices[1]]
                        conditions.append(gamma_score + table[row+1][col-1])
                    else:
                        conditions.append(table[row+1][col])
                        conditions.append(table[row][col-1])
                        conditions1 = []
                        for k in range(row, col):
                            conditions1.append(table[row][k] + table[k+1][col])
                        conditions.append(max(conditions1)) 
                    table[row][col] = max(conditions)
                else:
                    table[row][col] = 0
        return table


    def get_table_score(self, i, j):
        """Bypasses cyclic reading in Python to returns values from DP-Table.
        
        Arguments:
            i (int):
                Row in DP-Table
            j (int):
                Column in DP-Table

        Returns:
            (int):
                Value at position (i, j) in DP-Table
        """
        if j < 0:
            return 0
        else:
            return self.table[i][j]


    def trace_non_crossings(self, i, j):
        """Recursive function which traces the non-crossings in the RNA structure,
        denoting the left prediction of the pair with '{', and the right prediction
        of the base pair with '}'. This is an extension of the Nussinov Algorithm.

        Arguments:
            i (int):
                Left slice of the input sequence
            j (int):
                Right slice of the input sequence
        """
        if j <= i:
            return
        elif self.get_table_score(i, j) == self.get_table_score(i, j-1):
            self.trace_non_crossings(i, j-1)
            return
        else:
            for k in range(i, j):
                pairs = [self.n_base_pairs[k], self.n_base_pairs[j]]
                if self.is_legal(pairs[0]) and self.is_legal(pairs[1]):
                    score = self.legal_pair_matrix[self.pair_matrix_mapping[pairs[0]]][self.pair_matrix_mapping[pairs[0]]]
                    if self.get_table_score(i, j) == self.get_table_score(i, k-1) + self.get_table_score(k+1, j-1) + score:
                        self.base_pairs_braces.append(tuple([k, j]))
                        self.trace_non_crossings(i, k-1)
                        self.trace_non_crossings(k+1, j-1)
                        return


    def decode_output(self):
        """Decodes output from :class: ``PredictStableStructure`` attribute
        ``base_pair_braces``/

        Returns:
            output (str):
                Joined string of denotational braces. 
        """
        output = ['.' for i in range(len(self.n_base_pairs))]
        for braces in self.base_pairs_braces:
            output[braces[0]] = '{'
            if output[braces[0]+1] == '.':
                output[braces[0]+1] = '{'
            if output[braces[1]] == '.':
                output[braces[1]] = '}'
            output[braces[1]+1] = '}'
        return ''.join(output)
