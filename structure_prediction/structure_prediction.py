"""
This script contains the :class: ``PredictStructure`` which computes
the maximum number of permissible letter pairs from the input sequence
of 4 letters (A, C, G, U), of which A-U, U-A, G-C, C-G are legal pairs.
Furthermore, it denotes the nucleotide bases which are predicted to be
paired with '{' and '}', and those predicted to not be paired with any
other by '.'. 

Author:
-------
Aashish Yadavally
"""

LEGAL_PAIRS = {
    ('A','U') : 1,
    ('U','A') : 1,
    ('G','C') : 1,
    ('C','G') : 1,
}

LEGAL_PAIRS_ENERGY = {
    ('A','U'): 2,
    ('U','A'): 2,
    ('G','C'): 3,
    ('C','G'): 3,
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


class PredictStructure:
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
    def __init__(self, sequence, legal_pairs, energy_argument, gap):
        """Initializes :class: ``PredictStructure``.

        Arguments:
            sequence (str):
                RNA Sequence
            legal_pairs (list):
                Dictionary of legal pairs with their corresponding
                gamma-scores.
        """
        self.n_bases = list(sequence)
        self.legal_pairs = legal_pairs
        self.gap = gap
        initial_table = self.initialize_table()
        self.table = self.fill_table(initial_table)
        self.output = ['.' for i in range(len(self.n_bases))]
        self.trace_non_crossings(0, len(self.n_bases)-1)
        self.write_output(energy_argument)


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
        table = [[-1 for i in range(len(self.n_bases))] \
                    for j in range(len(self.n_bases))]
        # Setting up base cases
        for i in range(len(self.n_bases)):
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
        if pair in list(self.legal_pairs.keys()):
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
        n = len(self.n_bases)
        diagonal_elements = [[i, i] for i in range(n)]

        for distance_from_diagonal in range(1, n):
            for i in range(n - distance_from_diagonal):
                row, col = diagonal_elements[i][0], diagonal_elements[i][1] + distance_from_diagonal
                if col - row > self.gap:
                    pair = tuple([self.n_bases[row], self.n_bases[col]])
                    conditions = []

                    if self.is_legal(pair):             
                        gamma_score = self.legal_pairs[pair]
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
            for k in range(i, j-self.gap):
                pair = tuple([self.n_bases[k], self.n_bases[j]])
                if self.is_legal(pair):
                    if self.get_table_score(i, j) == self.get_table_score(i, k-1) + self.get_table_score(k+1, j-1) + self.legal_pairs[pair]:
                        self.output[k] = '{'
                        self.output[j] = '}'
                        self.trace_non_crossings(i, k-1)
                        self.trace_non_crossings(k+1, j-1)
                        return


    def write_output(self, energy_argument):
        """Writes output to the file 'output.txt'

        Arguments:
            energy_argument (bool):
                True, if LEGAL_PAIRS_ENERGY
                False, if LEGAL_PAIRS
        """
        if energy_argument:
            with open("output.txt", "w") as output_file:
                output_file.write(f'> total score\n')
                output_file.write(str(self.table[0][len(self.n_bases)-1]))
                output_file.close()
        else:
            with open("output.txt", "w") as output_file:
                output_file.write(f"> {''.join(self.n_bases)}\n")
                output_file.write(f"{''.join(self.output)}\n")
                output_file.write(f'> max count of pairs\n')
                output_file.write(str(self.table[0][len(self.n_bases)-1]))
                output_file.close()
