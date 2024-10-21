import sys
import csv

def matrix_builder(sequence_1, sequence_2):
    """
    Builds the scoring and traceback matrices using dynamic programming.

    Parameters:
    - sequence_1 (str): The first DNA sequence.
    - sequence_2 (str): The second DNA sequence.

    Returns:
    - traceback_matrix (list of lists): Matrix indicating the direction of movement (up, left, diagonally).
    - score (int): The final alignment score (number in the bottom right corner of the scoring matrix).
    
    """

    # Initializing variables
    m = len(sequence_1) + 1
    n = len(sequence_2) + 1
    d = -2 # Gap penalty
    MATCHING_SCORE = 1
    MISMATCH_PENALTY = -1
    scoring_matrix = []
    traceback_matrix = []

    # Initializing scoring matrix with all 0s and traceback matrix with all " "
    for i in range(m):
        row_scoring_matrix = []
        row_traceback_matrix = []
        for j in range(n):
            row_scoring_matrix.append(0)
            row_traceback_matrix.append(" ")
        scoring_matrix.append(row_scoring_matrix)
        traceback_matrix.append(row_traceback_matrix)
    
    # Initilization step for scoring and traceback matrix
    for i in range(1, m):
        scoring_matrix[i][0] = i * d
        traceback_matrix[i][0] = "U"
    for j in range(1, n):
        scoring_matrix[0][j] = j * d
        traceback_matrix[0][j] = "L"

    # Filling the scoring matrix and traceback matrix 
    for i in range(1, m):
        for j in range(1, n):
            if sequence_1[i-1] == sequence_2[j-1]:
                matching_score = MATCHING_SCORE
            else:
                matching_score = MISMATCH_PENALTY

            # Calculating score for a match, gap from above or gap from left 
            sequence_match = matching_score + scoring_matrix[i-1][j-1]
            up_gap_penalty = scoring_matrix[i-1][j] + d
            left_gap_penalty = scoring_matrix[i][j-1] + d

            # Getting the maximum score from the three
            scoring_matrix[i][j] = max(sequence_match, up_gap_penalty, left_gap_penalty)
            
            # Filling up traceback matrix according to the score received (from a match, mismatch or a gap penalty)
            # Tie-breaking priority: U < L < D
            if scoring_matrix[i][j] == up_gap_penalty:
                traceback_matrix[i][j] = "U"
            elif scoring_matrix[i][j] == left_gap_penalty:
                traceback_matrix[i][j] = "L"
            elif scoring_matrix[i][j] == sequence_match:
                traceback_matrix[i][j] = "D"
    
    # Save the final score for each pair of sequences       
    score = scoring_matrix[-1][-1]
    
    # Return traceback matrix and final score
    return traceback_matrix, score

def backtracking(traceback_matrix, sequence_1, sequence_2):
    """
    Backtracks the aligned sequences using the traceback matrix.

    Parameters:
    - traceback_matrix (list of lists): Matrix indicating the direction of movement (up, left, diagonally).
    - sequence_1 (str): The first DNA sequence.
    - sequence_2 (str): The second DNA sequence.

    Returns:
    - alignment1 (str): The alignment of sequence_1, including gaps ("-") where needed.
    - alignment2 (str): The alignment of sequence_2, including gaps ("-") where needed.

    """

    # Initialize alignments as empty strings
    alignment1 = ""
    alignment2 = ""

    # Sets pointers to the bottom-right corner of the traceback matrix, this is where the backtracking will begin
    i = len(sequence_1)
    j = len(sequence_2)

    # Loop through all the characters from the traceback matrix
    while i > 0 or j > 0:

        # If traceback matrix indicates "U" move up
        if traceback_matrix[i][j] == "U":
            alignment1 = sequence_1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i-=1

        # If traceback matrix indicates "L" move left
        elif traceback_matrix[i][j] == "L":
            alignment1 = "-" + alignment1
            alignment2 = sequence_2[j-1] + alignment2
            j-=1 

        # If traceback matrix indicates "D" move diagonally
        elif traceback_matrix[i][j] == "D":
            alignment1 = sequence_1[i-1] + alignment1
            alignment2 = sequence_2[j-1] + alignment2
            i-=1
            j-=1

    # Return the final alignments
    return alignment1, alignment2

def needleman_wunsch(sequence_1, sequence_2):
    ''' 
    Implements the Needleman-Wunsch algorithm for global sequence alignment of two DNA sequences, using a dynamic programming approach.
    It first calls the 'matrix_builder' function to build the scoring and traceback matrices, calculating the alignment score based on match, mismatch, and gap penalties.
    It then calls the 'backtracking' function that identifies the optimal alignment path using the traceback matrix.
    The function returns the aligned sequences and the final alignment score.

    Parameters:
    - sequence_1 (str): The first DNA sequence.
    - sequence_2 (str): The second DNA sequence.

    Returns:
    - alignment1 (str): The alignment of sequence_1, including gaps ("-") where needed.
    - alignment2 (str): The alignment of sequence_2, including gaps ("-") where needed.
    - score (int): The final alignment score (number in the bottom right corner of scoring matrix).

    '''

    # Build the matrices
    traceback_matrix, score = matrix_builder(sequence_1, sequence_2)
    
    # Backtracking to get the alignment
    alignment1, alignment2 = backtracking(traceback_matrix, sequence_1, sequence_2)

    # Return the final alignments and score
    return alignment1, alignment2, score

# Main code
if len(sys.argv) > 1:
    filename = sys.argv[1]

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        
        next(reader)  # skip header

        dna_sequences_list = []
        for row in reader:
            pair_of_sequences_list = []
            for sequence in row:
                pair_of_sequences_list.append(sequence)
            dna_sequences_list.append(pair_of_sequences_list)

    # Take each pair of sequences
    for sequences in dna_sequences_list:
        sequence_1 = sequences[0]
        sequence_2 = sequences[1]
        
        # Output the result for each pair of sequences
        alignment1, alignment2, score = needleman_wunsch(sequence_1, sequence_2)
        print(alignment1, alignment2, score)