import argparse

def check_sequence(seq):
    # Check if all characters in the sequence are valid amino acids
    valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    if not all(aa in valid_amino_acids for aa in seq):
        print("Sequence contains invalid characters; it should only include standard amino acids.")
        return 0

    # Check if the sequence length is less than 11
    if len(seq) < 11:
        print("Sequence is too short; it must be at least 12 amino acids long.")
        return 0

    # Check if positions 7, 8, 9 (indices 6, 7, 8) contain at least one E, D, or N
    critical_positions = seq[6:9]
    if not any(aa in critical_positions for aa in "ED"):
        print("No isopeptide found; the 7th/8th/9th amino acids should contain at least one E (Glu) or D (Asp).")
        return 0

    # Determine the first position within the critical positions where E, D, or N is found
    last_pos = min(critical_positions.rfind(aa) for aa in "ED")

    # Calculate the length from this position to the end of the sequence
    remaining_length = len(seq) - (6 + last_pos + 1)

    ## Check if the remaining length is over 50
    #if remaining_length > 50:
    #    print("Sequence's loop length may exceed the scaffold range.")
    #    return 0

    # Check if the remaining length is less than 5
    if remaining_length < 5:
        print("Sequence's loop length and tail length appear to be less than 5 (the shortest naturally discovered lasso peptide has a loop of 3 and a tail of 2).")
        return 0

    # If all conditions are satisfied
    return 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check the amino acid sequence for isopeptide conditions.")
    parser.add_argument("-seq", required=True, help="Amino acid sequence")

    args = parser.parse_args()
    result = check_sequence(args.seq)
    
    if result == 1:
        print(result)



