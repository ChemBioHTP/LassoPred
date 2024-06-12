###################
import os
import pickle
import numpy as np
import argparse
###################

def load_model_from_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)

def load_features_from_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)

def predict_iso_with_loaded_model(sequence, model_path, k):
    model = load_model_from_pickle(model_path)
    predictions = predict_sequence(sequence, model, k)
    # Find the top three k-mer pairs
    top_three_pairs_indices = find_consecutive_zero_kmers(sequence, model, k)
    penultimate_indices = []
    penultimate_amino_acid_list=[]
    for idx in top_three_pairs_indices:
        penultimate_index, penultimate_amino_acid = find_penultimate_amino_acid(idx, sequence, k)
        real_penultimate_index = penultimate_index + 1
        penultimate_indices.append(real_penultimate_index)
        penultimate_amino_acid_list.append(penultimate_amino_acid)
    #print(penultimate_amino_acid_list, penultimate_indices)
    pred_iso_pos = final_iso(sequence,penultimate_amino_acid_list, penultimate_indices)
    #print(pred_iso_pos)
    return pred_iso_pos



def predict_upper_plug_with_loaded_model(sequence, model_path, k, iso_position):
    model = load_model_from_pickle(model_path)
    predictions = predict_sequence(sequence, model, k)
    # Find the top three k-mer pairs
    top_three_pairs_indices = find_consecutive_zero_kmers_for_upper(sequence, model, k, iso_position)
    penultimate_indices = []
    penultimate_amino_acid_list=[]
    for idx in top_three_pairs_indices:
        penultimate_index, penultimate_amino_acid = find_penultimate_amino_acid(idx, sequence, k)
        real_penultimate_index = penultimate_index + 1
        penultimate_indices.append(real_penultimate_index)
        penultimate_amino_acid_list.append(penultimate_amino_acid)
    #print(penultimate_amino_acid_list, penultimate_indices)
    return penultimate_indices

def predict_sequence(sequence, model, k):
    # Padding to ensure coverage at the sequence end
    sequence += "B" * (k-1)
    kmers = [sequence[i:i+k] for i in range(len(sequence) + 1 - k)]
    # Encode and predict all kmers
    encoded_kmers = encode_kmers_with_features(kmers)
    prob_matrix = model.predict_proba(encoded_kmers)
    predictions = model.predict(encoded_kmers)
    return predictions

def encode_kmers_with_features(kmers):
    encoded_kmers = []
    for kmer in kmers:
        encoded_kmer = sum([extract_features(aa) for aa in kmer], [])
        encoded_kmers.append(encoded_kmer)
    return np.array(encoded_kmers)

# Initialize global variable
features = None
def extract_features(aa):
    global features
    # Load feature data on first call
    if features is None:
        # Define the path to the file in the same directory as this script
        current_path = os.path.dirname(__file__)
        feature_path = os.path.join(current_path, 'features.pkl')
        with open(feature_path, 'rb') as f:
            features = pickle.load(f)
    # Return the features for the requested amino acid, or default values if not found
    return features.get(aa, [0] * 20)

def find_consecutive_zero_kmers(sequence, model, k):
    # Padding to ensure coverage at the sequence end
    sequence += "B" * (k-1)
    full_length = len(sequence)
    kmers = [sequence[i:i+k] for i in range(len(sequence) + 1 - k)]
    # Define the valid range for k-mer consideration, excluding specified first and last k-mers
    start_index = 5
    end_index = 9  # Adjust for the actual end index exclusive of padding
    # Only encode and predict kmers within the valid range
    valid_kmers = kmers[start_index:end_index]
    if not valid_kmers:  # Check if there are any valid kmers to process
        return []
    encoded_kmers = encode_kmers_with_features(valid_kmers)
    prob_matrix = model.predict_proba(encoded_kmers)
    average_probs = []
    # Adjust the loop to consider only two consecutive kmers
    for i in range(len(prob_matrix) - 1):  # Adjust to -1 for two consecutive kmers
        avg_prob_zero = (prob_matrix[i][0] * prob_matrix[i+1][0])   # Calculate the average probability
        adjusted_index = i + start_index
        average_probs.append((adjusted_index, avg_prob_zero))
    # Sort the sequences by their average probability and select the top 3
    top_consecutive_sequences = sorted(average_probs, key=lambda x: x[1], reverse=True)[:3]
    return [seq[0] for seq in top_consecutive_sequences]


def find_consecutive_zero_kmers_for_upper(sequence, model, k, iso_position):
    # Padding to ensure coverage at the sequence end
    full_length = len(sequence)
    sequence += "B" * (k-1)
    kmers = [sequence[i:i+k] for i in range(len(sequence) + 1 - k)]
    # Define the valid range for k-mer consideration, excluding specified first and last k-mers
    start_index = iso_position + 1
    end_index = full_length - 2  # Adjust for the actual end index exclusive of padding
    # Only encode and predict kmers within the valid range
    valid_kmers = kmers[start_index:end_index]
    if not valid_kmers:  # Check if there are any valid kmers to process
        return []
    encoded_kmers = encode_kmers_with_features(valid_kmers)
    prob_matrix = model.predict_proba(encoded_kmers)
    average_probs = []
    # Adjust the loop to consider only two consecutive kmers
    for i in range(len(prob_matrix) - 1):  # Adjust to -1 for two consecutive kmers
        avg_prob_zero = (prob_matrix[i][0] * prob_matrix[i+1][0])   # Calculate the average probability
        adjusted_index = i + start_index
        average_probs.append((adjusted_index, avg_prob_zero))
    # Sort the sequences by their average probability and select the top 3
    top_consecutive_sequences = sorted(average_probs, key=lambda x: x[1], reverse=True)[:3]
    return [seq[0] for seq in top_consecutive_sequences]

def find_penultimate_amino_acid(kmer_index, sequence, k):
    # Calculate the index of the penultimate amino acid within the k-mer
    penultimate_amino_acid_index = kmer_index + k - 2
    # Find the amino acid at this index
    penultimate_amino_acid = sequence[penultimate_amino_acid_index]
    return penultimate_amino_acid_index, penultimate_amino_acid

def predict_sequence(sequence, model, k):
    # Padding to ensure coverage at the sequence end
    sequence += "B" * (k-1)
    kmers = [sequence[i:i+k] for i in range(len(sequence) + 1 - k)]
    
    # Encode and predict all kmers
    encoded_kmers = encode_kmers_with_features(kmers)
    prob_matrix = model.predict_proba(encoded_kmers)
    predictions = model.predict(encoded_kmers)
    return predictions

def final_iso(sequence, penultimate_amino_acid_list, penultimate_indices):
    # Pairs up penultimate_indices and penultimate_amino_acid_list using zip for easy comparison
    sequence_match = all(sequence[i-1] == aa for i, aa in zip(penultimate_indices, penultimate_amino_acid_list))
    # Checks if at least one of the specified positions in the sequence contains 'E' or 'D'
    contains_EDN = any(sequence[i-1] in ['E', 'D'] for i in penultimate_indices)
    # If the sequence at the specified positions matches the amino acid list and contains either 'E' or 'D'
    if sequence_match and contains_EDN:
        # Iterates over penultimate_indices again to find the specific positions where the residues are 'E' or 'D'
        for index in penultimate_indices:
            if sequence[index-1] in ['E', 'D']:
                return index  # Returns the first index where the condition is met as iso_position
    else:
        print("No isopeptide bond found.")
    return None  # If no suitable residue is found, returns None


def main():

    parser = argparse.ArgumentParser(description="Predict ring length and loop length using trained models.")
    parser.add_argument('-seq', '--sequence', type=str, required=True, help='The sequence for prediction.')
    args = parser.parse_args()

    #sequence="GGPLAGEEMGGITT"
    sequence = args.sequence
    # Get the directory path of the current script file
    current_path = os.path.dirname(__file__)
    # Define the paths for the model files
    iso_model_path = os.path.join(current_path, "iso_no_model.pkl")
    upper_plug_model_path = os.path.join(current_path, "upper_plup_model.pkl")
    k=3
    iso_pos = predict_iso_with_loaded_model(sequence=sequence, model_path=iso_model_path, k=k)
    upper_plug = predict_upper_plug_with_loaded_model(sequence=sequence, model_path=upper_plug_model_path, k=k, iso_position=iso_pos)
    ring_len = iso_pos
    loop_len = [x - iso_pos for x in upper_plug]
    print(sequence,",",ring_len,",",loop_len)
    

if __name__ == "__main__":
    main()
