import random
import torch
import esm
import numpy as np
import pandas as pd
import os
import pickle
import argparse

random.seed(42)
np.random.seed(42)

def load_scaler_and_pca(scaler_path, pca_path):
    with open(scaler_path, 'rb') as f:
        scaler = pickle.load(f)
    #print(f"Scaler loaded from {scaler_path}")
    if os.path.exists(pca_path):
        with open(pca_path, 'rb') as f:
            pca = pickle.load(f)
        #print(f"PCA loaded from {pca_path}")
    else:
        pca = None
        print(f"No PCA file found at {pca_path}. PCA will not be applied.")
    
    return scaler, pca

def load_model_from_pickle(file_path):
    with open(file_path, 'rb') as file:
        model = pickle.load(file)
    #print(f"Model loaded from {file_path}")
    return model

def load_model(device):
    """Load the ESM-2 model and return the model and alphabet."""
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device)
    model.eval()
    return model, alphabet

def get_pre_residue_representation(sequence, model, alphabet, device,truncation_seq_length=1022):
    """Compute per-residue representations for a list of protein sequences."""
    data = [("sequence", sequence)]
    batch_converter = alphabet.get_batch_converter()
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens[:, :truncation_seq_length].to(device)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33].cpu().numpy()
    representation_dict = {}
    for i, label in enumerate(batch_labels):
        original_sequence_length = len(batch_strs[i])
        representation_dict[label] = token_representations[i, 1:original_sequence_length + 1]
    representation_array = representation_dict['sequence']
    return representation_array

def get_features_from_esm(sequence,k):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, alphabet = load_model(device)
    feature_matrix = get_pre_residue_representation(sequence, model, alphabet, device)
    rows_to_add = k - 1
    num_features = feature_matrix.shape[1]
    padding = np.zeros((rows_to_add, num_features))
    feature_matrix_padded = np.vstack([feature_matrix, padding])

    return feature_matrix_padded

def get_selected_features(feature_matrix, sequence, aa_index):
    selected_features = feature_matrix[aa_index]
    corresponding_aa = sequence[aa_index]
    return feature_matrix, selected_features, aa_index, corresponding_aa

def pre_generate_kmers(sequence, k):
    """Generate k-mer sequences for a given sequence, appending 'B' for tail extension."""
    sequence += 'B' * (k-1)  # Append 'B' to ensure an additional character marked as 'T'
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]


def iso_pre_generate_truncated_kmers(sequence, k):
    """Generate k-mer sequences for a given sequence."""
    sequence = sequence[5:10]
    return [sequence[i:i+k] for i in range(len(sequence)-1)]

def plug_pre_generate_truncated_kmers(sequence, k, iso_position):
    """Generate k-mer sequences for a given sequence, appending 'B' for tail extension."""
    sequence = sequence[iso_position:]
    sequence += 'B' * (k-1)  # Append 'B' to ensure an additional character marked as 'T'
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def iso_pre_generate_part_labels(ring_len, loop_len, tail_len, k):
    """Generate labels for the parts of the sequence, including the appended 'B' marked as 'T'."""
    return ['R'] * ring_len + ['L'] * loop_len   # Include the 'B'

def plug_pre_generate_part_labels(ring_len, loop_len, tail_len, k):
    """Generate labels for the parts of the sequence, including the appended 'B' marked as 'T'."""
    return ['R'] * ring_len + ['L'] * loop_len + ['T'] * (tail_len + k-1)  # Include the 'B'


def pre_determine_label_for_iso(split_parts):
    """Determine the label based on the split parts."""
    if all(part == 'R' for part in split_parts):
        return 1
    elif 'R' in split_parts and 'L' in split_parts:
        return 0
    elif all(part == 'L' for part in split_parts):
        return 2
    elif 'L' in split_parts and 'T' in split_parts:
        return 2
    elif all(part == 'T' for part in split_parts):
        return 2
    else:
        return None

def pre_determine_label_for_upper(split_parts):
    """Determine the label based on the split parts."""
    if all(part == 'R' for part in split_parts):
        return 1
    elif 'R' in split_parts and 'L' in split_parts:
        return 1
    elif all(part == 'L' for part in split_parts):
        return 1
    elif 'L' in split_parts and 'T' in split_parts:
        return 0
    elif all(part == 'T' for part in split_parts):
        return 2
    else:
        return None



def encode_kmers_with_npy_features(kmers, sequence, k, feature_matrix):
    """
    Concatenate the features of all k amino acids in the k-mer.
    """
    encoded_kmers = []
    
    for i, kmer in enumerate(kmers):
        combined_features = []
        # Retrieve the features of each amino acid in the k-mer one by one and concatenate them
        for j in range(k):
            aa_index = i + j 
            _, selected_features, _, _ = get_selected_features(feature_matrix, sequence, aa_index)
            combined_features.extend(selected_features)  # Concatenate the feature vectors
        # Add to the list of encoded k-mers
        encoded_kmers.append(combined_features)
    return np.array(encoded_kmers)

def iso_find_consecutive_zero_kmers(sequence, model, k,  scaler, pca, feature_matrix):
    truncated_sequence = sequence[5:10]
    truncated_kmers = [truncated_sequence[i:i+k] for i in range(len(truncated_sequence)-1)]
    sequence += "B" * (k-1)
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

    #encode all k-mers in the entire sequence
    encoded_kmers = encode_kmers_with_npy_features(kmers, sequence, k, feature_matrix)
    truncated_encoded_kmers = encoded_kmers[5:9, ]

    # define the effective range
    start_index = 1
    end_index = 4 # Adjust for the actual end index exclusive of padding
    valid_kmers = truncated_kmers[start_index:end_index]
    valid_encoded_kmers = truncated_encoded_kmers[start_index:end_index]


    if not valid_encoded_kmers.any():
        return []

    scaled_features = scaler.transform(valid_encoded_kmers)
    if pca is not None:
        scaled_features = pca.transform(scaled_features)
    prob_matrix = model.predict_proba(scaled_features)
    #print(scaled_features.shape)
    average_probs = []

    for i in range(len(prob_matrix)):
        avg_prob_zero = prob_matrix[i][0]
        # Adjust the index and record the results
        adjusted_index = i + start_index
        average_probs.append((adjusted_index, avg_prob_zero))
    top_consecutive_sequences = sorted(average_probs, key=lambda x: x[1], reverse=True)[:3]

    return [seq[0] for seq in top_consecutive_sequences]


def iso_find_penultimate_amino_acid(kmer_index, sequence, k):
    """
    Finds the index and the value of the penultimate amino acid of a k-mer in a sequence.
    :param kmer_index: The index of the k-mer in the sequence.
    :param sequence: The amino acid sequence.
    :param k: The length of the k-mer.
    :return: A tuple containing the index of the penultimate amino acid and its value.
    """
    # Calculate the index of the penultimate amino acid within the k-mer
    penultimate_amino_acid_index = (kmer_index + 5) + k - 2
    #real_amino_acid_index = penultimate_amino_acid_index + 5
    # Find the amino acid at this index
    penultimate_amino_acid = sequence[penultimate_amino_acid_index]
    
    return penultimate_amino_acid_index, penultimate_amino_acid

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


def plug_find_consecutive_zero_kmers(sequence, model, k, iso_position, scaler, pca, feature_matrix):

    full_length = len(sequence[iso_position:])
    sequence += "B" * (k-1)
    truncated_sequence = sequence[iso_position:]
    kmers = [sequence[i:i+k] for i in range(len(sequence) + 1 - k)]
    truncated_kmers = [truncated_sequence[i:i+k] for i in range(len(truncated_sequence) - k + 1)]
    encoded_kmers = encode_kmers_with_npy_features(kmers, sequence, k, feature_matrix)
    truncated_encoded_kmers = encoded_kmers[iso_position + 1 - (k-1):, ]

    start_index = 4 - k
    end_index = full_length - 2  # Adjust for the actual end index exclusive of padding
    valid_kmers = truncated_kmers[start_index:end_index]
    valid_encoded_kmers = truncated_encoded_kmers[start_index:end_index]

    if not valid_encoded_kmers.any():
        return []

    scaled_features = scaler.transform(valid_encoded_kmers)
    if pca is not None:
        scaled_features = pca.transform(scaled_features)
    prob_matrix = model.predict_proba(scaled_features)
    average_probs = []


    for i in range(len(prob_matrix) - (k-1-1)):
        avg_prob_zero = prob_matrix[i][0]
  
        adjusted_index = i + start_index
        average_probs.append((adjusted_index, avg_prob_zero))
    top_consecutive_sequences = sorted(average_probs, key=lambda x: x[1], reverse=True)[:3]
    return [seq[0] for seq in top_consecutive_sequences]

def plug_find_penultimate_amino_acid(kmer_index, sequence, k, iso_position):
    """Finds the index and the value of the penultimate amino acid of a k-mer in a sequence."""
    penultimate_amino_acid_index = kmer_index + k - 2
    real_amino_acid_index = penultimate_amino_acid_index + iso_position
    penultimate_amino_acid = sequence[real_amino_acid_index]
    return real_amino_acid_index, penultimate_amino_acid

def predict_iso_with_loaded_model(sequence, model_path, pca_path, scaler_path, k, feature_matrix):
    iso_model = load_model_from_pickle(model_path)
    iso_scaler, iso_pca = load_scaler_and_pca(scaler_path, pca_path)
    iso_top_three_pairs_indices=iso_find_consecutive_zero_kmers(sequence,iso_model,k,iso_scaler,iso_pca,feature_matrix)
    iso_penultimate_indices = []
    iso_penultimate_amino_acid_list=[]
    for idx in iso_top_three_pairs_indices:
        penultimate_index, penultimate_amino_acid = iso_find_penultimate_amino_acid(idx, sequence, k)
        real_penultimate_index = penultimate_index + 1
        iso_penultimate_indices.append(real_penultimate_index)
        iso_penultimate_amino_acid_list.append(penultimate_amino_acid)
    pred_iso_pos = final_iso(sequence,iso_penultimate_amino_acid_list, iso_penultimate_indices)
    return pred_iso_pos

def predict_upper_plug_with_loaded_model(sequence, model_path, pca_path, scaler_path, k, feature_matrix, iso_position):
    plug_model = load_model_from_pickle(model_path)
    plug_scaler, plug_pca = load_scaler_and_pca(scaler_path, pca_path)
    plug_top_three_pairs_indices=plug_find_consecutive_zero_kmers(sequence,plug_model,k,iso_position,plug_scaler,plug_pca,feature_matrix)
    plug_penultimate_indices = []
    plug_penultimate_amino_acid_list=[]
    for idx in plug_top_three_pairs_indices:
        penultimate_index, penultimate_amino_acid = plug_find_penultimate_amino_acid(idx, sequence, k, iso_position)
        real_penultimate_index = penultimate_index + 1
        plug_penultimate_indices.append(real_penultimate_index)
        plug_penultimate_amino_acid_list.append(penultimate_amino_acid)
    #print(plug_penultimate_amino_acid_list, plug_penultimate_indices)
    return plug_penultimate_indices


def main():
    iso_model_path = "./iso_model.pkl"
    iso_scaler_path = "./iso_scaler.pkl"
    iso_pca_path = "./iso_pca.pkl"

    plug_model_path = "./plug_model.pkl"
    plug_scaler_path = "./plug_scaler.pkl"
    plug_pca_path = "./plug_pca.pkl"

    parser = argparse.ArgumentParser(description="Predict ring length and loop length using trained models.")
    parser.add_argument('-seq', '--sequence', type=str, required=True, help='The sequence for prediction.')
    args = parser.parse_args()
    sequence = args.sequence
    #sequence = "GGAGQYKEVEAGRWSDRIDSDDE"

    k=2

    feature_matrix = get_features_from_esm(sequence,k)

    iso_pos = predict_iso_with_loaded_model(sequence,iso_model_path,iso_pca_path,iso_scaler_path,k,feature_matrix)
    upper_plug = predict_upper_plug_with_loaded_model(sequence,plug_model_path,plug_pca_path,plug_scaler_path,k,feature_matrix,iso_position=iso_pos)
    ring_len = iso_pos
    loop_len = [x - iso_pos for x in upper_plug]
    print(sequence,",",ring_len,",",loop_len)


if __name__ == "__main__":
    main()
