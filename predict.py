import os
from datetime import datetime
import sys
import csv
import argparse
import shutil
import zipfile
# Default path
construc_package_path = './construc_package'
sys.path.append(construc_package_path) 
import construc_package.LassoPred_classifier as LassoPred_classifier
from construc_package.seq_check import check_sequence
from construc_package.Class_LassoConstructor import SequenceProcessor

def create_timed_folder(target_name, timestamp):
    folder_name = f"{target_name}_{timestamp}"
    path = os.path.join(os.getcwd(), folder_name)
    os.makedirs(path, exist_ok=True)
    return path

def create_result_timed_folder(target_name, timestamp):
    folder_name = f"{target_name}_result_{timestamp}"
    path = os.path.join(os.getcwd(), f"{target_name}_{timestamp}", folder_name)
    os.makedirs(path, exist_ok=True)
    return path

def predict_annotation(sequence):
    #sequence = args.sequence
    iso_model_path = os.path.join(construc_package_path, "iso_model.pkl")
    upper_plug_model_path = os.path.join(construc_package_path, "upper_plug_model.pkl")
    k=3
    iso_pos = LassoPred_classifier.predict_iso_with_loaded_model(sequence=sequence, model_path=iso_model_path, k=k)
    upper_plug = LassoPred_classifier.predict_upper_plug_with_loaded_model(sequence=sequence, model_path=upper_plug_model_path, k=k, iso_position=iso_pos)
    ring_len = iso_pos
    loop_lens = [x - iso_pos for x in upper_plug]
    print(sequence,",",ring_len,",",loop_lens)
    return ring_len, loop_lens

def compress_folder(folder_path, output_zip_file):
    """ Compress the folder into a zip file """
    with zipfile.ZipFile(output_zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                # Create the path relative to the folder being zipped
                file_path = os.path.join(root, file)
                zipf.write(file_path, os.path.relpath(file_path, os.path.join(folder_path, '..')))
    print(f"Folder {folder_path} has been compressed into {output_zip_file}")

def copy_files_to_result_folder(target_folder, result_folder, num_ranks=3):
    #target_file = os.path.join(os.getcwd(), f"{target_name}_{timestamp}")
    print(target_file)
    #result_folder = os.path.join(os.getcwd(), f"{target_name}_result_{timestamp}")
    print(result_folder)
    # Loop through each rank and copy necessary files
    for rank in range(1, num_ranks + 1):
        # Define source and destination paths for PDB files
        pdb_source = os.path.join(target_file, f"rank{rank}_wkdir/final_output/min1_clean.pdb")
        print(pdb_source)
        pdb_destination = os.path.join(result_folder, f"min{rank}.pdb")
        print(pdb_destination)
        
        # Check if source PDB file exists before copying
        if os.path.exists(pdb_source):
            shutil.copy(pdb_source, pdb_destination)
            print(f"Copied: {pdb_source} to {pdb_destination}")
        else:
            print(f"File not found: {pdb_source}, skipping copy.")
        
        # Define source and destination paths for MD production files
        md_source = os.path.join(target_file, f"rank{rank}_wkdir/raw_output/MDprod")
        print(pdb_source)
        md_destination = os.path.join(result_folder, f"MDfiles{rank}")
        print(md_destination)
        
        # Check if source MD directory exists before copying
        if os.path.exists(md_source):
            shutil.copytree(md_source, md_destination)
            print(f"Copied: {md_source} to {md_destination}")
        else:
            print(f"Directory not found: {md_source}, skipping copy.")
    

def create_summary_csv(filename, sequence, ring_len, loop_lens):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Rank', 'Sequence', 'Ring Length', 'Loop Length'])
        for i, loop_len in enumerate(sorted(loop_lens), start=1):
            writer.writerow([f'rank{i}', sequence, ring_len, loop_len])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process lasso peptide structure prediction")
    parser.add_argument('-seq', '--sequence', type=str, required=True, help='The sequence for prediction')
    parser.add_argument('-tname', '--target_name', type=str, required=True, help='The target directory name')
    parser.add_argument('-fdir', '--fold_direction', type=str, choices=['left', 'right'], default='right', help='Direction of lasso peptide ring folding (default: right)')
    args = parser.parse_args()

    args = parser.parse_args()

    sequence = args.sequence
    target_name = args.target_name
    fold_direction = args.fold_direction
    
    # Check the sequence and decide to continue or exit
    if check_sequence(sequence) == 0:
        print("Exiting program.")
        sys.exit()  

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    ring_len, loop_lens = predict_annotation(sequence)
    target_file = create_timed_folder(target_name, timestamp)
    result_folder = create_result_timed_folder(target_name, timestamp)
    print(result_folder)
    main_script_dir = os.path.dirname(os.path.abspath(__file__))

    workdirs_to_delete = []
    # Create working directories for each loop length
    for i, loop_len in enumerate(loop_lens, start=1):
        if loop_len > 50:
            #print(f"Skipping loop length {loop_len} as it is greater than 50.")
            continue
        workdir = f"rank{i}_wkdir"
        workdir_path = os.path.join(target_file, workdir)
        os.makedirs(workdir_path, exist_ok=True)
        workdirs_to_delete.append(workdir_path)
        # Initialize and process sequence with SequenceProcessor
        processor = SequenceProcessor(sequence=sequence, ring_len=ring_len, loop_len=loop_len, fold_direction=fold_direction, wk_dir=workdir_path, main_script_dir=main_script_dir)
        processor.process_sequence()
        print(f"Processing sequence with Ring Length: {ring_len} and Loop Length: {loop_len} in {workdir}")

    os.chdir(target_file)
    copy_files_to_result_folder(target_folder=target_file, result_folder=result_folder)
    summary_csv_path = os.path.join(result_folder, 'summary.csv')
    create_summary_csv(summary_csv_path, sequence, ring_len, loop_lens)
    # Compress the result folder into a ZIP file
    zip_file_name = f"{result_folder}.zip"
    #compress_folder(result_folder, zip_file_name)
    
    # Delete all work directories
    for workdir_path in workdirs_to_delete:
        shutil.rmtree(workdir_path)
        print(f"Deleted work directory {workdir_path}")