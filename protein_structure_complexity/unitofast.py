import csv
import requests
import os

# Function to download FASTA sequence from UniProt using the UniProt ID
def download_fasta(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to retrieve sequence for {uniprot_id}")
        return None

# Function to save FASTA sequence to a specific folder
def save_fasta(fasta_data, uniprot_id, folder_path):
    if fasta_data:
        
        os.makedirs(folder_path, exist_ok=True)
        
        
        file_name = f"{uniprot_id}.fasta"
        file_path = os.path.join(folder_path, file_name)
        
        # Save the FASTA file to the desired folder
        with open(file_path, "w") as fasta_file:
            fasta_file.write(fasta_data)
        print(f"FASTA sequence saved to {file_path}")

# Function to process the CSV file and download sequences
def process_csv_and_download(csv_file, folder_path):
    with open(csv_file, mode='r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            uniprot_id = row[0]  # Reading the UniProt ID from the first (and only) column
            print(f"Processing UniProt ID: {uniprot_id}")
            fasta_data = download_fasta(uniprot_id)
            save_fasta(fasta_data, uniprot_id, folder_path)

# Usage
csv_file = '/home/kuranage/Downloads/Protein Complexity - Sheet3 (3).csv'  # CSV file with uni_prot ids
folder_path = 'fasta_enzy'    # Output Folder
process_csv_and_download(csv_file, folder_path)
