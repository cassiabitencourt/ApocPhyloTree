from Bio import Entrez, SeqIO

# Define your email address for Entrez (required)
Entrez.email = "ca.biten@gmail.com"

# Define the family name and the gene of interest (rps16)
family_name = "Apocynaceae"
gene_of_interest = "rps16"

# Initialize a list to store the accession numbers
accession_numbers = []

# Set the maximum number of records to retrieve per query
retmax = 10000

# Search for sequences in GenBank and retrieve all matching records
query = f'{family_name} [ORGN] AND {gene_of_interest} [GENE]'
handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
record = Entrez.read(handle)
accession_numbers.extend(record["IdList"])

# Download and save the sequences in a FASTA file
output_file = f"{family_name}_{gene_of_interest}_sequences.fasta"
with open(output_file, "w") as fasta_file:
    for accession in accession_numbers:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        SeqIO.write(record, fasta_file, "fasta")

print(f"Downloaded {len(accession_numbers)} sequences and saved to {output_file}")
