from Bio import Entrez

# Replace with your email address
Entrez.email = "ca.biten@gmail.com"

def retrieve_species_info(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()
        
        # Parse GenBank record to get species information
        species = None
        for line in record.split("\n"):
            if "ORGANISM" in line:
                species = line.strip().split("ORGANISM")[1].strip()
                break
        
        return species
    except Exception as e:
        print(f"Error retrieving information for {accession}: {str(e)}")
        return None

def main(input_file):
    with open(input_file, "r") as f:
        accessions = [line.strip() for line in f]

    with open("results.txt", "w") as output_file:
        for accession in accessions:
            species = retrieve_species_info(accession)
            if species:
                output_line = f"Accession: {accession}\tSpecies: {species}\n"
            else:
                output_line = f"Accession: {accession}\tSpecies: Not found\n"
            
            # Write the result to the output file
            output_file.write(output_line)

if __name__ == "__main__":
    input_file = "genbank_accessions.txt"  # Replace with your input file name
    main(input_file)
