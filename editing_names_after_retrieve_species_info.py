# Create a dictionary to store the species information from the text file
species_info = {}

# Read species information from the text file
with open('species_info.txt', 'r') as species_file:
    for line in species_file:
        if line.strip():  # Ignore empty lines
            parts = line.strip().split('\t')
            accession = parts[0].split(': ')[1]
            species = parts[1].split(': ')[1].replace(' ', '_')
            species_info[accession] = species

# Open the input FASTA file for reading and the output FASTA file for writing
with open('input_to_retrieve.fasta', 'r') as input_file, open('output.fasta', 'w') as output_file:
    # Initialize variables to store the current sequence and accession
    sequence = ''
    accession = None

    # Iterate through the lines in the FASTA file
    for line in input_file:
        if line.startswith('>'):
            # If there's a previous sequence, write it to the output file with updated header
            if sequence and accession:
                formatted_header = f'{species_info[accession]} [{accession}]'
                output_file.write(f'>{formatted_header}\n')
                output_file.write(sequence)

            # Extract the accession from the FASTA header
            header = line.strip()[1:]
            accession = header.split()[0]
            sequence = ''  # Reset the sequence

        else:
            # Concatenate sequence lines
            sequence += line

    # Write the last sequence with updated header
    if sequence and accession:
        formatted_header = f'{species_info[accession]} [{accession}]'
        output_file.write(f'>{formatted_header}\n')
        output_file.write(sequence)
