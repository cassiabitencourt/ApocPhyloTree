# Open the input FASTA file for reading and the output FASTA file for writing
with open('Apocynaceae_rps16_sequences.fasta', 'r') as input_file, open('formatted_output.fasta', 'w') as output_file:
    # Iterate through the lines in the FASTA file
    for line in input_file:
        if line.startswith('>'):
            # Extract the Genbank accession and the first two words of the description
            header_parts = line.strip().split(' ')
            accession = header_parts[0][1:]
            species_name = '_'.join(header_parts[1:3])

            # Write the formatted header to the output file
            formatted_header = f'>{species_name} [{accession}]'
            output_file.write(formatted_header + '\n')
        else:
            # Write sequence lines as they are to the output file
            output_file.write(line)
