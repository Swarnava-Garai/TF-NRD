import os
import sys 

def generate_int_paths(ids_list, output_filename="generated_int_paths.txt"):
    """
    Generates .int file paths in the format 'Protein-NA/ID/ID.int'
    for a given list of IDs and saves them to an output file.

    Args:
        ids_list (list): A list of strings, where each string is an ID
                         in the format PDBID_Chain (e.g., '1A02_FA').
        output_filename (str): The name of the file to which the generated
                               paths will be written.
    """
    print(f"Generating .int file paths and saving to: {output_filename}")

    try:
        with open(output_filename, 'w') as outfile:
            for unique_id in ids_list:
                # Construct the path in the desired format
                # Example: Protein-NA/1A02_FA/1A02_FA.int
                generated_path = os.path.join("81_EM_NMR_TF_interface", unique_id, f"{unique_id}.int")
                outfile.write(f"{generated_path}\n")
        print("Successfully generated all paths.")
    except IOError as e:
        print(f"Error writing to file '{output_filename}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # --- Get input IDs from a user-specified text file ---
    input_ids_filename = input("Please enter the name of the text file containing the IDs (e.g., my_ids.txt): ").strip()

    input_ids = []
    try:
        with open(input_ids_filename, 'r') as infile:
            for line in infile:
                stripped_line = line.strip()
                if stripped_line: # Only add non-empty lines
                    input_ids.append(stripped_line)
        print(f"Successfully read {len(input_ids)} IDs from '{input_ids_filename}'.")
    except FileNotFoundError:
        print(f"Error: The file '{input_ids_filename}' was not found.")
        sys.exit(1) # Exit the script if the file is not found
    except Exception as e:
        print(f"An error occurred while reading '{input_ids_filename}': {e}")
        sys.exit(1) # Exit the script on other read errors

    if not input_ids:
        print("No IDs were found in the input file. Exiting.")
        sys.exit(0)

    # Call the function to generate and save the paths
    # The output file name remains "generated_int_paths.txt" by default, or you can change it here.
    generate_int_paths(input_ids, "generated_int_paths.txt")
