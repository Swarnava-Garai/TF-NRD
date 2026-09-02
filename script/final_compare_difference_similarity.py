import os
import sys

def parse_atoms(file_lines):
    atoms = set()
    for line in file_lines:
        if line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            res_num = line[22:26].strip()
            atoms.add((atom_name, res_num))
    return atoms

def compare_int_files(file1_lines, file2_lines):
    atoms1 = parse_atoms(file1_lines)
    atoms2 = parse_atoms(file2_lines)

    unique_to_file1 = atoms1 - atoms2
    unique_to_file2 = atoms2 - atoms1
    common_atoms = atoms1 & atoms2

    print(f"Atoms only in File 1: {len(unique_to_file1)}")
    print(f"Atoms only in File 2: {len(unique_to_file2)}")
    print(f"Common atoms in both files: {len(common_atoms)}")

    return len(common_atoms), unique_to_file1, unique_to_file2, common_atoms

def count_lines(file1_lines, file2_lines):
    len1 = len([l for l in file1_lines if l.strip()])
    len2 = len([l for l in file2_lines if l.strip()])
    print(f"Largest file has {max(len1, len2)} lines.")
    return max(len1, len2)

# === Configuration ===
base_directory = "/home/swarnava/PhD_CSB/Objective_1/Transcription_Factors/EM/"
input_list_filename = "generated_int_paths_NMR.txt"
output_results_filename = "atom_line_difference_similarity_summary.tsv"

def main():
    if not os.path.exists(input_list_filename):
        print(f"Error: Input list file '{input_list_filename}' not found.")
        sys.exit(1)

    print(f"Reading file paths from: {input_list_filename}")
    print(f"Results will be saved to: {output_results_filename}")

    with open(input_list_filename, 'r') as infile, open(output_results_filename, 'w') as outfile:
        outfile.write("File1\tFile2\tTotal_Common_Atoms\tUnique_to_File1\tUnique_to_File2\tMax_Lines\n")

        file_paths = infile.readlines()
        i = 0
        while i < len(file_paths):
            path1_raw = file_paths[i].strip()
            if i + 1 < len(file_paths):
                path2_raw = file_paths[i + 1].strip()
            else:
                print(f"Warning: Uneven number of paths. Skipping last path: {path1_raw}")
                break

            full_path1 = os.path.join(base_directory, path1_raw)
            full_path2 = os.path.join(base_directory, path2_raw)

            print(f"\nComparing files:\n  {full_path1}\n  {full_path2}")

            try:
                with open(full_path1, 'r') as f1, open(full_path2, 'r') as f2:
                    file1_lines = f1.readlines()
                    file2_lines = f2.readlines()

                    common_count, uniq1, uniq2, common_atoms = compare_int_files(file1_lines, file2_lines)
                    max_lines = count_lines(file1_lines, file2_lines)

                    outfile.write(f"{path1_raw}\t{path2_raw}\t{common_count}\t{len(uniq1)}\t{len(uniq2)}\t{max_lines}\n")

            except FileNotFoundError as e:
                print(f"Error: {e}")
                outfile.write(f"{path1_raw}\t{path2_raw}\tERROR\t-\t-\t-\n")
            except Exception as e:
                print(f"Unexpected error: {e}")
                outfile.write(f"{path1_raw}\t{path2_raw}\tERROR\t-\t-\t-\n")

            i += 2

    print(f"\nComparison complete. Results saved to '{output_results_filename}'.")

if __name__ == "__main__":
    main()
