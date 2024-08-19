import os
import csv

# for access to the data files from PMN by pathway:


def find_second_data_directories(root_path):
    second_data_directories = []

    # Define a function to search for second data directories
    def search_for_second_data(directory_path, depth):
        if depth == 3:
            second_data_directories.append(directory_path + "/data")
            return
        for entry in os.listdir(directory_path):
            full_path = os.path.join(directory_path, entry)
            if os.path.isdir(full_path):
                search_for_second_data(full_path, depth + 1)  # Recursively search within subdirectories


    # there are two different options, since there are directories with only the species name dir, an unknown dir and then the data
    # and one that has an additional directory inbetween with species_name.tar as name

    # Start the search from the root path
    search_for_second_data(root_path, 0)

    return second_data_directories


# Function to extract species name from path
def extract_species_name(path):
    return os.path.basename(os.path.dirname(path))


# Example usage:
root_directory = '/home/julia/Bachelor_Arbeit/PMN/data'

second_data_directories = find_second_data_directories(root_directory)

if second_data_directories:
    sorted_directories = sorted(second_data_directories, key=extract_species_name)

    print("Sorted species directory paths constructed.")
else:
    print("No second data directories found.")

# print(len(second_data_directories)) solution: #128 = |directories| -> correct!

##########################################################################

"""
# find all reactions by starting from list of lipids:
def create_dictionary_from_csv(csv_file):
    data_dict = {}
    with open(csv_file, 'r', newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:  # Ensure there are at least two columns in the row
                data_dict[row[1]] = row[0]
    return data_dict


csv_file = 'examples_for_data/standard_lipid_classes.csv'  # Replace 'example.csv' with the path to your CSV file
csv_data = create_dictionary_from_csv(csv_file) # dictionary of lipid names as key with class (e.g. fatty acyl) as value

"""









