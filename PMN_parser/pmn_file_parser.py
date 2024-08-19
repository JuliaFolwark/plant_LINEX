import json
import re
import create_pmn as species_reader
import pandas as pd
import os
import csv
import requests
from collections import OrderedDict
import tempfile
import shutil
import docker
import json
from pygoslin.parser.Parser import LipidParser

"""
- first step: test if the TYPES line contains either "Lipid-Biosynthesis" or "Fatty-acid-biosynthesis"
-> if it does:
    extract unique ID, types (pathway class), extract common-name, extract all reaction lists (there are multiple)
    (I think that predecessors and primaries are not needed)
"""


# --------------------------------------------------------------------------------------------------------------------#
# reading the standard_lipid_classes for classification

def read_standard_lipids():
    # reading the standard_lipid_classes >> will be used for classification
    df = pd.read_csv('examples_for_data/standard_lipid_classes.csv')
    lipid_stand_dict = OrderedDict((key, list(df[key])) for key in df.columns)

    lipid_stand_dict["Lipid class"] = [item.replace("\xa0", "").strip() for item in lipid_stand_dict["Lipid class"]]

    abbreviations_dict = {}
    for index, abbr in enumerate(lipid_stand_dict['Abbreviation']):
        abbr_data = {key: lipid_stand_dict[key][index] for key in lipid_stand_dict if key != 'Abbreviation'}
        abbreviations_dict[abbr] = abbr_data

    return lipid_stand_dict, abbreviations_dict


#
# --------------------------------------------------------------------------------------------------------------------
# filter classified compounds:

def find_unclassified_cmpds(file_path):
    filtered_rows = []

    with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        next(csvreader)
        for row in csvreader:
            if len(row) > 2 and row[2] == '[\'-\']':
                filtered_rows.append(row)

    return filtered_rows


# --------------------------------------------------------------------------------------------------------------------


def update_comp_file(file_path, dict_of_updates, dict_of_lines, dict_of_classifed_comps):
    """
    Modify the second column of the lines in a CSV file based on values from a dictionary.

    Parameters:
        file_path (str): Path to the CSV file.
        dict_of_updates (dict): Dictionary containing values to update in the CSV file.
        dict_of_lines (dict): Dictionary containing lines, where the compounds are saved in the CSV file.

    Returns:
        int: Number of lines modified.
    """

    # Create a temporary file to write the modified content
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, newline='', encoding='utf-8')
    temp_file_path = temp_file.name

    # find all the compounds in both dicts: (for super_comp_class)
    lines_to_modify_spc = []

    for compound in dict_of_lines.keys():
        if compound in dict_of_updates.keys():
            lines_to_modify_spc.append((dict_of_lines[compound], dict_of_updates[compound], compound))

    # Sort the list of tuples based on the integer value
    lines_to_modify_spc.sort(key=lambda x: x[0])

    cntr_tuple = 0
    cntr_line = 0

    print(lines_to_modify_spc)

    # Open the original CSV file for reading
    with open(file_path, 'r', newline='', encoding='utf-8') as file_in, temp_file:
        csv_reader = csv.reader(file_in, delimiter='\t')
        csv_writer = csv.writer(temp_file, delimiter='\t')

        for row in csv_reader:
            # first modify the super_comp_class
            if cntr_line == lines_to_modify_spc[cntr_tuple][0] and cntr_tuple < len(lines_to_modify_spc) - 1:

                if lines_to_modify_spc[cntr_tuple][1] != "COMP_NOT_FOUND":  # otherwise it will be used as supper class
                    if row[1] != "-":
                        new_string = row[1] + "; " + lines_to_modify_spc[cntr_tuple][1]
                    elif row[1] == lines_to_modify_spc[cntr_tuple][1]:
                        new_string = row[1]
                    else:
                        new_string = lines_to_modify_spc[cntr_tuple][1]
                    row[1] = new_string
                cntr_tuple += 1

            cntr_line += 1
            csv_writer.writerow(row)

    # Replace the original file with the temporary file
    shutil.move(temp_file_path, file_path)


# --------------------------------------------------------------------------------------------------------------------
# additional functions for parser:


def file_writer(content, file_path_for_saving):
    # filename = file_path_for_saving + "/" + filename
    try:
        # Try to open the file in append mode ('a')
        with open(file_path_for_saving, 'a', newline="") as file:
            csvwriter = csv.DictWriter(file, delimiter='\t', fieldnames=content.keys())
            # check if file is empty
            if file.tell() == 0:  # if empty write header
                csvwriter.writeheader()
            csvwriter.writerow(content)
    except FileNotFoundError:
        print(f"File '{file_path_for_saving}' not found.")
    except Exception as e:
        print(f"Error occurred: {str(e)}")


# --------------------------------------------------------------------------------------------------------------------

# for reaction ID comparison

def extract_first_words(filename):
    first_words = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:  # Exclude the first line
            first_word = line.strip().split()[0]
            first_words.append(first_word)
    return first_words


def compare_first_words(file1, file2):
    first_words_file1 = extract_first_words(file1)
    first_words_file2 = extract_first_words(file2)

    # Finding the common first words
    common_first_words = set(first_words_file1).intersection(set(first_words_file2))
    print("Common first reaction_IDs:", common_first_words)

    # Calculating the percentage of common first words
    total_unique_words = len(set(first_words_file1 + first_words_file2))
    print("Total unique reaction_IDs:", total_unique_words)

    if total_unique_words == 0:
        return 0

    percentage_common = (len(common_first_words) / total_unique_words) * 100
    return percentage_common


# --------------------------------------------------------------------------------------------------------------------

def find_super_compounds(compound, hierarchy_dict, hash_set_super_compounds):
    base_compounds = []

    for entry in hierarchy_dict:
        sub_comp_id = hierarchy_dict[entry]["UNIQUE-ID"][0]
        if compound in hierarchy_dict[entry]["TYPES"]:  # already in base_comp
            if sub_comp_id in base_compounds:
                continue
            else:
                if sub_comp_id not in hash_set_super_compounds:  # added to base_comp
                    base_compounds.append(sub_comp_id)
                else:
                    supper_comp = sub_comp_id
                    find_super_compounds(supper_comp, hierarchy_dict, hash_set_super_compounds)

    return base_compounds


# return list of all sub_compounds

# --------------------------------------------------------------------------------------------------------------------

# searches recursively for sub_pathways

def find_sub_pathway_reactions(set_of_pathways, set_of_reactions, hierarchy_dict, hash_set_super_pathways,
                               pathways_added, types_added):
    list_of_manual_pathways = (
        "PWY-401", "PWY-782", "PWY-6424", "PWY-6799", "NAGLIPASYN-PWY", "PWY-7590", "PWY-6468", "PWY-5973")
    if len(pathways_added) == 0:  # equivalent to first run of method

        # check if manually added pathways are in entry_dict -> if they are, add to set_of_pathways

        for pathway_in_list in list_of_manual_pathways:
            if pathway_in_list in hierarchy_dict.keys():
                set_of_pathways.add(pathway_in_list)

    for entry in hierarchy_dict.values():
        list_of_types = entry["TYPES"]

        # for hacking the manual pathways into the system

        if entry["UNIQUE-ID"][
            0] in list_of_manual_pathways:  # takes the ID and puts it in the types in order to simulate the correct pathway case
            list_of_types.append(entry["UNIQUE-ID"][0])

        temp_pathways = set()  # is used in order to update the set of pathways we search for in TYPES

        # case if in current entry are IDs of other pathways inside the links
        if "PATHWAY-LINKS" in entry:  # take the string and delete ( or ), split by " " and then take the ones that have PWY somewhere
            for link in entry[
                "PATHWAY-LINKS"]:  # -> these are then searched for with ID and added to the list of found pathways -> add reactions
                link = re.sub(r'[()]', '', link)
                pathways = link.split(" ")
                for pathway in pathways:
                    if pathway in pathways_added:
                        item = re.sub(r'[|()]+', '', pathway)  # delete "|" and "( or )"
                        if "PWY" in item and item in hierarchy_dict.keys():
                            if "SUB-PATHWAYS" not in hierarchy_dict[
                                item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                    hierarchy_dict[item]["TYPES"]:

                                # if PWY really is not a super-pathway -> need to take every reaction:
                                if "REACTION-LIST" in hierarchy_dict[item]:
                                    for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                        if reaction not in set_of_reactions and "RNX" in reaction:
                                            set_of_reactions.add(reaction)
                                if "PREDECESSORS" in hierarchy_dict[item]:
                                    for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                        pattern = r'"([^"]*)"'
                                        matches = re.findall(pattern, reaction)
                                        for match in matches:
                                            if match not in set_of_reactions and "RNX" in match:
                                                set_of_reactions.add(match)
                                pathways_added.add(entry["UNIQUE-ID"][0])
                            else:
                                set_of_pathways.add(entry["UNIQUE-ID"][0])

        if entry["UNIQUE-ID"][0] in pathways_added:
            if "PATHWAY-LINKS" in entry:  # take the string and delete ( or ), split by " " and then take the ones that have PWY somewhere
                for link in entry[
                    "PATHWAY-LINKS"]:  # -> these are then searched for with ID and added to the list of found pathways -> add reactions
                    link = re.sub(r'[()]', '', link)
                    pathways = link.split(" ")
                    for pathway in pathways:
                        item = re.sub(r'[|()]+', '', pathway)  # delete "|" and "( or )"
                        if "PWY" in item and item in hierarchy_dict.keys():
                            if "SUB-PATHWAYS" not in hierarchy_dict[
                                item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                    hierarchy_dict[item]["TYPES"]:

                                # if PWY really is not a super-pathway -> need to take every reaction:
                                if "REACTION-LIST" in hierarchy_dict[item]:
                                    for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                        if reaction not in set_of_reactions and "RNX" in reaction:
                                            set_of_reactions.add(reaction)
                                if "PREDECESSORS" in hierarchy_dict[item]:
                                    for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                        pattern = r'"([^"]*)"'
                                        matches = re.findall(pattern, reaction)
                                        for match in matches:
                                            if match not in set_of_reactions and "RNX" in match:
                                                set_of_reactions.add(match)
                                pathways_added.add(entry["UNIQUE-ID"][0])
                            else:
                                set_of_pathways.add(entry["UNIQUE-ID"][0])

        for pathway in set_of_pathways:
            if pathway not in list_of_types:  # no sub-pathway
                continue
            else:  # is sub-pathway
                # look if common-name has superpathway -> so can be taken as super-pathway
                col_name_super = False
                for name in entry["COMMON-NAME"]:
                    if "superpathway" in name:
                        col_name_super = True

                if "Super-Pathways" in list_of_types or col_name_super == True:
                    # take reaction-list
                    if "REACTION-LIST" in entry:
                        for item in entry["REACTION-LIST"]:
                            if "PWY" in item:
                                if item not in pathways_added and item in hierarchy_dict.keys():
                                    # add here to get every reaction, that the pathway has , when it is not Super-Pathway
                                    if "SUB-PATHWAYS" not in hierarchy_dict[
                                        item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                            hierarchy_dict[item]["TYPES"]:

                                        # if PWY really is not a super-pathway -> need to take every reaction:
                                        if "REACTION-LIST" in hierarchy_dict[item]:
                                            for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                                if reaction not in set_of_reactions and "RNX" in reaction:
                                                    set_of_reactions.add(reaction)
                                        if "PREDECESSORS" in hierarchy_dict[item]:
                                            for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                                pattern = r'"([^"]*)"'
                                                matches = re.findall(pattern, reaction)
                                                for match in matches:
                                                    if match not in set_of_reactions and "RNX" in match:
                                                        set_of_reactions.add(match)

                                    pathways_added.add(item)
                                    temp_pathways.add(
                                        item)  # add the pathway to the list of the ones we search, because it is super-pathway
                                    set_of_item = set()
                                    set_of_item.add(item)
                                    types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                    types_temp_set = set()

                                    for type_element in types_in_sub_pathway:
                                        if type_element not in set_of_pathways and type_element != "Super-Pathways" and type_element not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                            types_temp_set.add(type_element)
                                            types_added.add(type_element)

                                    # when should go into recursion? -> when new type is found we want to search for and Unique-ID that is a super-pathway
                                    temp_pathways.update(types_temp_set)
                                    set_of_item.update(
                                        types_temp_set)  # add the unique id of pathway and the unique type to the search query of TYPES
                                    list_new_reactions = find_sub_pathway_reactions(set_of_item,
                                                                                    set_of_reactions,
                                                                                    hierarchy_dict,
                                                                                    hash_set_super_pathways,
                                                                                    pathways_added, types_added)

                                    for reaction in list_new_reactions:
                                        if "RNX" in reaction:
                                            set_of_reactions.add(reaction)



                            else:  # if PWY not inside
                                if item not in set_of_reactions and "RNX" in item:
                                    set_of_reactions.add(item)

                    if "PREDECESSORS" in entry:
                        for item in entry["PREDECESSORS"]:
                            pattern = r'"([^"]*)"'
                            matches = re.findall(pattern, item)
                            for match in matches:
                                if "PWY" in match:
                                    if match not in pathways_added and match in hierarchy_dict.keys():

                                        if "SUB-PATHWAYS" not in hierarchy_dict[
                                            match] and match not in hash_set_super_pathways and "Super-Pathways" not in \
                                                hierarchy_dict[match]["TYPES"]:
                                            # if PWY really is not a super-pathway -> need to take every reaction:
                                            if "REACTION-LIST" in hierarchy_dict[match]:
                                                for reaction in hierarchy_dict[match]["REACTION-LIST"]:
                                                    if reaction not in set_of_reactions and "RNX" in reaction:
                                                        set_of_reactions.add(reaction)
                                            if "PREDECESSORS" in hierarchy_dict[match]:
                                                for reaction in hierarchy_dict[match]["PREDECESSORS"]:
                                                    pattern = r'"([^"]*)"'
                                                    matches = re.findall(pattern, reaction)
                                                    for match2 in matches:
                                                        if match2 not in set_of_reactions and "RNX" in match2:
                                                            set_of_reactions.add(match2)

                                        pathways_added.add(match)
                                        temp_pathways.add(match)
                                        set_of_match = set()
                                        set_of_match.add(match)
                                        types_in_sub_pathway = hierarchy_dict[match]["TYPES"]

                                        types_temp_set = set()

                                        for type_element in types_in_sub_pathway:
                                            if type_element not in set_of_pathways and type_element != "Super-Pathways" and type_element not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                                types_temp_set.add(type_element)
                                                types_added.add(type_element)

                                        temp_pathways.update(types_temp_set)
                                        set_of_match.update(
                                            types_temp_set)

                                        list_new_reactions = find_sub_pathway_reactions(set_of_match, set_of_reactions,
                                                                                        hierarchy_dict,
                                                                                        hash_set_super_pathways,
                                                                                        pathways_added, types_added)
                                        for reaction in list_new_reactions:
                                            if "RNX" in reaction:
                                                set_of_reactions.add(reaction)

                                else:  # if not PWY are found
                                    if item not in set_of_reactions and "RNX" in item:
                                        set_of_reactions.add(match)

                    if "SUB-PATHWAYS" in entry:
                        for item in entry["SUB-PATHWAYS"]:
                            if "PWY" in item:
                                if item not in pathways_added and item in hierarchy_dict.keys():

                                    if "SUB-PATHWAYS" not in hierarchy_dict[
                                        item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                            hierarchy_dict[item]["TYPES"]:
                                        # if PWY really is not a super-pathway -> need to take every reaction:

                                        if "REACTION-LIST" in hierarchy_dict[item]:
                                            for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                                if reaction not in set_of_reactions and "RNX" in reaction:
                                                    set_of_reactions.add(reaction)
                                        if "PREDECESSORS" in hierarchy_dict[item]:
                                            for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                                pattern = r'"([^"]*)"'
                                                matches = re.findall(pattern, reaction)
                                                for match in matches:
                                                    if match not in set_of_reactions and "RNX" in match:
                                                        set_of_reactions.add(match)

                                    pathways_added.add(item)
                                    temp_pathways.add(
                                        item)  # add the pathway to the list of the ones we search, because it is super-pathway
                                    set_of_item = set()
                                    set_of_item.add(item)
                                    types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                    types_temp_set = set()

                                    for type_element in types_in_sub_pathway:
                                        if type_element not in set_of_pathways and type_element != "Super-Pathways" and type_element not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                            types_temp_set.add(type_element)
                                            types_added.add(type_element)

                                    temp_pathways.update(types_temp_set)
                                    set_of_item.update(
                                        types_temp_set)

                                    list_new_reactions = find_sub_pathway_reactions(set_of_item, set_of_reactions,
                                                                                    hierarchy_dict,
                                                                                    hash_set_super_pathways,
                                                                                    pathways_added, types_added)
                                    for reaction in list_new_reactions:
                                        if "RNX" in reaction:
                                            set_of_reactions.add(reaction)



                            else:  # if PWY not inside
                                if item not in set_of_reactions and "RNX" in item:
                                    set_of_reactions.add(item)

                    if "PATHWAY-LINKS" in entry:  # take the string and delete ( or ), split by " " and then take the ones that have PWY somewhere
                        for link in entry[
                            "PATHWAY-LINKS"]:  # -> these are then searched for with ID and added to the list of found pathways -> add reactions
                            pathway = ()
                            link = re.sub(r'[()]', '', link)
                            pathway = link.split(" ")
                            for item in pathway:
                                item = re.sub(r'[|()]+', '', item)  # delete "|" and "( or )"
                                if "PWY" in item:
                                    if item not in pathways_added and item in hierarchy_dict.keys():

                                        if "SUB-PATHWAYS" not in hierarchy_dict[
                                            item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                                hierarchy_dict[item]["TYPES"]:

                                            # if PWY really is not a super-pathway -> need to take every reaction:
                                            if "REACTION-LIST" in hierarchy_dict[item]:
                                                for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                                    if reaction not in set_of_reactions and "RNX" in reaction:
                                                        set_of_reactions.add(reaction)
                                            if "PREDECESSORS" in hierarchy_dict[item]:
                                                for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                                    pattern = r'"([^"]*)"'
                                                    matches = re.findall(pattern, reaction)
                                                    for match in matches:
                                                        if match not in set_of_reactions and "RNX" in match:
                                                            set_of_reactions.add(match)

                                        pathways_added.add(item)
                                        temp_pathways.add(
                                            item)  # add the pathway to the list of the ones we search, because it is super-pathway
                                        set_of_item = set()
                                        set_of_item.add(item)
                                        types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                        types_temp_set = set()

                                        for type_element in types_in_sub_pathway:
                                            if type_element not in set_of_pathways and type_element != "Super-Pathways" and type_element not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                                types_temp_set.add(type_element)
                                                types_added.add(type_element)

                                        temp_pathways.update(types_temp_set)
                                        set_of_item.update(
                                            types_temp_set)

                                        list_new_reactions = find_sub_pathway_reactions(set_of_item, set_of_reactions,
                                                                                        hierarchy_dict,
                                                                                        hash_set_super_pathways,
                                                                                        pathways_added, types_added)
                                        for reaction in list_new_reactions:
                                            if "RNX" in reaction:
                                                set_of_reactions.add(reaction)


                                elif "RNX" in item:  # if PWY not inside
                                    if item not in set_of_reactions:
                                        set_of_reactions.add(item)

                    if "SUPER-PATHWAYS" in entry:  # strip the string and then take it and enter into the recursion
                        for super_pathway in entry["SUPER-PATHWAYS"]:
                            item = super_pathway.strip()  # do not have to check if is Pathway
                            if item not in pathways_added and item in hierarchy_dict.keys():  # is super-pathway!! -> put in recursion

                                if "SUB-PATHWAYS" not in hierarchy_dict[
                                    item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                        hierarchy_dict[item]["TYPES"]:

                                    # if PWY really is not a super-pathway -> need to take every reaction:
                                    if "REACTION-LIST" in hierarchy_dict[item]:
                                        for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                            if reaction not in set_of_reactions and "RNX" in reaction:
                                                set_of_reactions.add(reaction)
                                    if "PREDECESSORS" in hierarchy_dict[item]:
                                        for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                            pattern = r'"([^"]*)"'
                                            matches = re.findall(pattern, reaction)
                                            for match in matches:
                                                if match not in set_of_reactions and "RNX" in match:
                                                    set_of_reactions.add(match)
                                pathways_added.add(item)

                                temp_pathways.add(item)  # add the pathway to the list of the ones we search
                                set_of_item = set()
                                set_of_item.add(item)
                                types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                types_temp_set = set()

                                for type_element in types_in_sub_pathway:
                                    if type_element not in set_of_pathways and type_element != "Super-Pathways" and type_element not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                        types_temp_set.add(type_element)
                                        types_added.add(type_element)

                                temp_pathways.update(types_temp_set)
                                set_of_item.update(
                                    types_temp_set)

                                list_new_reactions = find_sub_pathway_reactions(set_of_item,
                                                                                set_of_reactions,
                                                                                hierarchy_dict,
                                                                                hash_set_super_pathways, pathways_added,
                                                                                types_added)
                                for reaction in list_new_reactions:
                                    if "RNX" in reaction:
                                        set_of_reactions.add(reaction)

                    # no else-statement needed, because it is a Super-pathway!

                    if "IN-PATHWAY" in entry:  # do the same as in super-pathways
                        for in_pathway in entry["IN-PATHWAY"]:
                            item = in_pathway.strip()
                            if in_pathway.strip() not in pathways_added:
                                if "PWY" in item:
                                    if item in hierarchy_dict.keys():

                                        if "SUB-PATHWAYS" not in hierarchy_dict[
                                            item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                                hierarchy_dict[item]["TYPES"]:

                                            # if PWY really is not a super-pathway -> need to take every reaction:
                                            if "REACTION-LIST" in hierarchy_dict[item]:
                                                for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                                    if reaction not in set_of_reactions and "RNX" in reaction:
                                                        set_of_reactions.add(reaction)
                                            if "PREDECESSORS" in hierarchy_dict[item]:
                                                for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                                    pattern = r'"([^"]*)"'
                                                    matches = re.findall(pattern, reaction)
                                                    for match in matches:
                                                        if match not in set_of_reactions and "RNX" in match:
                                                            set_of_reactions.add(match)

                                        pathways_added.add(item)
                                        temp_pathways.add(
                                            item)  # add the pathway to the list of the ones we search, because it is super-pathway
                                        set_of_item = set()
                                        set_of_item.add(item)
                                        types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                        types_temp_set = set()

                                        for type_element in types_in_sub_pathway:
                                            if type_element not in set_of_pathways and type_element != "Super-Pathways" and item not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                                types_temp_set.add(type_element)
                                                types_added.add(type_element)

                                        temp_pathways.update(types_temp_set)
                                        set_of_item.update(
                                            types_temp_set)

                                        list_new_reactions = find_sub_pathway_reactions(set_of_item, set_of_reactions,
                                                                                        hierarchy_dict,
                                                                                        hash_set_super_pathways,
                                                                                        pathways_added, types_added)
                                        for reaction in list_new_reactions:
                                            if "RNX" in reaction:
                                                set_of_reactions.add(reaction)


                                else:  # if PWY not inside
                                    if item not in set_of_reactions and "RNX" in item:
                                        set_of_reactions.add(item)

                    pathways_added.add(entry["UNIQUE-ID"][0])  # for the PWY that is in entry
                else:  # only has reactions -> take reaction list & predecessors -> new: still need to look at the pathway-links etc.

                    if "PATHWAY-LINKS" in entry:  # take the string and delete ( or ), split by " " and then take the ones that have PWY somewhere
                        for link in entry[
                            "PATHWAY-LINKS"]:  # -> these are then searched for with ID and added to the list of found pathways -> add reactions
                            pathway = ()
                            link = re.sub(r'[()]', '', link)
                            pathway = link.split(" ")
                            for item in pathway:
                                # delete "|" and "( or )"
                                item = re.sub(r'[|()]+', '', item)
                                if "PWY" in item:
                                    if item not in pathways_added and item in hierarchy_dict.keys():
                                        if "SUB-PATHWAYS" not in hierarchy_dict[
                                            item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                                hierarchy_dict[item]["TYPES"]:

                                            # if PWY really is not a super-pathway -> need to take every reaction:
                                            if "REACTION-LIST" in hierarchy_dict[item]:
                                                for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                                    if reaction not in set_of_reactions and "RNX" in reaction:
                                                        set_of_reactions.add(reaction)
                                            if "PREDECESSORS" in hierarchy_dict[item]:
                                                for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                                    pattern = r'"([^"]*)"'
                                                    matches = re.findall(pattern, reaction)
                                                    for match in matches:
                                                        if match not in set_of_reactions and "RNX" in match:
                                                            set_of_reactions.add(match)

                                        pathways_added.add(item)
                                        temp_pathways.add(
                                            item)  # add the pathway to the list of the ones we search
                                        set_of_item = set()
                                        set_of_item.add(
                                            item)  # in order to add to as pathways we search for in recursion
                                        types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                        types_temp_set = set()  # variable that decides if we should go into recursion or not -> so Unique ID is searched for in next iteration

                                        for type_element in types_in_sub_pathway:  # searches for types and puts it into search query
                                            if type_element not in set_of_pathways and type_element != "Super-Pathways" and item not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                                types_temp_set.add(
                                                    type_element)  # variable in order to save the types temporarily  and check then if already found
                                                types_added.add(type_element)

                                        temp_pathways.update(types_temp_set)
                                        set_of_item.update(
                                            types_temp_set)  # add the unique id of pathway and the unique type to the search query of TYPES
                                        list_new_reactions = find_sub_pathway_reactions(set_of_item,
                                                                                        set_of_reactions,
                                                                                        hierarchy_dict,
                                                                                        hash_set_super_pathways,
                                                                                        pathways_added, types_added)

                                        for reaction in list_new_reactions:
                                            if "RNX" in reaction:
                                                set_of_reactions.add(reaction)


                                elif "RNX" in item:  # if PWY not inside -> take reactions inside
                                    if item not in set_of_reactions:
                                        set_of_reactions.add(item)

                    if "SUPER-PATHWAYS" in entry:  # strip the string and then take it and enter into the recursion
                        for super_pathway in entry["SUPER-PATHWAYS"]:
                            item = super_pathway.strip()  # do not have to check if is Pathway
                            if item not in pathways_added and item in hierarchy_dict.keys():  # is super-pathway!! -> put in recursion

                                if "SUB-PATHWAYS" not in hierarchy_dict[
                                    item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                        hierarchy_dict[item]["TYPES"]:

                                    # if PWY really is not a super-pathway -> need to take every reaction:
                                    if "REACTION-LIST" in hierarchy_dict[item]:
                                        for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                            if reaction not in set_of_reactions and "RNX" in reaction:
                                                set_of_reactions.add(reaction)
                                    if "PREDECESSORS" in hierarchy_dict[item]:
                                        for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                            pattern = r'"([^"]*)"'
                                            matches = re.findall(pattern, reaction)
                                            for match in matches:
                                                if match not in set_of_reactions and "RNX" in match:
                                                    set_of_reactions.add(match)
                                pathways_added.add(item)
                                # super_pathways_added.add(item)
                                temp_pathways.add(item)  # add the pathway to the list of the ones we search
                                set_of_item = set()
                                set_of_item.add(item)
                                types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                types_temp_set = set()

                                for type_element in types_in_sub_pathway:
                                    if type_element not in set_of_pathways and type_element != "Super-Pathways" and item not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                        types_temp_set.add(type_element)
                                        types_added.add(type_element)

                                temp_pathways.update(types_temp_set)
                                set_of_item.update(
                                    types_temp_set)  # add the unique id of pathway and the unique type to the search query of TYPES

                                list_new_reactions = find_sub_pathway_reactions(set_of_item,
                                                                                set_of_reactions,
                                                                                hierarchy_dict,
                                                                                hash_set_super_pathways,
                                                                                pathways_added, types_added)
                                for reaction in list_new_reactions:
                                    if "RNX" in reaction:
                                        set_of_reactions.add(reaction)

                    # no else-statement needed, because it is a Super-pathway!

                    if "IN-PATHWAY" in entry:  # do the same as in super-pathways
                        for in_pathway in entry["IN-PATHWAY"]:
                            item = in_pathway.strip()
                            if item not in pathways_added and item in hierarchy_dict.keys():

                                if "SUB-PATHWAYS" not in hierarchy_dict[
                                    item] and item not in hash_set_super_pathways and "Super-Pathways" not in \
                                        hierarchy_dict[item]["TYPES"]:

                                    # if PWY really is not a super-pathway -> need to take every reaction:
                                    if "REACTION-LIST" in hierarchy_dict[item]:
                                        for reaction in hierarchy_dict[item]["REACTION-LIST"]:
                                            if reaction not in set_of_reactions and "RNX" in reaction:
                                                set_of_reactions.add(reaction)
                                    if "PREDECESSORS" in hierarchy_dict[item]:
                                        for reaction in hierarchy_dict[item]["PREDECESSORS"]:
                                            pattern = r'"([^"]*)"'
                                            matches = re.findall(pattern, reaction)
                                            for match in matches:
                                                if match not in set_of_reactions and "RNX" in match:
                                                    set_of_reactions.add(match)

                                pathways_added.add(item)
                                temp_pathways.add(item)  # add the pathway to the list of the ones we search
                                set_of_item = set()
                                set_of_item.add(item)
                                types_in_sub_pathway = hierarchy_dict[item]["TYPES"]

                                types_temp_set = set()

                                for type_element in types_in_sub_pathway:
                                    if type_element not in set_of_pathways and type_element != "Super-Pathways" and item not in types_added:  # in order to hopefully find every pathway that is missing -> in TYPES of already found
                                        types_temp_set.add(type_element)
                                        types_added.add(type_element)

                                temp_pathways.update(types_temp_set)
                                set_of_item.update(types_temp_set)

                                list_new_reactions = find_sub_pathway_reactions(set_of_item,
                                                                                set_of_reactions,
                                                                                hierarchy_dict,
                                                                                hash_set_super_pathways,
                                                                                pathways_added, types_added)
                                for reaction in list_new_reactions:
                                    if "RXN" in reaction:
                                        set_of_reactions.add(reaction)

                    pathways_added.add(entry["UNIQUE-ID"][0])  # for the entry element to be added
                    if "REACTION-LIST" in entry:
                        for item in entry["REACTION-LIST"]:
                            if item not in set_of_reactions and "PWY" not in item:
                                set_of_reactions.add(item)
                            elif "PWY" in item:
                                print(
                                    "There was an error, with saving the reactions of a sub-pathway, because PWY in there.")

                    if "PREDECESSORS" in entry:  # reaction ID is in parentheses e.g. ("RXN-5761" "RXN-5641")
                        for item in entry["PREDECESSORS"]:
                            pattern = r'"([^"]*)"'
                            matches = re.findall(pattern, item)
                            for match in matches:
                                if match not in set_of_reactions and "RNX" in match:
                                    set_of_reactions.add(match)
                                elif "PWY" in match:
                                    print(
                                        "There was an error, with saving the reactions of a sub-pathway, because PWY in there.")

        set_of_pathways.update(temp_pathways)
    return set_of_reactions


# returns a list of reactions, that has to be searched by the reaction parser


# --------------------------------------------------------------------------------------------------------------------

def find_unclassified_compounds(hierarchy_dict):
    unclassified_compounds = []
    for entry in hierarchy_dict:
        if "Unclassified-Compounds" in hierarchy_dict[entry]["TYPES"]:
            unclassified_compounds.append(hierarchy_dict[entry]["UNIQUE-ID"][0])
    return unclassified_compounds


# return list of all unclassified compounds

# --------------------------------------------------------------------------------------------------------------------
# Debugger with json file:

# with open("debug.json", 'w') as f:
# import json
# json.dump(entries_dict, f, indent=4)


# --------------------------------------------------------------------------------------------------------------------

# PARSER:


def parse_data_from_file(filename):
    start_reading = False
    dict_of_entries = {}
    line_count = 0

    #
    hash_set_super_types = set()
    #

    # count the lines and save it, in order to be able to stop the creation of empty entries because of //

    try:
        with open(filename, 'r', encoding="latin-1") as file:
            line_count = file.read().count('\n')
    except FileNotFoundError:
        print("There was an error with reading the line count of the file.")

    try:
        with open(filename, 'r', encoding="latin-1") as file:
            lines = file.readlines()
            counter_lines = 0

            for line in lines:
                # for first entry, since it does not start with //
                if start_reading == False and not line.startswith(
                        "#"):  # is always being skipped when first entry already read
                    start_reading = True

                    # Add the word separation here too, because the first line has to register the unique ID as well!

                    dict_of_lines = {}
                    words = line.split(" ")  # " " as partition symbol between category and content
                    del words[1]  # delete "-"
                    words[1:] = [''.join(words[1:]).strip()]  # join all words together after category and strip
                    unique_id = words[1]
                    dict_of_entries[unique_id] = dict_of_lines

                if line.startswith("//"):
                    counter_lines = counter_lines + 1  # to spot counter of lines correctly

                    # for stopping the creation of an empty dictionary at the end
                    if line_count == counter_lines:  # stops exactly, when line_count is reached
                        continue

                    start_reading = True
                    continue
                if start_reading:
                    if line.startswith("COMMENT") or (line.startswith('/') and line[
                        1] != '/') or line.startswith(
                        "\n"):  # Comments are not needed for data extraction and therefore not saved
                        counter_lines = counter_lines + 1
                        continue
                    words = line.split(" ")  # " " as partition symbol between category and content
                    del words[1]  # delete "-"

                    if line.startswith(
                            "PATHWAY-LINKS") and "|" not in line:  # for correct reading of pathway links -> when no | there is a problem with " "

                        words[1:] = [' '.join(words[1:]).strip()]
                    else:
                        words[1:] = [' '.join(words[1:]).strip()]  # join all words together after category and strip
                    category = words[0]
                    category_content = [words[1]]

                    if line.startswith("UNIQUE"):
                        dict_of_lines = {}
                        dict_of_entries[category_content[0]] = dict_of_lines

                    #
                    # added this part for creating a dict for all super-TYPES, e.g. for compounds, maybe even pathways
                    if line.startswith("TYPES"):
                        hash_set_super_types.update(category_content)
                    #

                    if category not in dict_of_lines:  # if is not multiple category lines
                        dict_of_lines[category] = category_content
                    else:
                        dict_of_lines[category].append(words[1])

                counter_lines = counter_lines + 1

        return dict_of_entries, hash_set_super_types

    except IndexError:
        print(str(counter_lines) + " lines were found before an error appeared.")

    except FileNotFoundError:
        print("File not found.")
    except UnicodeDecodeError:
        print("Error decoding file with US-ASCII encoding.")


# --------------------------------------------------------------------------------------------------------------------

# search for all pathways -> go through whole data and save as nested directories

def file_parser_pathways(file_path):
    try:
        # is hardcoded to be only for pathways
        entries_dict, hash_set_super_types = parse_data_from_file(file_path + "/" + "pathways.dat")

        set_of_pathways = set()
        set_of_reactions = set()
        set_of_pathways.add("Fatty-acid-biosynthesis")
        set_of_pathways.add("Lipid-Biosynthesis")

        pathways_added = set()
        types_added = set()

        list_of_fa_lipid_reaction_ids = find_sub_pathway_reactions(set_of_pathways, set_of_reactions, entries_dict,
                                                                   hash_set_super_types, pathways_added, types_added)

        return list_of_fa_lipid_reaction_ids


    except FileNotFoundError:
        print("File \"pathways.dat\" not found.")


def file_parser_reactions(reaction_ids_list, file_path, species, file_path_for_saving, standard_lipids_dict,
                          standard_lipids_key_dict, dict_of_classified_comps):
    file_path_for_react = file_path_for_saving + "/" + "PMN_fatty_acid_lipid_reactions_" + species + ".csv"

    print("species:" + species)

    # create a reaction file for every species

    try:
        entries_dict_reactions, unnecessary_hashset1 = parse_data_from_file(file_path + "/" + "reactions.dat")
        entries_dict_compounds, hash_set_super_compounds = parse_data_from_file(file_path + "/" + "compounds.dat")
        entries_dict_enzymes, unnecessary_hashset2 = parse_data_from_file(file_path + "/" + "enzrxns.dat")

        # defines which columns should be in the data file
        # old: column_categories = "ID\tReaction_type\tCommon_name\tEnzyme\tAll_Left_Compounds\tAll_Right_Compounds\tSynonyms\tConnected_Reactions"
        column_categories = ['ID', 'Reaction_type', 'Common_name', 'Enzyme', 'All_Left_Compounds',
                             'All_Right_Compounds', 'Synonyms', 'Connected_Reactions']

        # dicts of other information, so no duplicates are saved in the methods :
        dict_of_enzymes = {}
        dict_of_compounds = {}

        for reaction in reaction_ids_list:
            # take the reactions that are in id_list from entries and then extract the data
            # vars for text file creation

            unique_id = ["-"]
            types_list = ["-"]
            common_name_list = ["-"]
            enzymes = ["-"]
            left_comp_list = ["-"]
            right_comp_list = ["-"]
            synonyms = ["-"]
            connected_reactions = ["-"]

            dict_of_enzyme_line = {}  # dict to save the lines of enzymes as a dict for file_writer for all enzymes in reactions

            # in order to indicate, if there is an enzymatic reaction

            lines_dict = entries_dict_reactions[reaction]

            if "UNIQUE-ID" in lines_dict:
                unique_id = lines_dict["UNIQUE-ID"][0]

            if "TYPES" in lines_dict:
                types_list = lines_dict["TYPES"]

            if "COMMON-NAME" in lines_dict:
                common_name_list = lines_dict["COMMON-NAME"]
                cleaned_common_name_list = [re.sub('"', '', name) for name in common_name_list]
                common_name_list = cleaned_common_name_list

            if "ENZYMATIC-REACTION" in lines_dict:
                enzymes = lines_dict["ENZYMATIC-REACTION"]
                for enzyme in enzymes:
                    if enzyme not in dict_of_enzymes:
                        unique_id_enz = enzyme
                        types = ["-"]
                        common_name = ["-"]
                        cofactors = ["-"]
                        regulated_by = ["-"]
                        synonyms = ["-"]
                        if "TYPES" in entries_dict_enzymes[enzyme]:
                            types = entries_dict_enzymes[enzyme]["TYPES"]
                        if "COMMON-NAME" in entries_dict_enzymes[enzyme]:
                            common_name = entries_dict_enzymes[enzyme]["COMMON-NAME"]
                        if "COFACTORS" in entries_dict_enzymes[enzyme]:
                            cofactors = entries_dict_enzymes[enzyme]["COFACTORS"]
                        if "REGULATED-BY" in entries_dict_enzymes[enzyme]:
                            regulated_by = entries_dict_enzymes[enzyme]["REGULATED-BY"]
                        reaction = [unique_id]
                        if "SYNONYMS" in entries_dict_enzymes[enzyme]:
                            synonyms = entries_dict_enzymes[enzyme]["SYNONYMS"]

                        dict_of_enzyme_line["TYPES"] = types
                        dict_of_enzyme_line["COMMON-NAME"] = common_name
                        dict_of_enzyme_line["COFACTORS"] = cofactors

                        dict_of_enzyme_line["REGULATED-BY"] = regulated_by
                        dict_of_enzyme_line["REACTION-ID"] = reaction
                        dict_of_enzyme_line["SYNONYMS"] = synonyms

                        dict_of_enzymes[enzyme] = dict_of_enzyme_line

                    else:
                        dict_of_enzymes[enzyme]["REACTION-ID"].append(reaction)

            if "LEFT" in lines_dict:

                left_comp_list = lines_dict["LEFT"]
                final_left_comp_list = {}  # is a dict because the superclass has to be saved as well

                # first step: find the base of every molecule
                for compound in left_comp_list:

                    if compound in dict_of_compounds:  # if compound was already found
                        list_of_educts = dict_of_compounds[compound]
                        final_left_comp_list[compound] = list_of_educts

                    else:
                        if compound in hash_set_super_compounds:  # if compound is a supper-class:

                            list_of_educts = find_super_compounds(compound, entries_dict_compounds,
                                                                  hash_set_super_compounds)
                            dict_of_compounds[
                                compound] = list_of_educts  # save that in the dictionary, so it can be accessed if already went through
                            final_left_comp_list[compound] = list_of_educts

                        else:
                            dict_of_compounds[compound] = compound
                            final_left_comp_list[compound] = compound

                left_comp_list = final_left_comp_list

            if "RIGHT" in lines_dict:
                right_comp_list = lines_dict["RIGHT"]
                final_right_comp_list = {}

                for compound in right_comp_list:

                    if compound in dict_of_compounds:  # if compound was already found
                        list_of_products = dict_of_compounds[compound]
                        final_right_comp_list[compound] = list_of_products

                    else:
                        if compound in hash_set_super_compounds:  # if compound is a supper-class:

                            list_of_products = find_super_compounds(compound, entries_dict_compounds,
                                                                    hash_set_super_compounds)
                            dict_of_compounds[
                                compound] = list_of_products  # save that in the dictionary, so it can be accessed if already went through
                            final_right_comp_list[compound] = list_of_products

                        else:
                            dict_of_compounds[compound] = compound
                            final_right_comp_list[compound] = compound

                right_comp_list = final_right_comp_list

            if "SYNONYMS" in lines_dict:
                synonyms = lines_dict["SYNONYMS"]

            if "REACTION-LIST" in lines_dict:

                add_reactions = lines_dict["REACTION-LIST"]

                for reaction in add_reactions:  # for every reaction in the additional lists check for them as well
                    if reaction not in reaction_ids_list:
                        reaction_ids_list.append(reaction)

                connected_reactions = lines_dict["REACTION-LIST"]

            # combine all categories into one string and save into file

            combined_list_string = []

            combined_list_string.append("; ".join([unique_id]))
            combined_list_string.append("; ".join(types_list))
            combined_list_string.append("; ".join(common_name_list))
            combined_list_string.append("; ".join(enzymes))

            # writing educts and products as: superid = (sub_comp1; sub_comp2; ...)

            left_comp_content = "("

            if left_comp_list != ["-"]:  # added if statement, otherwise issues with getting the keys for dict
                for key in left_comp_list.keys():
                    if key == "-":
                        break
                    if key == left_comp_list[key]:  # if list has only one item
                        left_comp_content += str(key) + "; "
                        if key == list(left_comp_list.keys())[-1]:
                            left_comp_content = left_comp_content[:-2]
                    else:
                        left_comp_content += str(key) + "=("  # if list has multiple elements
                        for value in left_comp_list[key]:
                            left_comp_content += str(value) + "; "
                        left_comp_content = left_comp_content[:-2]
                        left_comp_content += "); "
                        if key == list(left_comp_list.keys())[-1]:
                            left_comp_content = left_comp_content[:-2]
                left_comp_content += ")"
            else:
                left_comp_content = "-"

            combined_list_string.append("; ".join([left_comp_content]))

            right_comp_content = "("

            if right_comp_list != ["-"]:
                for key in right_comp_list.keys():
                    if key == "-":
                        break
                    if key == right_comp_list[key]:
                        right_comp_content += str(key) + "; "
                        if key == list(right_comp_list.keys())[-1]:
                            right_comp_content = right_comp_content[:-2]
                    else:
                        right_comp_content += str(key) + "=("
                        for value in right_comp_list[key]:
                            right_comp_content += str(value) + "; "
                        right_comp_content = right_comp_content[:-2]
                        right_comp_content += "); "
                        if key == list(right_comp_list.keys())[-1]:
                            right_comp_content = right_comp_content[:-2]
                right_comp_content += ")"
            else:
                right_comp_content = "-"

            combined_list_string.append("; ".join([right_comp_content]))
            combined_list_string.append("; ".join(synonyms))
            combined_list_string.append("; ".join(connected_reactions))

            dict_to_write = dict(zip(column_categories, combined_list_string))

            file_writer(dict_to_write, file_path_for_react)

            # go through every enzyme in the beforehand saved dictionary and write every line
            # -> had to be done like this in order to add additional reaction IDs if needed and enzyme included in multiple reactions

        file_writer_enzymes(dict_of_enzymes, species, file_path_for_saving)

        file_writer_compounds(entries_dict_compounds, dict_of_compounds, species, file_path_for_saving,
                              standard_lipids_dict, standard_lipids_key_dict, dict_of_classified_comps)


    except FileNotFoundError:
        print("File \"reactions.dat\" not found.")


def file_writer_compounds(entries_of_compounds, dict_of_compounds, species, file_path_for_saving, standard_lipids_dict,
                          standard_lipids_key_dict, dict_of_classified_comps):
    # only take the compounds that where in the reactions -> done by using dict_of_compounds

    # save all noted compounds in a set, so there are no duplicates and in order to check that

    set_of_done_comps = set()

    # save the compounds in a dict with index and string, so they can be found and additional information added/modified, after the file was created >> e.g. super-compound

    dict_of_update_comps = OrderedDict()

    # dict in order to know where the lines are that we need to modify

    dict_of_comp_lines = dict()

    file_path_for_saving = file_path_for_saving + "/" + "PMN_fatty_acid_lipid_compounds_" + species + ".csv"
    counter_lines = 1  # because first line is column notation -> skip

    try:

        with open(file_path_for_saving, 'w') as file:

            column_categories = ['ID', 'Super_compound_class', 'Lipid_Category', 'Compound_type',
                                 'Main_Class', 'Common_name', 'Abbreviation',
                                 'Chemical_formula', 'DB_links',
                                 'SMILES', 'Synonyms', 'Systematic_name']

            for compound in dict_of_compounds:

                list_of_compounds = [compound]

                if type(dict_of_compounds[
                            compound]) == list:  # when is list then it means the compound is a super compound
                    list_of_compounds.extend(dict_of_compounds[compound])

                for comp in list_of_compounds:

                    lm_url = "https://www.lipidmaps.org/rest/compound/lm_id/"  # base link for REST usage
                    found_class = None

                    if comp not in set_of_done_comps:

                        if comp in entries_of_compounds.keys():

                            unique_id_comp = "-"
                            super_compound_class = "-"  # used to note down the super_class when not given in TYPES
                            lipid_category = "-"
                            types = ["-"]
                            lipid_main_class = "-"
                            common_name = ["-"]
                            abbrev = ["-"]
                            chemical_formula = ["-"]
                            DB_links = "-"
                            smiles = ["-"]
                            synonyms = ["-"]
                            systematic_name = ["-"]

                            if "UNIQUE-ID" in entries_of_compounds[comp]:
                                unique_id_comp = entries_of_compounds[comp]["UNIQUE-ID"][0]

                            if len(list_of_compounds) > 1 and comp != list_of_compounds[
                                0] and list_of_compounds[
                                0]:  # for instances of compound_class
                                super_compound_class = list_of_compounds[0]
                            elif len(list_of_compounds) > 1 and comp == list_of_compounds[
                                0]:  # for compounds without super class
                                super_compound_class = "SUPER"
                            else:
                                super_compound_class = "-"  # for compounds that are super-compounds of a class

                            if "ABBREV-NAME" in entries_of_compounds[comp]:
                                abbrev = entries_of_compounds[comp]["ABBREV-NAME"]

                            if "TYPES" in entries_of_compounds[comp]:
                                types = '; '.join(entries_of_compounds[comp]["TYPES"])

                            if "COMMON-NAME" in entries_of_compounds[comp]:
                                common_name = '; '.join(entries_of_compounds[comp]["COMMON-NAME"])
                                common_name = common_name.replace("<i>", "").replace("</i>", "").replace("<sup>",
                                                                                                         "").replace(
                                    "</sup>", "").replace("<sub>", "").replace("</sub>", "").replace("<SUB>",
                                                                                                     "").replace(
                                    "</SUB>", "").replace("<SUP>", "").replace("</SUP>", "").replace(";", "").replace("<I>", "").replace("</I>", "")
                                # for articles in the common name:

                                pattern = r'\b(a|an)\b'

                                common_name = re.sub(pattern, '', common_name)
                                common_name = re.sub(r'\s+', ' ', common_name).strip()

                            if "CHEMICAL-FORMULA" in entries_of_compounds[comp]:
                                chemical_formula = ';'.join(entries_of_compounds[comp]["CHEMICAL-FORMULA"])
                            if "DBLINKS" in entries_of_compounds[comp]:
                                DB_links = '; '.join(entries_of_compounds[comp]["DBLINKS"])

                            if "SYNONYMS" in entries_of_compounds[
                                comp]:  # add to abbrev because there are a lot of abbreviations in synonyms

                                for count, word in enumerate(entries_of_compounds[comp]["SYNONYMS"]):
                                    if ";" in word:
                                        entries_of_compounds[comp]["SYNONYMS"][count] = word.replace(";", "")

                                synonyms = '; '.join(entries_of_compounds[comp]["SYNONYMS"])
                                synonyms = synonyms.replace("<i>", "").replace("</i>", "").replace("<sup>",
                                                                                                   "").replace(
                                    "</sup>", "").replace("<sub>", "").replace("</sub>", "").replace("<SUB>",
                                                                                                     "").replace(
                                    "</SUB>", "").replace("<SUP>", "").replace("</SUP>", "").replace("<I>", "").replace("</I>", "")

                                pattern = r'\b(a|an)\b'
                                synonyms = re.sub(pattern, '', synonyms)
                                synonyms = re.sub(r'\s+', ' ', synonyms).strip()

                                if ":" in synonyms:
                                    for synonym in entries_of_compounds[comp]["SYNONYMS"]:
                                        if ":" in synonym and synonym not in abbrev:
                                            if abbrev != ["-"]:
                                                abbrev.append(synonym)
                                            else:
                                                abbrev = [synonym]


                            if "LIPID_MAPS" in DB_links and "COMP_NOT_FOUND" not in super_compound_class:
                                pattern_lm = r'(?<=LIPID_MAPS \")[A-Z0-9]+(?=\")'

                                match = re.findall(pattern_lm, DB_links)
                                lm_url += match[0] + "/all/" + "text"
                                r = requests.get(lm_url)

                                if r.status_code == 200:
                                    print("LM REST GET request was successful! ")

                                else:
                                    print("There was an issue with the LM REST GET: Error with code ",
                                          r.status_code)

                                lm_dict = json.loads(r.text)

                                if lm_dict != []:
                                    if "abbrev" in lm_dict.keys():
                                        if abbrev == ["-"]:
                                            abbrev = [lm_dict["abbrev"]]
                                        else:
                                            abbrev.append(lm_dict["abbrev"])

                                    # taking abbrev_chains:
                                    if "abbrev_chains" in lm_dict.keys():
                                        if abbrev == ["-"]:
                                            abbrev = [lm_dict["abbrev_chains"]]
                                        else:
                                            abbrev.append(lm_dict["abbrev_chains"])

                                    if "synonyms" in lm_dict.keys():
                                        if synonyms != ["-"]:
                                            synonyms.replace("\"", "")

                                        if ":" in synonyms:  # when there is an abbreviation in the synonyms category
                                            if abbrev == ["-"]:
                                                abbrev = [lm_dict["synonyms"]]
                                            else:
                                                abbrev.append(lm_dict["synonyms"])

                                        else:
                                            if synonyms != "-" and synonyms != ["-"]:
                                                synonyms = synonyms + "; " + lm_dict["synonyms"]
                                            else:
                                                synonyms = lm_dict["synonyms"]

                            if "SMILES" in entries_of_compounds[comp]:
                                smiles = '; '.join(entries_of_compounds[comp]["SMILES"])

                            if "SYSTEMATIC-NAME" in entries_of_compounds[comp]:
                                systematic_name = '; '.join(entries_of_compounds[comp]["SYSTEMATIC-NAME"])

                            if unique_id_comp in dict_of_classified_comps.keys() and found_class is None:
                                # set found_class on True, so we can distinguish if there is already a class found in the dict or not (when value > is list)
                                lipid_category = dict_of_classified_comps[unique_id_comp][0]
                                lipid_main_class = dict_of_classified_comps[unique_id_comp][1]
                                found_class = True
                            if types != "-":
                                for type_elem in types.split("; "):
                                    if type_elem in dict_of_classified_comps.keys() and found_class is None and type_elem != "Cofactors" and type_elem != "Compounds":
                                        lipid_category = dict_of_classified_comps[type_elem][0]
                                        lipid_main_class = dict_of_classified_comps[type_elem][1]
                                        found_class = True
                                        break

                            if super_compound_class in dict_of_classified_comps.keys() and found_class is None and type_elem != "Cofactors" and type_elem != "Compounds":
                                lipid_category = dict_of_classified_comps[super_compound_class][0]
                                lipid_main_class = dict_of_classified_comps[super_compound_class][1]
                                found_class = True

                            for item in abbrev:  # try to classify by the abbreviation
                                if abbrev != ["-"] and abbrev != "-":
                                    new_abbrev = ""
                                    for element in item.split("-"):
                                        stripped = element.strip()
                                        # in order to add manually CoA and FA/WE, otherwise can't be found and classified correctly
                                        if "CoA" in stripped and not "ACoA" in stripped:
                                            stripped = "ACoA"
                                            new_abbrev = item.strip().replace("CoA", "ACoA")
                                        if "COA" in stripped and not "ACoA" in stripped:
                                            stripped = "ACoA"
                                            new_abbrev = item.strip().replace("COA", "ACoA")
                                        if "FA" in stripped and not "NEFA" in stripped:
                                            stripped = "NEFA"
                                            new_abbrev = item.strip().replace("FA", "NEFA")
                                        if "WE" in stripped and not "NEFA" in stripped:     # here second clause of if not necessary, but still added.
                                            stripped = "NEFA"
                                            new_abbrev = item.strip().replace("WE", "NEFA")
                                        if "SPBP" in stripped and not "IPC" in stripped:
                                            stripped = "IPC"
                                            new_abbrev = item.strip().replace("SPBP", "IPC")
                                        if "SPB" in stripped and not "LCB" in stripped:
                                            stripped = "LCB"
                                            new_abbrev = item.strip().replace("SPB", "LCB")
                                        if "PtdGro" in stripped and not "PGP" in stripped:
                                            stripped = "PGP"
                                            new_abbrev = item.strip().replace("PtdGro", "PGP")

                                        if stripped in standard_lipids_dict["Abbreviation"]:
                                            found_class = stripped
                                            found_class.replace(";", "")
                                            if new_abbrev != "":
                                                if abbrev != ["-"]:
                                                    abbrev.append(new_abbrev)
                                                else:
                                                    abbrev = [new_abbrev]
                                            new_abbrev = ""
                                            break

                                        if stripped in standard_lipids_dict["AltAbb"]:
                                            stripped.replace(";", "")
                                            index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                            abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                            found_class = [abbrev_stand, stripped]
                                            if new_abbrev != "":
                                                if abbrev != ["-"]:
                                                    abbrev.append(new_abbrev)
                                                else:
                                                    abbrev = [new_abbrev]
                                            new_abbrev = ""
                                            break

                                    if found_class is None:
                                        new_abbrev = ""
                                        for element in item.split(" "):
                                            stripped = element.strip()
                                            if "CoA" in stripped and not "ACoA" in stripped:
                                                stripped = "ACoA"
                                                new_abbrev = item.replace("CoA", "ACoA")
                                            if "COA" in stripped and not "ACoA" in stripped:
                                                stripped = "ACoA"
                                                new_abbrev = item.strip().replace("COA", "ACoA")
                                            if "FA" in stripped and not "NEFA" in stripped:
                                                stripped = "NEFA"
                                                new_abbrev = item.strip().replace("FA", "NEFA")
                                            if "WE" in stripped and not "NEFA" in stripped:
                                                stripped = "NEFA"
                                                new_abbrev = item.strip().replace("WE", "")
                                            if "SPBP" in stripped and not "IPC" in stripped:
                                                stripped = "IPC"
                                                new_abbrev = item.strip().replace("SPBP", "IPC")
                                            if "SPB" in stripped and not "LCB" in stripped:
                                                stripped = "LCB"
                                                new_abbrev = item.strip().replace("SPB", "LCB")
                                            if "PtdGro" in stripped and not "PGP" in stripped:
                                                stripped = "PGP"
                                                new_abbrev = item.strip().replace("PtdGro", "PGP")

                                            if stripped in standard_lipids_dict["Abbreviation"]:
                                                found_class = stripped
                                                found_class.replace(";", "")
                                                if new_abbrev != "":
                                                    if abbrev != ["-"]:
                                                        abbrev.append(new_abbrev)
                                                    else:
                                                        abbrev = [new_abbrev]
                                                new_abbrev = ""
                                                break
                                            if stripped in standard_lipids_dict["AltAbb"]:
                                                stripped.replace(";", "")
                                                index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                found_class = [abbrev_stand, stripped]
                                                if new_abbrev != "":
                                                    if abbrev != ["-"]:
                                                        abbrev.append(new_abbrev)
                                                    else:
                                                        abbrev = [new_abbrev]
                                                new_abbrev = ""
                                                break

                                    if found_class is None:
                                        new_abbrev = ""
                                        for element in item.split("("):
                                            stripped = element.strip()
                                            if "CoA" in stripped and not "ACoA" in stripped:
                                                stripped = "ACoA"
                                                new_abbrev = item.replace("CoA", "ACoA")
                                            if "COA" in stripped and not "ACoA" in stripped:
                                                stripped = "ACoA"
                                                new_abbrev = item.strip().replace("COA", "ACoA")
                                            if "FA" in stripped and not "NEFA" in stripped:
                                                stripped = "NEFA"
                                                new_abbrev = item.strip().replace("FA", "NEFA")
                                            if "WE" in stripped and not "NEFA" in stripped:
                                                stripped = "NEFA"
                                                new_abbrev = item.strip().replace("WE", "NEFA")
                                            if "SPBP" in stripped and not "IPC" in stripped:
                                                stripped = "IPC"
                                                new_abbrev = item.strip().replace("SPBP", "IPC")
                                            if "SPB" in stripped and not "LCB" in stripped:
                                                stripped = "LCB"
                                                new_abbrev = item.strip().replace("SPB", "LCB")
                                            if "PtdGro" in stripped and not "PGP" in stripped:
                                                stripped = "PGP"
                                                new_abbrev = item.strip().replace("PtdGro", "PGP")

                                            if stripped in standard_lipids_dict["Abbreviation"]:
                                                found_class = stripped
                                                found_class.replace(";", "")
                                                if new_abbrev != "":
                                                    if abbrev != ["-"]:
                                                        abbrev.append(new_abbrev)
                                                    else:
                                                        abbrev = [new_abbrev]
                                                new_abbrev = ""
                                                break
                                            if stripped in standard_lipids_dict["AltAbb"]:
                                                stripped.replace(";", "")
                                                index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                found_class = [abbrev_stand, stripped]
                                                if new_abbrev != "":
                                                    if abbrev != ["-"]:
                                                        abbrev.append(new_abbrev)
                                                    else:
                                                        abbrev = [new_abbrev]
                                                new_abbrev = ""
                                                break

                            if found_class is None:
                                if synonyms != "-" and synonyms != ["-"] and "Gibberellins" not in types and type(synonyms) == str:
                                    for item in synonyms.split("; "):
                                        for element in item.split("-"):
                                            stripped = element.strip()
                                            if "CoA" in stripped and not "ACoA" in stripped:
                                                stripped = "ACoA"
                                            if "COA" in stripped and not "ACoA" in stripped:
                                                stripped = "ACoA"
                                            if "FA" in stripped and not "NEFA" in stripped:
                                                stripped = "NEFA"
                                            if "WE" in stripped:
                                                stripped = "NEFA"
                                            if "SPBP" in stripped:
                                                stripped = "IPC"
                                            if "SPB" in stripped:
                                                stripped = "LCB"
                                            if "PtdGro" in stripped:
                                                stripped = "PGP"

                                            if stripped in standard_lipids_dict["Abbreviation"]:
                                                found_class = stripped
                                                found_class.replace(";", "")
                                                break
                                            if stripped in standard_lipids_dict["AltAbb"]:
                                                stripped.replace(";", "")
                                                index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                found_class = [abbrev_stand, stripped]
                                                break
                                        if found_class is None:
                                            for element in item.split(" "):
                                                stripped = element.strip()
                                                if stripped in standard_lipids_dict["Abbreviation"]:
                                                    found_class = stripped
                                                    found_class.replace(";", "")
                                                    break
                                                if stripped in standard_lipids_dict["AltAbb"]:
                                                    stripped.replace(";", "")
                                                    index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                    abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                    found_class = [abbrev_stand, stripped]
                                                    break
                                        if found_class is None:
                                            for element in item.split("("):
                                                stripped = element.strip()
                                                if stripped in standard_lipids_dict["Abbreviation"]:
                                                    found_class = stripped
                                                    found_class.replace(";", "")
                                                    break
                                                if stripped in standard_lipids_dict["AltAbb"]:
                                                    stripped.replace(";", "")
                                                    index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                    abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                    found_class = [abbrev_stand, stripped]
                                                    break

                            if found_class is None:
                                if common_name != "-" and common_name != ["-"]:
                                    for item in common_name.split("; "):
                                        for element in item.split("-"):
                                            stripped = element.strip()
                                            if "CoA" in stripped:
                                                stripped = "ACoA"
                                            if "COA" in stripped:
                                                stripped = "ACoA"
                                            if "FA" in stripped and not "NEFA" in stripped:
                                                stripped = "NEFA"
                                            if "WE" in stripped:
                                                stripped = "NEFA"
                                            if "SPBP" in stripped:
                                                stripped = "IPC"
                                            if "SPB" in stripped:
                                                stripped = "LCB"
                                            if "PtdGro" in stripped:
                                                stripped = "PGP"
                                            if stripped in standard_lipids_dict["Abbreviation"]:
                                                found_class = stripped
                                                found_class.replace(";", "")
                                                break
                                            if stripped in standard_lipids_dict["AltAbb"]:
                                                stripped.replace(";", "")
                                                index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                found_class = [abbrev_stand, stripped]
                                                break
                                        if found_class is None:
                                            for element in item.split(" "):
                                                stripped = element.strip()
                                                if stripped in standard_lipids_dict["Abbreviation"]:
                                                    found_class = stripped
                                                    found_class.replace(";", "")
                                                    break
                                                if stripped in standard_lipids_dict["AltAbb"]:
                                                    stripped.replace(";", "")
                                                    index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                    abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                    found_class = [abbrev_stand, stripped]
                                                    break
                                        if found_class is None:
                                            for element in item.split("("):
                                                stripped = element.strip()
                                                if stripped in standard_lipids_dict["Abbreviation"]:
                                                    found_class = stripped
                                                    found_class.replace(";", "")
                                                    break
                                                if stripped in standard_lipids_dict["AltAbb"]:
                                                    stripped.replace(";", "")
                                                    index_of_abbr = standard_lipids_dict["AltAbb"].index(stripped)
                                                    abbrev_stand = standard_lipids_dict["Abbreviation"][index_of_abbr]
                                                    found_class = [abbrev_stand, stripped]
                                                    break

                            if found_class is None:

                                client = docker.from_env()
                                string = ""

                                for name in common_name.split("; "):
                                    output = {}
                                    try:
                                        if name != "-":
                                            test = client.containers.run("lipidlibrarian", [name], stdout=True,
                                                                         stderr=True)
                                            string = test.decode("utf-8")

                                    except:
                                        print(
                                            "LipidLibrarian failed to parse: " + unique_id_comp + " with common name: " + common_name)

                                    if 'WARNING:root:LipidQuery: The input could not be parsed.\n' not in string and "{" in string:
                                        try:
                                            json_start = string.find("{")

                                            # Extract the JSON part from that position onwards
                                            if json_start != -1:
                                                json_string = string[json_start:]

                                                pattern = r'"lipid_name":\s*"([^"]+)"'
                                                match = re.search(pattern, json_string)
                                                print(json_string + match.group(1))
                                                output = match.group(1)

                                                if output != None:
                                                    new_abbrev_elem = output
                                                    if "CoA" in output and not "ACoA" in output:
                                                        found_class = "ACoA"
                                                        new_abbrev_elem = output.strip().replace("CoA", "ACoA")
                                                    if "COA" in output and not "ACoA" in output:
                                                        found_class = "ACoA"
                                                        new_abbrev_elem = output.strip().replace("COA", "ACoA")
                                                    if "FA" in output and not "NEFA" in output:
                                                        found_class = "NEFA"
                                                        new_abbrev_elem = output.strip().replace("FA", "NEFA")
                                                    if "WE" in output:
                                                        found_class = "NEFA"
                                                        new_abbrev_elem = output.strip().replace("WE", "NEFA")
                                                    if "SPBP" in output:
                                                        found_class = "IPC"
                                                        new_abbrev_elem = output.strip().replace("SPBP", "IPC")
                                                    if "SPB" in output:
                                                        found_class = "LCB"
                                                        new_abbrev_elem = output.strip().replace("SPB", "LCB")
                                                    if "PtdGro" in output:
                                                        found_class = "PGP"
                                                        new_abbrev_elem = output.strip().replace("PtdGro", "PGP")

                                                    stripped = new_abbrev_elem.split(" ")
                                                    for item in stripped:
                                                        if item in standard_lipids_dict["Abbreviation"]:
                                                            found_class.replace(";", "")
                                                            found_class = item
                                                            if abbrev != ["-"]:
                                                                abbrev.append(new_abbrev_elem)
                                                            else:
                                                                abbrev = [new_abbrev_elem]
                                                            break
                                                        if item in standard_lipids_dict["AltAbb"]:
                                                            item.replace(";", "")
                                                            index_of_abbr = standard_lipids_dict["AltAbb"].index(
                                                                item)
                                                            abbrev_stand = standard_lipids_dict["Abbreviation"][
                                                                index_of_abbr]
                                                            found_class = [abbrev_stand, item]
                                                            if abbrev != ["-"]:
                                                                abbrev.append(new_abbrev_elem.replace(item,abbrev_stand))
                                                            else:
                                                                abbrev = [new_abbrev_elem.replace(item,abbrev_stand)]
                                                            break

                                                print(
                                                    "LipidLibrarian parsed: " + unique_id_comp + " with common_name: " + common_name)

                                            else:
                                                print(string)

                                        except:
                                            print(
                                                "LipidLibrarian failed to parse: " + unique_id_comp + " with common name: " + common_name)

                                    elif 'INFO:app_log:LipidLynxX Log started ...\nWARNING:root:name and ontology keys are both missing\n' in string and "{" not in string:  # for when the LipidLib doesn't work -> Goslin
                                        lipid_parser = LipidParser()
                                        pot_lipid = lipid_parser.parse(name)  # potential lipid
                                        output = pot_lipid.get_lipid_string()

                                        if output != None:
                                            new_abbrev_elem = output
                                            if "CoA" in output and not "ACoA" in output:
                                                found_class = "ACoA"
                                                new_abbrev_elem = output.strip().replace("CoA", "ACoA")
                                            if "COA" in output and not "ACoA" in output:
                                                found_class = "ACoA"
                                                new_abbrev_elem = output.strip().replace("COA", "ACoA")
                                            if "FA" in output and not "NEFA" in output:
                                                found_class = "NEFA"
                                                new_abbrev_elem = output.strip().replace("FA", "NEFA")
                                            if "WE" in output:
                                                found_class = "NEFA"
                                                new_abbrev_elem = output.strip().replace("WE", "NEFA")
                                            if "SPBP" in output:
                                                found_class = "IPC"
                                                new_abbrev_elem = output.strip().replace("SPBP", "IPC")
                                            if "SPB" in output:
                                                found_class = "LCB"
                                                new_abbrev_elem = output.strip().replace("SPB", "LCB")
                                            if "PtdGro" in output:
                                                found_class = "PGP"
                                                new_abbrev_elem = output.strip().replace("PtdGro", "PGP")
                                            if abbrev != ["-"]:
                                                abbrev.append(new_abbrev_elem)
                                            else:
                                                abbrev = [new_abbrev_elem]

                                        print("Goslin parsed: " + unique_id_comp + " with common_name: " + common_name)


                                    elif "Failed to decode Lipid:" in string:
                                        print(
                                            "LipidLibrarian failed to decode: " + unique_id_comp + " with common name: " + common_name)

                                    else:
                                        print(
                                            "LipidLibrarian failed to parse: " + unique_id_comp + " with common name: " + common_name)

                            # classification of FA and CoA at the end of the classification method, because I would rather have it being classified by Goslin first than by pattern matching, this only as last resort
                            if found_class is None:
                                # add when CoA or COA in ID in Types or in Common names >> had to be done manually, because it is often not found by any tool

                                if "CoA" in unique_id_comp or "COA" in unique_id_comp or "CO-A" in unique_id_comp or "CoA" in types or "COA" in types:
                                    if abbrev != ["-"]:
                                        abbrev.append("ACoA-miss")
                                    else:
                                        abbrev = ["ACoA-miss"]

                                    for element in abbrev:
                                        stripped = element.strip().split("-")
                                        for item in stripped:
                                            if item in standard_lipids_dict["Abbreviation"] and item != "-":
                                                found_class = stripped
                                                break
                                elif "Fatty-Acid" in unique_id_comp or "fatty-acid" in unique_id_comp or "Fatty-Acid" in types or "fatty-acid" in types or "FATTY-ACID" in unique_id_comp or "FATTY-ACID" in types:
                                    if abbrev != ["-"]:
                                        abbrev.append("NEFA-miss")
                                    else:
                                        abbrev = ["NEFA-miss"]

                                    for element in abbrev:
                                        stripped = element.strip().split("-")
                                        for item in stripped:
                                            if item in standard_lipids_dict["Abbreviation"] and item != "-":
                                                found_class = stripped
                                                break
                                # only part of the words for lyso taken, so LPC can also be matched
                                elif "Lysopho" in unique_id_comp or "LYSOPHO" in unique_id_comp or "Lysopho" in types:

                                    if abbrev != ["-"]:
                                        abbrev.append("LPC-miss")
                                    else:
                                        abbrev = ["LPC-miss"]

                                    for element in abbrev:
                                        stripped = element.strip().split("-")
                                        for item in stripped:
                                            if item in standard_lipids_dict["Abbreviation"] and item != "-":
                                                found_class = stripped
                                                break

                            # if "" in unique_id_comp or "COA" in unique_id_comp or "CO-A" in unique_id_comp or "CoA" in types or "COA" in types:

                            if abbrev != ["-"] and abbrev[-1:] == ";":
                                abbrev = abbrev[1:]

                            abbrev = '; '.join(abbrev)

                            pattern = r'\b(a|an)\b'
                            abbrev = re.sub(pattern, '', abbrev)
                            abbrev = re.sub(r'\s+', ' ', abbrev).strip()

                            DB_links = DB_links.replace("\"", "")

                            if found_class is not None and found_class != True:  # when there was some kind of classification found in the process.
                                if type(found_class) != list:
                                    lipid_category = standard_lipids_key_dict[found_class]['Lipid category']
                                    lipid_main_class = standard_lipids_key_dict[found_class]["Lipid class"]
                                    print(unique_id_comp+ " + " + common_name + " was parsed.")
                                else:
                                    lipid_category = standard_lipids_key_dict[found_class[0]]['Lipid category']
                                    lipid_main_class = standard_lipids_key_dict[found_class[0]]["Lipid class"]
                                    print(unique_id_comp + " + " + common_name + " was parsed.")


                            content = [unique_id_comp, super_compound_class, lipid_category, types, lipid_main_class,
                                       common_name, abbrev,
                                       chemical_formula, DB_links, smiles, synonyms,
                                       systematic_name]

                            dict_to_write = dict(
                                zip(column_categories, content))  # create dictionary for writing the file in csv format
                            file_writer(dict_to_write, file_path_for_saving)

                            set_of_done_comps.add(comp)  # add compound that was noted

                            # add to the set of classified
                            # unique-id, types and super_compound_class as ID of a dict, and then the category and main class as tuple

                            # add here that the dict is updated when there is new information about a compound regarding the classification


                            if unique_id_comp in dict_of_classified_comps.keys() and dict_of_classified_comps[unique_id_comp] == ('-','-') and (lipid_category != '-' and lipid_main_class != '-'):

                                dict_of_classified_comps.pop(unique_id_comp)
                                dict_of_classified_comps.update({unique_id_comp: (lipid_category, lipid_main_class)})  # here for the unique id
                            else:
                                dict_of_classified_comps.update(
                                    {unique_id_comp: (lipid_category, lipid_main_class)})  # here for the unique id

                            if types != "-":
                                for type_elem in types.split("; "):
                                    if type_elem in dict_of_classified_comps.keys() and dict_of_classified_comps[type_elem] == ("-","-") and (lipid_category != "-" and lipid_main_class != "-"):
                                        if "Compounds" not in type_elem and "Cofactors" not in type_elem:    # otherwise there is a false classification of lipids sometimes
                                            dict_of_classified_comps.pop(type_elem)
                                            dict_of_classified_comps.update({type_elem: (lipid_category, lipid_main_class)})  # here for the types
                                    else:
                                        if "Compounds" not in type_elem and "Cofactors" not in type_elem:    # otherwise there is a false classification of lipids sometimes
                                            dict_of_classified_comps.update({type_elem: (lipid_category, lipid_main_class)})

                            if super_compound_class in dict_of_classified_comps.keys() and dict_of_classified_comps[super_compound_class] == ("-","-") and (lipid_category != "-" and lipid_main_class != "-"):
                                if super_compound_class != "COMP_NOT_FOUND" and super_compound_class != "-" and super_compound_class != "SUPER":  # here for super_compound_class
                                    dict_of_classified_comps.pop(super_compound_class)
                                    dict_of_classified_comps.update({super_compound_class: (lipid_category, lipid_main_class)})
                            else:
                                if super_compound_class != "COMP_NOT_FOUND" and super_compound_class != "-" and super_compound_class != "SUPER":  # here for super_compound_class
                                    dict_of_classified_comps.update(
                                        {super_compound_class: (lipid_category, lipid_main_class)})

                        else:
                            unique_id_comp = comp  # because no ID in file
                            super_compound_class = "COMP_NOT_FOUND"  # used to note down the super_class when not given in TYPES
                            lipid_category = "-"
                            types = "-"
                            lipid_main_class = "-"
                            common_name = "-"
                            abbrev = "-"
                            chemical_formula = "-"
                            DB_links = "-"
                            smiles = "-"
                            synonyms = "-"
                            systematic_name = "-"
                            content = [unique_id_comp, super_compound_class, lipid_category, types, lipid_main_class,
                                       common_name, abbrev,
                                       chemical_formula, DB_links, smiles, synonyms,
                                       systematic_name]
                            dict_to_write = dict(
                                zip(column_categories, content))
                            file_writer(dict_to_write, file_path_for_saving)

                            set_of_done_comps.add(comp)
                            dict_of_classified_comps.update(
                                {unique_id_comp: (lipid_category, lipid_main_class)})  # here for the unique id

                        dict_of_comp_lines[comp] = counter_lines  # save the line where the compound is saved
                        counter_lines += 1

                    else:  # when in the set already noted, when "-" just skip

                        if super_compound_class != "-":
                            if comp not in dict_of_update_comps:  # when compound has until now the first element to add in super-comp-class
                                dict_of_update_comps[comp] = super_compound_class
                            else:  # if compound has multiple elements to add at least second additional super-comp-class
                                prev_string = dict_of_update_comps[comp]
                                new_string = prev_string + "; " + super_compound_class
                                dict_of_update_comps[comp] = new_string

    except FileNotFoundError:
        print("File \"compounds.dat\" not found.")

    update_comp_file(file_path_for_saving, dict_of_update_comps, dict_of_comp_lines, dict_of_classified_comps)

    """
    Description of strategy: 
    
    Run at least two times through the complete file again: 
    
    First run: add all the classifications by asserting the abbreviations (+ altabbr) to a main and sub class or with LipidLibrarian
    
    Second run: add all the classes to the other sub-compounds in the compound class. 
    
    In the end the file is being read and written again 3 times. ( first time for update_comp_file -> adjust hierarchy) 
    
    """



def file_writer_enzymes(dict_of_enzymes, species, file_path_for_saving):
    # enzymes is a dict of enzymes! It will be iterated through and then the file will be produced

    file_path_for_saving = file_path_for_saving + "/" + "PMN_fatty_acid_lipid_enzymes_" + species + ".csv"

    try:
        with open(file_path_for_saving, 'w') as file:
            # old: column_categories = "ID\tTypes\tEnzyme_name\tCommon_name\tCofactors\tRegulated-By\tReaction\tSynonyms\n"
            column_categories = ['ID', 'Types', 'Enzyme_name', 'Common_name', 'Cofactors', 'Regulated-By', 'Reaction',
                                 'Synonyms']

            for enzyme_entry in dict_of_enzymes:
                unique_id = str(enzyme_entry)
                types = '; '.join(dict_of_enzymes[enzyme_entry]["TYPES"])
                enzyme_name = '; '.join(dict_of_enzymes[enzyme_entry]["COMMON-NAME"])
                cofactors = '; '.join(dict_of_enzymes[enzyme_entry]["COFACTORS"])
                regulated_by = '; '.join(dict_of_enzymes[enzyme_entry]["REGULATED-BY"])
                reaction = '; '.join(dict_of_enzymes[enzyme_entry]["REACTION-ID"])
                synonyms = '; '.join(dict_of_enzymes[enzyme_entry]["SYNONYMS"])
                # old: content = f"{unique_id}\t{types}\t{enzyme_name}\t{cofactors}\t{regulated_by}\t{reaction}\t{synonyms}\n"
                content = [unique_id, types, enzyme_name, cofactors, regulated_by, reaction, synonyms]
                dict_to_write = dict(zip(column_categories, content))
                file_writer(dict_to_write, file_path_for_saving)

    except FileNotFoundError:
        print("File \"enzrxns.dat\" not found.")

    return


# --------------------------------------------------------------------------------------------------------------------


"""
Here is the main that is used to run this script! 
"""


def __main__():
    root_directory = '/home/julia/Bachelor_Arbeit/PMN/data'
    paths_species = species_reader.find_second_data_directories(root_directory)

    standard_lipids_dict, standard_lipids_key_dict = read_standard_lipids()

    dict_of_classified_comps = dict()  # have unique-id, types and super_compound_class as ID of a dict, and then the category and main class as tuple

    # create main directory
    directory_path = "output"
    try:
        os.makedirs(directory_path)
        print("Directory created successfully:", directory_path)
    except FileExistsError:
        print("Main directory already exists:", directory_path)
    except Exception as e:
        print("Error occurred while creating main directory:", e)

    ara_path = ""
    reaction_list = ()

    for one_species_path in paths_species:

        directory_path = "output"

        species = ""

        if "plantcyc" not in one_species_path:
            # search for species name in species.dat -> Species name is in common names
            try:
                with open(one_species_path + "/species.dat") as file:

                    for line in file:
                        if line.startswith("COMMON-NAME"):
                            species = line.strip().split(" ")
                            del species[0:2]  # only take the species name
                            species = "_".join(species)
                            break

            except FileNotFoundError as e:
                print("File \"species.dat\" not found.")

        # add special cases for plantcyc and lercyc, because they can not read the species names as usual

        if "plantcyc" in one_species_path:
            directory_path += "/Plantcyc"
            species = "plantcyc"
        else:
            directory_path += "/" + species

        if "lercyc" in one_species_path:
            directory_path += "_ler"

        # create directories for saving species

        try:
            os.makedirs(directory_path)
            print("Directory created successfully:", directory_path)
        except FileExistsError:
            print("Directory already exists:", directory_path)
        except Exception as e:
            print("Error occurred while creating directory:", e)

        reaction_ids = list(file_parser_pathways(one_species_path))  # searches for pathway entries we need!
        file_parser_reactions(reaction_ids, one_species_path, species, directory_path, standard_lipids_dict,
                              standard_lipids_key_dict, dict_of_classified_comps)


"""
        if '/aracyc.tar' in one_species_path:
            ara_path = one_species_path
            reaction_ids = list(file_parser_pathways(ara_path))  # searches for pathway entries we need!

            filename = ara_path
            
"""

__main__()

# only for testing:
"""
def main2():
    file1 = "/home/julia/Desktop/Bachelorarbeit_Stuff/Skripte_p-LINEX/pythonProject/pmn_creation/output/Arabidopsis_thaliana/PMN_fatty_acid_lipid_reactions_Arabidopsis_thaliana.csv"  # Replace with the actual filename
    file2 = "/home/julia/Desktop/Bachelorarbeit_Stuff/Skripte_p-LINEX/pythonProject/pmn_creation/output/Plantcyc/PMN_fatty_acid_lipid_reactions_plantcyc.csv"  # Replace with the actual filename
    percentage = compare_first_words(file1, file2)
    print(f"The percentage of common reaction_IDs (excluding the first line) in both files is: {percentage:.2f}%")
"""