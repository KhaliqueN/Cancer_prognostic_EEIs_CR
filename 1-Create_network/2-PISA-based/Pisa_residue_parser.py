#!/usr/bin/python3
"""
Code for parsing pisa xml files and saving data in a csv file.

  How to use
  ----------

First you need to have the interfacetable.xml and each hydrogenbond.xml and saltbridge.xml files
given by PDBePISA web server.

Then you can run the script with the following command :

    python Pisa_xml_parser.py interfacetable.xml

  Author
  ------
    Khalique Newaz

"""

import argparse
import pandas as pd
import re
import os.path


def residue_parser(xml_file):
    """
    The function to parse the interaction xml files.
    It haven't been tested on disulfide and covalent bonds.

    Parameters
    ----------
    xml_file : string
        the name of the xml file

    Returns
    -------
    list
    """
    value = []
    lst = []

    if os.path.exists(xml_file):
        with open(xml_file, "r") as f_xml:
            for line in f_xml :
                if line.startswith("<STRUCTURE>"):
                    value.append(re.split('<|>',line)[2])#print(re.split('<|>',line)[2])#
                elif line.startswith("<SOLVENTACCESSIBLEAREA>"):
                    value.append(float(re.split('<|>',line)[2]))
                elif line.startswith("<BURIEDSURFACEAREA>"):
                    value.append(float(re.split('<|>',line)[2]))
                elif line.startswith("<SOLVATIONENERGY>"):
                    value.append(float(re.split('<|>',line)[2]))
                    lst.append(value)
                    value = []
    else:
        print("No ",xml_file, "found")

    return(lst)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("xml_file", help="the pisa interface table xml file", type=str)
    PARSER.add_argument("store_dir", help="the folder to store results", type=str)
    PARSER.add_argument("pdbid", help="the pdb ID", type=str)

    ARGS = PARSER.parse_args()

    XML_FILE = ARGS.xml_file
    STR_DIR = ARGS.store_dir
    PDB_ID = ARGS.pdbid

    RP = residue_parser(XML_FILE)
    DF = pd.DataFrame(RP, columns = ['Chain', 'Solvent_accessible_area', 'Buried_area','Solvation_energy'])
    DF.to_csv(STR_DIR+"/"+PDB_ID+".csv")
