"""
Author: Jared Taylor
Date: 4/26/2021
Filename: tf_stringdb_analyzer.py

This program takes an input directory of csv files that contain lists of proteins
and compares the list of proteins to a list of transcription factors (tfs). The proteins
that match up with the list of tfs are then analyzed for interactions using the
Stringdb API. It will ouput a png file with the network map of interactions between the proteins
and the interaction scores for the interactions. It will complete this for every csv
file in the input directory. If a file doesn't have at least two tfs in it, it will be added
to an output txt file (files_without_networks.txt). The program will also output the network
that it is saving and wether it had interactions or not.

***** Please read the README.txt file for all instructions on how to use this program *****
"""


import pandas as pd
from pandas.core.common import flatten
import os
import requests
from time import sleep
from PIL import Image, ImageFont, ImageDraw

#requests a network image, from stringdb, of the known interactions between the proteins in gene_list
def get_network(input_file, gene_list):

    string_api_url = "https://version-11-0b.string-db.org/api"
    output_format = "image"
    method = "network"

    #construct url
    request_url = "/".join([string_api_url, output_format, method])

    #set parameters for api requests
    params = {

        "identifiers" : "%0d".join(gene_list), # protein list
        "species" : 9606, # species NCBI identifier for homo sapien
        "network_flavor": "evidence", # show evidence links
        "caller_identity" : "jared taylor" # your indentifier

    }

    #call stringdb
    response = requests.post(request_url, data=params)

    #Save the network to file
    new_file = input_file.split('.csv')[0]
    #new_file = input_file.split('./input2/')[1]

    file_name = "../output/%s_network.png" % new_file
    print("saving network for %s" % new_file)

    with open(file_name, 'wb') as fh:
        fh.write(response.content)

#request the interaction evidence scores, from stringdb, of the known interactions between the proteins in gene_list
#opens the network png file created by get_network, and writes this information on the png to prevent making two files_without_networks
#for each input file
def get_interactions(input_file, gene_list):

    string_api_url = "https://version-11-0b.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    ## construct url
    request_url = "/".join([string_api_url, output_format, method])

    ##
    ## Set parameters
    ##

    params = {

        "identifiers" : "%0d".join(gene_list), # your protein
        "species" : 9606, # species NCBI identifier
        "caller_identity" : "jared taylor" # your app name

    }

    ## call stringdb
    response = requests.post(request_url, data=params)

    x = int(0)
    print_indicator = 0

    for line in response.text.strip().split("\n"):

        l = line.strip().split("\t")

        if l == ['']:
            print('*** No known tf interactions ***\n')
            print_indicator = 1

        else:
            p1 = l[2]
            p2 = l[3]

            experimental_score = float(l[10])

            output_text = ("    ".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))

            new_file = input_file.split('.csv')[0]
            #new_file = input_file.split('./input2/')[1]

            file_name = "../output/%s_network.png" % new_file

            #opens the network image created by get_network
            my_image = Image.open(file_name)

            title_font = ImageFont.truetype('../Ubuntu-Medium.ttf', 12)

            #adds text from interaction scores request to the network image
            image_editable = ImageDraw.Draw(my_image)
            image_editable.text((0,x), output_text, (255, 255, 255), font=title_font)
            x = x + 15

            my_image.save(file_name)

    if print_indicator == 0:
            print('saving interactions\n')
            
    
#inputs tfs file, selects the subset of tfs df that are actually tfs
tfs = pd.read_csv('tfs.csv')
accession_num = pd.read_csv('accession.csv')
tfs = tfs[tfs['TF?'] == 'Yes']
tfs = tfs[['ID', 'Name']]
new_tfs = tfs.merge(accession_num, how='inner')

lookup_dict = {}

#create a dictonary for tfs
for entry in new_tfs.index:
    lookup_dict[new_tfs['Accession'][entry]] = [new_tfs['ID'][entry]]

#change directory to the input directory
cwd = os.getcwd()
new_dir = cwd + '/input'
os.chdir(new_dir)

#loops through each csv file in the input directory
for root,dirs,files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith('.csv'):
            my_genes = []
            input_csv = pd.read_csv(file)

            #compare input file to tfs dict and adds tfs to a new list
            for ind in input_csv.index:
                if input_csv['Accession Number'][ind] in lookup_dict.keys():
                    my_genes.append(lookup_dict.get(input_csv['Accession Number'][ind]))

            my_genes = list(flatten(my_genes))

            #outputs file name to txt file if there are less than two tfs
            if len(my_genes) < 2:
                output_file = "../output/files_without_networks.txt"

                with open(output_file, 'a') as of:
                    of.write(file)
                    of.write('\n')

            #runs get_network and get_interactions on files with more two or more tfs
            else:
                get_network(file, my_genes)
                sleep(1) #wait one second between requests
                get_interactions(file, my_genes)
                sleep(1) #wait one second between requests

os.chdir('..')
