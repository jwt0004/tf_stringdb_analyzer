Title / tf_stringdb_analyzer

Author / Jared Taylor

About / This program takes an input directory of csv files that contain lists of proteins
and compares the list of proteins to a list of transcription factors (tfs). The proteins
that match up with the list of tfs are then analyzed for interactions using the
Stringdb API. It will ouput a png file with the network map of interactions between the proteins
and the interaction scores for the interactions. It will complete this for every csv
file in the input directory. If a file doesn't have at least two tfs in it, it will be added
to an output txt file (files_without_networks.txt). The program will also output the network
that it is saving and wether it had interactions or not.

Files / tf_stringdb_analyzer.py, accession.csv, tfs.csv, Ubuntu-Medium.ttf, input, and output

Libraries / pandas, os, requests, time, and image

Requirements / os and time should be included with python, but pandas, image, and requests will need to be
downloaded before running

How to use / add csv files with proteins in the input directory

    *** Must be ran from the terminal ***

    macOS : open terminal, move to the tf_stringdb_analyzer directory, then input the command below

        python3 tf_stringdb_analyzer.py

    Windows : open powershell or bash, move to the tf_stringdb_analyzer directory, then input the command below

        python tf_stringdb_analyzer.py

Help / For any questions, email jared1taylor12@gmail.com
