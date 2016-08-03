#!/usr/bin/env python
#

# import modules used here -- sys is a very standard one
from __future__ import print_function
import argparse
import logging
from collections import OrderedDict
from glob import glob
import os
import sys

import nibabel
import json
import pandas as pd


# Gather our code in a main() function
from shutil import copy


def get_metadata_for_nifti(bids_root, path):

    sidecarJSON = path.replace(".nii.gz", ".json")

    pathComponents = os.path.split(sidecarJSON)
    filenameComponents = pathComponents[-1].split("_")
    sessionLevelComponentList = []
    subjectLevelComponentList = []
    topLevelComponentList = []
    ses = None;
    sub = None;

    for filenameComponent in filenameComponents:
        if filenameComponent[:3] != "run":
            sessionLevelComponentList.append(filenameComponent)
            if filenameComponent[:3] == "ses":
                ses = filenameComponent
            else:
                subjectLevelComponentList.append(filenameComponent)
                if filenameComponent[:3] == "sub":
                    sub = filenameComponent
                else:
                    topLevelComponentList.append(filenameComponent)

    topLevelJSON = os.path.join(bids_root, "_".join(topLevelComponentList))
    potentialJSONs = [topLevelJSON]

    subjectLevelJSON = os.path.join(bids_root, sub, "_".join(subjectLevelComponentList))
    potentialJSONs.append(subjectLevelJSON)

    if ses:
        sessionLevelJSON = os.path.join(bids_root, sub, ses, "_".join(sessionLevelComponentList))
        potentialJSONs.append(sessionLevelJSON)

    potentialJSONs.append(sidecarJSON)

    merged_param_dict = {}
    for json_file_path in potentialJSONs:
        if os.path.exists(json_file_path):
            param_dict = json.load(open(json_file_path, "r"))
            merged_param_dict.update(param_dict)

    return merged_param_dict


def run(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    if not os.path.exists(args.output_directory):
        logging.info("creating output directory at '{}'".format(args.output_directory))
        os.makedirs(args.output_directory)

    subject_ids = []
    study_dict = OrderedDict()
    for file in glob(os.path.join(args.bids_directory, "sub-*")):
        if os.path.isdir(file):
            subject_ids.append(os.path.split(file)[-1][4:])
    subject_ids.sort()
    study_dict["Source Name"] = subject_ids
    study_dict["Characteristics[organism]"] = "Homo sapiens"
    study_dict["Term Source REF"] = "NCBITAXON"
    study_dict["Term Accession Number"] = "NCBITaxon:9606"
    study_dict["Characteristics[organism part]"] = "brain"
    study_dict["Protocol REF"] = "Participant recruitment"
    study_dict["Sample Name"] = subject_ids
    df = pd.DataFrame(study_dict)

    participants_file = os.path.join(args.bids_directory, "participants.tsv")
    if os.path.exists(participants_file):
        participants_df = pd.read_csv(participants_file, sep="\t")
        participants_df.rename(columns={'participant_id': "Sample Name"}, inplace=True)
        participants_df["Sample Name"] = [s[4:] for s in list(participants_df["Sample Name"])]
        for col in participants_df.columns.tolist():
            if col != "Sample Name":
                participants_df.rename(columns={col:"Comment[%s]"%col}, inplace=True)
        df = pd.merge(df, participants_df, left_on="Sample Name", right_on="Sample Name")

    df.to_csv(os.path.join(args.output_directory, "s_study.txt"), sep="\t", index=False)

    assay_dict = OrderedDict()
    sample_names = []
    raw_file = []
    types = []
    resolutions = []
    resolutions_units = []
    rts = []
    rts_units = []
    assay_names = []
    other_fields = []
    mri_par_names = []

    for file in glob(os.path.join(args.bids_directory, "sub-*", "*", "sub-*.nii.gz")) + \
            glob(os.path.join(args.bids_directory, "sub-*", "ses-*", "*", "sub-*_ses-*.nii.gz")):
        sample_names.append(os.path.split(file)[-1].split("_")[0][4:])
        assay_names.append(os.path.split(file)[-1].split(".")[0])
        raw_file.append(file[len(args.bids_directory):])
        types.append(file.split("_")[-1].split(".")[0])
        header = nibabel.load(file).get_header()
        resolutions.append("x".join([str(i) for i in header.get_zooms()[:3]]))
        resolutions_units.append(header.get_xyzt_units()[0])
        if len(header.get_zooms()) > 3:
            rts.append(header.get_zooms()[3])
            rts_units.append(header.get_xyzt_units()[1])
        else:
            rts.append(None)
            rts_units.append(None)
        other_fields.append(get_metadata_for_nifti(args.bids_directory, file))

    assay_dict["Sample Name"] = sample_names
    assay_dict["Protocol REF"] = "Magnetic Resonance Imaging"
    assay_dict["Parameter Value[Modality]"] = types
    mri_par_names.append("Modality")
    assay_dict["Parameter Value[Resolution]"] = resolutions
    mri_par_names.append("Resolution")
    assay_dict["Unit"] = resolutions_units

    new_fields = set()
    for d in other_fields:
        new_fields = new_fields.union(set(d.keys()))

    for field in new_fields:
        assay_dict["Parameter Value[%s]"%field] = []

        for d in other_fields:
            if field in d:
                assay_dict["Parameter Value[%s]"%field].append(d[field])
            else:
                assay_dict["Parameter Value[%s]"%field].append(None)

    assay_dict["Assay Name"] = assay_names
    assay_dict["Raw Data File"] = raw_file

    df = pd.DataFrame(assay_dict)
    df.to_csv(os.path.join(args.output_directory, "a_assay.txt"), sep="\t", index=False)

    this_path = os.path.join(os.path.realpath(__file__))
    template_path = os.path.join(*(os.path.split(this_path)[:-1] + ("i_investigation_template.txt", )))
    investigation_template = open(template_path).read()

    title = os.path.split(args.bids_directory)[-1]

    if os.path.exists(os.path.join(args.bids_directory, "dataset_description.json")):
        with open(os.path.join(args.bids_directory, "dataset_description.json"), "r") as description_dict_fp:
            description_dict = json.load(description_dict_fp)
            if "Name" in description_dict:
                title = description_dict["Name"]

    investigation_template = investigation_template.replace("[TODO: TITLE]", title)
    investigation_template = investigation_template.replace("[TODO: MRI_PAR_NAMES]", ";".join(mri_par_names))

    with open(os.path.join(args.output_directory, "i_investigation.txt"), "w") as fp:
        fp.write(investigation_template)


def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(
        description="BIDS to ISA-Tab converter.",
        fromfile_prefix_chars='@')
    # TODO Specify your real parameters here.
    parser.add_argument(
        "bids_directory",
        help="Location of the root of your BIDS compatible directory",
        metavar="BIDS_DIRECTORY")
    parser.add_argument(
        "output_directory",
        help="Directory where ISA-TAB files will be stored",
        metavar="OUTPUT_DIRECTORY")
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true")
    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    run(args, loglevel)
    print("Metadata extraction complete.")


if __name__ == '__main__':
    main()
