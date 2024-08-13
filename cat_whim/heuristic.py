#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    Heuristic -- Version 1.2
Last edit:  2024/19/02
Author(s):  HeuDiConv devs (HD)
            Kobeleva, Xenia (XG)
            Geysen, Steven (SG)
            Leone, R
Notes:      - Heuristic script used for conversion to NIfTI
            - Release notes:
                * Convert fieldmaps
To do:      
Comments:   HD: Known BIDS labels:
                * anat - anatomical data.  Might also be collected multiple
                    times across runs (e.g. if subject is taken out of magnet
                    etc), so could (optionally) have "_run" definition
                    attached. For "standard anat" labels, please consult to
                    "8.3 Anatomy imaging data" but most common are 'T1w',
                    'T2w', 'angio', 'FLAIR'
                * func - functional (AKA task, including resting state) data.
                    Typically contains multiple runs, and might have multiple
                    different tasks different per each run
                    (e.g. _task-memory_run-01, _task-oddball_run-02)
                * fmap - field maps
                * dwi - diffusion weighted imaging (can as well have runs)
"""



#%% ~~ Imports ~~ %%#


import os



#%% ~~ Funcitons ~~ %%#
#######################


def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')

    return template, outtype, annotation_classes


def infotodict(seqinfo):
    """
    Heuristic evaluator for determining which runs belong where

    Allowed template fields - follow python string module:
        * item: index within category
        * subject: participant id
        * seqitem: run number during scanning
        * subindex: sub index within group
    """

    t1w = create_key(
        'sub-ADNI{subject}/{session}/anat/sub-ADNI{subject}_{session}_T1w'
        )
    flair = create_key(
        'sub-ADNI{subject}/{session}/anat/sub-ADNI{subject}_{session}_FLAIR'
        )
    keys = [t1w, flair]
    info = {keyi: [] for keyi in keys}
    
    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:
            * total_files_till_now
            * example_dcm_file
            * series_id
            * dcm_dir_name
            * unspecified2
            * unspecified3
            * dim1
            * dim2
            * dim3
            * dim4
            * TR
            * TE
            * protocol_name
            * is_motion_corrected
            * is_derived
            * patient_id
            * study_description
            * referring_physician_name
            * series_description
            * image_type
        """
        
        if (s.dim1 == 256) and (('spgr' in s.series_description.lower()) or ('rage' in s.series_description.lower()) or ('t1' in s.series_description.lower())):
            info[t1w].append(s.series_id)
        if (s.dim1 == 256) and ('flair' in s.series_description.lower()):
            info[flair].append(s.series_id)

    return info
