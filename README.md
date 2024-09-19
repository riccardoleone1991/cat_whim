# CAT-WHIM

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

In this project we want to assess the Combined effect of Amyloid, Tau and WHIte Matter hyperintensity-related disconnections on cortical thickness in different stages of Alzheimer's Diseases (CAT-WHIM).

## Project Organization

```
├── LICENSE            <- Open-source license if one is chosen
├── Makefile           <- Makefile with convenience commands like `make data` or `make train`
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── external       <- Data from third party sources.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump.
│
├── docs               <- A default mkdocs project; see mkdocs.org for details
│
├── models             <- Trained and serialized models, model predictions, or model summaries
│
├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
│
├── pyproject.toml     <- Project configuration file with package metadata for cat_whim
│                         and configuration for tools like black
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
│
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
│                         generated with `pip freeze > requirements.txt`
│
├── setup.cfg          <- Configuration file for flake8
│
└── cat_whim                <- Source code for use in this project.
    │
    ├── __init__.py    <- Makes cat_whim a Python module
    │
    ├── data           <- Scripts to download or generate data
    │   └── dataset.py
    │
    ├── features       <- Scripts to turn raw data into features for modeling
    │   └── build_features.py
    │
    ├── models         <- Scripts to train models and then use trained models to make
    │   │                 predictions
    │   ├── predict_model.py
    │   └── train_model.py
    │
    └── visualization  <- Scripts to create exploratory and results oriented visualizations
        └── visualize.py
```

--------


## How to use
NOTE: Thankds to the beauty of cookie-cutter data science (https://cookiecutter-data-science.drivendata.org/using-the-template/), this project aims at being almost fully automated. Nonetheless, it still needs a bit of human effort in order to select the subjects to download, since one needs to have access to ADNI in order to download their data. Also, please note that you will need a HPC cluster to perform all the preprocessing steps in a reasonable time frame, otherwise this could be the only project you do in your lifetime on your personal laptop. At the moment I am still looking to find the best way to run make data because several jobs that I run in parallel take days but I don't know how inform the main script when everything is finished, so I ran it in steps...

In this study we make use of amyloid and tau PET preprocessed data from UC Berkeley (e.g., UCBERKELEY_AMY_6MM_17Jul2024.csv). First, you need to download them from the ADNI website into the data/raw folder.

0. After downloading this repo, you should first create the conda environment that is needed to run all the scripts. This can easily be achieved by typing 'make create_environment' in a bash terminal. After this, install the requirements by typing 'make requirements'. Now you can move forward to the next steps!

1. You first need to run 'make select' to create a list of subjects with available PET data at the same timepoint (defined to be a 6-months interval between amyloid an tau PET or viceversa). This will create a .csv file (saved in data/utils) containing the names of the subjects to download. This command will also print on your screen a list of subjects separated by commas so that it is easier to just copy/paste them into the ADNI website under Download --> Image Collections --> Advanced Search --> Subject ID, select "T1" AND "T2" and also check the box "3D" and create a new collection of images called cat_whim. Please note that this will also add some sequences that are neither T1 nor FLAIR, but to avoid the tedious choice of selecting which sequence to download for each patient for all patient, we download some data that we are not going to use and then filter it later based on sequence names. If you have data storage concerns, you can select only the specific sequence names corresponding to T1 and FLAIR for each subject. Next, on the ADNI website, go to Data Collections --> cat_whim --> Not downloaded --> select all subjects --> click advanced download and download the .csv file containing the links (URL List, bottom right) saving it in data/utils/cat_whim_url_list.csv.
2. Now you should open a tmux terminal like 'tmux new -s cat_whim', then, in the new terminal, activate the cat_whim environment with 'source ~/miniconda3/etc/profile.d/conda.sh' (or different, depending on where you installed miniconda) followed by 'conda activate cat_whim' and then run 'python dataset.py'. Remember that to quit a tmux terminal you can hit "Ctrl+B" and then "D". If you hit "Ctrl+D" directly it closes the terminal directly. You can go again to your tmux by typing 'tmux attach -t cat_whim' This will automatically download data, check that downloaded subjects match those that you wanted to download and get rid of unused sequences.
3. For this step, you should have downloaded the bids-freesurfer singularity image from [here](https://hub.docker.com/r/bids/freesurfer/). We used freesurfer_7.4.1-202309. Change the settings of run_freesurfer_longitudinal.sh to your own folder/HPC settings, then type 'sbatch run_freesurfer_longitudinal.sh' in your terminal. This will perform the longitudinal pipeline for all subjects. The timeframes are roughly 6 hours for one session up to 1.5 days for 4 sessions.
4. For this step, you should have downloaded the smriprep singularity image. We used smriprep-0.16.0. Change the settings of run_smriprep.sh to your own folder/HPC setting and run 'sbatch run_smriprep.sh'. This will reorganize the output of freesurfer in order to run smriprep without having to recalculate everything, so it should be pretty fast, around 1-2 hours for subject.
5. After this, you can run 'sbatch run_wmh_segmentations.sh'. This will launch several instances of automatic segmentations of WMH for all subjects in parallel. It should take 6-10 minutes per subject. When all jobs are done, 'python wmh_segmentation_postprocessing.py' will create some empty masks for those who don't have WMH (so that further processing scripts don't throw errors) and create a new lesion_masks folder that is in the correct folder structure for LQT.
6. Note that we ran the Lesion Quantification Toolkit pipeline from R in Windows 11! The R version of LQT is available [here](https://github.com/jdwor/LQT). We have plans to implement a Python wrapper to LQT but this will take a little more time. So, for now, note that the lesion masks folder was copy/pasted into Windows where the 'run_LQT.R' script was run. Then, the results folder was recopied back into the main project folder in Linux. Running LQT took approximately 1 day on a HP Envy laptop with 12 cores.
7. The last step before starting the analyses is to run 'python make_final_dataset.py'. This will create the final dataset!

Now you are ready to go for your own analyses! If you want to follow ours, everything is in cat_whim/notebooks!
