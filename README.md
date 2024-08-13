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


# Current status: 2024-08-13
Running sMRIPREP longitudinal

TODO: 

--------

- Run sMRIPREP longitudinal
- Collate results from sMRIPrep Freesurfer
- Start analyzing!


## How to use
NOTE: This project is almost fully automated, but still needs a bit of human effort in order to select the subjects to download, since one needs to have access to ADNI in order to download their data.

In this study we make use of amyloid and tau PET preprocessed data from UC Berkeley (e.g., UCBERKELEY_AMY_6MM_17Jul2024.csv). First, you need to download them from the ADNI website into the data/raw folder.

1. You first need to run 'make select' to create a list of subjects with available PET data at the same timepoint (defined to be 6 months interval between amyloid an tau PET or viceversa). This will create a .csv file (saved in data/utils) containing the names of the subjects to download. Copy/paste it into the ADNI website under Download --> Image Collections --> Advanced Search --> Subject ID, select T1 and T2 and create a new collection of images called cat_whim. Please note that this will also add some sequences that are neither T1 nor FLAIR, but to avoid the tedious choice of selecting which sequence to download for each patient for all patient, we download some data that we are not going to use and then filter it later based on sequence names. If you have data storage concerns, you can select only the specific sequence names corresponding to T1 and FLAIR for each subject. Next go to Data Collections --> cat_whim --> Not downloaded --> select all subjects --> click advanced download and download the .csv file containing the links (URL List, bottom right) saving it in data/utils/cat_whim_url_list.csv
2. Change the settings of run_smriprep.sh to your own folder/HPC settings --> my cat_whim is under /home/leoner/Projects/cat_whim, you should change to yours...
3. Now you can just run 'make data'. This will automatically download data, check that downloaded subjects match those that you wanted to download and get rid of unused sequences.