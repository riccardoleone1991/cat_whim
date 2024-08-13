from pathlib import Path

from dotenv import load_dotenv
from loguru import logger

# Load environment variables from .env file if it exists
load_dotenv()

# Paths
PROJ_ROOT = Path(__file__).resolve().parents[1]
logger.info(f"PROJ_ROOT path is: {PROJ_ROOT}")

SCRIPTS_DIR = PROJ_ROOT / "cat_whim"

DATA_DIR = PROJ_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"
BIDS_DATA_DIR = INTERIM_DATA_DIR / "bids"
UTILS_DATA_DIR = DATA_DIR / "utils"

MODELS_DIR = PROJ_ROOT / "models"

DOCS_DIR = PROJ_ROOT / "docs"

REPORTS_DIR = PROJ_ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
INCL_EXCL_FLW_DIR = REPORTS_DIR / "incl_excl_flowchart"

SINGIMS_DIR = PROJ_ROOT / "singims"


# If tqdm is installed, configure loguru with tqdm.write
# https://github.com/Delgan/loguru/issues/135
try:
    from tqdm import tqdm

    logger.remove(0)
    logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True)
except ModuleNotFoundError:
    pass
