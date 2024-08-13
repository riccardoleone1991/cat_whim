#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = cat_whim
PYTHON_VERSION = 3.11
PYTHON_INTERPRETER = python
SHELL=/bin/bash

#################################################################################
# COMMANDS                                                                      #
#################################################################################


## Install Python Dependencies
.PHONY: requirements
requirements:
	$(PYTHON_INTERPRETER) -m pip install -U pip
	$(PYTHON_INTERPRETER) -m pip install -r requirements.txt
	



## Delete all compiled Python files
.PHONY: clean
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Lint using flake8 and black (use `make format` to do formatting)
.PHONY: lint
lint:
	flake8 cat_whim
	isort --check --diff --profile black cat_whim
	black --check --config pyproject.toml cat_whim

## Format source code with black
.PHONY: format
format:
	black --config pyproject.toml cat_whim




## Set up python interpreter environment
.PHONY: create_environment
create_environment:
	
	conda create --name $(PROJECT_NAME) python=$(PYTHON_VERSION) -y
	
	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"
	



#################################################################################
# PROJECT RULES                                                                 #
#################################################################################
## Select Subjects
.PHONY: select
data: requirements
	$(PYTHON_INTERPRETER) cat_whim/select_subjects_to_download.py


## Make Dataset
.PHONY: data
data: requirements
	module load singularity
	$(PYTHON_INTERPRETER) cat_whim/dataset.py
	sbatch cat_whim/run_smriprep.sh
	


#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys; \
lines = '\n'.join([line for line in sys.stdin]); \
matches = re.findall(r'\n## (.*)\n[\s\S]+?\n([a-zA-Z_-]+):', lines); \
print('Available rules:\n'); \
print('\n'.join(['{:25}{}'.format(*reversed(match)) for match in matches]))
endef
export PRINT_HELP_PYSCRIPT

help:
	@$(PYTHON_INTERPRETER) -c "${PRINT_HELP_PYSCRIPT}" < $(MAKEFILE_LIST)
