# Abaqus Stress Intensity Factor Calculator

![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

This repository contains a command-line tool that can be used to post-process the output of finite element analysis performed on a cracked model in Abaqus, in order to calculate the J-integral and the stress intensity factor. The tool operates on either 2D or 3D models, meshed using either CPS6 triangular elements or C3D10 tetrahedral elements. The main novelty of the work is to calculate the J-integral for unstructured tetrahedral meshes in three dimensions, which cannot be performed in Abaqus.

---

## Getting Started

Follow these instructions to set up the environment and run the analysis.

### Prerequisites

- **Abaqus CAE:** This tool interfaces with Abaqus, so a local installation of Abaqus CAE and accompanying license is required to run the model creation and data export scripts. The analysis script can be run stand-alone, using previously exported results.
- **Python:** Python 3.9 or higher.
- **pip:** Python package installer.

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/allabright/abaqus-j-integral-sif.git
    cd https://github.com/allabright/abaqus-j-integral-sif.git
    ```

2.  **Create and activate a virtual environment (recommended):**
    ```bash
    # For macOS/Linux
    python3 -m venv venv
    source venv/bin/activate

    # For Windows
    python -m venv venv
    .\venv\Scripts\activate
    ```

3.  **Install the required dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

---

## Documentation

This work was created as part of an MSc thesis.  The full report is available in `report/full_report.pdf` (), and the accompanying LaTeX files are also available in `report/latex`. The results of the parametric studies are available in the `results` directory, along with a Jupyter Notebook used to plot the necessary graphs.

**[➡️ Read the full thesis here (PDF)](./report/final_report.pdf)**


---

## Usage

### Arguments

#### Mandatory

`--dimension` - The dimension of the model being analysed, either `2d` or `3d`.

`--model` - The name of the model being analysed. The only model supported at present is `Edge_Crack`. If you have your own ready-made model that you want to analyse, it needs to be dragged into the `data/models` directory, and the name of the directory needs to be passed as the value for this argument.

#### Optional

`--input` -- By default, the number of the last directory is used. For example, if the last directory in the `data/models/Edge_Crack_2D` directory is `Model_008`, then that will be used for the export by default. If you want to use another model - for example `Model_005`, then `005` needs to be passed as an argument here.

`--config` -- By default the `config.ini` file in the `config` directory is used to specify the configuration of the analysis. If you want to use a different file, then drag it into the `config` directory and pass the filename as the value for this argument.

### Functions

#### Model Creation

The model must first be created in Abaqus, using the `model.py` file in the `src` directory. This can be done manually as long as the correct approach is followed, but it is automated here for the edge-crack case study. Usage is as follows:

```
python model.py --dimension X --model Y
```

This creates a model with dimension `X` using the model `Y`, using the parameters specified in the configuration file, and runs the analysis. A new directory is created within the `data/models/Edge_Crack_{dimension}` directory for this specific model. The configuration file is copied to this directory, and the results are saved here. This script must be run on a machine which has Abaqus CAE installed, along with the appropriate license.

#### Results Export

Before the analysis can be performed, it is necessary to dump the raw data (e.g. stresses, strains) from Abaqus. This is done using the `export.py` file contained in the `src` directory. This is a command line script, and is run using:

```
python export.py --dimension X --model Y
```

This script creates a new directory for the export under `data/exports/Edge_Crack_{dimension}`, and copies over the model data, and the configuration data used for that model. The data is then dumped using the Abaqus API into a JSON file called `export.json`. This script must be run on a machine which has Abaqus CAE installed, along with the appropriate license.

#### Results Analysis

After exporting the raw data from Abaqus, the analysis can be performed. This is done using the `analysis.py` file contained in the `src` directory. This is also a command line script, and is run using the same arguments as the previous script:

```
python export.py --dimension X --model Y
```

This script creates a new for the analysis under `data/analysis/Edge_Crack_{dimension}`, and copies over the model data, export data, and configuration data. The analysis is then performed, and the results are stored in a JSON file and a CSV file within the run directory.

### Example Use

#### Creating a New Model

No model has been created yet, so you call:

```python model.py --dimension 3d --model edge_crack```

This creates the directory `/data/models/Edge_Crack_3D`, copies over the `config/config.ini` configuration file, then creates and runs a 3D analysis within Abaqus, for a 3D edge-cracked part.

#### Performing an Export

Now, for example, you have run four models. However, you want to perform an export using the third. You call:

```python export.py --dimension 3d --model edge_crack --input 003```

A new directory is created within `data/exports/Edge_Crack_3D`. The latest model located within `data/models/Edge_Crack_3D/Model_004` is ignored because the `--input` argument was passed. Instead, the data is copied from `data/models/Edge_Crack_3D/Model_003` to the run directory. The results of the analysis are exported, and saved in the run directory.

#### Performing an Analysis

Now you have a set of exports within the `data/exports` directory. However, you want to change the weight function used for the analysis, without having to re-run your model and export. You therefore alter the `config/config.ini` file to change the weight function, and save it as `config/config_edited.ini`. You then call:

```python analyse.py --dimension 3d --model edge_crack --config config_edited.ini```

Since no `--input` argument was specified, the most recent export in the `data/exports/Edge_Crack_3D` directory is used. A run directory is created for the analysis in `data/analysis/Edge_Crack_3D`, and the model, export, and configuration data is copied over. Since a `--config` argument was specified, the customized configuration file is used for the analysis. The analysis is performed, and the results are saved in the run directory.
