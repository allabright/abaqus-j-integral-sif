import argparse
import json
import os
import pathlib
import shutil
import subprocess
from analysis.analyse_results import *

class Runner:

    def __init__(self, run_type):
        self.run_type = run_type
        self.base_dir = pathlib.Path(__file__).resolve().parent.parent.parent
        
    def print_spacers(self, title=None, spacers=False, spacer_line = "-" * 120):
            if spacers:
                print(spacer_line)
            if title:
                print(title)
                if spacers:
                    print(spacer_line)
            return 0

    def parse_arguments(self):

        parser = argparse.ArgumentParser(
            description="Run export and/or analysis for an Abaqus job."
        )
        parser.add_argument(
            "--dimensions",
            help="Perform the export step (process the Abaqus output file).")

        parser.add_argument(
            "--model",
            help="The model to be used for the export or analysis.")

        parser.add_argument(
            "--input",
            help="Directory number to use for existing results.")
        
        parser.add_argument(
            "--config",
            help="Filename of the config file to use"
        )
    
        args = parser.parse_args()

        self.dimensions = args.dimensions
        self.model_name = "_".join(word.capitalize() for word in args.model.split("_")) + f"_{self.dimensions.upper()}"

        if args.input:
            self.input = args.input
        else:
            self.input = None

        if args.config:
            self.config_path = args.config
        else:
            self.config_path = os.path.join(self.base_dir, "config", "config.ini")

        return self
    
    def get_last_directory_number(self, dir, prefix):
        """
        Get the number of the last directory. Directories are incremented progressively for each analysis,
        i.e. 001, 002 etc. We need the last number to create the new directory.
        """
    
        existing_dirs = [
                d for d in os.listdir(dir)
                if d.startswith(prefix.title()) and d[len(prefix.title())+1:].isdigit()
            ]
        
        existing_numbers = [int(d[-3:]) for d in existing_dirs] if existing_dirs else []
        last_number = max(existing_numbers, default=0)  # Increment highest number

        # Format the new directory name with zero padding (e.g., "001", "002", etc.)
        return last_number
    
    def create_run_directories(self):
        """
        Create a directory for this specific run of the tool, and populate with subdirectories.
        """

        self.print_spacers(f"Creating Directories for Job: {self.model_name}", spacers=True)

        self.data_dir = os.path.join(self.base_dir, "data")

        if self.run_type == "analysis":
            self.analysis_dir = os.path.join(self.data_dir, "analysis")
        elif self.run_type == "export":
            self.analysis_dir = os.path.join(self.data_dir, "exports")
        elif self.run_type == "model":
            self.analysis_dir = os.path.join(self.data_dir, "models")

        self.model_dir = os.path.join(self.analysis_dir, self.model_name)

        dirs = [
            self.data_dir,
            self.analysis_dir,
            self.model_dir,      
        ]
        
        for dir in dirs:
            try:
                os.makedirs(dir, exist_ok=False)
                print(f"Created directory -> {dir}")                
            except OSError:
                print(f"Directory exists  -> {dir}")
                os.makedirs(dir, exist_ok=True)

        next_number = self.get_last_directory_number(self.model_dir, self.run_type) + 1

        run_dir_name = f"{self.run_type.title()}_{next_number:03d}"
        self.run_dir = os.path.join(self.model_dir, run_dir_name)

        try:
            os.makedirs(self.run_dir, exist_ok=False)
            print(f"Created directory -> {self.run_dir}")                
        except OSError:
            print(f"Directory exists  -> {self.run_dir}")
            os.makedirs(self.run_dir, exist_ok=True)
               
        return self
    
    def copy_config_data(self):
        """
        Copy the configuration data from the old directory to the new. We want to keep track of which configuration was used for each analysis.
        """

        src_config_file = self.config_path
        dest_config_file = os.path.join(self.run_dir, "config.ini")

        if os.path.isfile(src_config_file):  # Ensure it's a file, not a directory
                shutil.copy2(src_config_file, dest_config_file)
                print(f"Copied: {src_config_file} -> {dest_config_file}")

        return self


    def copy_and_read_input_data(self):
        """
        If the run is being performed with some already existing exported data as an input, this function copies that data
        to the run directory, to be used as the input to the analysis.
        """

        self.print_spacers("Copying Exported Data to the Run Directory", spacers=True)
        
        if self.run_type == "analysis":
            input_data_dir = "exports"
            prefix = "Export"
        elif self.run_type == "export":
            input_data_dir = "models"
            prefix = "Model"

        if self.input:
            last_number = str(self.input).zfill(3)
        else:
            last_number = str(self.get_last_directory_number(os.path.join(self.data_dir, input_data_dir, self.model_name), prefix + "_")).zfill(3)
            print(last_number)

        data_dir_name = prefix + "_" + last_number

        input_data_path = os.path.join(
            self.data_dir,
            input_data_dir,
            self.model_name,
            data_dir_name,
        )

        output_data_path = os.path.join(
            self.run_dir,
            data_dir_name
        )

        try:
            os.makedirs(output_data_path, exist_ok=False)
            print(f"Created directory -> {output_data_path}")                
        except OSError:
            print(f"Directory exists  -> {output_data_path}")
            os.makedirs(output_data_path, exist_ok=True)

        print("Copying exported data to the run directory...")
        for file in os.listdir(input_data_path):

            src_file = os.path.join(input_data_path, file)

            if file == "config.ini":
                dest_file = os.path.join(self.run_dir, file)
            else:
                dest_file = os.path.join(output_data_path, file)

            if os.path.isfile(src_file):  # Ensure it's a file, not a directory
                shutil.copy2(src_file, dest_file)
                print(f"Copied: {src_file} -> {dest_file}")

            elif os.path.isdir(src_file):  # Copy directories properly
                dest_file = os.path.join(self.run_dir, file)
                shutil.copytree(src_file, dest_file, dirs_exist_ok=True)
                print(f"Copied directory: {src_file} -> {dest_file}")

        if self.run_type == "analysis":
            print("Reading the exported data into memory...")
            json_file = os.path.join(output_data_path, prefix + "_" + self.model_name + ".json")
            with open(json_file, "r") as f:
                self.input_data = json.load(f)
        elif self.run_type == "export":
            self.odb_path = os.path.join(output_data_path, self.model_name + ".odb")

        return self
        
    def submit_abaqus_command(self):
        """
        This method submits the command to call the other script and run Abaqus. The export script is called with no GUI, the modelling script is called
        with the GUI.
        """

        self.print_spacers(f"Connecting to Abaqus", spacers=True)

        self.src_dir = os.path.join(self.base_dir, "src")

        if self.run_type == "export":
            script_path = os.path.join(self.src_dir, "export", "export_results.py")
            self.abaqus_command = f"cd {self.base_dir} && abaqus cae noGUI={script_path} -- {self.dimensions} {self.model_name} {self.run_dir} {self.odb_path}"

        elif self.run_type == "model":
            script_path = os.path.join(self.src_dir, "model", self.model_name.lower() + ".py")
            self.model_path = os.path.join(self.run_dir, self.model_name + ".cae")
            self.abaqus_command = f"cd {self.base_dir} && abaqus cae script={script_path} -- {self.dimensions} {self.model_name} {self.run_dir} {self.model_path}"
        
        subprocess.run(self.abaqus_command, shell=True, check=True)
        
        return self
    
    def run_analysis(self):

        self.print_spacers(f"Running Analysis for {self.dimensions.upper()} Scenario", spacers=True)
        self.analysis = AnalyseResults(self.dimensions, self.input_data, self.run_dir).analyse_results()

        return self
    
    def save_output_data_to_json(self):
        """
        Save the processed results data to a JSON file inside the run directory.
        """

        self.results_filename = "Analysis_" + self.model_name + ".json"
        self.results_filename_csv = "Analysis_" + self.model_name + ".csv"

        output_data_path = os.path.join(self.run_dir, self.results_filename)
        output_data_path_csv = os.path.join(self.run_dir, self.results_filename_csv)
        
        with open(output_data_path, "w") as f:
            json.dump(self.analysis.output_data_dict, f, indent=4)

        import pandas as pd

        self.analysis.output_data_df.to_csv(output_data_path_csv)

        self.print_spacers(f"Results data saved to -> {output_data_path}")
        self.print_spacers(spacers=True)

        return self
