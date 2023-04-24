from utils import is_valid_pdb
from utils import rel_path
from utils import run_bash_command
from utils import make_list_of_pdbs
from utils import get_lines
from utils import write_lines
from rna_tools.rna_tools_lib import *
import os
import shutil
import subprocess
import glob
import json
import csv
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# TODO: Consider using decorations instead
class Metric:
    
    project_path = os.path.abspath(os.getcwd())


    def __init__(self, name):
        self.name = name
        self.binary_path = None

        os.makedirs("scores", exist_ok=True)


    def chdir_run_folder_ref_mod(self, run_name, r, m):
        run_folder_ref_mod = os.path.abspath(f"runs/{run_name}/{r[0]}/{r[0]}.{m[0]}/{self.name}")

        # Mkdir
        os.makedirs(run_folder_ref_mod, exist_ok=True)
        os.chdir(run_folder_ref_mod)

    def chdir_run_folder_ref(self, run_name, r):
        run_folder_ref = os.path.abspath(f"runs/{run_name}/{r[0]}/{self.name}")

        # Mkdir
        os.makedirs(run_folder_ref, exist_ok=True)
        os.chdir(run_folder_ref)

    def calculate(self, reference, model, run_name="default"):
        reference = os.path.abspath(reference)
        model = os.path.abspath(model)

        assert is_valid_pdb(reference) and is_valid_pdb(model), "reference and model PDB must be valid."

        self.cleanup_stack = []

        r = os.path.splitext(os.path.basename(reference))
        m = os.path.splitext(os.path.basename(model))

        self.chdir_run_folder_ref_mod(run_name, r, m)

        # LGA does not take paths outside the current working directory, so workaround needed
        # Create a folder called references and models
        
        shutil.copy(rel_path(reference), "".join(r))
        shutil.copy(rel_path(model), "".join(m))
        self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), ''.join(r))))
        self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), ''.join(m))))
    
    def cleanup(self):
        while self.cleanup_stack:
            os.remove(self.cleanup_stack[0])
            self.cleanup_stack.pop(0)

        os.chdir(self.project_path)      

    def calc_bulk(self, reference, models, force=False, run_name="default"):
        reference = os.path.abspath(reference)
        r = os.path.splitext(os.path.basename(reference))

        models_abs = [os.path.abspath(model) for model in models]

        self.chdir_run_folder_ref(run_name, r)

        assert is_valid_pdb(reference), f"reference PDB {reference} must be valid."

        self.cleanup_stack = []
        r = os.path.splitext(os.path.basename(reference))
        shutil.copy(rel_path(reference), "".join(r))
        self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), "".join(r))))


        for model in models_abs:
            
            m = os.path.splitext(os.path.basename(model))

            shutil.copy(rel_path(model), "".join(m))
            self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), ''.join(m))))

    def consol_bulk(self, target_name):
        raise NotImplementedError("This functionality has not yet been implemented.")

    def consolidate(self, target_name):
        raise NotImplementedError("This functionality has not yet been implemented.")
    
    def process_out_file(self, out_file):
        raise NotImplementedError("This functionality has not yet been implemented.")

class GDT(Metric):

    def __init__(self):
        super().__init__("gdt")  

        lga_path = os.path.abspath("bins/lga/LGA_package_src")
        path = os.path.abspath(lga_path)
        self.binary_path = path
        os.makedirs("figures", exist_ok=True)

    def consolidate(self, target_name):
        subprocess.run(["echo", "reference,model,GDT_HA,GDT_TS"], stdout=open(f"scores/{self.name}.{target_name}.csv", "w")) # run echo command and write output to file
        for txt in glob.glob(f"runs/*/*/*/{self.name}/{self.name}.*.txt"): # loop through files in runs directory
            print(txt)
            bn = os.path.basename(txt)
            noext = os.path.splitext(bn)[0] # get file name without extension
            print(noext)
            gdt, reference, model = noext.split(".") # split file name by dot
            ha = subprocess.check_output(["grep", "GDT_HA =", txt], encoding='utf-8') # run grep command and get output
            ha = ha.split()[2] # get last element of output
            print(f"ha: {ha}")
            ts = subprocess.check_output(["grep", "GDT_TS =", txt], encoding='utf-8') # run grep command and get output
            ts = ts.split()[5] # get last element of output
            print(f"ts: {ts}")

            p1 = subprocess.Popen(["echo", f"{reference}.pdb,{model}.pdb,{ha},{ts}"], stdout=subprocess.PIPE, encoding='utf-8') # create first process with echo command and pipe output
            p2 = subprocess.Popen(["tee", "-a", f"scores/{self.name}.{target_name}.csv"], stdin=p1.stdout, encoding='utf-8') # create second process with tee command and pipe input from first process

            p1.stdout.close() # close first process output stream
            p2.communicate() # wait for second process to finish

    def calculate(self, reference, model, force=False, run_name="default"):
        super().calculate(reference, model, run_name=run_name)

        r = os.path.splitext(os.path.basename(reference))
        m = os.path.splitext(os.path.basename(model))
        output = f"{self.name}.{r[0]}.{m[0]}.txt"

        if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
            print(f"{output} exists. Skipping calculations...")
        
        
        else:    
            ## Run commands here
            
            run_lga = os.path.join(self.binary_path, "runlga.mol_mol.pl")
            c1 = f"ulimit -s unlimited && {run_lga} {''.join(m)} {''.join(r)} -4 -d:4 -atom:C4, -stral -o2"
            run_bash_command(c1)

            run_gdt = os.path.join(self.binary_path, "run_GDT_for_structures_with_unknown_residue_residue_correspondences.sh")
            c2 = f"ulimit -s unlimited && {run_gdt} {''.join(m)} {''.join(r)} 2 > {output}"
            run_bash_command(c2)

        ## Cleanup
        self.cleanup()
        
    # Extract "GDT DIST_CUTOFF" and "GDT PERCENT_AT" from runs/rna_only_r1107_processed/rna_only_r1107_processed.R1107TS029_1_processed/gdt/RESULTS/GDT.R1107TS029_1_processed.pdb.rna_only_r1107_processed.pdb.gdt_res and store as a numpy array where "GDT PERCENT_AT" is the x-axis and "GDT DIST_CUTOFF" is the y-axis
    def extract_gdt_percent(self, gdt_res_file):
        # Initialize variables to store the values
        dist_cutoff = []
        percent_at = []

        # Open the file and read the lines
        with open(gdt_res_file, 'r') as f:
            lines = f.readlines()

            # Loop through the lines
            for line in lines:
                # Check if the line contains the required values
                if 'GDT DIST_CUTOFF' in line:
                    dist_cutoff = [float(val) for val in line.split()[2:]]
                elif 'GDT PERCENT_AT' in line:
                    percent_at = [float(val) for val in line.split()[2:]]

        # Convert the lists to numpy arrays
        dist_cutoff = np.array(dist_cutoff)
        percent_at = np.array(percent_at)

        return percent_at, dist_cutoff
    
    # Write a command that says hello world!
    
    
    # runs/rna_only_r1107_processed/rna_only_r1107_processed.R1107TS029_1_processed/gdt/RESULTS/GDT.R1107TS029_1_processed.pdb.rna_only_r1107_processed.pdb.gdt_res
    # Get the files in runs/*/*/*/gdt/RESULTS/GDT.*.pdb.rna_only_r1107_processed.pdb.gdt_res and run extract_gdt_percent(). Take the results and make a dataframe. Set dist_cut_off as the index and percent_at as the columns. The column name is name of the directory at runs/*/*/*
    def consol_gdt_percent(self, target_name):
        # Get the files in the runs directory
        files = glob.glob(f"runs/{target_name}/*/*/gdt/RESULTS/GDT.*.gdt_res")
        merged = pd.DataFrame()
        for file in files:
            # Extract the "GDT PERCENT_AT" and "GDT DIST_CUTOFF" from the file
            percent_at, dist_cutoff = self.extract_gdt_percent(file)
            
            # Get the directory name
            match = re.search(r"(\w+.\w+.)\/gdt\/RESULTS", file)
            model_ref = ""
            
            model_ref = match.group(1)
            print(f"model_ref: {model_ref}")  # Output: abc

            df = pd.DataFrame(percent_at, index=pd.Index(data=dist_cutoff, name="dist_cutoff"), columns=[model_ref])

            if merged.empty:
                merged = df
            else:
                merged = merged.join(df)
                    
        out = f"scores/{self.name}_percent.{target_name}.csv"
        merged.to_csv(out, header=True, sep=',')
        return merged
    
    def save_consol_gdt_percent_fig(self, target_name, df):
        sns.lineplot(data=df)
        plt.title('Multiple curves in a single line plot')
        plt.show()
        # save figure
        plt.savefig(f"figures/{self.name}_percent.{target_name}.png")
    
    # def consol_gdt_percent(self, target_name):
    #     # Get the files in the runs directory
    #     files = glob.glob(f"runs/*/*/*/{self.name}/{self.name}.*.txt")

    #     # Initialize the numpy array
    #     xy = np.empty((0, 2))

    #     # Loop through the files
    #     for file in files:
    #         # Extract the "GDT PERCENT_AT" and "GDT DIST_CUTOFF" from the file
    #         xy = np.append(xy, self.extract_gdt_percent(file), axis=0)

    #     # Sort the array
    #     xy = xy[xy[:,0].argsort()]

    #     # Write the array to a file
    #     np.savetxt(f"scores/{self.name}.{target_name}.csv", xy, delimiter=",", header="GDT PERCENT_AT,GDT DIST_CUTOFF", comments="")


class INF(Metric):

    def __init__(self):
        super().__init__("inf")  
        self.binary_path = "rna_calc_inf.py"
                

    def consolidate(self, target_name):
        second_lines = get_lines(f"runs/*/*/*/{self.name}/{self.name}.*.txt", 2, 2)
        write_lines(f"scores/{self.name}.{target_name}.csv", second_lines)


    def calculate(self, reference, model, force=False, run_name="default"):
        super().calculate(reference, model, run_name=run_name)

        r = os.path.splitext(os.path.basename(reference))
        m = os.path.splitext(os.path.basename(model))
        output = f"{self.name}.{r[0]}.{m[0]}.txt"
        print(f"output: {output}")

        if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
            print(f"{output} exists. Skipping calculations...")
        
        else:    
            ## Run commands here
            # rna_calc_inf.py -f -pr -t $tmp_ground ${tmp_model}/normalized*.pdb

            c = f"{self.binary_path} -f -pr -o {output} -t {''.join(r)} {''.join(m)}"
            run_bash_command(c)        

        ## Cleanup
        self.cleanup()

    def calc_bulk(self, reference, models, force=False):
        super().calc_bulk(reference, models)

        r = os.path.splitext(os.path.basename(reference))
        # m = [''.join(os.path.splitext(os.path.basename(m))) for m in models]
        
        # Intentionally left a space in ' '.join(m)
        models = [os.path.basename(model) for model in models]
        m = ' '.join(models)
        
        output = f"{self.name}.{r[0]}.txt"
        print(f"output: {output}")

        if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
            print(f"{output} exists. Skipping calculations...")
        else:    
            c = f"{self.binary_path} -f -pr -o {output} -t {''.join(r)} {m}" 
            run_bash_command(c)        

        ## Cleanup
        self.cleanup()

    def consol_bulk(self, target_name):
        # Loop through all the csv files that match the path pattern
        print(os.getcwd())
        exclude_header = get_lines(f"runs/*/*/{self.name}/{self.name}.*.txt", 2, -1)
        write_lines(f"scores/{self.name}.{target_name}.csv", exclude_header)

class Clashscores(Metric):

    def __init__(self):
        super().__init__("clashscores")
        
        phenix_path = os.path.abspath("bins/phenix/phenix.clashscore")
        self.binary_path  = os.path.abspath(phenix_path)

    def consolidate(self, target_name):
        print()

        lines = []
        
        # Read the files using glob.glob(f"runs/*/*/*/{self.name}/{self.name}.*.txt")
        for filename in glob.glob(f"runs/*/*/*/{self.name}/{self.name}.*.txt"):
            print(filename)
            with open(filename) as f:
                for line in f:
                    lines.append(line)

        with open(f"scores/{self.name}.{target_name}.csv", "w") as outfile:
            outfile.flush()

            outfile.write("reference,model,clashscore\n")

            for line in lines:
                outfile.write(line + "\n")

    
    def calculate(self, reference, model, force=False, run_name="default"):
        super().calculate(reference, model)

        r = os.path.splitext(os.path.basename(reference))
        m = os.path.splitext(os.path.basename(model))
        output = f"{self.name}.{r[0]}.{m[0]}.txt"

        if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
            print(f"{output} exists. Skipping calculations...")
        
        else:    
            m_name = ''.join(m)
            c = f"{self.binary_path} {m_name} nuclear=True keep_hydrogens=True"
            phenix_out = str(run_bash_command(c))

            print(phenix_out)

            match = re.search(r'clashscore = (\d+\.\d+)', phenix_out)
            assert match, "Clashscore was not found in phenix_out. Clashscore must be a number and phenix must be ran successfully."

            clashscore = float(match.group(1))

            with open(output, 'a') as f:
                f.write(''.join(r) + ',' + ''.join(m) + ',' + str(clashscore))

        ## Cleanup
        self.cleanup()

class TMScore(Metric):
    
        def __init__(self):
            super().__init__("tm_score")

            us_align_path = os.path.abspath("bins/us-align/USalign")
            path = os.path.abspath(us_align_path)
            self.binary_path = path

        def consolidate(self, target_name):
            # Raise notyetimplemented error
            raise NotImplementedError("TMScore consolidation is not yet implemented")
        
        def calculate(self, reference, model, force=False):
            raise NotYetImplementedError("TMScore calculation is not yet implemented")
        
        def calc_bulk(self, reference, models, force=False):
            # USalign -dir1 ${models_path}/ ${models_path}/list -suffix .pdb $grounds_file -outfmt 2 > output_table/usalign_${name}.${grounds_name}.csv
            super().calc_bulk(reference, models)
            make_list_of_pdbs(os.getcwd())
            # print("somettest")


            r = os.path.splitext(os.path.basename(reference))
            
            # Intentionally left a space in ' '.join(m)
            models = [os.path.basename(model) for model in models]
            m = ' '.join(models)
            
            output = f"{self.name}.{r[0]}.txt"

            if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
                print(f"{output} exists. Skipping calculations...")
            else:    
                c = f"{self.binary_path} -dir1 ./ list -suffix .pdb {''.join(r)} -outfmt 2 > {output}" 
                run_bash_command(c)        

            ## Cleanup
            self.cleanup()

        def consol_bulk(self, target_name):
            # Define the columns to rename
            column_map = {'#PDBchain1': 'model', 'PDBchain2': 'reference', 'TM1': 'tm_score'}

            # Create an empty DataFrame to hold the combined data from all input files
            combined_df = pd.DataFrame()

            # Loop through each file in the input directory
            for filename in glob.glob(f"runs/*/*/{self.name}/{self.name}.*.txt"):
                # Read the input CSV file into a DataFrame
                input_df = pd.read_csv(filename, delimiter='\t')
                input_df['#PDBchain1'] = input_df['#PDBchain1'].str.split(':').str[0]
                input_df['PDBchain2'] = input_df['PDBchain2'].str.split(':').str[0]
                input_df = input_df[['#PDBchain1', 'PDBchain2', 'TM1']]


                # Rename the columns
                input_df = input_df.rename(columns=column_map)

                # Append the renamed data to the combined DataFrame
                combined_df = combined_df.append(input_df, ignore_index=True)

            # Write the combined data to the output CSV file
            combined_df.to_csv(f"scores/{self.name}.{target_name}.csv", index=False)

class LDDT(Metric):
    
        def __init__(self):
            super().__init__("lddt")
            
            lddt_path = os.path.abspath("bins/lddt")
            self.binary_path  = os.path.abspath(lddt_path)
    
        def consolidate(self, target_name):
            # Define the headers for the CSV file
            headers = ['reference', 'model', 'lddt']

            # Open the CSV file for writing
            with open(f"scores/{self.name}.{target_name}.csv", 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)

                # Iterate through all the JSON files in the directory
                for json_file in glob.glob(f"runs/*/*/*/{self.name}/{self.name}.*.txt"):
                    # Open the JSON file and extract the relevant data
                    with open(json_file) as jf:
                        data = json.load(jf)
                        reference = data['trg_file']
                        model = data['mdl_file']
                        lddt = data['lDDT']

                        # Write the data to the CSV file
                        writer.writerow([reference, model, lddt])
        
        def calculate(self, reference, model, force=False):
            super().calculate(reference, model)
    
            r = os.path.splitext(os.path.basename(reference))
            m = os.path.splitext(os.path.basename(model))
            output = f"{self.name}.{r[0]}.{m[0]}.txt"
    
            if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
                print(f"{output} exists. Skipping calculations...")
            
            else:    
                c = f"{self.binary_path}/ema --mount {os.getcwd()} {self.binary_path}/scoring/monomer_lddt_no_stereocheck.py {''.join(m)} {''.join(r)} {output}"
                run_bash_command(c)
    
            ## Cleanup
            self.cleanup()


