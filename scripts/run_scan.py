from glob import glob
import os
import json
import yaml
import argparse
import subprocess

### USAGE EXAMPLE
# ```
# python3 scripts/run_scan.py 
# ```

cfg = 'configs/harvestDatacards.yml'

with open(cfg, 'r') as file:
   setup = yaml.safe_load(file)

folder = setup['output_folder']

def get_args():
    parser = argparse.ArgumentParser(description="Run combine workflow for CP analysis")
    return parser.parse_args()

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
        print("\n")
        print("*"*140)
        print(f"\033[1;32m\nCommand '{command}' ran successfully!\n\033[0m")
        print("*"*140)
        print("\n")
    except subprocess.CalledProcessError as e:
        print("\n")
        print("*"*140)
        print(f"\033[1;31m\nCommand '{command}' failed with exit code: {e.returncode}\n\033[0m")
        print("*"*140)
        print("\n")
        raise

def run_scan(): #vsjet, name):
    # HARVEST DATACARDS
    run_command(f"python3 scripts/harvestDatacards.py") 
    # CREATE COMBINE WORKSPACE
    run_command(f"combineTool.py -m 125 -M T2W -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays -i {folder}/cmb -o ws.root --parallel 8")
    # RUN MAX LIKELIHOOD FIT
    run_command(f"combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d {folder}/cmb/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1")
    # PLOT 1D SCAN OF ALPHA
    run_command(f"python3 scripts/plot1DScan.py --main={folder}/cmb/higgsCombine.alpha.MultiDimFit.mH125.root --POI=alpha --output={folder}/alpha_cmb --no-numbers --no-box --x-min=-90 --x-max=90 --y-max=8")


def main():
    # args = get_args()
    try:
        run_scan() #args.vsjet, args.name)
    except subprocess.CalledProcessError:
        raise RuntimeError("Production failed!")

if __name__ == "__main__":
    main()