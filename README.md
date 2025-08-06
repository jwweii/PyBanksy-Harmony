## 2025-08-06 Simon
Modified from Evan's Banksy pipeline to speed up anndata generation and streamline the process
https://github.com/jwweii/PyBanksy-Harmony/tree/main?tab=readme-ov-file
- Allow multitasking in the first step for faster AnnData generation
- Streamlined workflow using 1 master script ``Run_All.sh`` run 3 workflows in ``workflow2_subworkflow`` automatically. 

## To run
- Copy 2 files in ``workflow2_master`` folder. 
1. Modify paths and parameters in the ``Run_All.sh`` file
2. Modify the ``Xenium_tracking.csv`` to include sample info to run.

The run usually takes a long time, so it's recommended to start a tmux session and then  
 ``bash Run_All.sh`` to start. 
- Check log files for any initial error messages

## Note 
Currently referencing script as absolute path to the katmai folder - aka will not work out-of-box if copied to different places
