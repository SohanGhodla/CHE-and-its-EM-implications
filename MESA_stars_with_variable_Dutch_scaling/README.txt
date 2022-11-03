This is to create evolved (till core carbon depletion) CHE models.

The main file that creates the MESA model grid is located in the directory named **Symbolic_link** and is called **Make_models.py**.

**XX-clean_script** is the directory that is cloned everytime (with a variable input for mass, metallicity, rotation rate etc) to create a grid.
The resulting data is stored in the directory named **Makemodels**.

The **pick_LOGS_files.py** file can be used to remove unnecessary files from the **Makemodels** directory (one's the directory is created). This would 
only leave the LOGS folder (containing the history and profile data) and remove the rest.

The **Symbolic_link** directory contains Symbolic links that point to the models in the **Makemodels** directory.

To schedule a job on a cluster, one would need to execute the **run_mesaarray.sh** file. This would run the models in the **Makemodels** directory.
The **Saved_ZAMS_models** file contains pre-evoled zero age main sequence models that would be required by the **Make_models.py** file. Note that this needs uncompressed before scheduling the job.


The directory **create main sequence models** creates a larger grid of models (including CHE and non CHE) with variable rotation rate, mass and metallicity. This helps to generate the data required for Fig. 4. It would also need the ZAMS models stored in The **Saved_ZAMS_models** file, so please copy them into the 'create main sequence models** directory.