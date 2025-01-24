# SLiMulations_STRs
## Getting Started


## 1. **Install Dependencies and Packages**
This simulation works by implementing **gwf**, a tool for building workflows and automating pipelines.  
To install [gwf](https://gwf.app/), use **conda** as described in the official documentation:

```bash
conda install -c conda-forge gwf or
conda install gwf
```
After installed, if you’re running gwf on a cluster you may want to use a backend that can submit targets to your clusters’ queueing system, run the command:

```bash
gwf config set backend slurm
gwf config set backend sge

```

Once configured, within the **Workflow.py** script, the  **simulate** function must be modified to adjust the **option** parameters according to the requirements of the cluster you are using, either **SGE** or **SLURM**.
Detailed information  configuring the option parameters based on the cluster type can be found here: [GWF Backends](https://gwf.app/reference/backends/)

Install [SLiM](https://messerlab.org/slim/) in the cluster or in your local machine
## 2. **Code Overview**

This repository contains the following scripts:

1. **`Workflow.py`**  
   - A Python script implementing `gwf` to optimize the number of simulations.  
   - This script is designed to efficiently manage and execute simulation workflows.

2. **`paralel.slim`**  
   - A SLiM script to perform simulations under neutrality based on tandem repeat mutation models.  
   - It simulates mutation processes for short tandem repeats (STRs).

3. **`transitionMatrixMSprimeMethod.py`**  
   - A Python script based on the stepwise mutation model to simulate expansions and contractions of STRs.  
   - This script generates a **transition matrix** filled with probabilities of stepwise mutations, capturing the evolutionary dynamics of STRs.

  
## 2. **How to run the simulations**

Once gwf is installed and configured with its specific settings, run the following commands:

```bash
gwf status
```
then
```bash
gwf run
```

