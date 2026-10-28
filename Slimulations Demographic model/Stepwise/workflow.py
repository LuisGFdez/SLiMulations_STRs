
import os
from gwf import Workflow, AnonymousTarget
import pandas as pd
import numpy as np

gwf = Workflow()


def simulate(input,trees):
    inputs = []
    outputs = trees
    options = {
    	'cores': 1,
        'memory': '60g',
        'walltime': '40:00:00',
        'queue': 'memorylong',
        'account': 'rohlfslab',
    }
    spec = f'''
    

    module load python3/3.11.4
    make -f {input} -j 3

    
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
def newick_simplify_vcf_remove(tree,vcf,mutation,genotype,genotype_sorted,str_positions,popsize,species):
    inputs = [tree]
    outputs = [vcf,genotype,genotype_sorted,str_positions]
    options = {
        'cores': 1,
        'memory': '5g',
        'walltime': '4:00:00',
        'queue': 'memorylong',
        'account': 'rohlfslab',
    }
    spec = f'''
    module load python3/3.11.4
    python3 modify_trees_get_vcf_infer_arg.py {tree} {mutation} {popsize} {species}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def new_treva_smm(mutation,str_genotype,samples,overall_likelihood):
   
    inputs = [str_genotype]
    outputs = [overall_likelihood]
    options = {
        'cores': 1,
        'memory': '30g',
        'walltime': '50:00:00',
        'queue': 'memorylong',
        'account': 'rohlfslab',
    }
    spec = f'''
    export LC_ALL=C
    export LANG=C
    module load R
    Rscript new_treva_SMM.R  {mutation} {str_genotype} {samples}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
def add_truth_values(truth_mu,samples,species,overall_file):
    mu=float(truth_mu)
    samples=int(popsize)
    species=str(species)
    #output_file=f"{truth_mu}/samples_{samples}/Overall_best_likelihood_2.tsv"
    base, ext = os.path.splitext(overall_file)
    output_file = f'{base}_2{ext}'
    inputs = overall_file
    #print(inputs)
    outputs = [output_file]
    #print(outputs)
    options = {
        'cores': 1,
        'memory': '1g',
        'walltime': '10:00',
        'queue': 'memorylong',
        'account': 'rohlfslab',
    }
    spec = f'''
     awk -v mu="{mu}" -v samples="{samples}" -v species="{species}" 'BEGIN {{OFS="\\t"}} NR==1 {{print $0, "truth_mu",  "samples"; next}} {{print $0, mu, samples}}' {inputs} > {output_file}
    
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
def cat_likelihood_file(merge_files_smm):
 
    sim_likelihoods_smmm = gwf.glob(f"0.0*/*/Overall_best_likelihood_2*")

    inputs = [merge_files_smm]
    outputs = ["sim_likelihoods_smm.tsv" ]
    options = {
        'cores': 1,
        'memory': '1g',
        'walltime': '10:00',
        'queue': 'memory',
        'account': 'rohlfslab',
    }
    ##awk 'FNR==1 && NR!=1 {{next}} {{print}}' {' '.join(sim_likelihoods_tpm)} > sim_likelihoods_tpm.tsv
    ####awk 'FNR==1 && NR!=1 {{next}} {{print}}' {' '.join(inf_likelihoods_smm)} > inf_likelihoods_smm.tsv   
    spec = f'''
    	 awk 'FNR==1 && NR!=1 {{next}} {{print}}' {' '.join(sim_likelihoods_smmm)} > sim_likelihoods_smm.tsv
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
df = pd.read_csv("/home/luigf/rohlfslab/great_apes_str_simulations/metadata/edges_meta.tsv", sep="\t")
df = df.fillna('')
df["infile"] = df.parent + ".trees"
df["outfile"] = df.edge + ".trees"
df.loc[df["infile"]==".trees", "infile"] = ""
df=df.rename({"N":"popsize","edge":"child"},axis=1)
df["is_leaf"] = ~df.child.isin(df.parent)
df=df.loc[8:18,:]

dad= " "
df.loc[:8,["parent",]]=dad
df.loc[:8,["infile",]]=dad
##rename columns
df=df.rename({"N":"popsize","edge":"child"},axis=1)
###### Testing values for popsize and generations
##pop_size=[500,60,150,40,90,90,90,90,90,90,90]
##edges=[2000,250,200,50,50,50,50,50,50,50,50]

##set size testing values 
##old_values=df.loc[:18,"popsize"].tolist()
##df.popsize=df.loc[:,["popsize",]].replace(old_values,pop_size)

##generations testing values
##old_time_values=df.loc[:18,"gens"].tolist()
##df.gens=df.loc[:,["gens",]].replace(old_time_values,edges)

###########IMPORTANT:Pay attention here###
####MODIFY "l" values based on the number of decimal places in the mutation rate
l=4
merge_files_smm=[]
for mutation in [1e-5,1e-6]:
    l+=1
    mutation=format(float(mutation), f'.{l}f')
    file=f'{mutation}.make'
    f = open(file, "w")
    tree_outputs = [f"{mutation}/{outfile}" for outfile in df.outfile]
    print(f"all: {' '.join(tree_outputs)}\n", file=f)
    for i, row in df.iterrows():
        outfile = f"{mutation}/{row.outfile}"
        infile = f"{mutation}/{row.infile}" if row.infile.strip() else row.infile
        print(f"{outfile}:{infile} paralel.slim", file=f)
        print(f"\t~/.conda/envs/simulations/bin/slim -d \"infile='{infile}'\" -d popsize={row.popsize} "
              f"-d \"popname=\'{row.child}\'\" "
              f"-d num_gens={row.gens} " f"-d mutation_rate={mutation} "  f"-d \"folder='{mutation}'\" "  f"-d \"outfile='{outfile}'\" "  "paralel.slim\n",file=f)
    f.close()

    jobname = f'neutrality_{mutation}_sim'
    os.makedirs(f"{mutation}", exist_ok=True)
    gwf.target_from_template(jobname, simulate(file,tree_outputs))
    for idx,species in df['child'].items():
            popsize = df.loc[idx, 'popsize']
            tree=f"{mutation}/{species}.trees"
        ###create specific dirs
            os.makedirs(f"{mutation}/{species}", exist_ok=True)
            os.makedirs(f"{mutation}/{species}/vcf", exist_ok=True)
            os.makedirs(f"{mutation}/{species}/infered_likelihoods_smm", exist_ok=True)

        ##define job names
            
            jobname2= f'newick{mutation}_{species}'
            jobname5=f'new_treva_{mutation}_{species}'
            jobname8=f'add_col_names_{mutation}_{species}'
        ##defining outputs
            vcf_file=f"{mutation}/{species}/vcf/{mutation}_new.vcf"
            str_genotype=f"{mutation}/{species}/{species}.csv"
            str_genotype_sorted=f"{mutation}/{species}/{species}_sorted.csv"
            str_positions=f"{mutation}/{species}/positions.txt"

            best_likelihood_file_new=f"{mutation}/{species}/Best_likelihoods_smm_{mutation}.tsv"
            overall_likelihood=f"{mutation}/{species}/Overall_best_likelihood.tsv"
            merge_infered_likelihoods_smm=f"{mutation}/{species}/infered_likelihoods_merge_smm_{mutation}.tsv"

        ###targets
            gwf.target_from_template(jobname2, newick_simplify_vcf_remove(tree,vcf_file,mutation,str_genotype,str_genotype_sorted,str_positions,popsize,species))
            ###treva target smm
            new_model=gwf.target_from_template(jobname5, new_treva_smm(mutation,str_genotype,popsize,overall_likelihood))
            like_file_smm=gwf.target_from_template(jobname8, add_truth_values(mutation,popsize,species,overall_likelihood))
            merge_files_smm.extend(like_file_smm.outputs)


gwf.target_from_template(
    name="cat_likelihoods_file",
    template=cat_likelihood_file(merge_files_smm))
    





