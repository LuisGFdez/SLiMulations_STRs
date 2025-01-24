
from gwf import Workflow, AnonymousTarget
import pandas as pd
import numpy as np

gwf = Workflow()


def simulate(input,output):
    inputs = []
    outputs = [output]
    ##CHANGE the option bracket according to your cluster requirements
    options = {
    	'cores': 1,
        'memory': '60g',
        'walltime': '40:00:00',
        'queue': 'rohlfs',
        'account': 'rohlfslab',
    }
    spec = f'''
    
    mkdir {output}
    module load python3/3.11.4
    make -f {input} 
    
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




popsize=100
popname="Humans"
num_gens=50000
l=1
for mutation in [1e-2,1e-3,1e-4,1e-5]:
    l+=1
    for STR in [5,15,35,55,75,100]:
        mutation=format(float(mutation), f'.{l}f')
        file=f'{mutation}_{STR}.make'
        f = open(file, "w")
        print(f"all:{mutation}_{STR}.trees\n", file=f)
        print(f"{mutation}_{STR}.trees:  paralel.slim", file=f)
        print(f"\t~/.conda/envs/simulations/bin/slim -d \"infile=' '\" -d popsize={popsize} "
                f"-d \"popname=\'{popname}\'\" "
                f"-d num_gens={num_gens} " f"-d mutation_rate={mutation} " f"-d str_count={STR} " f"-d \"folder='{mutation}_{STR}'\" "  f"-d \"outfile='{mutation}_{STR}/{mutation}_{STR}.trees'\" "  
                "paralel.slim\n",
                file=f)
        f.close()
        output=f'{mutation}_{STR}'
        jobname = f'simulate_{mutation}_{STR}'
        gwf.target_from_template(jobname, simulate(file,output))





