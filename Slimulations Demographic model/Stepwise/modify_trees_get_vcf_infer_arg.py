import tskit
import msprime
import numpy as np
import pandas as pd
import pyslim
import os
import subprocess
import sys
import argparse
import io
import msprime
import gc
import pyslim
import glob

tree_file = str(sys.argv[1]) 
print(tree_file)
mutation = str(sys.argv[2]) 
print(mutation)
samples = int(sys.argv[3]) 
print(samples)
species = str(sys.argv[4])
print(species)

tree_nodes=samples*2
print(tree_nodes)


def get_m2_positions(ts_original):
    mutation_table = ts_original.tables.mutations
    site_table = ts_original.tables.sites

    # Extract positions where mutation_type == 2
    positions = []

    for mutation in mutation_table:
     ###get the metadat
        metadata = mutation.metadata
        if len(metadata["mutation_list"])>0:
           if metadata['mutation_list'][-1]["mutation_type"] == 2:
              site_index = mutation.site  # Get site index
              position = site_table.position[site_index]  # Get genomic position
              positions.append(position)
            
    pos=np.unique(positions)
    return pos
    
def remove_str_mutations(ts):
    
    ###remove str mutations
    ts = pyslim.update(ts)
    tables = ts.dump_tables()
    tables.mutations.clear()
    for mut in ts.mutations():
        if mut.derived_state !='':
            if mut.metadata['mutation_list'][-1]["mutation_type"] != 2:
                 tables.mutations.append(mut)
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    tables.sort()
    new_ts = tables.tree_sequence()
    assert ts.num_mutations-new_ts.num_mutations, "No mutations were removed!"
    print(f"# of mutations before: {ts_sim.num_mutations}, # of muts after: {new_ts.num_mutations}, # of removed: {ts_sim.num_mutations-new_ts.num_mutations}", flush=True)
    return new_ts

def recapitate_simplify_tree(ts_original):
    rts = pyslim.recapitate(ts_original,
            recombination_rate=1e-8,
            ancestral_Ne=10000, random_seed=3)
    rng = np.random.default_rng(seed=3)
    alive_inds = pyslim.individuals_alive_at(rts, 0)
    keep_indivs = rng.choice(alive_inds, samples, replace=False)
    keep_nodes = []
    for i in keep_indivs:
    	keep_nodes.extend(rts.individual(i).nodes)
    sts = rts.simplify(keep_nodes, keep_input_roots=True)
    ids=sts.metadata["SLiM"]["user_metadata"]["ids"]
    strs=sts.metadata["SLiM"]["user_metadata"]["tags"]
    df = pd.DataFrame(list(zip(ids, strs)), columns=["id", "allele"])
    return sts,df
 
def get_vcf_file(sts,mutation):
   
    nts=pyslim.update(sts)
    nts = pyslim.generate_nucleotides(nts)
    nts = pyslim.convert_alleles(nts)
    with open(f"{mutation}/{species}/vcf/{mutation}_new.vcf", "w") as vcffile:
        nts.write_vcf(vcffile)

def create_newick(positions,sts,df):
    allele_dict = {}
    allele_dict_sorted={}
    cont=0
    variant=tskit.Variant(sts)
    individuals=sts.tables.nodes[0:tree_nodes].individual 
    sorted_indices = np.argsort(individuals,kind="stable")
    for pos in positions:
    	cont+=1
    	key = f"STR{cont}"
    	tree=sts.at(pos)
    	output_file=f"{mutation}/{species}/STR{cont}_newick.txt" 
    	with open(output_file, "w") as file_newick:
            file_newick.write(tree.as_newick(precision=10,include_branch_lengths=True) + "\n")
    	str_alleles=[]
    	str_alleles_sorted=[]
    	variant.decode(sts.site(position=pos).id)

    	genotype=variant.genotypes
    	sorted_genotypes=genotype[sorted_indices]
    	alleles=variant.alleles

    	for gt in genotype:
            mutation_id=alleles[gt]
            #print(type(mutation_id))
            if len(mutation_id.split(","))>1:
            	mutation_id=mutation_id.split(",")[0]
            #print(mutation_id)
            str_alleles.append(df.loc[df["id"] == int(mutation_id), "allele"].values)
    	allele_dict[key] =str_alleles
    	#####sorted
    	for gt in sorted_genotypes:
            mutation_id=alleles[gt]
            #print(type(mutation_id))
            if len(mutation_id.split(","))>1:
            	mutation_id=mutation_id.split(",")[0]
            #print(mutation_id)
            str_alleles_sorted.append(df.loc[df["id"] == int(mutation_id), "allele"].values)
    	allele_dict_sorted[key] =str_alleles_sorted
        #print(var)
    df_strs = pd.DataFrame.from_dict(allele_dict)
    df_strs.to_csv(f"{mutation}/{species}/{species}.csv",index=False)
    df_strs_sorted = pd.DataFrame.from_dict(allele_dict_sorted)
    df_strs_sorted.to_csv(f"{mutation}/{species}/{species}_sorted.csv",index=False)


##first get the str positions to create the newick format and get the str genotype of each genome 
ts_sim=tskit.load(tree_file)
positions=get_m2_positions(ts_sim)
print(positions)
np.savetxt(f"{mutation}/{species}/positions.txt", positions, delimiter=" ",fmt="%.6f")
##then simplify the tree to run treva with small sample size
new_ts,df=recapitate_simplify_tree(ts_sim)
##third step: create newick format and get the str genotype of every genome:
create_newick(positions,new_ts,df)

###now get vcf file without multiallelic loci for arg inference
str_tree=remove_str_mutations(new_ts)

get_vcf_file(str_tree,mutation)

str_tree.dump(f"{mutation}/{species}/{mutation}_simplify.trees")



