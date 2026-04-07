import pandas
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='Directory for model')
    parser.add_argument('gene_path', metavar='gene-path', type=str,
                        help='Genes to KO')
    parser.add_argument('--output', 
                        type=str,
                        help='Name of file with disabled genes',
                        default=None)
    parser.add_argument('--prevent_reaction_removal', action='store_true',
                        help='Prevent disabling reactions from metabolic network')

    args = parser.parse_args()

    #model_dir=args.model_dir
    #model_dir="../Athal_model__unspecific"
    #gene_path=args.gene_path
    #gene_path="../Genes_to_remove_RNAseq/RBAgene_not_in_rnaseq_with_threshold_6.7942.csv"

    genes= list(pandas.read_csv(args.gene_path,sep="\t",header=None).iloc[:,0])
    protein_curation_file=pandas.read_csv("{}/data/helper_files/protein_curation.tsv".format(args.model_dir),sep="\t")

    if not args.prevent_reaction_removal:
        for gene in genes:
            if gene in list(protein_curation_file["GENE"]):
                protein_curation_file.loc[protein_curation_file["GENE"]==gene,"DISABLED"]=1
    else:
        for gene in genes:
            if gene in list(protein_curation_file["GENE"]):
                rxns_associated=list(protein_curation_file.loc[protein_curation_file["GENE"]==gene,"REACTION"])
                for rxn in rxns_associated:
                    subunit=protein_curation_file.loc[(protein_curation_file["GENE"]==gene)&(protein_curation_file["REACTION"]==rxn),"SUBUNIT"].values[0]
                    disable=False
                    for su_alternative in list(set(list(protein_curation_file.loc[(protein_curation_file["REACTION"]==rxn)&(protein_curation_file["SUBUNIT"]==subunit),"GENE"]))):
                        if su_alternative not in genes:
                            disable=True
                    if disable:
                        protein_curation_file.loc[(protein_curation_file["GENE"]==gene)&(protein_curation_file["REACTION"]==rxn),"DISABLED"]=1
                    else:
                        print("WARNING: gene-reaction association {}-{} not disabled, due to avoiding removal of all isenzymes".format(gene,rxn))

    if args.output is not None:
        protein_curation_file.to_csv("{}/data/helper_files/{}".format(args.model_dir,args.output),sep="\t",index=None)
    else:
        protein_curation_file.to_csv("{}/data/helper_files/protein_curation.tsv".format(args.model_dir),sep="\t",index=None)


if __name__ == '__main__':
    main()

#python disable_genes.py ../Athal_model__unspecific ../Genes_to_remove_RNAseq/RBAgene_not_in_rnaseq_with_threshold_6.7942.csv protein_curation_KO3.tsv