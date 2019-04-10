import pandas as pd
import argparse
import copy
import scipy.stats as stats
import math
import sys
import multiprocessing
from datetime import datetime

shared_columns = [
        'Hugo_Symbol',
        'Chromosome',
        'Start_position',
        'End_position',
        'Reference_Allele',
        'Tumor_Seq_Allele1',
        'Tumor_Seq_Allele2',
        'Variant_Classification',
        'Variant_Type',
        'Genome_Change',
        'Codon_Change',
        'Protein_Change',
        'Tumor_Sample_Barcode',
        'COSMIC_tissue_types_affected',
        'MUTSIG_Published_Results',
        'Transcript_Exon',
        'Other_Transcripts',
        'Entrez_Gene_Id',
        'DNARepairGenes_Role',
        'Transcript_Strand',
        'Strand',
        'cDNA_Change',
        'Annotation_Transcript',
        'GO_Cellular_Component',
        'Transcript_Position',
        'GO_Biological_Process',
        'GO_Molecular_Function',
        'COSMIC_total_alterations_in_gene',
        'gc_content'
]

# CONSTANTS
LOW_DNA_TUMOR_COVERAGE = "low_dna_tumor_coverage"
LOW_DNA_GERMLINE_COVERAGE = "low_dna_germline_coverage"
LOW_RNA_TUMOR_COVERAGE = "low_rna_tumor_coverage"
LOW_RNA_GERMLINE_COVERAGE = "low_rna_germline_coverage"
LOW_COVERAGE = "low_coverage"
ARTIFACT = "artifact"

POTENTIAL_LOH = "potential_loh"
TUMOR_ONLY_VSL = "tumor_only_vsl"
TUMOR_ONLY_VSE = "tumor_only_vse"
SOMATIC = "somatic"
GENERAL_VSE = "general_vse"
GENERAL_RNA_EDITING = "general_rna_editing"
TUMOR_SPECIFIC_RNA_EDITING = "tumor_specific_rna_editing"
GENERAL_VSL = "general_vsl"
TUMOR_SPECIFIC_VSL = "tumor_specific_vsl"
TUMOR_SPECIFIC_VSE = "tumor_specific_vse"
GERMLINE = "germline"
GERMLINE_SPECIFIC_VSL = "germline_specific_vsl"
GERMLINE_VSL_TUMOR_VSE = "germline_vsl_tumor_vse"
GERMLINE_VSE_TUMOR_VSL = "germline_vse_tumor_vsl"
GERMLINE_SPECIFIC_VSE = "germline_specific_vse"



all_categories = [
    POTENTIAL_LOH,
    TUMOR_ONLY_VSL,
    TUMOR_ONLY_VSE,
    SOMATIC,
    GENERAL_VSE,
    GENERAL_RNA_EDITING,
    TUMOR_SPECIFIC_RNA_EDITING,
    GENERAL_VSL,
    TUMOR_SPECIFIC_VSL,
    TUMOR_SPECIFIC_VSE,
    GERMLINE,
    GERMLINE_SPECIFIC_VSL,
    GERMLINE_VSL_TUMOR_VSE,
    GERMLINE_VSE_TUMOR_VSL,
    GERMLINE_SPECIFIC_VSE
]


def make_index(r):
    """Generate a unique index for a variant from a few key column values."""
    index = '{}:{}:{}-{}({}/{})'.format(r.Hugo_Symbol,
                                        r.Chromosome,
                                        int(r.Start_position),
                                        int(r.End_position),
                                        r.Reference_Allele,
                                        r.Tumor_Seq_Allele2)
    return index


def make_clinvar_index(r):
    """Generate a unique index for a clinvar variant that can be used to join with our variants"""
    index = '{}:{}:{}-{}({}/{})'.format(r.GeneSymbol,
                                        r.Chromosome,
                                        int(r.Start),
                                        int(r.Stop),
                                        r.ReferenceAllele,
                                        r.AlternateAllele)
    return index


def fishers_exact_two_sided_dt_dg(r):
    """Use a two-sided Fisher's Exact test to check whether the MAFs for DNA tumor and DNA normal at each site come
    from the same underlying distribution or not."""
    for v in [r.t_alt_count_dt, r.t_ref_count_dt, r.t_alt_count_dg, r.t_ref_count_dg]:
        if math.isnan(v):
            return None

    oddsratio, pvalue = stats.fisher_exact([[r.t_alt_count_dt, r.t_ref_count_dt],
                                            [r.t_alt_count_dg, r.t_ref_count_dg]],
                                           alternative='two-sided')
    return pvalue


def fishers_exact_two_sided_rt_rg(r):
    """Use a two-sided Fisher's Exact test to check whether the MAFs for RNA tumor and RNA normal at each site come
    from the same underlying distribution or not."""
    for v in [r.t_alt_count_rt, r.t_ref_count_rt, r.t_alt_count_rg, r.t_ref_count_rg]:
        if math.isnan(v):
            return None

    oddsratio, pvalue = stats.fisher_exact([[r.t_alt_count_rt, r.t_ref_count_rt],
                                            [r.t_alt_count_rg, r.t_ref_count_rg]],
                                           alternative='two-sided')
    return pvalue


def fishers_exact_one_sided_dg_rg(r, alternative):
    """Use a one-sided Fisher's Exact DNA germline vs. RNA germline"""
    for v in [r.t_alt_count_dg, r.t_ref_count_dg, r.t_alt_count_rg, r.t_ref_count_rg]:
        if math.isnan(v):
            return None

    oddsratio, pvalue = stats.fisher_exact([[r.t_alt_count_dg, r.t_ref_count_dg],
                                            [r.t_alt_count_rg, r.t_ref_count_rg]],
                                           alternative=alternative)

    return pvalue


def fishers_exact_one_sided_dt_rt(r, alternative):
    """Use a one-sided Fisher's Exact DNA tumor vs. RNA tumor"""
    for v in [r.t_alt_count_dt, r.t_ref_count_dt, r.t_alt_count_rt, r.t_ref_count_rt]:
        if math.isnan(v):
            return None

    oddsratio, pvalue = stats.fisher_exact([[r.t_alt_count_dt, r.t_ref_count_dt],
                                                            [r.t_alt_count_rt, r.t_ref_count_rt]],
                                           alternative=alternative)
    return pvalue


def fishers_exact_one_sided_rg_rt(r, alternative):
    """Use a one-sided Fisher's Exact RNA tumor vs RNA normal"""
    for v in [r.t_alt_count_rg, r.t_ref_count_rg, r.t_alt_count_rt, r.t_ref_count_rt]:
        if math.isnan(v):
            return None

    oddsratio, pvalue = stats.fisher_exact([[r.t_alt_count_rg, r.t_ref_count_rg],
                                            [r.t_alt_count_rt, r.t_ref_count_rt]],
                                           alternative=alternative)
    return pvalue


def fishers_exact_one_sided_dt_dg(r, alternative):
    """Use a one-sided Fisher's Exact DNA tumor vs. DNA normal"""
    for v in [r.t_alt_count_dt, r.t_ref_count_dt, r.t_alt_count_dg, r.t_ref_count_dg]:
        if math.isnan(v):
            return None

    oddsratio, pvalue = stats.fisher_exact([[r.t_alt_count_dt, r.t_ref_count_dt],
                                            [r.t_alt_count_dg, r.t_ref_count_dg]],
                                           alternative=alternative)
    return pvalue


def combine_dfs(dg=None, dt=None, rg=None, rt=None):
    """Join all the allele counts from each sample type into one dataframe."""
    altered_dfs = {}
    all_columns = shared_columns + ['dt_i_judgement', 'rt_i_judgement']
    column_suffixes = ['dg', 'dt', 'rg', 'rt']
    for dataframe, column_suffix in zip([dg, dt, rg, rt], column_suffixes):
        if dataframe is not None:
            dataframe.index = dataframe.apply(make_index, axis=1)
            dataframe = dataframe.rename(index=str,
                                         columns=
                                         {
                                            't_alt_count': 't_alt_count_{}'.format(column_suffix),
                                            't_ref_count': 't_ref_count_{}'.format(column_suffix)
                                         })

            altered_dfs[column_suffix] = copy.copy(dataframe)

    joined_df = altered_dfs['dg']\
        .join(altered_dfs['dt'], rsuffix='_dt', how='inner')\
        .join(altered_dfs['rt'], rsuffix='_rt', how='inner')\
        .join(altered_dfs['rg'], rsuffix='_rg', how='inner')

    joined_df['dna_germline_af'] = joined_df.t_alt_count_dg.divide(joined_df.t_alt_count_dg + joined_df.t_ref_count_dg)
    joined_df['dna_tumor_af'] = joined_df.t_alt_count_dt.divide(joined_df.t_alt_count_dt + joined_df.t_ref_count_dt)
    joined_df['rna_germline_af'] = joined_df.t_alt_count_rg.divide(joined_df.t_alt_count_rg + joined_df.t_ref_count_rg)
    joined_df['rna_tumor_af'] = joined_df.t_alt_count_rt.divide(joined_df.t_alt_count_rt + joined_df.t_ref_count_rt)

    final_cols = joined_df[all_columns + ['t_alt_count_dg', 't_ref_count_dg', 'dna_germline_af',
                                          't_alt_count_dt', 't_ref_count_dt', 'dna_tumor_af',
                                          't_alt_count_rg', 't_ref_count_rg', 'rna_germline_af',
                                          't_alt_count_rt', 't_ref_count_rt', 'rna_tumor_af']]
    return final_cols


def run_fishers_exact_tests(combined_df):
    # Two-sided Fisher's Exact Test on tumor/normal DNA comparisons
    combined_df['dt_dg_test'] = combined_df.apply(fishers_exact_two_sided_dt_dg, axis=1)
    combined_df['rt_rg_test'] = combined_df.apply(fishers_exact_two_sided_rt_rg, axis=1)

    # One-sided Fisher's Exact Tests on DNA/RNA comparisons
    combined_df['dt_dg_test_greater'] = combined_df.apply(fishers_exact_one_sided_dt_dg, args=('greater',), axis=1)
    combined_df['rg_rt_test_greater'] = combined_df.apply(fishers_exact_one_sided_rg_rt, args=('greater',), axis=1)
    combined_df['rg_rt_test_less'] = combined_df.apply(fishers_exact_one_sided_rg_rt, args=('less',), axis=1)
    combined_df['dg_rg_test_less'] = combined_df.apply(fishers_exact_one_sided_dg_rg, args=('less',), axis=1)
    combined_df['dg_rg_test_greater'] = combined_df.apply(fishers_exact_one_sided_dg_rg, args=('greater',), axis=1)
    combined_df['dt_rt_test_less'] = combined_df.apply(fishers_exact_one_sided_dt_rt, args=('less',), axis=1)
    combined_df['dt_rt_test_greater'] = combined_df.apply(fishers_exact_one_sided_dt_rt, args=('greater',), axis=1)

    return combined_df


def categorize(combined_df):
    """Based on a set of comparisons of rna and dna tumor and normal read counts, categorize each variant"""
    sys.stdout.write("Categorizing {} variants\n".format(len(combined_df)))

    combined_df = run_fishers_exact_tests(combined_df)

    # Alt count is 0
    dna_alt_count_zero = (combined_df['t_alt_count_dg'] == 0) & (combined_df['t_alt_count_dt'] == 0)
    dna_germline_alt_count_zero = (combined_df['t_alt_count_dg'] == 0)
    rna_germline_alt_count_zero = (combined_df['t_alt_count_rg'] == 0)

    # Artifact flag if alt count is 0 across all four samples
    artifact_flag = (combined_df['t_alt_count_dg'] == 0) \
                    & (combined_df['t_alt_count_dg'] == 0) \
                    & (combined_df['t_alt_count_rg'] == 0) \
                    & (combined_df['t_alt_count_rt'] == 0)

    dna_tumor_same_dna_germline = combined_df['dt_dg_test'] >= 0.05
    dna_tumor_different_dna_germline = combined_df['dt_dg_test'] < 0.05
    dna_tumor_gt_dna_germline = combined_df['dt_dg_test_greater'] < 0.05
    rna_germline_gt_rna_tumor = combined_df['rg_rt_test_greater'] < 0.05
    rna_germline_lt_rna_tumor = combined_df['rg_rt_test_less'] < 0.05

    combined_df['rna_germline_gt_rna_tumor'] = rna_germline_gt_rna_tumor
    combined_df['rna_germline_lt_rna_tumor'] = rna_germline_lt_rna_tumor
    combined_df['dna_tumor_greater_dna_germline'] = dna_tumor_gt_dna_germline
    combined_df['dna_tumor_different_dna_germline'] = dna_tumor_different_dna_germline
    combined_df['dna_tumor_same_dna_germline'] = dna_tumor_same_dna_germline

    flagged_as_somatic = (combined_df['dt_i_judgement'] == 'KEEP')

    dg_rg_less_passed = (combined_df['dg_rg_test_less'] < 0.05)
    dg_rg_less_failed = (combined_df['dg_rg_test_less'] >= 0.05)

    dg_rg_greater_passed = (combined_df['dg_rg_test_greater'] < 0.05)
    dg_rg_greater_failed = (combined_df['dg_rg_test_greater'] >= 0.05)

    dt_rt_less_passed = (combined_df['dt_rt_test_less'] < 0.05)
    dt_rt_less_failed = (combined_df['dt_rt_test_less'] >= 0.05)

    dt_rt_greater_passed = (combined_df['dt_rt_test_greater'] < 0.05)
    dt_rt_greater_failed = (combined_df['dt_rt_test_greater'] >= 0.05)

    # Compare dna tumor and dna germline
    combined_df['dna_tumor_gt_dna_germline'] = dna_tumor_gt_dna_germline

    # Compare germline dna and germline rna
    dna_germline_lt_rna_germline = dg_rg_less_passed
    dna_germline_gt_rna_germline = dg_rg_greater_passed
    dna_germline_same_rna_germline = dg_rg_less_failed & dg_rg_greater_failed

    combined_df['dna_germline_lt_rna_germline'] = dna_germline_lt_rna_germline
    combined_df['dna_germline_gt_rna_germline'] = dna_germline_gt_rna_germline
    combined_df['dna_germline_same_rna_germline'] = dna_germline_same_rna_germline

    # Compare tumor dna and tumor rna
    dna_tumor_lt_rna_tumor = dt_rt_less_passed
    dna_tumor_gt_rna_tumor = dt_rt_greater_passed
    dna_tumor_same_rna_tumor = dt_rt_less_failed & dt_rt_greater_failed
    combined_df['dna_tumor_same_rna_tumor'] = dna_tumor_same_rna_tumor
    combined_df['dna_tumor_gt_rna_tumor'] = dna_tumor_gt_rna_tumor
    combined_df['dna_tumor_lt_rna_tumor'] = dna_tumor_lt_rna_tumor

    # Left side of graph

    combined_df[TUMOR_ONLY_VSL] = dna_tumor_different_dna_germline \
                                  & flagged_as_somatic \
                                  & dna_tumor_gt_rna_tumor

    combined_df[TUMOR_ONLY_VSE] = dna_tumor_different_dna_germline \
                                  & flagged_as_somatic \
                                  & dna_tumor_lt_rna_tumor

    combined_df[SOMATIC] = dna_tumor_gt_dna_germline & flagged_as_somatic & dna_tumor_same_rna_tumor

    # Right side of graph
    combined_df[GENERAL_VSE] = dna_tumor_same_dna_germline \
                               & dna_germline_lt_rna_germline \
                               & dna_tumor_lt_rna_tumor \
                               & ~dna_alt_count_zero

    combined_df[GENERAL_RNA_EDITING] = dna_germline_lt_rna_germline \
                                       & dna_tumor_lt_rna_tumor\
                                       & dna_alt_count_zero

    combined_df[TUMOR_SPECIFIC_RNA_EDITING] = dna_tumor_lt_rna_tumor \
                                              & dna_alt_count_zero \
                                              & rna_germline_alt_count_zero

    combined_df[GENERAL_VSL] = dna_tumor_same_dna_germline \
                               & dna_germline_gt_rna_germline\
                               & dna_tumor_gt_rna_tumor

    combined_df[TUMOR_SPECIFIC_VSL] = (dna_tumor_same_dna_germline | dna_tumor_gt_dna_germline) \
                                      & dna_germline_same_rna_germline \
                                      & dna_tumor_gt_rna_tumor \
                                      & rna_germline_gt_rna_tumor

    combined_df[TUMOR_SPECIFIC_VSE] = ~dna_tumor_gt_dna_germline \
                                      & dna_germline_same_rna_germline \
                                      & dna_tumor_lt_rna_tumor \
                                      & ~dna_alt_count_zero \
                                      & rna_germline_lt_rna_tumor

    # potential_loh is on left side of graph but will be default when tumor_specific_vsl/vse are not true
    combined_df[POTENTIAL_LOH] = dna_tumor_different_dna_germline & ~flagged_as_somatic & ~combined_df[TUMOR_SPECIFIC_VSE] & ~combined_df[TUMOR_SPECIFIC_VSL]

    combined_df[GERMLINE] = ~dna_germline_alt_count_zero \
                            & dna_tumor_same_rna_tumor \
                            & dna_germline_same_rna_germline \
                            & dna_tumor_same_dna_germline \
                            & ~flagged_as_somatic

    # Less probable categories
    combined_df[GERMLINE_SPECIFIC_VSL] = dna_tumor_same_dna_germline\
                                         & dna_germline_gt_rna_germline \
                                         & dna_tumor_same_rna_tumor \
                                         & rna_germline_lt_rna_tumor

    combined_df[GERMLINE_VSL_TUMOR_VSE] = dna_tumor_same_dna_germline & dna_germline_gt_rna_germline \
                                            & dna_tumor_lt_rna_tumor & ~dna_tumor_same_rna_tumor \
                                          & ~dna_germline_same_rna_germline
    combined_df[GERMLINE_VSE_TUMOR_VSL] = dna_tumor_same_dna_germline & dna_germline_lt_rna_germline \
                                            & dna_tumor_gt_rna_tumor & ~dna_tumor_same_rna_tumor \
                                          & ~ dna_germline_same_rna_germline

    combined_df[GERMLINE_SPECIFIC_VSE] = dna_tumor_same_dna_germline \
                                         & dna_germline_lt_rna_germline \
                                         & ~dna_alt_count_zero \
                                         & dna_tumor_same_rna_tumor\
                                         & rna_germline_gt_rna_tumor

    # Add a new category for variants that just have very low (<5) dna tumor alt read count
    combined_df[LOW_DNA_TUMOR_COVERAGE] = combined_df.t_alt_count_dt.isin([0, 1, 2, 3, 4])
    combined_df[LOW_DNA_GERMLINE_COVERAGE] = combined_df.t_alt_count_dg.isin([0, 1, 2, 3, 4])
    combined_df[LOW_RNA_TUMOR_COVERAGE] = combined_df.t_alt_count_rt.isin([0, 1, 2, 3, 4])
    combined_df[LOW_RNA_GERMLINE_COVERAGE] = combined_df.t_alt_count_rg.isin([0, 1, 2, 3, 4])
    combined_df[LOW_COVERAGE] = (combined_df[LOW_DNA_TUMOR_COVERAGE]) \
                                | combined_df[LOW_DNA_GERMLINE_COVERAGE] \
                                | combined_df[LOW_RNA_TUMOR_COVERAGE] \
                                | combined_df[LOW_RNA_GERMLINE_COVERAGE]

    combined_df[ARTIFACT] = artifact_flag

    return combined_df


def number_of_trues(r):
    """Return the number of True values in the row"""
    return r.sum()


def first_valid_column_name(r):
    """In each row, which will contain several Falses and one True, return the first True column name"""
    idx = r.idxmax()
    # If there were no Trues at all in the given columns, there is no assignment so return None
    if not r[idx]:
        return None
    else:
        return idx


def infer_category(categorized_df):
    category_columns = categorized_df[all_categories]
    assignment = category_columns.apply(first_valid_column_name, axis=1)
    categorized_df['assignment'] = assignment
    categorized_df['num_categories'] = category_columns.apply(number_of_trues, axis=1)
    return categorized_df


def return_new_category_if_unassigned(r):
    """If assignment column is unassigned, figure out an alternative category for it. Otherwise just return its
    current assignment"""
    assignment = r.assignment
    if assignment:
        return assignment
    else:
        if r.t_alt_count_dg >= 5:
            return GERMLINE
        else:
            return 'uncategorized'


def reassign_unassigned(categorized_df):
    categorized_df['assignment'] = categorized_df.apply(return_new_category_if_unassigned, axis=1)
    return categorized_df


def summarize(categorized_df):
    """Cursory overview of how many variants were assigned to each category"""
    counts_for_categories = {}
    for category in all_categories:
        num_in_category = len(categorized_df[categorized_df[category] == True])
        counts_for_categories[category] = num_in_category
    sys.stdout.write("~~~~~~~~~~~~~~\n")
    for k,v in counts_for_categories.items():
        sys.stdout.write('{}: {}\n'.format(k,v))

    # Also figure out how many were assigned to multiple categories...
    sys.stdout.write("Sum of all values: {}\n".format(sum(counts_for_categories.values())))


def main():
    parser = argparse.ArgumentParser(description='Run R2D2')
    parser.add_argument('dna_germline', metavar='dna_germline', type=str)
    parser.add_argument('dna_tumor', metavar='dna_tumor', type=str)
    parser.add_argument('rna_germline', metavar='rna_germline', type=str)
    parser.add_argument('rna_tumor', metavar='rna_tumor', type=str)
    parser.add_argument('pair_id', metavar='pair_id', type=str)
    parser.add_argument('--num_threads', metavar='num_threads', type=int, default=8)
    args = parser.parse_args()

    dna_germline = args.dna_germline
    dna_tumor = args.dna_tumor
    rna_germline = args.rna_germline
    rna_tumor = args.rna_tumor
    pair_id = args.pair_id
    num_threads = args.num_threads

    dg = pd.read_csv(dna_germline, comment='#', sep='\t', encoding="ISO-8859-1")
    dt = pd.read_csv(dna_tumor, comment='#', sep='\t', encoding="ISO-8859-1")
    rg = pd.read_csv(rna_germline, comment='#', sep='\t', encoding="ISO-8859-1")
    rt = pd.read_csv(rna_tumor, comment='#', sep='\t', encoding="ISO-8859-1")

    combined_df = combine_dfs(dg=dg, dt=dt, rg=rg, rt=rt)

    # in docker image this should be changed to "/clinvar_variant_summary.lite.txt" (extra forward slash at front)
    clinvar = pd.read_csv('/clinvar_variant_summary.lite.txt', sep='\t', header='infer')
    clinvar.index = clinvar.apply(make_clinvar_index, axis=1)
    clinvar = clinvar.drop_duplicates()
    """
    Multi-processing
    """
    num_rows = len(combined_df)
    chunk_size = 800
    num_chunks = math.ceil(num_rows/chunk_size)
    chunks = []
    # Split up the dataframe into chunks of rows
    indices = combined_df.index
    i_range = range(0, num_chunks, 1)
    sys.stdout.write("Num chunks: {}\n".format(num_chunks))
    sys.stdout.write("Num rows in original dataframe: {}\n".format(num_rows))
    for i in i_range:
        sys.stdout.write("Subsetting chunk {}\n".format(i))
        chunk_indices = indices[i*chunk_size:(i+1)*chunk_size]
        chunk = combined_df.loc[chunk_indices]
        chunks.append(chunk)

    # Sanity check
    sys.stdout.write("Num total rows among all chunks being processed: {}\n".format(sum([len(c) for c in chunks])))

    sys.stdout.write("Using {} threads\n".format(num_threads))
    pool = multiprocessing.Pool(num_threads)

    start = datetime.now()

    outputs = pool.map(categorize, chunks)

    merged_categorized_df = pd.concat(outputs)

    merged_categorized_df = infer_category(merged_categorized_df)

    # For uncategorized ones, assign a backup category (mostly will be germline probably)
    merged_categorized_df = reassign_unassigned(merged_categorized_df)

    summarize(merged_categorized_df)

    # Join with clinvar
    merged_categorized_df = merged_categorized_df.join(clinvar, how='left', rsuffix='_clinvar')
    # Drop the unnecessary clinvar columns
    merged_categorized_df = merged_categorized_df.drop(['Chromosome_clinvar',
                                                        'GeneSymbol',
                                                        'Start',
                                                        'Stop',
                                                        'ReferenceAllele',
                                                        'AlternateAllele'],
                                                       axis=1)

    merged_categorized_df['pair_id'] = pair_id
    merged_categorized_df.to_csv('{}.categorized.maf'.format(pair_id), sep='\t', index=False)

    sys.stdout.write('Total runtime: {}'.format((datetime.now() - start)))


if __name__ == '__main__':
    main()
