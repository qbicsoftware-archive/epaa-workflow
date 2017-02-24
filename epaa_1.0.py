#! /usr/bin/python
import os
import sys
import logging
import csv
import re
import vcf
import argparse
import urllib2
import itertools
import pandas as pd
import numpy as np
import Fred2.Core.Generator as generator

from collections import defaultdict
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.Core.Variant import Variant, VariationType, MutationSyntax
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.ADBAdapter import EIdentifierTypes
from Fred2.IO.UniProtAdapter import UniProtDB
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.IO import FileReader

from datetime import datetime
from string import Template

__author__ = 'mohr, walzer'
VERSION = "1.0"

ID_SYSTEM_USED = EIdentifierTypes.ENSEMBL
transcriptProteinMap = {}

# Adapt depending on system
REF_PATH = ['/lustre_cfc/qbic/reference_genomes/IRMA_references/uniprot_all_swissprot_UP000005640.fasta', '/lustre_cfc/qbic/reference_genomes/IRMA_references/Homo_sapiens.GRCh38.pep.all.fa']

REPORT_TEMPLATE = """
###################################################################

        EPAA - EPITOPE PREDICTION REPORT

###################################################################

Persons in Charge: Christopher Mohr, Mathias Walzer

Date: $date

Pipeline Version: 1.0
Workflow Version: 1.0

Sample ID: $sample

Alleles
-------------
$alleles

Used Prediction Methods
-------------
$methods

Used Reference
-------------
$reference

Binding Assessment Criteria
-------------
Syfpeithi predictions: prediction score > half max score of corresponding allele
netMHC/netMHCpan predictions: affinity (as IC50 value in nM) <= 500

Additional Steps
-------------
NO filtering for peptide input.

Filtering of self-peptides (Reviewed (Swiss-Prot) UP000005640 uniprot-all.fasta.gz - 29/02/16, ENSEMBL release 84 Homo_sapiens.GRCh38.pep.all.fa.gz - 27/04/2016)
When personalized protein sequences are provided, peptides will be filtered against those as well.

Stats
-------------
Number of Variants: $variants
Number of Peptides: $peptides
Number of Peptides after Filtering: $filter
Number of Predictions: $predictions
Number of Predicted Binders: $binders
Number of Predicted Non-Binders: $nonbinders
Number of Binding Peptides: $uniquebinders
Number of Non-Binding Peptides: $uniquenonbinders

Contacts
-------------
mohr@informatik.uni-tuebingen.de
walzer@informatik.uni-tuebingen.de
University of Tuebingen, Applied Bioinformatics,
Center for Bioinformatics, Quantitative Biology Center,
and Dept. of Computer Science,
Sand 14, 72076 Tuebingen, Germany
"""


def get_fred2_annotation(vt, p, r, alt):
    if vt == VariationType.SNP:
        return p, r, alt
    elif vt == VariationType.DEL or vt == VariationType.FSDEL:
        # more than one observed ?
        if alt != '-':
            alternative = '-'
            reference = r[len(alt):]
            position = p + len(alt)
        else:
            return p, r, alt
    elif vt == VariationType.INS or vt == VariationType.FSINS:
        if r != '-':
            position = p
            reference = '-'
            if alt != '-':
                alt_new = alt[len(r):]
                alternative = alt_new
            else:
                alternative = str(alt)
        else:
            return p, r, alt

    return position, reference, alternative


def check_min_req_GSvar(row):
    """
    checking the presence of mandatory columns
    :param row: dictionary of a GSvar row
    :return: boolean, True if min req met
    """
    if ("#chr" in row.keys() and "start" in row.keys() and "end" in row.keys() and "ref" in row.keys() and "obs" in row.keys() and ("coding_and_splicing_details" in row.keys() or "coding" in row.keys() or "coding_and_splicing" in row.keys())):
        return True
    return False


def read_GSvar(filename, pass_only=True):
    """
    reads GSvar and tsv files (tab sep files in context of genetic variants), omitting and warning about rows missing
    mandatory columns
    :param filename: /path/to/file
    :return: list FRED2 variants
    """
    global ID_SYSTEM_USED
    RE = re.compile("(\w+):([\w.]+):([&\w]+):\w*:exon(\d+)\D*\d*:(c.\D*([_\d]+)\D*):(p.\D*(\d+)\w*)")

    list_vars = list()
    lines = list()
    transcript_ids = []

    with open(filename, 'rb') as tsvfile:
        tsvreader = csv.DictReader((row for row in tsvfile if not row.startswith('##')), delimiter='\t')
        for row in tsvreader:
            if not check_min_req_GSvar(row):
                logging.warning("read_GSvar: Omitted row! Mandatory columns not present in: \n"+str(row))
                continue
            lines.append(row)
    for mut_id, line in enumerate(lines):
        if pass_only and line["filter"].strip():
            continue
        genome_start = int(line["start"]) - 1
        genome_stop = int(line["end"]) - 1
        chrom = line["#chr"]
        ref = line["ref"]
        alt = line["obs"]

        # metadata
        variation_dbid = line.get("dbSNP", '')
        norm_depth = line.get("normal_dp", '')
        tum_depth = line.get("tumor_dp", '')
        tum_af = line.get("tumor_af", '')
        normal_af = line.get("normal_af", '')
        rna_tum_freq = line.get("rna_tum_freq", '')
        rna_tum_dp = line.get("rna_tum_depth", '')

        gene = line.get("gene", '')

        isHomozygous = True if (('tumour_genotype' in line) and (line['tumour_genotype'].split('/')[0] == line['tumour_genotype'].split('/')[1])) else False
        # old GSvar version
        if "coding_and_splicing_details" in line:
            mut_type = line.get("variant_details", '')
            annots = RE.findall(line["coding_and_splicing_details"])
        else:
            mut_type = line.get("variant_type", '')
            annots = RE.findall(line["coding_and_splicing"])
        isyn = mut_type == "synonymous_variant"

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
        """
        vt = VariationType.UNKNOWN
        if mut_type == 'missense_variant' or 'missense_variant' in mut_type:
            vt = VariationType.SNP
        elif mut_type == 'frameshift_variant':
            if (ref == '-') or (len(ref) < len(alt)):
                vt = VariationType.FSINS
            else:
                vt = VariationType.FSDEL
        elif mut_type == "inframe_deletion":
            vt = VariationType.DEL
        elif mut_type == "inframe_insertion":
            vt = VariationType.INS

        coding = dict()

        for annot in annots:
            a_gene, nm_id, a_mut_type, exon, trans_coding, trans_pos, prot_coding, prot_start = annot
            if 'NM' in nm_id:
                ID_SYSTEM_USED = EIdentifierTypes.REFSEQ
            if "stop_gained" not in mut_type:
                if not gene:
                    gene = a_gene
                if not mut_type:
                    mut_type = a_mut_type
                nm_id = nm_id.split(".")[0]

                coding[nm_id] = MutationSyntax(nm_id, int(trans_pos.split('_')[0])-1, int(prot_start)-1, trans_coding, prot_coding)
                transcript_ids.append(nm_id)  
        if coding:
            var = Variant(mut_id, vt, chrom.strip('chr'), int(genome_start),
                        ref.upper(), alt.upper(), coding, isHomozygous, isSynonymous=isyn)
            var.gene = gene
            var.log_metadata("vardbid", variation_dbid)
            var.log_metadata("normal_dp", norm_depth)
            var.log_metadata("tumor_dp", tum_depth)
            var.log_metadata("tumor_af", tum_af)
            var.log_metadata("normal_af", normal_af)
            var.log_metadata("rna_tum_freq", rna_tum_freq)
            var.log_metadata("rna_tum_depth", rna_tum_dp)
            list_vars.append(var)

    return list_vars, transcript_ids, ["vardbid", "normal_dp", "tumor_dp", "tumor_af", "normal_af", "rna_tum_freq", "rna_tum_depth"]


def read_vcf(filename, pass_only=True):
    """
    reads vcf files
    returns a list of FRED2 variants
    :param filename: /path/to/file
    :return: list of FRED2 variants
    """
    global ID_SYSTEM_USED

    vl = list()
    with open(filename, 'rb') as tsvfile:
        vcf_reader = vcf.Reader(tsvfile)
        vl = [r for r in vcf_reader]

    list_vars = []
    transcript_ids = []
    genotye_dict = {"het": False, "hom": True, "ref": True}

    for num, record in enumerate(vl):
        c = record.CHROM.strip('chr')
        p = record.POS - 1
        variation_dbid = record.ID
        r = str(record.REF)
        v_list = record.ALT
        f = record.FILTER
        if pass_only and f:
            continue

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
        """
        vt = VariationType.UNKNOWN
        if record.is_snp:
            vt = VariationType.SNP
        elif record.is_indel:
            if len(v_list) % 3 == 0:  # no frameshift
                if record.is_deletion:
                    vt = VariationType.DEL
                else:
                    vt = VariationType.INS
            else:  # frameshift
                if record.is_deletion:
                    vt = VariationType.FSDEL
                else:
                    vt = VariationType.FSINS
        gene = ''

        for alt in v_list:
            isHomozygous = False
            if 'HOM' in record.INFO:
                isHomozygous = record.INFO['HOM'] == 1
            elif 'SGT' in record.INFO:
                zygosity = record.INFO['SGT'].split("->")[1]
                if zygosity in genotye_dict:
                    isHomozygous = genotye_dict[zygosity]
                else:
                    if zygosity[0] == zygosity[1]:
                        isHomozygous = True
                    else:
                        isHomozygous = False
            else:
                for sample in record.samples:
                    if 'GT' in sample.data:
                        isHomozygous = sample.data['GT'] == '1/1'

            if record.INFO['ANN']:
                isSynonymous = False
                coding = dict()
                types = []
                for annraw in record.INFO['ANN']:  # for each ANN only add a new coding! see GSvar
                    annots = annraw.split('|')
                    obs, a_mut_type, impact, a_gene, a_gene_id, feature_type, transcript_id, exon, tot_exon, trans_coding, prot_coding, cdna, cds, aa, distance, warnings = annots
                    types.append(a_mut_type)

                    tpos = 0
                    ppos = 0
                    positions = ''

                    # get cds/protein positions and convert mutation syntax to FRED2 format
                    if trans_coding != '':
                        positions = re.findall(r'\d+', trans_coding)
                        ppos = int(positions[0]) - 1

                    if prot_coding != '':
                        positions = re.findall(r'\d+', prot_coding)
                        tpos = int(positions[0]) - 1

                    isSynonymous = (a_mut_type == "synonymous_variant")

                    gene = a_gene_id
                    # there are no isoforms in biomart
                    transcript_id = transcript_id.split(".")[0]

                    if 'NM' in transcript_id:
                        ID_SYSTEM_USED = EIdentifierTypes.REFSEQ

                    #take online coding variants into account, FRED2 cannot deal with stopgain variants right now
                    if not prot_coding or 'stop_gained' in a_mut_type:
                        continue

                    coding[transcript_id] = MutationSyntax(transcript_id, ppos, tpos, trans_coding, prot_coding)
                    transcript_ids.append(transcript_id)

                if coding:
                    pos, reference, alternative = get_fred2_annotation(vt, p, r, str(alt))
                    var = Variant("line" + str(num), vt, c, pos, reference, alternative, coding, isHomozygous, isSynonymous)
                    var.gene = gene
                    var.log_metadata("vardbid", variation_dbid)
                    list_vars.append(var)

    return list_vars, transcript_ids


def read_peptide_input(filename):
    peptides = []
    metadata = []

    '''expected columns (min required): id sequence'''
    with open(filename, 'r') as peptide_input:
        reader = csv.DictReader(peptide_input, delimiter='\t')
        for row in reader:
            pep = Peptide(row['sequence'])

            for col in row:
                if col != 'sequence':
                    pep.log_metadata(col, row[col])
                    metadata.append(col)
            peptides.append(pep)

    metadata = set(metadata)
    return peptides, metadata


# parse protein_groups of MaxQuant output to get protein intensitiy values
def read_protein_quant(filename):
    # protein id: sample1: intensity, sample2: instensity:
    intensities = {}

    with open(filename, 'r') as inp:
        inpreader = csv.DictReader(inp, delimiter='\t')
        for row in inpreader:
            if 'CON' in row['Protein IDs'] or 'REV' in row['Protein IDs']:
                pass
            else:
                valuedict = {}
                for key, val in row.iteritems():
                    if 'LFQ intensity' in key:
                        valuedict[key.replace('LFQ intensity ', '')] = val
                # do we have to split the ID somehow ?
                intensities[row['Protein IDs']] = valuedict
    return intensities


# parse different expression analysis results (DESeq2), link log2fold changes to transcripts/genes
def read_diff_expression_values(filename):
    # feature id: log2fold changes
    fold_changes = {}

    with open(filename, 'r') as inp:
        inp.readline()
        for row in inp:
            values = row.strip().split('\t')
            fold_changes[values[0]] = values[1]

    return fold_changes


def create_length_column_value(pep):
    return int(len(pep[0]))


def create_protein_column_value(pep):
    all_proteins = [transcriptProteinMap[x.transcript_id.split(':')[0]] for x in set(pep[0].get_all_transcripts())]
    return ','.join(set([item for sublist in all_proteins for item in sublist]))


def create_transcript_column_value(pep):
    return ','.join(set([x.transcript_id.split(':')[0] for x in set(pep[0].get_all_transcripts())]))


def create_mutationsyntax_column_value(pep):
    transcript_ids = [x.transcript_id for x in set(pep[0].get_all_transcripts())]
    variants = []
    syntaxes = []
    for t in transcript_ids:
        variants.extend([v for v in pep[0].get_variants_by_protein(t)])
    transcript_ids = set([t.split(':')[0] for t in transcript_ids])
    for v in set(variants):
        for c in v.coding:
            if c in transcript_ids:
                syntaxes.append(v.coding[c])
    return ','.join(set([y.aaMutationSyntax for y in syntaxes]))


def create_mutationsyntax_genome_column_value(pep):
    transcript_ids = [x.transcript_id for x in set(pep[0].get_all_transcripts())]
    variants = []
    syntaxes = []
    for t in transcript_ids:
        variants.extend([v for v in pep[0].get_variants_by_protein(t)])
    transcript_ids = set([t.split(':')[0] for t in transcript_ids])
    for v in set(variants):
        for c in v.coding:
            if c in transcript_ids:
                syntaxes.append(v.coding[c])
    return ','.join(set([y.cdsMutationSyntax for y in syntaxes]))


def create_variationfilelinenumber_column_value(pep):
    v = [x.vars.values() for x in pep[0].get_all_transcripts()]
    vf = list(itertools.chain.from_iterable(v))
    return ','.join([str(int(y.id.replace('line', ''))+1) for y in vf])


def create_gene_column_value(pep):
    transcript_ids = [x.transcript_id for x in set(pep[0].get_all_transcripts())]
    variants = []
    for t in transcript_ids:
        variants.extend([v for v in pep[0].get_variants_by_protein(t)])
    return ','.join(set([y.gene for y in set(variants)]))


def create_variant_pos_column_value(pep):
    transcript_ids = [x.transcript_id for x in set(pep[0].get_all_transcripts())]
    variants = []
    for t in transcript_ids:
        variants.extend([v for v in pep[0].get_variants_by_protein(t)])
    return ','.join(set(['{}'.format(y.genomePos) for y in set(variants)]))


def create_variant_chr_column_value(pep):
    transcript_ids = [x.transcript_id for x in set(pep[0].get_all_transcripts())]
    variants = []
    for t in transcript_ids:
        variants.extend([v for v in pep[0].get_variants_by_protein(t)])
    return ','.join(set(['{}'.format(y.chrom) for y in set(variants)]))


def create_variant_type_column_value(pep):
    types = {0: 'SNP', 1: 'DEL', 2: 'INS', 3: 'FSDEL', 4: 'FSINS', 5: 'UNKNOWN'}

    transcript_ids = [x.transcript_id for x in set(pep[0].get_all_transcripts())]
    variants = []
    for t in transcript_ids:
        variants.extend([v for v in pep[0].get_variants_by_protein(t)])
    return ','.join(set([types[y.type] for y in set(variants)]))


def create_variant_syn_column_value(pep):

    v = [x.vars.values() for x in pep[0].get_all_transcripts()]
    vf = list(itertools.chain.from_iterable(v))
    return ','.join(set([str(y.isSynonymous) for y in set(vf)]))


def create_variant_hom_column_value(pep):
    v = [x.vars.values() for x in pep[0].get_all_transcripts()]
    vf = list(itertools.chain.from_iterable(v))
    return ','.join(set([str(y.isHomozygous) for y in set(vf)]))


def create_coding_column_value(pep):
    v = [x.vars.values() for x in pep[0].get_all_transcripts()]
    vf = list(itertools.chain.from_iterable(v))
    return ','.join(set([str(y.coding) for y in set(vf)]))


def create_metadata_column_value(pep, c):
    v = [x.vars.values() for x in pep[0].get_all_transcripts()]
    vf = list(itertools.chain.from_iterable(v))
    return ','.join(set([str(y.get_metadata(c)[0]) for y in set(vf)]))


def create_quant_column_value(row, dict):
    if row[1] in dict:
        value = dict[row[1]]
    else:
        value = np.nan
    return value


def create_quant_column_value_for_result(row, dict):
    values = []
    for p in row['proteins'].split(','):
        if p in dict:
            values.append(dict[p])
    return ','.join(values)


#defined as : RPKM = (10^9 * C)/(N * L)
# L = exon length in base-pairs for a gene
# C = Number of reads mapped to a gene in a single sample
# N = total (unique)mapped reads in the sample
def create_expression_column_value_for_result(row, dict, deseq, transcript_objects):
    ts = row['transcripts'].split(',')
    if deseq:
        for t in ts:
            if t in dict:
                value = dict[t]
            else:
                value = np.nan
    else:
        for t in ts:
            if t in dict:
                value = (10**9 * dict[t]) / (len(transcript_objects[t]) * sum([dict[k] for k in dict.keys() if not (k.startswith('__'))]))
            else:
                value = np.nan
    return value


def create_quant_column_value_for_result(row, dict):
    values = []
    for p in row['proteins'].split(','):
        if p in dict:
            values.append(dict[p])
    return ','.join(values)


def write_prediction_report(values):
    s = Template(REPORT_TEMPLATE)
    return s.substitute(values)


def get_protein_ids_for_transcripts(idtype, transcripts, ensembl_url):
    result = {}

    biomart_url = "{}/biomart/martservice?query=".format(ensembl_url)
    biomart_head = """
    <?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query client="true" processor="TSV" limit="-1" header="1" uniqueRows = "1" >
            <Dataset name="%s" config="%s">
    """.strip()
    biomart_tail = """
            </Dataset>
        </Query>
    """.strip()
    biomart_filter = """<Filter name="%s" value="%s" filter_list=""/>"""
    biomart_attribute = """<Attribute name="%s"/>"""

    ENSEMBL = False
    if idtype == EIdentifierTypes.ENSEMBL:
        idname = "ensembl_transcript_id"
        ENSEMBL = True
    elif idtype == EIdentifierTypes.REFSEQ:
        idname = "refseq_mrna"

    input_lists = []

    # too long requests will fail
    if len(transcripts) > 800:
        input_lists = [transcripts[i:i + 3] for i in xrange(0, len(transcripts), 3)]

    else:
        input_lists += [transcripts]

    tsvselect = []
    for l in input_lists:
        rq_n = biomart_head % ('hsapiens_gene_ensembl', 'default') \
             + biomart_filter % (idname, ','.join(l)) \
             + biomart_attribute % ("ensembl_peptide_id") \
             + biomart_attribute % ("uniprot_swissprot") \
             + biomart_attribute % ("refseq_peptide") \
             + biomart_attribute % (idname) \
             + biomart_tail

        tsvreader = csv.DictReader(urllib2.urlopen(biomart_url + urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')

        tsvselect += [x for x in tsvreader]

    if(ENSEMBL):
        key = 'Ensembl Transcript ID'
        for dic in tsvselect:
            if dic[key] in result:
                #merged = result[dic[key]] + [dic['Ensembl Protein ID']]
                merged = result[dic[key]] + [dic['Protein ID']]
                result[dic[key]] = merged
            else:
                #result[dic[key]] = [dic['Ensembl Protein ID']]
                result[dic[key]] = [dic['Protein ID']]
    else:
        key = 'RefSeq mRNA [e.g. NM_001195597]'
        for dic in tsvselect:
            if dic[key] in result:
                merged = result[dic[key]] + [dic['RefSeq Protein ID [e.g. NP_001005353]']]
                result[dic[key]] = merged
            else:
                result[dic[key]] = [dic['RefSeq Protein ID [e.g. NP_001005353]']]

    return result


def get_matrix_max_score(allele, length):
    allele_model = "%s_%i" % (allele, length)
    try:
        pssm = getattr(__import__("Fred2.Data.pssms.syfpeithi.mat."+allele_model, fromlist=[allele_model]),
                           allele_model)
        return sum([max(scrs.values()) for pos, scrs in pssm.iteritems()])
    except:
        return np.nan


def create_affinity_values(allele, length, j, method, max_scores, allele_strings):
    if not pd.isnull(j):
        if 'syf' in method:
            return round(((100.0 / float(max_scores[allele_strings[('%s_%s' % (str(allele), length))]]) * float(j)) / 100.0) * 100, 2)
        else:
            return round((50000**(1.0-float(j))), 2)
    else:
        return np.nan


def create_binder_values(aff, method):
    if not pd.isnull(aff):
        if 'syf' in method:
            return True if aff > 50.0 else False
        else:
            return True if aff <= 500.0 else False
    else:
        return np.nan

def make_predictions_from_variants(variants_all, methods, alleles, minlength, maxlength, martsadapter, protein_db, identifier, metadata, transcriptProteinMap):
    # list for all peptides and filtered peptides
    all_peptides = []
    all_peptides_filtered = []

    # dictionaries for syfpeithi matrices max values and allele mapping
    max_values_matrices = {}
    allele_string_map = {}

    # list to hold dataframes for all predictions
    pred_dataframes = []

    # list to hold peptide transcript connections
    pep_transcripts = []

    # list to hold peptide protein connections
    pep_proteins = []

    transcripts = {}

    # get variants and group them by chromosome
    variants_dict = defaultdict(list)
    for v in variants_all:
        variants_dict[v.chrom].append(v)

    for c in variants_dict:
        variants = variants_dict[c]

        prots = [p for p in generator.generate_proteins_from_transcripts(generator.generate_transcripts_from_variants(variants, martsadapter, ID_SYSTEM_USED))]

        for peplen in range(minlength, maxlength):
            logging.info("Prediction of length %i" % (peplen))

            peptide_gen = generator.generate_peptides_from_proteins(prots, peplen)

            peptides_var = [x for x in peptide_gen]

            # remove peptides which are not 'variant relevant'
            peptides = [x for x in peptides_var if any(x.get_variants_by_protein(y) for y in x.proteins.keys())]

            # filter out self peptides
            selfies = [str(p) for p in peptides if protein_db.exists(str(p))]
            filtered_peptides = [p for p in peptides if str(p) not in selfies]

            all_peptides = all_peptides + peptides
            all_peptides_filtered = all_peptides_filtered + filtered_peptides

            for pep in filtered_peptides:
                for t in pep.get_all_transcripts():
                    transcripts[t.transcript_id.split(':')[0]] = t
                    pep_transcripts.append([str(pep), t.transcript_id.split(':')[0]])

            for pep in filtered_peptides:
                for t in pep.get_all_transcripts():
                    #protein IDs will be filtered already according to ID system
                        pep_proteins.append([str(pep), set(transcriptProteinMap[t.transcript_id.split(':')[0]])])

            # run prediction for peptides of length peplen for all alleles and all methods
            results = []

            for m in methods:
                try:
                    results.extend([EpitopePredictorFactory(m.split('-')[0], version=m.split('-')[1]).predict(filtered_peptides, alleles=alleles)])
                except:
                    logging.warning("Prediction for length {length} and allele {allele} not possible with {method}. No model available.".format(length=peplen, allele=','.join([str(a) for a in alleles]), method=m))

            if(len(results) == 0):
                continue

            df = results[0].merge_results(results[1:])

            for a in alleles:
                conv_allele = "%s_%s%s" % (a.locus, a.supertype, a.subtype)
                allele_string_map['%s_%s' % (a, peplen)] = '%s_%i' % (conv_allele, peplen)
                max_values_matrices['%s_%i' % (conv_allele, peplen)] = get_matrix_max_score(conv_allele, peplen)

            df.insert(0, 'length', df.index.map(create_length_column_value))
            df['chr'] = df.index.map(create_variant_chr_column_value)
            df['pos'] = df.index.map(create_variant_pos_column_value)
            df['gene'] = df.index.map(create_gene_column_value)
            df['transcripts'] = df.index.map(create_transcript_column_value)
            df['proteins'] = df.index.map(create_protein_column_value)
            df['variant type'] = df.index.map(create_variant_type_column_value)
            df['synonymous'] = df.index.map(create_variant_syn_column_value)
            df['homozygous'] = df.index.map(create_variant_hom_column_value)
            df['variant details (genomic)'] = df.index.map(create_mutationsyntax_genome_column_value)
            df['variant details (protein)'] = df.index.map(create_mutationsyntax_column_value)

            # reset index to have index as columns
            df.reset_index(inplace=True)

            for c in df.columns:
                if '*' in str(c):
                    idx = df.columns.get_loc(c)
                    df.insert(idx + 1, '%s affinity' % c, df.apply(lambda x: create_affinity_values(str(c), int(x['length']), float(x[c]), x['Method'], max_values_matrices, allele_string_map), axis=1))
                    df.insert(idx + 2, '%s binder' % c, df.apply(lambda x: create_binder_values(float(x['%s affinity' % c]), x['Method']), axis=1))
                    df = df.rename(columns={c: '%s score' % c})
                    df['%s score' % c] = df['%s score' % c].map(lambda x: round(x, 4))

            for c in metadata:
                df[c] = df.apply(lambda row: create_metadata_column_value(row, c), axis=1)

            df = df.rename(columns={'Seq': 'sequence'})
            df = df.rename(columns={'Method': 'method'})
            pred_dataframes.append(df)

    statistics = {'date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"), 'sample': identifier, 'alleles': '\n'.join([str(a) for a in alleles]),
        'methods': '\n'.join(methods), 'variants': len(variants_all), 'peptides': len(all_peptides), 'filter': len(all_peptides_filtered)}

    return pred_dataframes, statistics, pep_transcripts, transcripts, pep_proteins


def make_predictions_from_peptides(peptides, methods, alleles, protein_db, identifier, metadata):
    # dictionaries for syfpeithi matrices max values and allele mapping
    max_values_matrices = {}
    allele_string_map = {}

    # list to hold dataframes for all predictions
    pred_dataframes = []

    # filter out self peptides if specified
    selfies = [str(p) for p in peptides if protein_db.exists(str(p))]
    peptides_filtered = [p for p in peptides if str(p) not in selfies]

    # sort peptides by length (for predictions)
    sorted_peptides = {}

    for p in peptides_filtered:
        length = len(str(p))
        if length in sorted_peptides:
            sorted_peptides[length].append(p)
        else:
            sorted_peptides[length] = [p]

    for peplen in sorted_peptides:
        all_peptides_filtered = sorted_peptides[peplen]
        results = []
        for m in methods:
            try:
                results.extend([EpitopePredictorFactory(m.split('-')[0], version=m.split('-')[1]).predict(all_peptides_filtered, alleles=alleles)])
            except:
                logging.warning("Prediction for length {length} and allele {allele} not possible with {method}. No model available.".format(length=peplen, allele=','.join([str(a) for a in alleles]), method=m))

        # merge dataframes of the performed predictions
        if(len(results) == 0):
            continue;
        df = results[0].merge_results(results[1:])

        df.insert(0, 'length', df.index.map(create_length_column_value))

        for a in alleles:
            conv_allele = "%s_%s%s" % (a.locus, a.supertype, a.subtype)
            allele_string_map['%s_%s' % (a, peplen)] = '%s_%i' % (conv_allele, peplen)
            max_values_matrices['%s_%i' % (conv_allele, peplen)] = get_matrix_max_score(conv_allele,peplen)

        # reset index to have index as columns
        df.reset_index(inplace=True)

        mandatory_columns = ['chr', 'pos', 'gene', 'transcripts', 'proteins', 'variant type', 'synonymous', 'homozygous', 'variant details (genomic)', 'variant details (protein)']

        for header in mandatory_columns:
            if header not in metadata:
                df[header] = np.nan
            else:
                df[header] = df.apply(lambda row: row[0].get_metadata(header)[0], axis=1)

        for c in list(set(metadata) - set(mandatory_columns)):
            df[c] = df.apply(lambda row: row[0].get_metadata(c)[0], axis=1)

        for c in df.columns:
            if '*' in str(c):
                idx = df.columns.get_loc(c)
                df.insert(idx + 1, '%s affinity' % c, df.apply(lambda x: create_affinity_values(str(c), int(x['length']), float(x[c]), x['Method'], max_values_matrices, allele_string_map), axis=1))
                df.insert(idx + 2, '%s binder' % c, df.apply(lambda x: create_binder_values(float(x['%s affinity' % c]), x['Method']), axis=1))
                df = df.rename(columns={c: '%s score' % c})

        df = df.rename(columns={'Seq': 'sequence'})
        df = df.rename(columns={'Method': 'method'})
        pred_dataframes.append(df)

    # write prediction statistics
    statistics = {'date': str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")), 'sample': identifier, 'alleles': '\n'.join([str(a) for a in alleles]), 'methods': '\n'.join(methods),
    'variants': '-', 'peptides': len(peptides), 'filter': len(peptides_filtered), 'reference': '-'}

    return pred_dataframes, statistics


def __main__():
    parser = argparse.ArgumentParser(description="""EPAA 1.0 \n Pipeline for prediction of MHC class I and II epitopes from variants or peptides for a list of specified alleles. 
        Additionally predicted epitopes can be annotated with protein quantification values for the corresponding proteins or differential expression values for the corresponding transcripts.""", version=VERSION)
    parser.add_argument('-s', "--somatic_mutations", help='Somatic variants')
    parser.add_argument('-g', "--germline_mutations", help="Germline variants")
    parser.add_argument('-p', "--peptides", help="File with one peptide per line")
    parser.add_argument('-i', "--identifier", help="<Required> Predictions will be written with this name prefix", required=True)
    parser.add_argument('-c', "--mhcclass", default="I", help="MHC class I or II")
    parser.add_argument('-a', "--alleles", help="<Required> MHC Alleles", required=True)
    parser.add_argument('-r', "--reference", help="Reference, retrieved information will be based on this ensembl version", required=False, default='GRCh37', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('-f', "--filter_self", help="Filter peptides against human proteom", required=False, action='store_true')
    parser.add_argument('-pq', "--protein_quantification", help="File with protein quantification values")
    parser.add_argument('-de', "--differential_expression", help="File with differential expression analysis results (DESeq2 Output)")
    parser.add_argument('-ip', "--individualized_proteins", help="Fasta with individualized protein sequences")
    parser.add_argument('-o', "--output_dir", help="All files written will be put in this directory")

    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if((args.peptides is not None) & (args.somatic_mutations is not None or args.germline_mutations is not None)):
        parser.error("Please specify peptide file OR variant file.")
    elif ((args.peptides is None) & (args.somatic_mutations is None)):
        parser.error("Please specify at least one of the following options: --somatic_mutations, --peptides")

    logging.basicConfig(filename=os.path.join(args.output_dir,'{}_prediction.log'.format(args.identifier)), filemode='w+',
                        level=logging.DEBUG)
    logging.info("Starting predictions at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    if args.output_dir is not None:
        try:
            os.chdir(args.output_dir)
            logging.info("Using provided data directory: {}".format(str(args.output_dir)))
        except:
            logging.info("No such directory, using current.")
    else:
        logging.info("Using current data directory.")

    '''start the actual IRMA functions'''
    metadata = []
    references = {'GRCh37': 'http://grch37.ensembl.org', 'GRCh38': 'http://ensembl.org'}
    global transcriptProteinMap

    '''read in variants or peptides'''
    if args.peptides is not None:
        peptides, metadata = read_peptide_input(args.peptides)
    else:
        if args.somatic_mutations.endswith('.GSvar'):
            vl, transcripts, metadata = read_GSvar(args.somatic_mutations)
        elif args.somatic_mutations.endswith('.vcf'):
            vl, transcripts = read_vcf(args.somatic_mutations)

        if args.germline_mutations is not None:
            if args.germline_mutations.endswith('.GSvar'):
                vl_normal, transcripts_germline, metadata = read_GSvar(args.germline_mutations)
            elif args.germline_mutations.endswith('.vcf'):
                vl_normal, transcripts_germline = read_vcf(args.germline_mutations)

            # combine germline and somatic variants
            vl = vl + vl_normal
            transcripts = transcripts_germline + transcripts

        transcripts = list(set(transcripts))
        transcriptProteinMap = get_protein_ids_for_transcripts(ID_SYSTEM_USED, transcripts, references[args.reference])

    # get the alleles
    alleles = FileReader.read_lines(args.alleles, in_type=Allele)

    # initialize MartsAdapter, GRCh37 or GRCh38 based
    ma = MartsAdapter(biomart=references[args.reference])

    # create protein db instance for filtering self-peptides
    up_db = UniProtDB('sp')
    if args.filter_self:
        logging.info('Reading human proteome')
        up_db.read_seqs(REF_PATH[0])
        up_db.read_seqs(REF_PATH[1])

    if args.individualized_proteins is not None:
        up_db.read_seqs(args.individualized_proteins)

    # MHC class I or II predictions
    if args.mhcclass == "I":
        methods = ['netmhc-4.0', 'syfpeithi-1.0', 'netmhcpan-3.0']
        if args.peptides:
            pred_dataframes, statistics = make_predictions_from_peptides(peptides, methods, alleles, up_db, args.identifier, metadata)
        else:
            pred_dataframes, statistics, pep_transcripts, transcripts_objs, pep_proteins = make_predictions_from_variants(vl, methods, alleles, 8, 12, ma, up_db, args.identifier, metadata, transcriptProteinMap)
    else:
        methods = ['netmhcII-2.2', 'syfpeithi-1.0', 'netmhcIIpan-3.1']
        if args.peptides:
            pred_dataframes, statistics = make_predictions_from_peptides(peptides, methods, alleles, up_db, args.identifier, metadata)
        else:
            pred_dataframes, statistics, pep_transcripts, pep_proteins = make_predictions_from_variants(vl, methods, alleles, 15, 17, ma, up_db, args.identifier, metadata, transcriptProteinMap)

    # concat dataframes for all peptide lengths
    try:
        complete_df = pd.concat(pred_dataframes)
    except:
        complete_df = pd.DataFrame()
        logging.error("NO PREDICTIONS AVAILABLE.")

    # store version of used methods
    method_map = {}
    for m in methods:
        method_map[m.split('-')[0]] = m

    # replace method names with method names with version
    complete_df.replace({'method': method_map}, inplace=True)

    # Change the order (the index) of the columns
    columns_tiles = ['sequence', 'length', 'chr', 'pos', 'gene', 'transcripts', 'proteins', 'variant type', 'method']
    for c in complete_df.columns:
        if c not in columns_tiles:
            columns_tiles.append(c)
    complete_df = complete_df.reindex(columns=columns_tiles)

    binder_cols = [col for col in complete_df.columns if 'binder' in col]

    binders = []
    non_binders = []
    pos_predictions = []
    neg_predictions = []

    for i, r in complete_df.iterrows():
        binder = False
        for c in binder_cols:
            if r[c] is True:
                binder = True
                continue
        if binder:
            binders.append(str(r['sequence']))
            pos_predictions.append(str(r['sequence']))
        else:
            neg_predictions.append(str(r['sequence']))
            if str(r['sequence']) not in binders:
                non_binders.append(str(r['sequence']))

    # parse protein quantification results, annotate proteins for samples
    if args.protein_quantification is not None:
        protein_quant = read_protein_quant(args.protein_quantification)
        unique_data = [list(x) for x in set(tuple(x) for x in pep_proteins)]
        df = pd.DataFrame(unique_data, columns=['peptide', 'protein'])
        df['log2ratio'] = df.apply(lambda row: create_quant_column_value(row, protein_quant), axis=1)
        # write dataframe to tsv
        df.fillna('')
        df.to_csv("{}_protein_values.tsv".format(args.identifier), '\t', index=False)

        # add column to result dataframe
        complete_df['protein log2ratio'] = complete_df.apply(lambda row: create_quant_column_value_for_result(row, protein_quant), axis=1)


    # parse differential expression analysis results (DESe2), annotate features (genes/transcripts)
    if args.differential_expression is not None:
        fold_changes = read_diff_expression_values(args.differential_expression)
        unique_data = [list(x) for x in set(tuple(x) for x in pep_transcripts)]
        df = pd.DataFrame(unique_data, columns=['peptide', 'transcript'])

        if 'HTSeq' in args.differential_expression:
            col_name = 'RNA expression (rkpm)'
        else:
            col_name = 'RNA normal_vs_tumor.log2FoldChange'
        df[col_name] = df.apply(lambda row: create_quant_column_value(row, fold_changes), axis=1)
        # write dataframe to tsv
        df.fillna('')
        df.to_csv("{}_transcript_values.tsv".format(args.identifier), '\t', index=False)

        # add column to result dataframe
        complete_df[col_name] = complete_df.apply(lambda row: create_quant_column_value_for_result(row, fold_changes, transcripts_objs), axis=1)

    # write dataframe to tsv
    complete_df.fillna('')
    complete_df.to_csv("{}_prediction_results.tsv".format(args.identifier), '\t', index=False)

    statistics['predictions'] = complete_df.shape[0]
    statistics['binders'] = len(pos_predictions)
    statistics['nonbinders'] = len(neg_predictions)
    statistics['uniquebinders'] = len(set(binders))
    statistics['uniquenonbinders'] = len(set(non_binders) - set(binders))

    if 'reference' not in statistics:
        statistics['reference'] = args.reference

    with open('{}_prediction_statistics.txt'.format(args.identifier), 'w') as stats:
        stats.write(write_prediction_report(statistics))
    logging.info("Finished predictions at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))


if __name__ == "__main__":
    __main__()
