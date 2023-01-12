import re
import os
import sys
import itertools
import unidecode
import subprocess
sys.path.append("/mnt/chaelab/rachelle/src")
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "scripts"))

from data_manip import make_custom_get, splitlines
from datetime import datetime

## TODO: replace grep_bedmerge with whatever I did for MINORg

chrom_pref="Chr"
data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "data", "1135acc.csv")
bed_path = "/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed"
ref_fasta_pref = "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_"
ref_fasta = {chrom_pref + '1': ref_fasta_pref + "chr01.fasta",
             chrom_pref + '2': ref_fasta_pref + "chr02.fasta",
             chrom_pref + '3': ref_fasta_pref + "chr03.fasta",
             chrom_pref + '4': ref_fasta_pref + "chr04.fasta",
             chrom_pref + '5': ref_fasta_pref + "chr05.fasta",
             chrom_pref + 'M': ref_fasta_pref + "mt.fasta",
             chrom_pref + 'C': ref_fasta_pref + "pltd.fasta"}

# fields={"gff3": {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
#                          "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
#                          "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
#                          "Is_circular": "Is_circular"}},
#         "gtf2": {"all": {"ID": "transcript_id", "Parent": "gene_id"}}}
fields={"gff3": {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                         "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                         "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                         "Is_circular": "Is_circular"}},
        "gtf2": {"all": {"ID": "transcript_id", "Parent": "gene_id"}}}

def parse_genetic_code(c):
    try:
        return int(c)
    except ValueError:
        return c

def has_overlap(r1, r2):
    r1 = sorted(r1)
    r2 = sorted(r2)
    return (r2[0] <= r1[0] <=r2[1]) or (r1[0] <= r2[0] <= r1[1])

def has_any_overlap(l1, l2):
    any_overlap = [has_overlap(r1, r2) for r1 in l1 for r2 in l2]
    return True in any_overlap

def has_cat_overlap(l1, l2):
    return set(l1).intersection(set(l2)) != set()

def merge_ranges(*l):
    # print("l:", l)
    all_ranges = list(sorted(set(itertools.chain(*l))))
    # print("all_ranges:", all_ranges)
    if len(all_ranges) <= 1:
        return(all_ranges)
    final_ranges = [tuple(all_ranges.pop(0))]
    while all_ranges:
        if has_overlap(final_ranges[-1], all_ranges[0]):
            r1, r2 = tuple(final_ranges.pop(-1)), tuple(all_ranges.pop(0))
            # final_ranges.append(tuple(min(*r1, *r2), max(*r1, *r2)))
            final_ranges.append((min(*r1, *r2), max(*r1, *r2)))
        else:
            final_ranges.append(all_ranges.pop(0))
    # print("final_ranges:", final_ranges)
    return(final_ranges)

def store_fname(fout, fname):
    if fout:
        open(fout, "a+").write(fname + '\n')
    return

get_gffbed = make_custom_get(["chrom", "start", "end", "ID", "score", "strand", "source", "type", "phase", "attributes"])

def make_get_attribute(attribute_fields, attribute_mod = {}, **kwargs):
    ## update attribute_fields w/ attribute_mod
    from copy import deepcopy
    attribute_fields = deepcopy(attribute_fields)
    for feature, mods in attribute_mod.items():
        if feature in attribute_fields:
            for canon, atypical in mods.items():
                attribute_fields[feature][canon] = atypical
        else:
            attribute_fields[feature] = mods
    def parse_gff_attributes(s, field_sep_inter = ';', field_sep_intra = ','):
        ## if multiple separate entries for same field (e.g. 'Parent=abc.1;Parent=abc.2'), parse properly
        attributes_l = [x.split('=') for x in re.findall(f"[^{field_sep_intra}{field_sep_inter}=]+=.+?(?=[^{field_sep_intra}{field_sep_inter}=]+=.+?|$)", s)]
        attributes = {}
        for attribute in attributes_l:
            attributes[attribute[0]] = attributes.get(attribute[0], []) + \
                                       re.search(f'^.+?(?=[{field_sep_intra}{field_sep_inter}]?$)',
                                                 attribute[1]).group(0).split(field_sep_intra)
        return attributes
    def get_attribute(entry, a, feature = ''):
        ## map the field name to what's used in the file
        if not isinstance(entry, str):
            attributes = parse_gff_attributes(get_gffbed(entry, "attributes"), **kwargs)
            feature = get_gffbed(entry, "type")
        else:
            attributes = parse_gff_attributes(entry, **kwargs)
        ## get attributes
        # if feature in attribute_mod and a in attribute_mod[feature]:
        #     mapped_field = attribute_mod[feature][a]
        # elif feature in attribute_fields and a in attribute_fields[feature]:
        if feature in attribute_fields and a in attribute_fields[feature]:
            mapped_field = attribute_fields[feature][a]
        elif a in attribute_fields["all"]:
            mapped_field = attribute_fields["all"][a]
        else:
            mapped_field = a
        return attributes.get(mapped_field, [])
    return get_attribute

def grep_bedmerge(gene, bed, feature, out_dir, merge = False, encoding = "utf-8",
                  attribute_fields = fields["gff3"], field_sep_inter = ';', field_sep_intra = ',',
                  attribute_mod = {}, store_bed = None):
    ## get chrom, start, end of gene's features
    data = [x.split('\t') for x in splitlines(bed) if len(x) > 1]
    # data = [x.replace('\n', '').split('\t') for x in open(bed, 'r').readlines() if len(x) > 1]
    ## make local get_attribute function based on gff/gtf format provided
    get_attribute = make_get_attribute(attribute_fields, attribute_mod = attribute_mod)
    ## get relevant feature entries
    if feature == "gene":
        bed_raw = [x for x in data if get_gffbed(x, "type") == "gene" and gene in get_attribute(x, "ID")]
    elif feature == "mRNA":
        bed_raw = [x for x in data if get_gffbed(x, "type") == "mRNA" and (gene in get_attribute(x, "ID") or
                                                                           gene in get_attribute(x, "Parent"))]
    elif feature == "CDS":
        tmp_mrna = set(itertools.chain(*[get_attribute(x, "ID") for x in data
                                         if get_gffbed(x, "type") == "mRNA"
                                         and (gene in get_attribute(x, "ID") or
                                              gene in get_attribute(x, "Parent"))]))
        bed_raw = [x for x in data if (get_gffbed(x, "type") == "CDS" and
                                       has_cat_overlap(get_attribute(x, "Parent"), tmp_mrna))]
    else:
        bed_raw = [x for x in data if gene in get_gffbed(x, "attributes")]
    ## coerce into final columns (and merge overlapping entries if so required)
    if merge:
        printf = subprocess.Popen(("printf", '\n'.join(['\t'.join(x) for x in bed_raw])), stdin=None, stdout=subprocess.PIPE)
        output = [entry.split('\t') for entry in subprocess.check_output(("bedtools","merge","-c","6,8,10","-o","collapse,collapse,collapse"), stdin=printf.stdout).decode(encoding).split('\n') if entry]
    else:
        output = [[x[i] for i in (0,1,2,5,7,9)] for x in bed_raw]
    ## write bed file of regions used (0-indexed)
    fout = os.path.join(out_dir, "bed", gene + '_' + feature + ".bed")
    os.makedirs(os.path.dirname(fout), exist_ok=True)
    f = open(fout, "w+")
    f.write('\n'.join(['\t'.join(x) for x in output]) + '\n')
    f.close()
    store_fname(store_bed, fout)
    return {"fout": fout, "data": output}

# ## only accepts .tsv and .csv file if 'delim' is not specified.
# ## accepts a tuple (or list) of accession names
# ## returns [(acc_num:acc_name)] list, all strings
# def get_acc_num(acc_names, acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
#     ## read data
#     if not delim:
#         delim = ',' if acc_ref_file[-4:] == ".csv" else '\t' if acc_ref_file[-4:] == ".tsv" else ''
#     if not delim:
#         print("Unable to detect delimiter. Please specify delimiter with 'delim' keyword.")
#         sys.exit(1)
#     with open(acc_ref_file, 'r') as f:
#         ## convert all names to lower case, only keep alphanumeric characters (data)
#         data = {(re.sub(r'[^a-zA-Z0-9]', '',
#                         unidecode.unidecode(x.split(delim)[name_col]))).lower():x.split(delim)[num_col]\
#                 for x in f.readlines()[0 if not header else 1:]}
#     ## convert all names to lower case, only keep alphanumeric characters (query), and get acc num
#     return [(data[(re.sub(r'[^a-zA-Z0-9]', '', unidecode.unidecode(x))).lower()], x) for x in acc_names]

# ## only accepts .tsv and .csv file if 'delim' is not specified.
# ## accepts a tuple (or list) of accession names
# ## returns [(acc_num:acc_name)] list, all strings
# def get_acc_name(acc_nums, acc_ref_file=data_path, num_col=1, name_col=0, delim='', header=True):
#     ## read data
#     with open(acc_ref_file, 'r') as f:
#         if not delim:
#             delim = ',' if (acc_ref_file[-4:] == ".csv" or acc_ref_file[-4:] == ".txt") else \
#                     '\t' if acc_ref_file[-4:] == ".tsv" else ''
#         if not delim:
#             print("Unable to detect delimiter. Please specify delimiter with 'delim' keyword.")
#             sys.exit(1)
#         data = {x.split(delim)[num_col]:x.split(delim)[name_col]
#                 for x in f.readlines()[0 if not header else 1:]}
#     return [(str(x), data[str(x)]) for x in acc_nums]

# ## returns address from which to retrieve sequence data from 1001genomes
# def make_1001pseudogenome_api(chrom, start, end, accs, add_chr="Chr", remove_chr=True):
#     chrom = add_chr + (chrom if not remove_chr else re.search("[Cc]hr(.+)", chrom).group(1))
#     return "http://tools.1001genomes.org/api/v1/pseudogenomes/strains/" + ','.join(map(str, accs)) +\
#         "/regions/" + chrom + ':' + str(start) + ".." + str(end)

# ## parse accession num/name
# def raw_accs_to_id(accs, acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
#     accs = [acc for acc in accs if acc]
#     if str(accs[0]).isdigit():
#         accs = dict(get_acc_name(accs, acc_ref_file=acc_ref_file, num_col=num_col, name_col=name_col,
#                                  delim=delim, header=header))
#     else:
#         accs = dict(get_acc_num(accs, acc_ref_file=acc_ref_file, num_col=num_col, name_col=name_col,
#                                 delim=delim, header=header))
#     return accs

# def get_1001pseudogenome_raw(accs, chrom, start, end, fasta_out,
#                              acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
#     accs = raw_accs_to_id(accs, acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True)
#     subprocess.run(("curl", "-o", fasta_out,
#                     make_1001pseudogenome_api(chrom, start, end, accs.keys())))
#     return fasta_out

# ## get sequence info using genomic range
# def get_1001pseudogenome_by_range(accs, chrom, start, end, out_dir,
#                                   acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
#     fasta_out = os.path.join(out_dir, "chr{}_p{}-{}.fasta".format(chrom, start, end))
#     get_1001pseudogenome_raw(accs, chrom_pref + chrom, start, end, fasta_out,
#                              acc_ref_file=acc_ref_file, num_col=num_col, name_col=name_col, delim=delim,
#                              header=header)
#     ## ensure nice line width
#     from fasta_manip import fasta_to_dict
#     tmp = fasta_to_dict(fasta_out)
#     from fasta_manip import dict_to_fasta
#     dict_to_fasta(tmp, fasta_out)
#     print("Sequences were successfully written to: {}".format(fasta_out))
    
# ## get sequence info of genes
# def get_1001pseudogenome_by_gene(accs, gene, feature, out_dir, bed=bed_path, encoding="utf-8",
#                                  acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True,
#                                  complete=False, domain="", domain_f="", merge=False, adj_dir=False,
#                                  **for_get_domain_in_genome_coords):
#     data = grep_bedmerge(gene, bed, feature, out_dir, encoding = encoding, merge = merge)["data"]
#     ## extract sequences
#     chrom = data[0][0]
#     start = min([int(x[1]) for x in data])
#     end = max([int(x[2]) for x in data])
#     strand = data[0][3]
#     print(data[0], chrom, start, end)
#     if feature == "gene":
#         isoforms = {gene: data}
#     else:
#         isoforms = {re.search("=(.+\.\d+)(?=[,;]|$)", isoform).group(1): [x for x in data if x[-1] == isoform]
#                     for isoform in set([x[-1] for x in data])}
#     fasta_out_l = []
#     ## iterate through isoforms
#     for isoform, isoform_dat in isoforms.items():
#         ## if extracting only domain-specific range
#         if domain_f:
#             d_start, d_end = get_domain_in_genome_coords(gene, domain, domain_f, out_dir,
#                                                          bed=bed, encoding=encoding,
#                                                          isoform = isoform,
#                                                          start_inc = start_inc, end_inc = end_inc,
#                                                          **{k: v for k, v
#                                                             in for_get_domain_in_genome_coords.items()
#                                                             if k in ["qname_dname", "qstart_qend"]})
#             if (not d_start) or (not d_end):
#                 continue
#         ## get fasta file of sequence data
#         fasta_out = os.path.join(out_dir, isoform + '_' + feature + ("_complete" if complete else '') + \
#                                  (('_' + ("domain" if not domain else domain))if domain_f else '') + ".fasta")
#         get_1001pseudogenome_raw(accs, chrom, start + 1, end, fasta_out)
#         accs = raw_accs_to_id(accs, acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True)
#         ## rename seqs + trim if required
#         from fasta_manip import fasta_to_dict
#         tmp_seqs = fasta_to_dict(fasta_out)
#         seqs = {"{}|{}|{}|{}|{}".format(k.split('|')[3], accs[k.split('|')[3]], gene, feature, isoform) +\
#                 (('|' + ("domain" if not domain else domain)) if domain_f else '' +
#                  ("|complete" if complete else '')): seq \
#                 for k, seq in tmp_seqs.items()}
#         if (not complete) or domain_f:
#             if complete and domain_f and d_start and d_end:
#                 ranges = [(max(start, d_start) - start, min(end, d_end) - start)]
#             elif domain_f and d_start and d_end:
#                 ranges = [(max(int(x[1]), d_start) - start, min(int(x[2]), d_end) - start) \
#                           for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
#             else:
#                 ranges = [(int(x[1])-start, int(x[2])-start) for x in isoform_dat]
#             from fasta_manip import extract_ranges
#             seqs = {k: extract_ranges(seq, ranges) for k, seq in seqs.items()}
#         if adj_dir and strand == '-':
#             seqs = {k + "|revcomp": seq.reverse_complement() for k, seq in seqs.items()}
#         from fasta_manip import dict_to_fasta
#         dict_to_fasta(seqs, fasta_out)
#         fasta_out_l.append(fasta_out)
#     fasta_out_final = os.path.join(out_dir, gene + '_' + feature + ("_complete" if complete else '') + \
#                                    (('_' + ("domain" if not domain else domain)) if domain_f else '') + ".fasta")
#     if fasta_out_l and fasta_out_l[0] != fasta_out_final:
#         from file_manip import cat_files
#         cat_files(sorted(fasta_out_l), fasta_out_final)
#         for fasta_out in fasta_out_l:
#             os.remove(fasta_out)
#     print("Sequences were successfully written to: {}".format(fasta_out_final))

###################
## get reference ##
###################

def get_ref_raw(chrom, start, end, fasta_out, encoding="utf-8",
                ref_fasta_files=ref_fasta, mode="col0", ref_pref="Col-0_ref",
                store_fasta = None, **kwargs):
    # from fasta_manip import fasta_to_dict, dict_to_fasta
    from fasta_manip import dict_to_fasta
    from index_fasta import IndexedFasta
    if isinstance(ref_fasta_files, dict):
        # ref_seq = list(fasta_to_dict(ref_fasta_files[chrom]).values())[0][start:end] ## 0-indexed
        ref_seq = IndexedFasta(ref_fasta_files[chrom])[chrom][start:end]
        if mode == "col0":
            dict_to_fasta({"Col-0_ref|{}:{}..{}".format(chrom, start + 1, end): ref_seq}, fasta_out)
        else:
            dict_to_fasta({"Reference|{}:{}..{}".format(chrom, start + 1, end): ref_seq}, fasta_out)
    else:
        # ref_seq = fasta_to_dict(ref_fasta_files)[chrom][start:end]
        ref_seq = IndexedFasta(ref_fasta_files)[chrom][start:end]
        dict_to_fasta({"Reference|{}:{}..{}".format(chrom, start + 1, end): ref_seq}, fasta_out)
    store_fname(store_fasta, fasta_out)

## note: setting by_gene to True will collapse identical entries from all isoforms
def get_ref_by_gene(gene, feature, out_dir, bed=bed_path, encoding="utf-8", mode="col0",
                    ref_fasta_files=ref_fasta, complete=False, domain="", domain_f="", ref_pref="Col-0_ref",
                    start_inc=True, end_inc=True, merge=False, translate=False, adj_dir=False,
                    attribute_fields = fields["gff3"], by_gene=False, attribute_mod={}, genetic_code=1,
                    store_fasta = None, store_bed = None, **kwargs):
    ## strip quotes from 'gene' if present
    if gene:
        gene = [y for y in [x for x in gene.split('"') if len(x) > 0][0].split("'") if len(y) > 0][0]
    else:
        gene = gene
    ## get relevant gffbed entries
    original_feature = feature
    if feature == "protein":
        feature = "CDS"
        translate = True
        adj_dir = True
        complete = False
    if domain_f:
        feature = "CDS"
    import re
    attribute_mod = attribute_mod if isinstance(attribute_mod, dict) else {} if not re.search("(?<=').*(?=')", attribute_mod) else eval(re.search("(?<=').*(?=')", attribute_mod).group(0))
    data = grep_bedmerge(gene, bed, feature, out_dir, encoding = encoding, merge = merge, attribute_fields = attribute_fields, attribute_mod = attribute_mod, store_bed = store_bed)["data"]
    if not data:
        return
    ## get chrom, strand, and max boundaries
    chrom = data[0][0]
    start = min([int(x[1]) for x in data])
    end = max([int(x[2]) for x in data])
    strand = data[0][3]
    ## extract sequences from fasta file
    get_attribute = make_get_attribute(attribute_fields)
    if merge or (feature == "gene"):
        isoforms = {gene: data}
    elif feature == "mRNA" or feature == "protein":
        isoforms = {isoform: [x for x in data if ((isoform in get_attribute(x[-1], "Parent"))
                                                  or (isoform in get_attribute(x[-1], "ID")))]
                    for isoform in set(itertools.chain(*[get_attribute(x[-1], "ID") for x in data]))}
    else:
        isoforms = {isoform: [x for x in data if isoform in get_attribute(x[-1], "Parent")]
                    for isoform in set(itertools.chain(*[get_attribute(x[-1], "Parent") for x in data]))}
    isoforms = {k: sorted(v, key = lambda x: get_gffbed(x, "start"))
                for k, v in isoforms.items()}
    fasta_out_l = []
    seq_ranges = {}
    ## get fasta file of sequence data
    fasta_out = os.path.join(out_dir, (gene + "_ref_" + original_feature + ("_complete" if complete else '') + \
                                       (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                       ("_protein" if (translate and (feature=="CDS")) else '') + \
                                       ".fasta").replace('|', '_'))
    get_ref_raw(chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files,
                mode=mode, ref_pref=ref_pref)
    from fasta_manip import fasta_to_dict
    ref_seq_original = list(fasta_to_dict(fasta_out).values())[0]
    ## iterate through isoforms
    for isoform, isoform_dat in isoforms.items():
        ref_seq = ref_seq_original
        ## get fasta file of sequence data
        fasta_out = os.path.join(out_dir, (isoform + "_ref_" + original_feature + \
                                           ("_complete" if complete else '') + \
                                           (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                           ("_protein" if (translate and (feature=="CDS")) else '') + \
                                           ".fasta").replace('|', '_'))
        ## if extracting only domain-specific range
        if domain_f:
            domain_data = get_domain_in_genome_coords(gene, domain, domain_f, out_dir,
                                                      bed=bed, encoding=encoding,
                                                      isoform = isoform, attribute_fields = attribute_fields,
                                                      attribute_mod = attribute_mod,
                                                      start_inc = start_inc, end_inc = end_inc,
                                                      store_bed = store_bed,
                                                      **{k: v for k, v
                                                         in kwargs.items()
                                                         if k in ["qname_dname", "qstart_qend"]})
            if (not domain_data):
                continue
        else:
            domain_data = [(start, end)]
        seqs_to_write = {}
        for i, domain_range in enumerate(domain_data):
            d_start, d_end = domain_range
            ranges = [(d_start, d_end)]
            ## trim sequence if complete flag not raised or if domain required
            if (not complete) or domain_f:
                if complete and domain_f and d_start and d_end:
                    ranges = [(max(start, d_start) - start, min(end, d_end) - start)]
                elif domain_f and d_start and d_end:
                    ranges = [(max(int(x[1]), d_start) - start, min(int(x[2]), d_end) - start) \
                              for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
                else:
                    ranges = [(int(x[1])-start, int(x[2])-start) for x in isoform_dat]
                from fasta_manip import extract_ranges
                ref_seq = extract_ranges(ref_seq_original, ranges)
            if (adj_dir or translate) and strand == '-':
                ref_seq = ref_seq.reverse_complement()
            ## translate sequence if translate flag raised AND feature is CDS
            if translate:
                if feature == "CDS" and not complete:
                    ref_seq = ref_seq.translate(to_stop = True, table = parse_genetic_code(genetic_code))
                else:
                    print("Translation is only possible when the selected feature is 'CDS' and the flag 'complete' is not raised.")
            # seq_name = "{}|{}|{}|{}".format("Col-0_ref" if mode=="col0" else "Reference",
            #                                 gene, original_feature, isoform) + \
            #            (('|' + ("domain" if not domain else domain) + f"|{i+1}") if domain_f else '') + \
            #            ("|complete" if complete else '') + \
            #            ("|revcomp" if adj_dir and strand == '-' else '')
            seq_name = ("Col-0_ref" if mode=="col0" else "Reference", gene, original_feature, isoform) + \
                       ( (("domain" if not domain else domain), f"{i+1}") if domain_f else () ) + \
                       ( ("complete",) if complete else () ) + \
                       ( ("revcomp",) if adj_dir and strand == '-' else () )
            seqs_to_write[seq_name] = ref_seq
            ## for by_gene
            if by_gene:
                overlap_ranges = []
                overlap_seq_names = []
                ## iterate through processed ranges
                for logged_ranges, logged_seq_names in seq_ranges.items():
                    ## if new range overlaps with already processed ranges
                    if has_any_overlap(ranges, logged_ranges):
                        ## note processed ranges that overlap with new range
                        overlap_ranges.append(logged_ranges)
                        overlap_seq_names.extend(logged_seq_names)
                ## if the new range overlaps w/ any of the processed ranges
                if overlap_ranges:
                    ## replace overlapped processed ranges w/ new merged range
                    for logged_ranges in overlap_ranges:
                        del(seq_ranges[logged_ranges])
                    ranges = merge_ranges(*overlap_ranges, ranges)
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name] + overlap_seq_names
            else:
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name]
        if seqs_to_write:
            from fasta_manip import dict_to_fasta
            dict_to_fasta(seqs_to_write, fasta_out)
            fasta_out_l.append(fasta_out)
    fasta_out_final = os.path.join(out_dir, (gene + "_ref_" + original_feature + \
                                             ("_complete" if complete else '') + \
                                             (('_' + ("domain" if not domain else domain))
                                              if domain_f else '') + \
                                             ("_protein" if (translate and (feature=="CDS") and
                                                             not complete) else '')+\
                                             ".fasta").replace('|', '_'))
    if fasta_out_l:
        if (not by_gene) and (fasta_out_l[0] != fasta_out_final):
            from file_manip import cat_files
            cat_files(sorted(fasta_out_l), fasta_out_final)
            # for fasta_out in fasta_out_l:
            #     os.remove(fasta_out)
        if by_gene:
            isoform_seqs = fasta_to_dict(fasta_out_final)
            final_seqs = {}
            i = 0
            ## get fasta file of sequence data
            fasta_out = os.path.join(out_dir, (isoform + "_ref_" + original_feature +
                                               ("_complete" if complete else '') + \
                                               (('_' + ("domain" if not domain else domain))
                                                if domain_f else '') + \
                                               ("_protein" if (translate and (feature=="CDS")) else '') + \
                                               ".fasta").replace('|', '_'))
            get_ref_raw(chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files,
                        mode=mode, ref_pref=ref_pref)
            for ranges, seq_names in sorted(seq_ranges.items()):
                # seq_name_l = seq_names[0].split('|')
                seq_name_l = list(seq_names[0])
                # seq_name_l_prefgid, seq_name_l_post = seq_names[0].split(f"|{original_feature}|")
                # seq_name_l_tid = seq_name_l_post.split("|revcomp")
                # seq_name_l = [seq_name_l_pre.split('|')[0], '|'.join(seq_name_l_pre.split('|')[1:]),
                #               original_feature, seq_name_l_tid] + \
                #               (["revcomp"] if re.search("revcomp$", seq_names[0]) else [])
                seq_name_l[3] = ','.join(f'{r[0]}-{r[1]}' for r in ranges)
                if domain_f:
                    seq_name_l[5] = str(i + 1)
                seq_name = '|'.join(seq_name_l)
                # final_seqs[seq_name] = isoform_seqs[seq_names[0]]
                seq_out = ''.join([str(ref_seq_original[r[0]:r[1]]) for r in sorted(ranges)])
                if (adj_dir or translate) and strand == '-':
                    from Bio import Seq
                    seq_out = str(Seq.Seq(seq_out).reverse_complement())
                final_seqs[seq_name] = seq_out
                i += 1
            dict_to_fasta(final_seqs, fasta_out_final)
        ## remove intermediate files
        for fasta_out in fasta_out_l:
            if fasta_out != fasta_out_final:
                os.remove(fasta_out)
        print("{}\tSequences were successfully written to: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), fasta_out_final))
    elif not fasta_out_l:
        f = open(fasta_out_final, "w+")
        f.write('')
        f.close()
        print("{}\t{} is an empty file".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), fasta_out_final))
    store_fname(store_fasta, fasta_out_final)

def get_ref_by_range(chrom, start, end, out_dir, encoding="utf-8",
                     ref_fasta_files=ref_fasta, mode="col0", ref_pref="Col-0_ref",
                     store_fasta = None):
    fasta_out = os.path.join(out_dir, "chr{}_p{}-{}_ref.fasta".format(chrom, start + 1, end))
    get_ref_raw(chrom_pref + chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files,
                mode=mode, ref_pref=ref_pref)
    store_fname(store_fasta, fasta_out)
    print("Sequences were successfully written to: {}".format(fasta_out))

## output is...0-indexed, I think...it seems like it follows .bed conventions...
def get_domain_in_genome_coords(gene, domain, domain_f, out_dir, pos_type="aa", isoform='',
                                bed=bed_path, encoding="utf-8", start_inc = True, end_inc = True,
                                attribute_fields = fields["gff3"], attribute_mod = {}, store_bed = None,
                                qname_dname=("name", "domain name"), qstart_qend=("start", "end")):
    
    qname, dname = qname_dname
    qstart, qend = qstart_qend
    domain_raw = [x.split('\t') for x in splitlines(domain_f)]
    domain_header = domain_raw[0]
    domain_data = [x for x in domain_raw[1:] if len(domain_raw) > 1 and len(x) == len(domain_header) and \
                   (((not domain) or domain == x[domain_header.index(dname)]) and \
                    re.search(f"\|{re.escape(gene if not isoform else isoform)}(\||$)",
                              x[domain_header.index(qname)]))]
                    # re.search(gene if not isoform else isoform) in x[domain_header.index(qname)].split('|'))]
    output = []
    for domain_dat in domain_data:
        ## convert to 0-index, start inclusive, stop exclusive + adjust for unit (e.g. aa or nt)
        def adj_pos(start, end):
            account_unit = lambda x: x * (3 if pos_type == "aa" else 1)
            return (account_unit(int(start) - (1 if start_inc else 0)),
                    account_unit(int(end) - (0 if end_inc else 1)))
        domain_start, domain_end = adj_pos(int(domain_dat[domain_header.index(qstart)]),
                                           int(domain_dat[domain_header.index(qend)]))
        ## get CDS to define boundaries of translated nucleotides
        cds_dat = grep_bedmerge((gene if not isoform else isoform), bed, "CDS", out_dir, encoding = encoding,
                                attribute_fields = attribute_fields, attribute_mod = attribute_mod,
                                store_bed = store_bed)["data"]
        if len(cds_dat) == 0:
            continue
        ## extract boundaries
        chrom = cds_dat[0][0]
        plus_strand = (cds_dat[0][3] == '+')
        bounds = sorted([(int(x[1]), int(x[2])) for x in cds_dat],
                        key = lambda x: x[0], reverse = (not plus_strand))
        last_end = 0
        genome_start, genome_end = None, None
        for i, v in enumerate(bounds):
            curr_cds_start = last_end
            curr_cds_end = last_end + (v[1] - v[0])
            get_g_pos = lambda x: ((v[0] + (x - curr_cds_start)) if plus_strand else\
                                   (v[1] - (x - curr_cds_start)))
            if curr_cds_start <= domain_start < curr_cds_end:
                genome_start = get_g_pos(domain_start)
            if curr_cds_start < domain_end <= curr_cds_end: ## recall that end is exclusive
                genome_end = get_g_pos(domain_end)
            last_end = curr_cds_end
        # if genome_start and genome_end:
        #     output.append((min(genome_start, genome_end), max(genome_start, genome_end)))
        output.append((min(genome_start, genome_end), max(genome_start, genome_end)))
    return output

# ## get sequence info of genes
# def get_1001pseudogenome(accs, gene, feature, out_dir, bed=bed_path, encoding="utf-8",
#                          acc_ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
#     accs = [acc for acc in accs if acc]
#     if str(accs[0]).isdigit():
#         accs = dict(get_acc_name(accs, acc_ref_file=acc_ref_file, num_col=num_col, name_col=name_col,
#                                  delim=delim, header=header))
#     else:
#         accs = dict(get_acc_num(accs, acc_ref_file=acc_ref_file, num_col=num_col, name_col=name_col,
#                                 delim=delim, header=header))
#     ## get chrom, start, end of gene's features
#     grep = subprocess.Popen(("grep", gene, bed), stdout=subprocess.PIPE)
#     awk = subprocess.Popen(("awk", "$8 == \"" + feature + "\""), stdin=grep.stdout, stdout=subprocess.PIPE)
#     merged = [entry.split('\t') for entry in
#               subprocess.check_output(("bedtools", "merge"),
#                                       stdin=awk.stdout).decode(encoding).split('\n') if entry]
#     ## write bed file of regions used (0-indexed)
#     f = open(os.path.join(out_dir, gene + '_' + feature + ".bed"), "w+")
#     f.write('\n'.join(['\t'.join(x) for x in merged]))
#     f.close()
#     ## extract sequences
#     chrom = merged[0][0]
#     start = min([int(x[1]) for x in merged]) + 1
#     end = max([int(x[2]) for x in merged])
#     ## get fasta file of sequence data
#     fasta_out = os.path.join(out_dir, gene + '_' + feature + ".fasta")
#     subprocess.run(("curl", "-o", fasta_out,
#                     make_1001pseudogenome_api(chrom, start, end, accs.keys())))
#     ## trim sequences in file
#     from fasta_manip import extract_feature
#     seqs_trimmed = extract_feature(fasta_out, bed, tuple(accs.keys())[0], gene, start - 1, fasta_out, feature,
#                                    encoding = encoding)
#     from fasta_manip import fasta_to_dict
#     seqs = fasta_to_dict(fasta_out)
#     seqs_fin = {k.split('|')[3] + '|' + accs[k.split('|')[3]] + '|' + gene + '|' + feature:v \
#                 for k,v in seqs.items()}
#     from fasta_manip import dict_to_fasta
#     dict_to_fasta(seqs_fin, fasta_out)
#     print("Sequences were successfully written to: {}".format(fasta_out))

    
# def get_ref(gene, feature, out_dir, bed=bed_path, encoding="utf-8",
#             ref_fasta_files=ref_fasta, include_within_features=False):
#     ## get chrom, start, end of gene's features
#     grep = subprocess.Popen(("grep", gene, bed), stdout=subprocess.PIPE)
#     awk = subprocess.Popen(("awk", "$8 == \"" + feature + "\""), stdin=grep.stdout, stdout=subprocess.PIPE)
#     merged = [entry.split('\t') for entry in
#               subprocess.check_output(("bedtools", "merge"),
#                                       stdin=awk.stdout).decode(encoding).split('\n') if entry]
#     ## write bed file of regions used (0-indexed)
#     f = open(os.path.join(out_dir, gene + '_' + feature + ".bed"), "w+")
#     f.write('\n'.join(['\t'.join(x) for x in merged]))
#     f.close()
#     ## extract sequences from fasta file
#     chrom = merged[0][0]
#     start = min([int(x[1]) for x in merged])
#     end = max([int(x[2]) for x in merged])
#     ## get fasta file of sequence data
#     fasta_out = os.path.join(out_dir, gene + '_' + feature + ("_includewithinfeatures" if include_within_features else '') + ".fasta")
#     from fasta_manip import fasta_to_dict
#     ref_seq = list(fasta_to_dict(ref_fasta_files[chrom]).values())[0][start:end]
#     if not include_within_features:
#         ranges = [(int(x[1])-start, int(x[2])-start) for x in merged]
#         from fasta_manip import extract_ranges
#         ref_seq = extract_ranges(ref_seq, ranges)
#     from fasta_manip import dict_to_fasta
#     dict_to_fasta({"Col-0_ref|{}|{}".format(gene, feature): ref_seq}, fasta_out)
#     print("Sequences were successfully written to: {}".format(fasta_out))    
