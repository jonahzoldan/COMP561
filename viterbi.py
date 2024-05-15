# Import packages


import pprint
# from BCBio.GFF import GFFExaminer
from gff3 import Gff3
import gffutils
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=np.inf)

vibrio_cholerae_annot = "C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\COMP 561\\Vibrio_cholerae.GFC_11.37.gff3"
fn = gffutils.example_filename(vibrio_cholerae_annot)
print('Running program')
print('----------------------------------------------------------------------------------- \n \
----------------------------------------------------------------------------------- \n \
----------------------------------------------------------------------------------- ')
# print(open(fn).read())

# Create database for all features
vibrio_cholerae_db = gffutils.create_db(fn, dbfn='vc.db', force=True, keep_order=True,
                                        merge_strategy="warning", sort_attribute_values=True)

# Database of all features
vc_db = gffutils.FeatureDB('vc.db')
vc_db_all = vc_db.all_features()

# Create database of CDS features
vc_db_cds = vc_db.features_of_type('CDS')

# Create list of CDS features on positive strand
vc_cds_feats = []
for feat in vc_db_cds:
    if feat.strand == '+':
        vc_cds_feats.append(feat)


# Now we look at the fasta file
vibrio_cholerae_fa = "C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\COMP 561\\Vibrio_cholerae.GFC_11.dna.nonchromosomal.fa"

vc_fa = open(vibrio_cholerae_fa)
# print(vc_fa.read())
vc_fa_parsed = SeqIO.parse(vc_fa, 'fasta')
vc_fa_dict = SeqIO.to_dict(vc_fa_parsed)


# Create dictionary of CDS features
# key = contig ID and value = list of features

vc_cds_dict = {}
for vc_feat in vc_cds_feats:
    if vc_feat.seqid not in vc_cds_dict.keys():
        vc_cds_dict[vc_feat.seqid] = [vc_feat]
    else:
        vc_cds_dict[vc_feat.seqid].append(vc_feat)


inter_start = 0
nb_inters = 0
len_inters = 0
nb_genics = 0
len_genics = 0
nb_a = 0
nb_c = 0
nb_g = 0
nb_t = 0
len_seqs = 0

vc_codon_freqs = {}
vc_start_freqs = {}
vc_stop_freqs = {}

def nucleotide_frequency(sequence):
    a = 0
    c = 0
    g = 0
    t = 0
    for nuc in sequence:
        if nuc == 'A':
            a += 1
        elif nuc == 'C':
            c += 1
        elif nuc == 'G':
            g += 1
        else:
            t += 1
    return a, c, g, t



for contig in vc_fa_dict:
    inter_start = 0
    # contig = 'DN38.contig00084'
    fa_seq = vc_fa_dict[contig].seq
    len_seq = len(fa_seq)
    len_seqs += len_seq

    # not all contigs fit the requirements of being CDS and on the positive strand
    if contig not in vc_cds_dict.keys():
        len_inters += len_seq
        nb_inters += 1
        new_a, new_c, new_g, new_t = nucleotide_frequency(fa_seq)
        nb_a += new_a
        nb_c += new_c
        nb_g += new_g
        nb_t += new_t
        continue

    for vc_feat in vc_cds_dict[contig]:

        len_genics = len_genics + (vc_feat.end - vc_feat.start + 1)
        nb_genics += 1

        # Codon frequencies
        for g in range(vc_feat.start-1, vc_feat.end, 3):
            cur_codon = fa_seq[g:g+3]

            if g == vc_feat.start-1:
                if cur_codon not in vc_start_freqs.keys():
                    vc_start_freqs[cur_codon] = 1
                else:
                    vc_start_freqs[cur_codon] += 1
            if g == vc_feat.end - 3:
                if cur_codon not in vc_stop_freqs.keys():
                    vc_stop_freqs[cur_codon] = 1
                else:
                    vc_stop_freqs[cur_codon] += 1

            if cur_codon not in vc_codon_freqs.keys():
                vc_codon_freqs[cur_codon] = 1
            else:
                vc_codon_freqs[cur_codon] += 1


        # if current feature is not last in config
        last_feat = vc_cds_dict[contig][len(vc_cds_dict[contig])-1]
        if vc_feat != last_feat:

            len_inters = len_inters + (vc_feat.start - 1) - inter_start
            if vc_feat.start - 1 - inter_start > 0:
                nb_inters += 1

            new_a, new_c, new_g, new_t = nucleotide_frequency(fa_seq[inter_start:vc_feat.start-1])
            nb_a += new_a
            nb_c += new_c
            nb_g += new_g
            nb_t += new_t

            cur_index = vc_cds_dict[contig].index(vc_feat)
            after_index = cur_index + 1

            if vc_feat.end - 1 >= vc_cds_dict[contig][after_index].start - 1:
                inter_start = vc_cds_dict[contig][after_index].start - 1
            else:
                inter_start = vc_feat.end

        # if current feature IS last in config
        else:
            # if current feature extends until end of sequence
            if vc_feat.end == len_seq:
                len_inters = len_inters + (vc_feat.start - 1) - inter_start
                if vc_feat.start - 1 - inter_start > 0:
                    nb_inters += 1

                new_a, new_c, new_g, new_t = nucleotide_frequency(fa_seq[inter_start:vc_feat.start-1])
                nb_a += new_a
                nb_c += new_c
                nb_g += new_g
                nb_t += new_t

            # if current feature stops before end of sequence
            else:
                len_inters = len_inters + (vc_feat.start - 1) - inter_start
                if vc_feat.start - 1 - inter_start > 0:
                    nb_inters += 1
                len_inters = len_inters + len_seq - vc_feat.end - 1
                nb_inters += 1

                # nucleotide frequencies before gene
                new_a, new_c, new_g, new_t = nucleotide_frequency(fa_seq[inter_start:vc_feat.start-1])
                nb_a += new_a
                nb_c += new_c
                nb_g += new_g
                nb_t += new_t

                # nucleotide frequencies after gene until end
                new_a, new_c, new_g, new_t = nucleotide_frequency(fa_seq[vc_feat.end:len_seq])
                nb_a += new_a
                nb_c += new_c
                nb_g += new_g
                nb_t += new_t

# statistics:
# ---------------------------------------------------------------------------
# average length of intergenic:
avg_inter = len_inters / float(nb_inters)
print('Average intergenic length (in nucleotides) = ', avg_inter, '\n')
#
# # average length of genic:
avg_genic = len_genics / float(nb_genics)
print('Average genic length (in nucleotides) = ', avg_genic, '\n')
print('Average genic length (in codons) = ', avg_genic/3, '\n')
#
# # frequency of nucleotides in intergenic:
a_freq = nb_a / float(len_inters)
c_freq = nb_c / float(len_inters)
g_freq = nb_g / float(len_inters)
t_freq = nb_t / float(len_inters)
print(f'Nucleotide frequences: \n A: {a_freq} \n C: {c_freq} \n G: {g_freq} \n T: {t_freq}')

# Turn codon total into codon frequency
for codon in vc_codon_freqs:
    vc_codon_freqs[codon] /= (len_genics/3.0)

# Turn start and stop total into frequencies
for start in vc_start_freqs:
    vc_start_freqs[start] /= float(nb_genics)

for stop in vc_stop_freqs:
    vc_stop_freqs[stop] /= float(nb_genics)

# There is no chance of a stop codon within a genic region
for codon in vc_codon_freqs:
    if codon in vc_stop_freqs.keys():
        vc_codon_freqs[codon] = 0

# The start dictionary should have a 0 set for all non-start codons:
for codon in vc_codon_freqs:
    if codon not in vc_start_freqs.keys():
        vc_start_freqs[codon] = 0

# The stop dictionary should have 0 set for all non-stop codons:
for codon in vc_codon_freqs:
    if codon not in vc_stop_freqs.keys():
        vc_stop_freqs[codon] = 0

print('All codons: ', vc_codon_freqs)
print('Start codons: ', vc_start_freqs)
print('Stop codons: ', vc_stop_freqs)


# Turn nucleotide frequencies into dictionary
nuc_freqs = {}
nuc_freqs['A'] = a_freq
nuc_freqs['C'] = c_freq
nuc_freqs['G'] = g_freq
nuc_freqs['T'] = t_freq

codon_genic_freq_table = np.empty((len(vc_codon_freqs), 2))

transition_mat = np.zeros((4,4))
i = 0
for codon in vc_codon_freqs:
    print(codon, ' : ', vc_codon_freqs[codon])

# from intergenic
transition_mat[0,1] = 1/avg_inter # to start
transition_mat[0,0] = 1 - transition_mat[0,1] # to itself
transition_mat[0,2] = 0 # to genic
transition_mat[0,3] = 0 # to stop codon

# from start
transition_mat[1,1] = 0
transition_mat[1,2] = 1
transition_mat[1,3] = 0
transition_mat[1,0] = 0

# from genic
transition_mat[2,0] = 0
transition_mat[2,1] = 0
transition_mat[2,3] = 1/(avg_genic / 3)
transition_mat[2,2] = 1 - transition_mat[2,3]

# from stop
transition_mat[3,0] = 1
transition_mat[3,1] = 0
transition_mat[3,2] = 0
transition_mat[3,3] = 0

print('Intergenic -> Start: ', transition_mat[0,1])
print('Intergenic -> Intergenic: ', transition_mat[0,0])
print('Start -> Genic: ', transition_mat[1,2])
print('Genic -> Genic: ', transition_mat[2,2])
print('Genic -> Stop: ', transition_mat[2,3])
print('Stop -> Intergenic: ', transition_mat[3,0])


def score(prev_prob, emission, transition):
    log_emission = 0
    log_transition = 0
    if emission != 0:
        log_emission = np.log(emission)
    if transition != 0:
        log_transition = np.log(transition)

    if emission != 0 and transition != 0:
        return prev_prob + log_emission + log_transition
    else:
        return np.NINF


def unify_sequences(fasta_dict):
    unified = ""

    for seq in fasta_dict:
        cur_seq = fasta_dict[seq].seq
        unified = unified + cur_seq

    return unified

def get_sequence_indices(fasta_dict):
    seq_indices = {}
    running_index = 0

    for seq in fasta_dict:
        seq_indices[seq] = (running_index + 1, len(fasta_dict[seq].seq) + running_index)
        running_index = running_index + len(fasta_dict[seq].seq)

    return seq_indices


# Try looking backwards in this case

def viterbi2(fasta_dict):

    genome = unify_sequences(fasta_dict)

    indices_dict = get_sequence_indices(fasta_dict)

    total_predictions = ""

    for c in range(1):

        print('Now filling out viterbi matrix...')
        # Initialize viterbi matrix
        total_len = len(genome)
        # 0 -> Intergenic
        # 1 -> Start
        # 2 -> Gene
        # 3 -> Stop
        value = np.empty((), dtype=object)
        value[()] = (0, 0)
        v = np.full((4, total_len), value, dtype=object)

        # Initialize first column
        v[0, 0] = np.log(1), (-1, -1)
        v[1, 0] = np.NINF, (-1, -1)
        v[2, 0] = np.NINF, (-1, -1)
        v[3, 0] = np.NINF, (-1, -1)
        # We always start in intergenic state

        for i in range(1, total_len):
            cur_nuc = genome[i]

            if i == total_len/4:
                print('25% done filling viterbi matrix!')

            if i == total_len/2:
                print('50% done filling viterbi maxtrix!')

            if i == total_len - (total_len/4):
                print('75% done filling viterbi matrix')

            # Calculate for intergenic
            inter_best_score = np.NINF
            inter_best_dir = 0,0

            for prev_state in range(4):
                inter_cur_score = score(v[prev_state, i-1][0], nuc_freqs[cur_nuc], transition_mat[prev_state, 0])

                if inter_cur_score > inter_best_score:
                    inter_best_score = inter_cur_score
                    inter_best_dir = prev_state, i-1

            v[0, i] = inter_best_score, inter_best_dir

            # Calculate for start
            start_best_score = np.NINF
            start_best_dir = 0,0


            for prev_state in range(4):
                if i > 2:
                    start_cur_score = score(v[prev_state, i-3][0], vc_start_freqs[genome[i-2:i+1]], transition_mat[prev_state, 1])
                else:
                    start_best_score = np.NINF
                    break

                if start_cur_score > start_best_score:
                    start_best_score = start_cur_score
                    start_best_dir = prev_state, i-3

            v[1, i] = start_best_score, start_best_dir

            # Calculate for genic
            genic_best_score = np.NINF
            genic_best_dir = 0,0


            for prev_state in range(4):
                if i > 2:
                    genic_cur_score = score(v[prev_state, i-3][0], vc_codon_freqs[genome[i-2:i+1]], transition_mat[prev_state, 2])
                else:
                    genic_best_score = np.NINF
                    break

                if genic_cur_score > genic_best_score:
                    genic_best_score = genic_cur_score
                    genic_best_dir = prev_state, i-3

            v[2, i] = genic_best_score, genic_best_dir

            # Calculate for stop
            stop_best_score = np.NINF
            stop_best_dir = 0,0

            for prev_state in range(4):
                if i > 2:
                    stop_cur_score = score(v[prev_state, i-3][0], vc_stop_freqs[genome[i-2:i+1]], transition_mat[prev_state, 3])
                else:
                    stop_best_score = np.NINF
                    break

                if stop_cur_score > stop_best_score:
                    stop_best_score = stop_cur_score
                    stop_best_dir = prev_state, i-3

            v[3, i] = stop_best_score, stop_best_dir


        # Now we backtrack
        best_final_score = np.NINF
        best_final_state = 0,0
        for state in range(4):
            if v[state, total_len-1][0] > best_final_score:
                best_final_score = v[state, total_len-1][0]
                best_final_state = state, total_len - 1

        print('100% done filling viterbi matrix!')
        print("Now backtracking...")

        indices_list = []
        cur_row, cur_col = best_final_state
        cur_stop = None

        while True:
            if cur_row == 3:
                cur_stop = cur_col + 1
            if cur_row == 1:
                if cur_stop is not None:
                    indices_list.insert(0, (cur_col - 1, cur_stop))
                else:
                    indices_list.insert(0, (cur_col - 1, total_len + 5))

            cur_row, cur_col = v[cur_row, cur_col][1]
            if cur_row == -1 and cur_col == -1:
                break

            halfway = False
            if cur_col < total_len/2 and halfway:
                print('Halfway done backtracking!')
                halfway = True

        print('Done backtracking! Producing predictions...')
        predictions = ""

        for gene in indices_list:
            pred_start, pred_end = gene
            for this_contig in indices_dict:
                this_start, this_stop = indices_dict[this_contig]
                if pred_start >= this_start and pred_end <= this_stop:
                    predictions += f"{this_contig}\tena\tCDS\t{pred_start - this_start + 1}\t{pred_end - this_start + 1}\t.\t+\t0\t.\n"
                    break
                if pred_start >= this_start and pred_start < this_stop:
                    predictions += f"{this_contig}\tena\tCDS\t{pred_start - this_start + 1}\tnext contig\t.\t+\t0\t.\n"
                    break

        total_predictions += predictions

    return total_predictions


### The below code runs viterbi on vibrio cholerae and outputs the predictions to a new file

vc_output_preds = open('C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\\COMP 561\\vc_predictions.gff', 'w')
vc_predictions = viterbi2(vc_fa_dict)
vc_output_preds.write(vc_predictions)
vc_output_preds.close()

# Now we run viterbi on vibrio vulnificus
vibrio_vulnificus = "C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\COMP 561\\Vibrio_vulnificus.ASM74310v1.dna.nonchromosomal.fa"
vv_fa = open(vibrio_vulnificus)
vv_fa_parsed = SeqIO.parse(vv_fa, 'fasta')
vv_fa_dict = SeqIO.to_dict(vv_fa_parsed)

### The below code runs viterbi on vibrio vulnificus and outputs the predictions to a new file
vv_output_preds = open('C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\\COMP 561\\vv_predictions.gff', 'w')
vv_predictions = viterbi2(vv_fa_dict)
vv_output_preds.write(vv_predictions)
vv_output_preds.close()




# Format vv GFF in same way as vc GFF
vibrio_vulnificus_annot = "C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\COMP 561\\Vibrio_vulnificus.ASM74310v1.37.gff3"
vv = gffutils.example_filename(vibrio_vulnificus_annot)


# Create database for all features
vibrio_cholerae_db = gffutils.create_db(vv, dbfn='vv.db', force=True, keep_order=True,
                                        merge_strategy="warning", sort_attribute_values=True)

vv_db = gffutils.FeatureDB('vv.db')
vv_db_cds = vv_db.features_of_type('CDS')

# Generate list of positive strand genes
vv_cds_feats = []
for feat in vv_db_cds:
    if feat.strand == '+':
        vv_cds_feats.append(feat)

# Split each line of the output file into a comma separated list

vv_output_preds_lists = open('C:\\Users\\jonah\\OneDrive\\Desktop\\F2022\\\COMP 561\\vv_predictions.gff', 'r')

def grade_viterbis():
    '''
    No arguments necessary. The function only evaluates on the vibrio vulnificus dataset
    :return: a list of missed or partially missed annotated genes and a list of false positive predictions
    '''
    annot_both_ends = 0
    annot_starts = 0
    annot_ends = 0
    annot_neither = 0
    # Fraction of annotated genes on positive strand that:
    vv_output_lines = vv_output_preds_lists.readlines()
    vv_missed_genes = []

    for annot_gene in vv_cds_feats:
        miss = True
        for line in vv_output_lines:
            line = line.split('\t')
            if line[0] == str(annot_gene.seqid):

                if line[3] == str(annot_gene.start) and line[4] == str(annot_gene.end):
                    annot_both_ends += 1
                    miss = False
                    break
                if line[3] == str(annot_gene.start) and line[4] != str(annot_gene.end):
                    annot_starts += 1
                    miss = False
                    vv_missed_genes.append(annot_gene)
                    break
                if line[3] != str(annot_gene.start) and line[4] == str(annot_gene.end):
                    annot_ends += 1
                    miss = False
                    vv_missed_genes.append(annot_gene)
                    break
        if miss:
            vv_missed_genes.append(annot_gene)
            annot_neither += 1

    pred_both_ends = 0
    pred_starts = 0
    pred_ends = 0
    pred_neither = 0

    vv_false_positives = []

    nb_vv_preds = 0
    for pred_gene in vv_output_lines:
        pred_gene = pred_gene.split('\t')
        nb_vv_preds += 1
        miss = True
        for annot_gene in vv_cds_feats:
            if pred_gene[0] == str(annot_gene.seqid):

                if str(annot_gene.start) == pred_gene[3] and str(annot_gene.end) == pred_gene[4]:
                    pred_both_ends += 1
                    miss = False
                    break

                if str(annot_gene.start) == pred_gene[3] and str(annot_gene.end) != pred_gene[4]:
                    pred_starts += 1
                    miss = False
                    if pred_gene[4] == 'next contig':
                        vv_false_positives.append(pred_gene)
                    break

                if str(annot_gene.start) != pred_gene[3] and str(annot_gene.end) == pred_gene[4]:
                    pred_ends += 1
                    miss = False
                    break
        if miss:
            vv_false_positives.append(pred_gene)
            pred_neither += 1

    annot_both_ends /= float(len(vv_cds_feats))
    annot_starts /= float(len(vv_cds_feats))
    annot_ends /= float(len(vv_cds_feats))
    annot_neither /= float(len(vv_cds_feats))

    print('Fraction of annotated genes that:')
    print('\t match both ends of a predicted gene', annot_both_ends)
    print('\t match only the start of a predicted gene', annot_starts)
    print('\t match only the end of a predicted gene', annot_ends)
    print('\t are completely missed by a predicted gene', annot_neither)

    pred_both_ends /= float(nb_vv_preds)
    pred_starts /= float(nb_vv_preds)
    pred_ends /= float(nb_vv_preds)
    pred_neither /= float(nb_vv_preds)

    print('Fraction of predicted genes that:')
    print('\t match both ends of an annotated gene', pred_both_ends)
    print('\t match only the start of an annotated gene', pred_starts)
    print('\t match only the end of an annotated gene', pred_ends)
    print('\t completely miss an annotated gene', pred_neither)

    return vv_missed_genes, vv_false_positives

vv_missed_genes, vv_false_positives = grade_viterbis()

def evaluate_mistakes(missed_genes, false_positives):
    missed_dict = {}
    for annot in missed_genes:
        if annot.seqid not in missed_dict.keys():
            missed_dict[annot.seqid] = 1
        else:
            missed_dict[annot.seqid] += 1

    missed_list = sorted(missed_dict.items(), key=lambda x: x[1], reverse=True)
    contigs1 = []
    frequencies1 = []
    for gene in missed_list:
        contigs1.append(gene[0])
        frequencies1.append(gene[1])

    fig1, ax1 = plt.subplots()
    bars = ax1.barh(contigs1[0:20], frequencies1[0:20])
    ax1.bar_label(bars)
    plt.xlabel('Frequency')
    plt.title('Frequency of contig for missed genes (Top 20)')
    # plt.show()

    false_pos_dict = {}
    for pred in false_positives:
        if pred[0] not in false_pos_dict.keys():
            false_pos_dict[pred[0]] = 1
        else:
            false_pos_dict[pred[0]] += 1

    false_list = sorted(false_pos_dict.items(), key=lambda x: x[1], reverse=True)
    contigs2 = []
    frequencies2 = []
    for pred in false_list:
        contigs2.append(pred[0])
        frequencies2.append(pred[1])

    fig2, ax2 = plt.subplots()
    bars = ax2.barh(contigs2[0:20], frequencies2[0:20], color='Orange')
    ax2.bar_label(bars)
    plt.xlabel('Frequency')
    plt.title('Frequency of contig for false positive predictions (Top 20)')

    # Analyze frequency of start and stop codons
    start_codons = {'ATG', 'GTG', 'TTG'}
    stop_codons = {'TAA', 'TGA', 'TAG'}
    nb_starts_missed = 0
    nb_stops_missed = 0
    genic_length_in_codons = 0
    for missed_gene in missed_genes:
        real_start = missed_gene.start
        real_end = missed_gene.end
        seq = vv_fa_dict[missed_gene.seqid].seq
        genic_length_in_codons += len(seq) / float(3)
        for n in range(real_start-1, real_end, 3):
            if seq[n:n+3] in start_codons:
                nb_starts_missed += 1
            else:
                if seq[n:n+3] in stop_codons:
                    nb_stops_missed += 1

    freq_start_in_missed = nb_starts_missed / float(genic_length_in_codons)
    freq_stop_in_missed = nb_stops_missed / float(genic_length_in_codons)
    print('Frequency of start codons in missed genes: ', freq_start_in_missed)
    print('Frequency of stop codons in missed genes: ', freq_stop_in_missed)

    nb_starts_false = 0
    nb_stops_false = 0
    pred_genic_length_in_codons = 0
    for false_pos in false_positives:
        seq = vv_fa_dict[false_pos[0]].seq
        pred_start = int(false_pos[3])
        if false_pos[4] != 'next contig':
            pred_end = int(false_pos[4])
        else:
            pred_end = len(seq)
        pred_genic_length_in_codons += len(seq) / float(3)
        for n in range(pred_start-1, pred_end, 3):
            if seq[n:n+3] in start_codons:
                nb_starts_false += 1
            else:
                if seq[n:n+3] in stop_codons:
                    nb_stops_false += 1

    freq_start_in_false = nb_starts_false / float(pred_genic_length_in_codons)
    freq_stops_in_false = nb_stops_false / float(pred_genic_length_in_codons)

    print('Frequency of start codons in false positive predictions: ', freq_start_in_false)
    print('Frequency of stop codons in false positive predictions: ', freq_stops_in_false)

    plt.show()


    return

evaluate_mistakes(vv_missed_genes, vv_false_positives)

vv_output_preds_lists.close()
vc_fa.close()
vv_fa.close()