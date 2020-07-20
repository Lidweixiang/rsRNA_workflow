import pandas as pd
import numpy as np
import os
import collections
from itertools import chain
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as mg
from mpl_toolkits.axes_grid1 import ImageGrid
import seaborn as sns
from sklearn.manifold import TSNE


def fastq_dump(sra_file_path):
    '''
    Transfer format from .sra to .fastq use fastq-dump
    :return:NaN
    '''
    print('Now fastq-dump your sra file...')
    fastq_dump_cmd = 'fastq-dump --split-3 ' + sra_file_path
    os.system(fastq_dump_cmd)


def fastqc(fastq_file_path):
    '''
    Quality control of smallRNAseq data use FastQC
    :return:NaN
    '''
    print('Now FastQC your fastq file...')
    fastqc_cmd = 'FastQC --extract ' + fastq_file_path
    os.system(fastqc_cmd)


def recognize_adaptor(fastqc_data_path):
    '''
    Recognize reads adaptor through FastQC result file:fastqc_data.txt
    :return:A list such as ['-a adaptor1','-a adaptor2']
    '''
    print('Now recognize sequencing adaptors...')
    adaptor_list = {
        'Illumina_Universal_Adapter': 'AGATCGGAAGAGCACACGTCT',
        'Illumina_Small_RNA_3_Adapter': 'TGGAATTCTCGGGTGCCAAGG',
        'Illumina_Small_RNA_5_Adapter': 'GATCGTCGGACTGTAGAACTCTGAAC',
        'Nextera_Transposase_Sequence': 'CTGTCTCTTATACACATCT',
        'SOLID_Small_RNA_Adapter': 'CGCCTTGGCCGTACAGCAG',
        'Illumina_3_old_RNA_Adapter': 'TCGTATGCCGTCTTCTGCTTGT',
        'Illumina_v1_5_Small_RNA_3_Adapter': 'ATCTCGTATGCCGTCTTCTGCTTG'
    }
    adaptors = []
    line = 0
    # extract adaptor data from fastqc_data.txt and open with dataframe format
    with open(fastqc_data_path) as data:
        fastqc_datas = data.readlines()
    for data_line in fastqc_datas:
        line += 1
        if data_line[:17] == '>>Adapter Content':
            break
    adaptor_datas = fastqc_datas[line:-1]
    adaptors_file = fastqc_data_path[:-4] + '2.txt'
    with open(adaptors_file, 'w') as af:
        af.writelines(adaptor_datas)
    adaptor_dataframe = pd.read_csv(adaptors_file, sep=r'\s+', index_col='0')
    # set the threshold to 2% ,if the max content of adaptor > 2%, we think
    # existing the adaptor
    for adaptor_name, adaptor_sequence in adaptor_list.items():
        if adaptor_dataframe[adaptor_name].max() > 2:
            adaptor = '-a ' + adaptor_sequence
            adaptors.append(adaptor)
    return adaptors


def cutadaptor(fastq_file_path, cut_fastq_file_path, adaptors):
    """
    Use software cutadapt to remove sequence adaptor
    :param adaptors:A list such as ['-a adaptor1','-a adaptor2']
    :return:NaN
    """
    print('Now remove adaptors from your fastq file...')
    adaptor_num = len(adaptors)
    adaptors_parameter = ' '.join(adaptors)
    # adaptors exist
    # set -n str(adaptor_num + 1) to remove adaptors completely
    if adaptors:
        cutadapt_cmd = ('cutadapt -n ' + str(adaptor_num + 1) +
                        ' -m 16 -M 45 -O 5 -e 0.1 --max-n 0 -j 2 ' +
                        adaptors_parameter + ' -o ' + cut_fastq_file_path + ' ' +
                        fastq_file_path)
        os.system(cutadapt_cmd)
    else:
        print("The fastq file has been processed or owns alternative "
              "adaptors，therefore we can't detect adaptors!")


def fastq_transfer_fasta(fastq_file_path):
    """
    Transfer fastq to fasta
    :return:NaN
    """
    print('Now transfer fastq to fasta...')
    line_num = 0
    sequence_num = 0
    fasta_datas = []
    fasta_path = fastq_file_path[:-5] + 'fa'
    # extract sequence
    with open(fastq_file_path) as cffd:
        for line in cffd:
            line_num += 1
            if line_num % 4 == 1:
                sequence_num += 1
                fasta_line = '>sequence_' + str(sequence_num) + '\n'
                fasta_datas.append(fasta_line)
            elif line_num % 4 == 2:
                fasta_datas.append(line)
    # write sequence to .fa format file
    with open(fasta_path, 'w') as fd:
        fd.writelines(fasta_datas)


def fasta_sequence_count(fasta_path, count_fasta_path):
    '''
    Count the reads number of every sequence from fasta file
    :return:NaN
    '''
    print('Now count the sequence reads number from fasta file')
    line_num = 0
    fasta_dictionary = {}
    sorted_fasta_dictionary = collections.OrderedDict()
    with open(fasta_path) as fad:
        for fasta_line in fad:
            line_num += 1
            if line_num % 2 == 0:
                if fasta_line in fasta_dictionary.keys():
                    fasta_dictionary[fasta_line] += 1
                else:
                    fasta_dictionary[fasta_line] = 1
    # .items() return tuple list
    # key = sequecne_num will reorder tuple list in numerical descending order
    # Finally return a sorted tuple list
    fasta_list = sorted(fasta_dictionary.items(),
                        key=lambda sequence_num: sequence_num[1], reverse=True)
    for sequence_reads in fasta_list:
        sorted_fasta_dictionary[sequence_reads[0]] = sequence_reads[1]
    # counted fasta file format：>s1  sequence_num
    #                           sequence
    with open(count_fasta_path, 'w') as cfdy:
        sequence_rank = 0
        for sequence, sequence_num in sorted_fasta_dictionary.items():
            sequence_rank += 1
            count_fasta_line = '>s%08d_%d\n' % (sequence_rank, sequence_num)
            cfdy.write(count_fasta_line)
            cfdy.write(sequence)


def fastq_sequence_count(fastq_file_path):
    """
    Count the reads number of every sequence from fastq file and transfer to fasta
    :return:NaN
    """
    print('Now count the sequence number from fastq file and transfer to fasta...')
    line_num = 0
    fasta_dictionary = {}
    sorted_fasta_dictionary = collections.OrderedDict()
    with open(fastq_file_path) as fqd:
        for fastq_line in fqd:
            line_num += 1
            if line_num % 4 == 2:
                if fastq_line in fasta_dictionary.keys():
                    fasta_dictionary[fastq_line] += 1
                else:
                    fasta_dictionary[fastq_line] = 1
    fasta_list = sorted(fasta_dictionary.items(),
                        key=lambda sequence_num: sequence_num[1], reverse=True)
    for sequence_reads in fasta_list:
        sorted_fasta_dictionary[sequence_reads[0]] = sequence_reads[1]
    count_fastq_path = fastq_file_path[:-5] + 'fa'
    with open(count_fastq_path, 'w') as cqdy:
        sequence_rank = 0
        for sequence, sequence_num in sorted_fasta_dictionary.items():
            sequence_rank += 1
            count_fastq_line = '>s%08d_%d\n' % (sequence_rank, sequence_num)
            cqdy.write(count_fastq_line)
            cqdy.write(sequence)


def bowtie(reference_index_path, fastq_or_fasta_path):
    """
    Bowtie sequence align
    :return: NaN
    """
    print('Now use bowtie to get sequence map...')
    # -v 1 -a -p 2 result:sam format,fasta represent unique,fastq represent redundance
    sam_filename = fastq_or_fasta_path.split('.')[0] + '.sam'
    if fastq_or_fasta_path[-1] == 'a':
        bowtie_cmd = (
            'bowtie -v 1 -f -a -p 2 --sam-nohead ' +
            reference_index_path +
            ' ' +
            fastq_or_fasta_path +
            ' --sam ' +
            sam_filename)
    elif fastq_or_fasta_path[-1] == 'q':
        bowtie_cmd = (
            'bowtie -v 1 -q -a -p 2 --sam-nohead ' +
            reference_index_path +
            ' ' +
            fastq_or_fasta_path +
            ' --sam ' +
            sam_filename)
    os.system(bowtie_cmd)


def get_and_process_sam_dataframe(rRNA_sam_path):
    """
    Open initial sam file,alter and add columns
    :return: A dataframe of modified sam file
    """
    colname = [
        'QNAME',
        'FLAG',
        'RNAME',
        'POS_START',
        'LENGTH',
        'SEQUENCE',
        'MISMATCH']

    def get_rank(qname):
        rank = int(qname.split('_')[0][1:])
        return rank

    def get_sequence_count(qname):
        sequence_number = int(qname.split('_')[1])
        return sequence_number

    def get_length(length):
        if length == '*':
            length = 1
        else:
            length = int(length[:-1])
        return length

    def get_position(position_start):
        if position_start == '*':
            position_start = 0
        else:
            position_start = int(position_start)
        return position_start

    def add_45s_internal_region_row(modified_sam_dataframe):
        region_of_45s = modified_sam_dataframe[modified_sam_dataframe['RNAME'] == '45S']
        internal_region_of_45s = region_of_45s[
            ((region_of_45s['POS_START'] < 3655) & (
                region_of_45s['POS_END'] > 3655)) | (
                (region_of_45s['POS_START'] < 5523) & (
                    region_of_45s['POS_END'] > 5523)) | (
                    (region_of_45s['POS_START'] < 6601) & (
                        region_of_45s['POS_END'] > 6601)) | (
                            (region_of_45s['POS_START'] < 6757) & (
                                region_of_45s['POS_END'] > 6757)) | (
                                    (region_of_45s['POS_START'] < 7925) & (
                                        region_of_45s['POS_END'] > 7925)) | (
                                            (region_of_45s['POS_START'] < 12990) & (
                                                region_of_45s['POS_END'] > 12990))]
        internal_region_of_45s.loc[:, 'RNAME'] = '45S_across_region'
        modified_sam_dataframe = pd.concat(
            [modified_sam_dataframe, internal_region_of_45s])
        return modified_sam_dataframe
    sam_dataframe = pd.read_csv(rRNA_sam_path, sep='\t', header=None,
                                usecols=[0, 1, 2, 3, 5, 9, 11], names=colname)
    ranks = list(map(get_rank, sam_dataframe['QNAME']))
    sequence_numbers = list(map(get_sequence_count, sam_dataframe['QNAME']))
    lengths = list(map(get_length, sam_dataframe['LENGTH']))
    position_starts = list(map(get_position, sam_dataframe['POS_START']))
    modified_sam_dataframe = sam_dataframe.copy()
    modified_sam_dataframe.drop(
        ['QNAME', 'LENGTH', 'POS_START'], axis=1, inplace=True)
    modified_sam_dataframe['RANK'] = ranks
    modified_sam_dataframe['COUNT'] = sequence_numbers
    modified_sam_dataframe['LENGTH'] = lengths
    modified_sam_dataframe['POS_START'] = position_starts
    modified_sam_dataframe['POS_END'] = (modified_sam_dataframe['POS_START'] +
                                         modified_sam_dataframe['LENGTH'] - 1)
    modified_sam_dataframe.set_index('RANK', inplace=True)
    modified_sam_dataframe = add_45s_internal_region_row(
        modified_sam_dataframe)
    modified_sam_dataframe.sort_index(inplace=True)
    return modified_sam_dataframe


class SamFile():
    """
    Process sam dataframe
    """

    def __init__(self, rRNA_sam_path, species='human'):
        """
        initial_attribution,default species:human
        """
        self.sam_file_name = rRNA_sam_path.split('/')[-1]
        self.species = species
        self.sam = self.get_and_process_sam_dataframe(rRNA_sam_path, species)
        self.total_unique_read_number = max(self.sam.index)
        self.total_redundant_read_number = sum(
            self.sam[(1 - self.sam.index.duplicated('first')) .astype(np.bool)]['COUNT'])
        self.total_unique_length_distribution = pd.DataFrame(
            (self.sam[~self.sam.index .duplicated('last')]['SEQUENCE']).apply(len)).groupby('SEQUENCE').size()
        # When counldn't align to reference genome,length present 0,so I write such suffocating
        # code to get length from sequence and at last get the all sequence
        # length distribution
        self.total_redundant_length_distribution = (pd.DataFrame(
            self.sam[~self.sam.index.duplicated(keep='last')]['SEQUENCE'].apply(len))).merge(
            pd.DataFrame(self.sam[~self.sam.index.duplicated(keep='last')]['COUNT']), left_index=True,
            right_index=True).groupby('SEQUENCE').apply(lambda x: x['COUNT'].sum())

    @staticmethod
    def get_and_process_sam_dataframe(rRNA_sam_path, species):
        """
           default species:human
           Open initial sam file,alter and add columns
           :return: A dataframe of modified sam file
           """
        colname = [
            'QNAME',
            'FLAG',
            'RNAME',
            'POS_START',
            'LENGTH',
            'SEQUENCE',
            'MISMATCH']

        def get_rank(qname):
            rank = int(qname.split('_')[0][1:])
            return rank

        def get_sequence_count(qname):
            sequence_number = int(qname.split('_')[1])
            return sequence_number

        def get_length(length):
            if length == '*':
                length = 1
            else:
                length = int(length[:-1])
            return length

        def get_position(position_start):
            if position_start == '*':
                position_start = 0
            else:
                position_start = int(position_start)
            return position_start

        def add_45s_internal_region_row(modified_sam_dataframe):
            region_of_45s = modified_sam_dataframe[modified_sam_dataframe['RNAME'] == '45S']
            if species == 'human':
                internal_region_of_45s = (
                    region_of_45s[
                        ((region_of_45s['POS_START'] < 3655) & (
                            region_of_45s['POS_END'] >= 3655)) | (
                            (region_of_45s['POS_START'] <= 5523) & (
                                region_of_45s['POS_END'] > 5523)) | (
                            (region_of_45s['POS_START'] < 6601) & (
                                region_of_45s['POS_END'] >= 6601)) | (
                                (region_of_45s['POS_START'] <= 6757) & (
                                    region_of_45s['POS_END'] > 6757)) | (
                                        (region_of_45s['POS_START'] < 7925) & (
                                            region_of_45s['POS_END'] >= 7925)) | (
                                                (region_of_45s['POS_START'] <= 12990) & (
                                                    region_of_45s['POS_END'] > 12990))]).copy()
            elif species == 'mouse':
                internal_region_of_45s = (
                    region_of_45s[
                        ((region_of_45s['POS_START'] <= 4007) & (
                            region_of_45s['POS_END'] > 4007)) | (
                            (region_of_45s['POS_START'] < 5878) & (
                                region_of_45s['POS_END'] >= 5878)) | (
                            (region_of_45s['POS_START'] <= 6877) & (
                                region_of_45s['POS_END'] > 6877)) | (
                                (region_of_45s['POS_START'] < 7035) & (
                                    region_of_45s['POS_END'] >= 7035)) | (
                                        (region_of_45s['POS_START'] <= 8122) & (
                                            region_of_45s['POS_END'] > 8122)) | (
                                                (region_of_45s['POS_START'] < 12850) & (
                                                    region_of_45s['POS_END'] >= 12850))]).copy()
            internal_region_of_45s['RNAME'] = '45S_across_region'
            modified_sam_dataframe = pd.concat(
                [modified_sam_dataframe, internal_region_of_45s])
            return modified_sam_dataframe
        sam_dataframe = pd.read_csv(
            rRNA_sam_path, sep='\t', header=None, usecols=[
                0, 1, 2, 3, 5, 9, 11], names=colname)
        ranks = list(map(get_rank, sam_dataframe['QNAME']))
        sequence_numbers = list(
            map(get_sequence_count, sam_dataframe['QNAME']))
        lengths = list(map(get_length, sam_dataframe['LENGTH']))
        position_starts = list(map(get_position, sam_dataframe['POS_START']))
        modified_sam_dataframe = sam_dataframe.copy()
        modified_sam_dataframe.drop(
            ['QNAME', 'LENGTH', 'POS_START'], axis=1, inplace=True)
        modified_sam_dataframe['FLAG'] = modified_sam_dataframe['FLAG'].apply(
            str)
        modified_sam_dataframe['RANK'] = ranks
        modified_sam_dataframe['COUNT'] = sequence_numbers
        modified_sam_dataframe['LENGTH'] = lengths
        modified_sam_dataframe['POS_START'] = position_starts
        modified_sam_dataframe['POS_END'] = (
            modified_sam_dataframe['POS_START'] +
            modified_sam_dataframe['LENGTH'] -
            1)
        modified_sam_dataframe.set_index('RANK', inplace=True)
        modified_sam_dataframe = add_45s_internal_region_row(
            modified_sam_dataframe)
        modified_sam_dataframe.sort_index(inplace=True)
        return modified_sam_dataframe

    @staticmethod
    def get_unique_sequence_dataframe(select_dataframe):
        """
        replication code
        :param select_dataframe:redundant dataframe(mean a sequence align multiple genes or sites)
        :return:unique dataframe(mean a sequence align just one gene one site)
        """
        bool_query_of_unique = (
            1 -
            select_dataframe.index.duplicated('last')).astype(
            np.bool)
        unique_sequence_dataframe = select_dataframe[bool_query_of_unique]
        return unique_sequence_dataframe

    @staticmethod
    def length_distribution(length_series):
        """
        transfer and format the length distribution series to a sorted ordered dictionary,
        when the needed key isn't existed,assigned its value 0
        when key > 45 ，calculate the sum of the values >45
        :param length_series:length_distribution_series or array(with a index)
        :return:formatted ordered dictionary
        """
        length_series = dict(length_series)
        for length in range(16, 46):
            if length not in length_series.keys():
                length_series[length] = 0
        sum_more_45 = sum([x for k, x in length_series.items() if k > 45])
        length_series = {k: v for k, v in length_series.items() if k < 46}
        length_series = collections.OrderedDict(sorted(length_series.items()))
        length_series['>45'] = sum_more_45
        return length_series

    @staticmethod
    def plot_stack_bar(
            ax,
            referencename,
            bar_location=0.,
            width=0.9,
            legend=True,
            ylabel='Read counts'):
        """
        plot length distribution stack bar
        :param ax: subplot name
        :param referencename: data dictionary，key：reference gene value：length distribution
        :param bar_location: defualt 0 single column status
        :param width: columns width
        :return: stack bar chart
        """
        ind = np.arange(16, 47)
        cmap = plt.get_cmap('Set3')
        colornumber = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1]
        i = 0
        for name, length_distribution in referencename.items():
            if i == 0:
                ax.bar(
                    ind + bar_location,
                    length_distribution,
                    width=width,
                    label=name,
                    color=cmap(
                        colornumber[i]))
                i += 1
                k = length_distribution.copy()
            else:
                ax.bar(
                    ind + bar_location,
                    length_distribution,
                    bottom=k,
                    width=width,
                    label=name,
                    color=cmap(
                        colornumber[i]))
                i += 1
                k += length_distribution
        ax.set_xlabel('Length', fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.set_xticks([15, 20, 25, 30, 35, 40, 45, 46])
        ax.set_xticklabels(['15', '20', '25', '30', '35', '40', '45', '>45'])
        if legend:
            ax.legend(
                bbox_to_anchor=(
                    0., 1.02, 1., .102), loc=3, fontsize=10, ncol=len(
                    referencename.keys()), mode="expand", borderaxespad=0.)

    def recursion_savefig(self):
        """
        replication code:save plot
        """
        try:
            save_directory = input(
                'type the directory and figure name:(eg:/Users/lid/Desktop/foo)')
            plt.savefig(
                save_directory,
                dpi=200,
                bbox_inches='tight',
                format='pdf')
        except BaseException:
            print('please type correct directory and figure name!')
            self.recursion_savefig()

    def rpm_normalization(self, fastq_or_sam, ref_path=''):
        """
        use rpm to normalization for further multiple samples analysis such as plot heatmap
        two method :
        one deliver a sam file path
        second deliver a fastq or fasta file path and genome index path
        :param genome_sam_file:to get map genome reads number,
        input a genome sam path or clean fastq path or clean fasta path
        :param ref_path :reference genome index path if fastq or fasta given
        :return:None but add rpm value as colname 'RPM'
        you had better not use rpm normalization data for SamFile class function analysis because the y label
        is read counts but not rpm and for single sample,read counts is enough,if use rpm,when plot,add parameter
        ylabel='RPM' to amend ylabel and a fatal question is only redundant is revised to rpm value and though
        unique ylabel may show rpm label，but real value is Read counts,so this function aim to extract rpm value
        of every sample
        """

        def get_sequence_count(qname):
            """
            function used map
            """
            sequence_number = int(qname.split('_')[1])
            return sequence_number
        if 'sam' in fastq_or_sam:
            colname = [
                'QNAME',
                'FLAG',
                'RNAME',
                'POS_START',
                'LENGTH',
                'SEQUENCE',
                'MISMATCH']
            genome_sam = pd.read_csv(
                fastq_or_sam, sep='\t', header=None, usecols=[
                    0, 1, 2, 3, 5, 9, 11], names=colname)
            sequence_numbers = list(
                map(get_sequence_count, genome_sam['QNAME']))
            genome_sam['COUNT'] = sequence_numbers
            self.map_genome_number = sum(
                genome_sam[genome_sam['FLAG'] == 0 | genome_sam['FLAG'] == 16]['COUNT'])
        elif 'fastq' or 'fq' in fastq_or_sam:
            bowtie = 'bowtie -q -v 1 -p 4 ' + ref_path + \
                ' ' + fastq_or_sam + ' --sam sam.tmp 2>&1'
            map_info = os.popen(bowtie).readlines()
            os.system('rm -rf sam.tmp')
            self.map_genome_number = int(
                map_info[1].split(':')[1].split(' ')[1])
        elif 'fasta' or 'fa' in fastq_or_sam:
            bowtie = 'bowtie -f -v 1 -p 4 ' + ref_path + \
                ' ' + fastq_or_sam + ' --sam sam.tmp 2>&1'
            map_info = os.popen(bowtie).readlines()
            os.system('rm -rf sam.tmp')
            self.map_genome_number = int(
                map_info[1].split(':')[1].split(' ')[1])
        if self.map_genome_number == 0:
            print("Your sequencing result can't match the genome")
        else:
            self.sam['RPM'] = self.sam['COUNT'] / \
                self.map_genome_number * 1000000

    def count_different_map(
            self,
            flag='',
            rname='',
            threshold=0,
            sequence_dataframe=''):
        """
        count different source alignment number
        :param flag: positive or negative chain
        :param rname:5S 5.8S 18S 28S 12S 16S 45S ITS1 ITS2 5ETS 3ETS
        :return: (unique_sequence_counts,redundant_sequence_counts)
        """
        def deduplicate_and_count(select_dataframe):
            """
            deduplicate the multialignment of single sequence to avoid counting duplicately
            :param select_dataframe
            :return: (unique_sequence_counts,redundant_sequence_counts)
            """
            unique_sequence_dataframe = self.get_unique_sequence_dataframe(
                select_dataframe)
            unique_sequence_dataframe = unique_sequence_dataframe[
                unique_sequence_dataframe['COUNT'] > threshold]
            unique_sequence_numbers = len(unique_sequence_dataframe.index)
            redundant_sequence_numbers = sum(
                unique_sequence_dataframe['COUNT'])
            return unique_sequence_numbers, redundant_sequence_numbers
        if sequence_dataframe:
            sequence_numbers = deduplicate_and_count(sequence_dataframe)
        else:
            if flag and rname:
                select_dataframe = self.sam[(self.sam['FLAG'] == flag) & (
                    self.sam['RNAME'] == rname)]
                sequence_numbers = deduplicate_and_count(select_dataframe)
            elif rname and not flag:
                select_dataframe = self.sam[self.sam['RNAME'] == rname]
                sequence_numbers = deduplicate_and_count(select_dataframe)
            elif flag and not rname:
                select_dataframe = self.sam[self.sam['FLAG'] == flag]
                sequence_numbers = deduplicate_and_count(select_dataframe)
            else:
                select_dataframe = self.sam[self.sam['FLAG'] != '4']
                sequence_numbers = deduplicate_and_count(select_dataframe)
        return sequence_numbers

    def different_map_dataframe(
            self,
            flag='',
            rname='',
            threshold=0,
            sequence_dataframe=''):
        """
        Get designed dataframe
        :param flag：positive or negative chain
        :param rname:5S 5.8S 18S 28S 12S 16S 45S ITS1 ITS2 5ETS 3ETS
        :return:design dataframe
        """
        def get_chunk_dataframe(initial_dataframe):
            """
            Get designed dataframe
            """
            if flag and rname:
                map_chunk_dataframe = (initial_dataframe[(
                    self.sam['FLAG'] == flag) & (initial_dataframe['RNAME'] == rname)])
                map_chunk_dataframe = map_chunk_dataframe[map_chunk_dataframe['COUNT'] > threshold]
            elif flag and not rname:
                map_chunk_dataframe = initial_dataframe[initial_dataframe['FLAG'] == flag]
                map_chunk_dataframe = map_chunk_dataframe[map_chunk_dataframe['COUNT'] > threshold]
            elif rname and not flag:
                map_chunk_dataframe = initial_dataframe[initial_dataframe['RNAME'] == rname]
                map_chunk_dataframe = map_chunk_dataframe[map_chunk_dataframe['COUNT'] > threshold]
            else:
                map_chunk_dataframe = initial_dataframe[initial_dataframe['FLAG'] != '4']
                map_chunk_dataframe = map_chunk_dataframe[map_chunk_dataframe['COUNT'] > threshold]
            return map_chunk_dataframe
        if sequence_dataframe:
            map_chunk_dataframe = get_chunk_dataframe(sequence_dataframe)
        else:
            map_chunk_dataframe = get_chunk_dataframe(self.sam)
        return map_chunk_dataframe

    def get_most_abundant_hits(
            self,
            number,
            flag='',
            rname='',
            threshold=0,
            sequence_dataframe=''):
        """
        top number hits
        :param number: hits number
        :param flag:positive or negative chain
        :param rname:5S 5.8S 18S 28S 12S 16S 45S ITS1 ITS2 5ETS 3ETS
        :return:top hits dataframe
        """
        def judge_hit_numbers(select_dataframe):
            """
            judge the designated hit numbers whether more than current dataframe，
            if more ，return current abundant hits dataframe
            """
            select_dataframe = select_dataframe[select_dataframe['COUNT'] > threshold]
            if len(select_dataframe.index) >= number:
                abundant_hits = select_dataframe[0:number][[
                    'SEQUENCE', 'COUNT', 'LENGTH', 'POS_START', 'POS_END']]
            else:
                abundant_hits = select_dataframe[0:len(select_dataframe.index)][
                    ['SEQUENCE', 'COUNT', 'LENGTH', 'POS_START', 'POS_END']]
            return abundant_hits
        if sequence_dataframe:
            abundant_hits = judge_hit_numbers(sequence_dataframe)
        else:
            if flag and rname:
                select_dataframe = self.sam[(self.sam['FLAG'] == flag) & (
                    self.sam['RNAME'] == rname)]
                abundant_hits = judge_hit_numbers(select_dataframe)
            elif rname and not flag:
                select_dataframe = self.sam[self.sam['RNAME'] == rname]
                abundant_hits = judge_hit_numbers(select_dataframe)
            elif flag and not rname:
                select_dataframe = self.sam[self.sam['FLAG'] == flag]
                abundant_hits = judge_hit_numbers(select_dataframe)
            else:
                select_dataframe = self.sam[self.sam['FLAG'] != '4']
                abundant_hits = judge_hit_numbers(select_dataframe)
        return abundant_hits

    def get_length_distribution(self, flag='', rname='', type='', threshold=0,
                                sequence_dataframe='', count_format='COUNT'):
        """
        get length distribution series
        :param flag: 0：positive 16：negative
        :param rname: reference gene
        :param type: unique or redundant or else(both unique and redundant)
        :param threshold: filter sequence whose read counts are more than threshold
        :param sequence_dataframe: when given dataframe，you don't need flag and rname
        :return:length distribution series
        """
        def deduplicate_and_length_distribution(select_dataframe):
            """
            deduplicate and return unique or redundant or both length distribution series
            """
            unique_sequence_dataframe = select_dataframe[~select_dataframe.index.duplicated(
                keep='last')]
            unique_sequence_dataframe = unique_sequence_dataframe[
                unique_sequence_dataframe[count_format] > threshold]
            if type == 'unique':
                unique_length_distribution = unique_sequence_dataframe.LENGTH.value_counts()
                return unique_length_distribution
            elif type == 'redundant':
                redundant_length_distribution = (
                    unique_sequence_dataframe.groupby(
                        by='LENGTH') .apply(
                        lambda x: x[count_format].sum()))
                return redundant_length_distribution
            else:
                unique_length_distribution = unique_sequence_dataframe.LENGTH.value_counts()
                redundant_length_distribution = (
                    unique_sequence_dataframe.groupby(
                        by='LENGTH') .apply(
                        lambda x: x[count_format].sum()))
                return unique_length_distribution, redundant_length_distribution
        if sequence_dataframe:
            length_distributions = deduplicate_and_length_distribution(
                sequence_dataframe)
        else:
            if flag and rname:
                select_dataframe = self.sam[(self.sam['FLAG'] == flag) & (
                    self.sam['RNAME'] == rname)]
                length_distributions = deduplicate_and_length_distribution(
                    select_dataframe)
            elif flag and not rname:
                select_dataframe = self.sam[self.sam['FLAG'] == flag]
                length_distributions = deduplicate_and_length_distribution(
                    select_dataframe)
            elif rname and not flag:
                select_dataframe = self.sam[self.sam['RNAME'] == rname]
                length_distributions = deduplicate_and_length_distribution(
                    select_dataframe)
            else:
                select_dataframe = self.sam[self.sam['FLAG'] != '4']
                length_distributions = deduplicate_and_length_distribution(
                    select_dataframe)
        return length_distributions

    def get_location_distribution_list(
            self,
            flag='',
            rname='',
            threshold=0,
            sequence_dataframe=''):
        """
        get location distribution list
        :param flag: 0：positive 16：negative
        :param rname: reference gene
        :param threshold: filter sequence whose read counts are more than threshold
        :param sequence_dataframe: when given dataframe，you don't need flag and rname
        :return: a large list own duplicate digital which is base site and can count read base number，
        positive chain，digital is positive value，negative chain，digital is negative value
        througn counting different digital number，we can get location distribution dictionary and histgram
        return unique_location_distribution_list,redundant_location_distribution_list
        """
        def get_location_list(select_dataframe):
            """
            get total base site list
            :return unique and redundant base location list
            """
            def get_sequence_location_list(
                    location_start,
                    location_end,
                    location_flag,
                    location_counts=1):
                """
                for map parameter funcation to get total base site list
                """
                if location_flag == '0':
                    single_sequence_location_list = (
                        [i for i in range(location_start, location_end + 1)] * location_counts)
                elif location_flag == '16':
                    single_sequence_location_list = (
                        [i for i in range(location_start, location_end + 1)] * location_counts)
                    single_sequence_location_list = list(
                        map(lambda i: i * (-1), single_sequence_location_list))
                return single_sequence_location_list
            select_dataframe = select_dataframe[select_dataframe['COUNT'] > threshold]
            unique_location_distribution = list(
                map(
                    get_sequence_location_list,
                    select_dataframe['POS_START'],
                    select_dataframe['POS_END'],
                    select_dataframe['FLAG']))
            unique_location_distribution = list(
                chain.from_iterable(unique_location_distribution))
            redundant_location_distribution = list(
                map(
                    get_sequence_location_list,
                    select_dataframe['POS_START'],
                    select_dataframe['POS_END'],
                    select_dataframe['FLAG'],
                    select_dataframe['COUNT']))
            redundant_location_distribution = list(
                chain.from_iterable(redundant_location_distribution))
            return unique_location_distribution, redundant_location_distribution
        if sequence_dataframe:
            location_distributions = get_location_list(sequence_dataframe)
        else:
            if flag and rname:
                select_dataframe = self.sam[(self.sam['FLAG'] == flag) & (
                    self.sam['RNAME'] == rname)]
                location_distributions = get_location_list(select_dataframe)
            elif rname and not flag:
                select_dataframe = self.sam[self.sam['RNAME'] == rname]
                location_distributions = get_location_list(select_dataframe)
        return location_distributions

    def get_location_distribution_dic(
            self,
            rname,
            unique_location_distribution_list,
            redundant_location_distribution_list):
        """
        get location distribution dictionary
        :param rname: for assigning value 0 to specified reference gene base site which don't be mapped
        :param unique_location_distribution_list: get from function get_location_distribution_list,
        we use function map and abs to transfer negative value in list to positive value
        :param redundant_location_distribution_list: get from function get_location_distribution_list,
        we use function map and abs to transfer negative value in list to positive value
        :return: unique_location_distribution_dic,redundant_location_distribution_dic
        """
        if self.species == 'human':
            rname_length = {
                '45S': 13351,
                '28S': 5066,
                '18S': 1869,
                '5.8S': 157,
                '5S': 121,
                '12S': 954,
                '16S': 1559,
                'ITS1': 1077,
                'ITS2': 1167,
                '5ETS': 3654,
                '3ETS': 361,
                '45S_across_region': 13351}
        elif self.species == 'mouse':
            rname_length = {
                '45S': 13400,
                '28S': 4730,
                '18S': 1870,
                '5.8S': 157,
                '5S': 121,
                '4.5S': 174,
                '16S': 110,
                'ITS1': 1000,
                'ITS2': 1088,
                '5ETS': 4007,
                '3ETS': 551,
                '45S_across_region': 13400}

        def location_distribution_dictionary(location_distribution_list):
            """
            count base location number, assign value 0 to base site which don't be mapped,
            and sort the base location distribution dictionary
            :param location_distribution_list:get from function get_location_distribution_list
            :return:base location distribution dic
            """
            location_distribution_dic = dict(
                collections.Counter(location_distribution_list))
            for i in range(1, rname_length[rname] + 1):
                if i not in location_distribution_dic.keys():
                    location_distribution_dic[i] = 0
            location_distribution_dic = collections.OrderedDict(
                sorted(location_distribution_dic.items()))
            return location_distribution_dic
        unique_location_distribution_list = list(
            map(abs, unique_location_distribution_list))
        redundant_location_distribution_list = list(
            map(abs, redundant_location_distribution_list))
        unique_location_distribution_dic = location_distribution_dictionary(
            unique_location_distribution_list)
        redundant_location_distribution_dic = location_distribution_dictionary(
            redundant_location_distribution_list)
        return unique_location_distribution_dic, redundant_location_distribution_dic

    def quick_location_dic(
            self,
            flag='',
            rname='',
            type='',
            threshold=0,
            sequence_dataframe='',
            count_format='COUNT'):
        """
        quick get location distribution dictionary
        dictionary type:location site is positive value but may represent positive or negative chain
        positive:support production from self.sam
        negative:support production from self.sam
        positive+negative:same site,sum(value from positive and negative chain),
        support production from self.sam and sequence dataframe
        :param flag: positive or negative
        :param rname: reference gene
        :param threshold: rpm or Read counts threshold
        :param sequence_dataframe: you can add a eligible dataframe to get positive+negative type location
        distribution
        :param count_format:'RPM' or 'COUNT'，when use RPM，you need use function rpm_normalization firstly
        :return:unique and redundant location distribution dictionary
        """
        if self.species == 'human':
            rname_length = {
                '45S': 13351,
                '28S': 5066,
                '18S': 1869,
                '5.8S': 157,
                '5S': 121,
                '12S': 954,
                '16S': 1559,
                'ITS1': 1077,
                'ITS2': 1167,
                '5ETS': 3654,
                '3ETS': 361,
                '45S_across_region': 13351}
        elif self.species == 'mouse':
            rname_length = {
                '45S': 13400,
                '28S': 4730,
                '18S': 1870,
                '5.8S': 157,
                '5S': 121,
                '4.5S': 174,
                '16S': 110,
                'ITS1': 1000,
                'ITS2': 1088,
                '5ETS': 4007,
                '3ETS': 551,
                '45S_across_region': 13400}

        def get_location_dic(select_dataframe):
            """
            get location distribution dictionary
            """
            select_dataframe = select_dataframe[select_dataframe[count_format] > threshold]
            location_dic_key = np.array(
                [i for i in range(1, rname_length[rname] + 1)])
            if type == 'unique':
                unique_location_dic_value = np.array(
                    [0.] * rname_length[rname])
                for i, p in zip(
                        select_dataframe['POS_START'], select_dataframe['POS_END']):
                    unique_location_dic_value[i - 1:p] += 1
                unique_location_dic = collections.OrderedDict(
                    zip(location_dic_key, unique_location_dic_value))
                return unique_location_dic
            elif type == 'redundant':
                redundant_location_dic_value = np.array(
                    [0.] * rname_length[rname])
                for i, p, j in zip(select_dataframe['POS_START'], select_dataframe['POS_END'],
                                   select_dataframe[count_format]):
                    redundant_location_dic_value[i - 1:p] += j
                redundant_location_dic = collections.OrderedDict(
                    zip(location_dic_key, redundant_location_dic_value))
                return redundant_location_dic
            else:
                unique_location_dic_value = np.array(
                    [0.] * rname_length[rname])
                redundant_location_dic_value = unique_location_dic_value.copy()
                for i, p, j in zip(select_dataframe['POS_START'], select_dataframe['POS_END'],
                                   select_dataframe[count_format]):
                    unique_location_dic_value[i - 1:p] += 1
                    redundant_location_dic_value[i - 1:p] += j
                unique_location_dic = collections.OrderedDict(
                    zip(location_dic_key, unique_location_dic_value))
                redundant_location_dic = collections.OrderedDict(
                    zip(location_dic_key, redundant_location_dic_value))
                return unique_location_dic, redundant_location_dic
        if sequence_dataframe:
            location_distributions = get_location_dic(sequence_dataframe)
        else:
            if flag and rname:
                select_dataframe = self.sam[(self.sam['FLAG'] == flag) & (
                    self.sam['RNAME'] == rname)]
                location_distributions = get_location_dic(select_dataframe)
            elif rname and not flag:
                select_dataframe = self.sam[self.sam['RNAME'] == rname]
                location_distributions = get_location_dic(select_dataframe)
        return location_distributions

    def plot_location_distribution_area(
            self,
            rname='',
            chain='',
            unique_location_distribution_dic='',
            redundant_location_distribution_dic='',
            save_directory='',
            n_unique_location_distribution_dic='',
            n_redundant_location_distribution_dic='',
            p_unique_location_distribution_dic='',
            p_redundant_location_distribution_dic=''):
        """
        plot location distribution area figure
        two type :first is positive + negative unique and redudant in the same axes
                  the other is positive:unique and redundant negative:unique and redundant in different axes
        :param rname: reference gene
        :param chain: 'all':first type else:the other type
        :param unique_location_distribution_dic: first type unique_location_distribution_dic
        :param redundant_location_distribution_dic: first type redundant_location_distribution_dic
        :param save_directory: storing path
        :param n_unique_location_distribution_dic:second type negative unique
        :param n_redundant_location_distribution_dic: second type negative redundant
        :param p_unique_location_distribution_dic: second type positive unique
        :param p_redundant_location_distribution_dic: second type positive redundant
        :return: location distribution line area chart
        """
        def transfer_numpy(location_distribution_dic):
            """
            repication code:transfer dic to numpy ndarray
            :param location_distribution_dic
            :return: 2d array 1:location 2:read counts
            """
            location_distribution_dic = location_distribution_dic.items()
            location_distribution_dic = list(zip(*location_distribution_dic))
            location_distribution_dic = np.array(location_distribution_dic)
            return location_distribution_dic

        def area_plot(
                ax_1,
                ax_2,
                label1,
                label2,
                unique_location_distribution_array,
                redundant_location_distribution_array,
                ylabel='Read counts'):
            """
            plot location distribution area figure
            :param ax_1: same axes subplot1
            :param ax_2: same axes subplot2
            :param label1:legend label1
            :param label2: legend label2
            :param unique_location_distribution_array
            :param redundant_location_distribution_array
            :return: area figure
            """
            ax_1.plot(
                unique_location_distribution_array[0],
                unique_location_distribution_array[1],
                color='purple',
                label=label1)
            ax_1.fill_between(
                unique_location_distribution_array[0],
                unique_location_distribution_array[1],
                where=unique_location_distribution_array[1] > 0,
                facecolor='purple',
                alpha=0.6)
            ax_2.plot(
                redundant_location_distribution_array[0],
                redundant_location_distribution_array[1],
                color='blue',
                label=label2)
            ax_2.fill_between(
                redundant_location_distribution_array[0],
                redundant_location_distribution_array[1],
                where=redundant_location_distribution_array[1] > 0,
                facecolor='blue',
                alpha=0.2)
            ax_1.legend(loc=2, fontsize=13)
            ax_2.legend(loc=1, fontsize=13)
            ax_1.set_xlabel('Location', fontsize=13)
            ax_1.set_ylabel('Unique Read counts', fontsize=13)
            ax_2.set_ylabel('Redundant ' + ylabel, fontsize=13)
        if chain == 'all':
            unique_location_distribution_array = transfer_numpy(
                unique_location_distribution_dic)
            redundant_location_distribution_array = transfer_numpy(
                redundant_location_distribution_dic)
            plt.figure(figsize=(20, 6))
            ax_1 = plt.subplot(111)
            ax_2 = ax_1.twinx()
            area_plot(
                ax_1,
                ax_2,
                'unique',
                'redundant',
                unique_location_distribution_array,
                redundant_location_distribution_array)
            ax_1.set_title(
                rname +
                'rsRNA Location Distribution:Unique vs Redundant',
                fontsize=16)
        else:
            p_unique_location_distribution_array = transfer_numpy(
                p_unique_location_distribution_dic)
            p_redundant_location_distribution_array = transfer_numpy(
                p_redundant_location_distribution_dic)
            n_unique_location_distribution_array = transfer_numpy(
                n_unique_location_distribution_dic)
            n_redundant_location_distribution_array = transfer_numpy(
                n_redundant_location_distribution_dic)
            plt.figure(figsize=(20, 12))
            ax_1 = plt.subplot(211)
            ax_2 = ax_1.twinx()
            area_plot(
                ax_1,
                ax_2,
                'positive unique',
                'positive redundant',
                p_unique_location_distribution_array,
                p_redundant_location_distribution_array)
            ax_1.set_title(
                rname +
                'rsRNA Positive Location Distribution:Unique vs Redundant',
                fontsize=16)
            ax_3 = plt.subplot(212)
            ax_4 = ax_3.twinx()
            area_plot(
                ax_3,
                ax_4,
                'negative unique',
                'negative redundant',
                n_unique_location_distribution_array,
                n_redundant_location_distribution_array)
            ax_3.set_title(
                rname +
                'rsRNA Negative Location Distribution:Unique vs Redundant',
                fontsize=16)
        if save_directory:
            plt.savefig(save_directory + '.tiff', dpi=200)
        else:
            ask_save_fig = input('Do you want to save the figure?(y/n)')
            if ask_save_fig == 'y':
                self.recursion_savefig()

    def plot_location_distribution_hist(
            self,
            rname,
            unique_location_distribution_list,
            redundant_location_distribution_list,
            save_directory=''):
        """
        according to function get_location_distribution_list return value
        to draw location distribution hist
        three type figure：first positive second negative third positive + negative
        rname is a certainly needed value in function get_location_distribution_list
        and draw location distribution hist
        :param rname:reference gene--for determining bins number
        :param unique_location_distribution_list from function get_location_distribution_list
        :param redundant_location_distribution_list from function get_location_distribution_list
        :param save_directory:figure storing path
        :return: location_distribution_hist
        """
        def plot_hist(ax, location_distribution_list, ylabel='Read counts'):
            """
            plot location distribution hist
            :param ax:subplot
            :param location_distribution_list
            :return: location distribution hist
            """
            ax.set_xlabel('Location', fontsize=13)
            ax.set_ylabel(ylabel, fontsize=13)
            ax.set_xlim(1, rname_length[rname])
            N, bins, patches = ax.hist(location_distribution_list,
                                       bins=max(location_distribution_list))
            fracs = N / \
                (N.max() + 0.00000000000000000000000000000000000000000000000000000000000000001)
            norm = colors.Normalize(fracs.min(), fracs.max())
            for thisfrac, thispatch in zip(fracs, patches):
                color = plt.cm.viridis(norm(thisfrac))
                thispatch.set_facecolor(color)
        if self.species == 'human':
            rname_length = {
                '45S': 13351,
                '28S': 5066,
                '18S': 1869,
                '5.8S': 157,
                '5S': 121,
                '12S': 954,
                '16S': 1559,
                'ITS1': 1077,
                'ITS2': 1167,
                '5ETS': 3654,
                '3ETS': 361,
                '45S_across_region': 13351}
        elif self.species == 'mouse':
            rname_length = {
                '45S': 13400,
                '28S': 4730,
                '18S': 1870,
                '5.8S': 157,
                '5S': 121,
                '4.5S': 174,
                '16S': 110,
                'ITS1': 1000,
                'ITS2': 1088,
                '5ETS': 4007,
                '3ETS': 551,
                '45S_across_region': 13400}
        if min(unique_location_distribution_list) < 0 and max(
                unique_location_distribution_list) > 0:
            unique_location_distribution_list = list(
                map(abs, unique_location_distribution_list))
            redundant_location_distribution_list = list(
                map(abs, redundant_location_distribution_list))
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(211)
            ax1.set_title(
                rname +
                'rRNA location distribution-unique sequence',
                fontsize=15)
            plot_hist(ax1, unique_location_distribution_list)
            ax2 = fig.add_subplot(212)
            ax2.set_title(
                rname +
                'rRNA location distribution-redundant sequence',
                fontsize=15)
            plot_hist(ax2, redundant_location_distribution_list)
            plt.tight_layout()
        elif min(unique_location_distribution_list) > 0:
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(211)
            ax1.set_title(
                rname +
                'rRNA location distribution-unique positive sequence',
                fontsize=15)
            plot_hist(ax1, unique_location_distribution_list)
            ax2 = fig.add_subplot(212)
            ax2.set_title(
                rname +
                'rRNA location distribution-redundant positive sequence',
                fontsize=15)
            plot_hist(ax2, redundant_location_distribution_list)
            plt.tight_layout()
        elif max(unique_location_distribution_list) < 0:
            unique_location_distribution_list = np.array(
                unique_location_distribution_list) * (-1)
            redundant_location_distribution_list = np.array(
                redundant_location_distribution_list) * (-1)
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(211)
            ax1.set_title(
                rname +
                'rRNA location distribution-unique negative sequence',
                fontsize=15)
            plot_hist(ax1, unique_location_distribution_list)
            ax2 = fig.add_subplot(212)
            ax2.set_title(
                rname +
                'rRNA location distribution-redundant negative sequence',
                fontsize=15)
            plot_hist(ax2, redundant_location_distribution_list)
            plt.tight_layout()
        if save_directory:
            plt.savefig(save_directory + '.tiff', dpi=200)
        else:
            ask_save_fig = input('Do you want to save the figure?(y/n)')
            if ask_save_fig == 'y':
                self.recursion_savefig()

    def plot_pie_number_distribution(
            self,
            chain='',
            threshold=0,
            save_directory=False):
        """
        plot pie figure about sequence number distribution
        :param chain:'All':distinct positive and negative pie,
                     'positive':distinct positive pie,
                     'negative':distinct negative pie,
                     '':all distinct positive+negative pie(default)
        :param threshold: filter sequences
        :param save_directory: eg:/Users/lid/Desktop/counts_picture
        :return:reads number distribution pie figure
        """
        species = self.species

        def plot_pie(ax, unique_data, redundant_data, label):
            """
            plot pie figure
            :param ax: subplot name
            :param unique_data: unique sequence reads number list
            :param redundant_data: redundant sequence reads number list
            :param label: label list
            :return: pie figure
            """
            def percentage_distribution(data):
                """
                make percentage of every element reads number for legend
                :param data:list of raw data
                :return: a array which every element reads number takes percentage of sum(data)
                """
                data = np.array(data)
                data = data / \
                    (sum(data) + 0.000000000000000000000000000000000000000000000000000000000001)
                return data
            size = 0.3
            unique_datas = percentage_distribution(unique_data)
            redundant_datas = percentage_distribution(redundant_data)
            ax.pie(
                unique_data,
                labels=label,
                radius=1,
                colors=colors,
                wedgeprops=dict(
                    width=size,
                    edgecolor='w'),
                textprops={
                    'size': 'large'})
            wedges, texts = ax.pie(redundant_data, radius=1 - size, colors=colors,
                                   wedgeprops=dict(width=size, edgecolor='w'))
            if chain == 'all':
                if species == 'human':
                    percentage_legend = ['28SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[0], redundant_datas[0]),
                                         '18SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[1], redundant_datas[1]),
                                         '5.8SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[2], redundant_datas[2]),
                                         '45S_rRNA_otherp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[3], redundant_datas[3]),
                                         '12SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[4], redundant_datas[4]),
                                         '16SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[5], redundant_datas[5]),
                                         "5'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[6], redundant_datas[6]),
                                         "3'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[7], redundant_datas[7]),
                                         'ITS1p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[8], redundant_datas[8]),
                                         'ITS2p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[9], redundant_datas[9]),
                                         '5SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[10], redundant_datas[10]),
                                         '28SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[11], redundant_datas[11]),
                                         '18SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[12], redundant_datas[12]),
                                         '5.8SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[13], redundant_datas[13]),
                                         '45S_rRNA_othern\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[14], redundant_datas[14]),
                                         '12SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[15], redundant_datas[15]),
                                         '16SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[16], redundant_datas[16]),
                                         "5'ETSn\nO:{:.1%} ,I:{:.1%}".format(unique_datas[17], redundant_datas[17]),
                                         "3'ETSn\nO:{:.1%} ,I:{:.1%}".format(unique_datas[18], redundant_datas[18]),
                                         'ITS1n\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[19], redundant_datas[19]),
                                         'ITS2n\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[20], redundant_datas[20]),
                                         '5SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[21], redundant_datas[21])]
                elif species == 'mouse':
                    percentage_legend = ['28SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[0], redundant_datas[0]),
                                         '18SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[1], redundant_datas[1]),
                                         '5.8SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[2], redundant_datas[2]),
                                         '45S_rRNA_otherp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[3],
                                                                                      redundant_datas[3]),
                                         '4.5SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[4], redundant_datas[4]),
                                         '16SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[5], redundant_datas[5]),
                                         "5'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[6], redundant_datas[6]),
                                         "3'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[7], redundant_datas[7]),
                                         'ITS1p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[8], redundant_datas[8]),
                                         'ITS2p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[9], redundant_datas[9]),
                                         '5SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[10], redundant_datas[10]),
                                         '28SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[11], redundant_datas[11]),
                                         '18SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[12], redundant_datas[12]),
                                         '5.8SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[13], redundant_datas[13]),
                                         '45S_rRNA_othern\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[14],
                                                                                      redundant_datas[14]),
                                         '4.5SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[15], redundant_datas[15]),
                                         '16SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[16], redundant_datas[16]),
                                         "5'ETSn\nO:{:.1%} ,I:{:.1%}".format(unique_datas[17], redundant_datas[17]),
                                         "3'ETSn\nO:{:.1%} ,I:{:.1%}".format(unique_datas[18], redundant_datas[18]),
                                         'ITS1n\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[19], redundant_datas[19]),
                                         'ITS2n\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[20], redundant_datas[20]),
                                         '5SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[21], redundant_datas[21])]
            elif chain == 'positive':
                if species == 'human':
                    percentage_legend = ['28SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[0], redundant_datas[0]),
                                         '18SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[1], redundant_datas[1]),
                                         '5.8SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[2], redundant_datas[2]),
                                         '45S_rRNA_otherp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[3],
                                                                                      redundant_datas[3]),
                                         '12SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[4], redundant_datas[4]),
                                         '16SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[5], redundant_datas[5]),
                                         "5'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[6], redundant_datas[6]),
                                         "3'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[7], redundant_datas[7]),
                                         'ITS1p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[8], redundant_datas[8]),
                                         'ITS2p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[9], redundant_datas[9]),
                                         '5SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[10], redundant_datas[10])]
                elif species == 'mouse':
                    percentage_legend = ['28SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[0], redundant_datas[0]),
                                         '18SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[1], redundant_datas[1]),
                                         '5.8SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[2], redundant_datas[2]),
                                         '45S_rRNA_otherp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[3],
                                                                                      redundant_datas[3]),
                                         '4.5SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[4], redundant_datas[4]),
                                         '16SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[5], redundant_datas[5]),
                                         "5'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[6], redundant_datas[6]),
                                         "3'ETSp\nO:{:.1%} ,I:{:.1%}".format(unique_datas[7], redundant_datas[7]),
                                         'ITS1p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[8], redundant_datas[8]),
                                         'ITS2p\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[9], redundant_datas[9]),
                                         '5SrRNAp\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[10], redundant_datas[10])]
            elif chain == 'negative':
                if species == 'human':
                    percentage_legend = [
                        '28SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[0], redundant_datas[0]), '18SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[1], redundant_datas[1]), '5.8SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[2], redundant_datas[2]), '45S_rRNA_othern\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[3], redundant_datas[3]), '12SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[4], redundant_datas[4]), '16SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[5], redundant_datas[5]), "5'ETSn\nO:{:.1%} ,I:{:.1%}".format(
                                unique_datas[6], redundant_datas[6]), "3'ETSn\nO:{:.1%} ,I:{:.1%}".format(
                                    unique_datas[7], redundant_datas[7]), 'ITS1n\nO:{:.1%} ,I:{:.1%}'.format(
                                        unique_datas[8], redundant_datas[8]), 'ITS2n\nO:{:.1%} ,I:{:.1%}'.format(
                                            unique_datas[9], redundant_datas[9]), '5SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(
                                                unique_datas[10], redundant_datas[10])]
                elif species == 'mouse':
                    percentage_legend = ['28SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[0], redundant_datas[0]),
                                         '18SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[1], redundant_datas[1]),
                                         '5.8SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[2], redundant_datas[2]),
                                         '45S_rRNA_othern\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[3],
                                                                                      redundant_datas[3]),
                                         '4.5SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[4], redundant_datas[4]),
                                         '16SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[5], redundant_datas[5]),
                                         "5'ETSn\nO:{:.1%} ,I:{:.1%}".format(unique_datas[6], redundant_datas[6]),
                                         "3'ETSn\nO:{:.1%} ,I:{:.1%}".format(unique_datas[7], redundant_datas[7]),
                                         'ITS1n\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[8], redundant_datas[8]),
                                         'ITS2n\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[9], redundant_datas[9]),
                                         '5SrRNAn\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[10], redundant_datas[10])]
            else:
                if species == 'human':
                    percentage_legend = [
                        '28SrRNA\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[0], redundant_datas[0]), '18SrRNA\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[1], redundant_datas[1]), '5.8SrRNA\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[2], redundant_datas[2]), '45S_rRNA_other\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[3], redundant_datas[3]), '12SrRNA\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[4], redundant_datas[4]), '16SrRNA\nO:{:.1%} ,I:{:.1%}'.format(
                            unique_datas[5], redundant_datas[5]), "5'ETS\nO:{:.1%} ,I:{:.1%}".format(
                                unique_datas[6], redundant_datas[6]), "3'ETS\nO:{:.1%} ,I:{:.1%}".format(
                                    unique_datas[7], redundant_datas[7]), 'ITS1\nO:{:.1%} ,I:{:.1%}'.format(
                                        unique_datas[8], redundant_datas[8]), 'ITS2\nO:{:.1%} ,I:{:.1%}'.format(
                                            unique_datas[9], redundant_datas[9]), '5SrRNA\nO:{:.1%} ,I:{:.1%}'.format(
                                                unique_datas[10], redundant_datas[10])]
                elif species == 'mouse':
                    percentage_legend = ['28SrRNA\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[0], redundant_datas[0]),
                                         '18SrRNA\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[1], redundant_datas[1]),
                                         '5.8SrRNA\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[2], redundant_datas[2]),
                                         '45S_rRNA_other\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[3],
                                                                                     redundant_datas[3]),
                                         '4.5SrRNA\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[4], redundant_datas[4]),
                                         '16SrRNA\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[5], redundant_datas[5]),
                                         "5'ETS\nO:{:.1%} ,I:{:.1%}".format(unique_datas[6], redundant_datas[6]),
                                         "3'ETS\nO:{:.1%} ,I:{:.1%}".format(unique_datas[7], redundant_datas[7]),
                                         'ITS1\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[8], redundant_datas[8]),
                                         'ITS2\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[9], redundant_datas[9]),
                                         '5SrRNA\nO:{:.1%} ,I:{:.1%}'.format(unique_datas[10], redundant_datas[10])]
            ax.legend(
                wedges,
                percentage_legend,
                title="         rRNA\n   Outer:Uniqe\n   Inner:redundant",
                loc="center left",
                bbox_to_anchor=(
                    1.2,
                    0,
                    1,
                    1),
                fontsize=12)
        if species == 'human':
            referencename = collections.OrderedDict(
                [
                    ('28S',
                     'rRNA_28S'),
                    ('18S',
                     'rRNA_18S'),
                    ('5.8S',
                     'rRNA_5_8S'),
                    ('45S_across_region',
                     'junction_45S'),
                    ('12S',
                     'rRNA_12S'),
                    ('16S',
                     'rRNA_16S'),
                    ('5ETS',
                     "ETS5"),
                    ('3ETS',
                     "ETS3"),
                    ('ITS1',
                     'ITS_1'),
                    ('ITS2',
                     'ITS_2'),
                    ('5S',
                     'rRNA_5S')])
            flag_references = collections.OrderedDict(
                [('28S/0', 'rRNA_28S'), ('18S/0', 'rRNA_18S'), ('5.8S/0', 'rRNA_5_8S'),
                 ('45S_across_region/0', 'junction_45S'), ('12S/0', 'rRNA_12S'), ('16S/0', 'rRNA_16S'),
                 ('5ETS/0', "ETS5"), ('3ETS/0', "ETS3"), ('ITS1/0', 'ITS_1'), ('ITS2/0', 'ITS_2'),
                 ('5S/0', 'rRNA_5S'), ('28S/16', 'rRNA_28S'), ('18S/16', 'rRNA_18S'), ('5.8S/16', 'rRNA_5_8S'),
                 ('45S_across_region/16', 'junction_45S'), ('12S/16', 'rRNA_12S'), ('16S/16', 'rRNA_16S'),
                 ('5ETS/16', "ETS5"), ('3ETS/16', "ETS3"), ('ITS1/16', 'ITS_1'), ('ITS2/16', 'ITS_2'),
                 ('5S/16', 'rRNA_5S')])
            flag_refs = [
                '28SrRNA+',
                '18SrRNA+',
                '5.8SrRNA+',
                '45S_rRNA_other+',
                '12SrRNA+',
                '16SrRNA+',
                '',
                '',
                '',
                '',
                '5SrRNA+',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                'negative chain']
            positive_refs = [
                '28SrRNA+',
                '18SrRNA+',
                '5.8SrRNA+',
                '45S_rRNA_other+',
                '12SrRNA+',
                '16SrRNA+',
                '',
                '',
                '',
                '',
                '5SrRNA+']
            negative_refs = [
                '28SrRNA-',
                '18SrRNA-',
                '5.8SrRNA-',
                '45S_rRNA_other-',
                '12SrRNA-',
                '16SrRNA-',
                '',
                '',
                '',
                '',
                '5SrRNA-']
            refs = [
                '28SrRNA',
                '18SrRNA',
                '5.8SrRNA',
                '45S_rRNA_other',
                '12SrRNA',
                '16SrRNA',
                '',
                '',
                '',
                '',
                '5SrRNA']
            gene_length = [
                5066,
                1869,
                157,
                540,
                954,
                1559,
                3654,
                361,
                1077,
                1167,
                121,
                5066,
                1869,
                157,
                540,
                954,
                1559,
                3654,
                361,
                1077,
                1167,
                121]
        elif species == 'mouse':
            referencename = collections.OrderedDict(
                [
                    ('28S',
                     'rRNA_28S'),
                    ('18S',
                     'rRNA_18S'),
                    ('5.8S',
                     'rRNA_5_8S'),
                    ('45S_across_region',
                     'junction_45S'),
                    ('4.5S',
                     'rRNA_12S'),
                    ('16S',
                     'rRNA_16S'),
                    ('5ETS',
                     "ETS5"),
                    ('3ETS',
                     "ETS3"),
                    ('ITS1',
                     'ITS_1'),
                    ('ITS2',
                     'ITS_2'),
                    ('5S',
                     'rRNA_5S')])
            flag_references = collections.OrderedDict(
                [('28S/0', 'rRNA_28S'), ('18S/0', 'rRNA_18S'), ('5.8S/0', 'rRNA_5_8S'),
                 ('45S_across_region/0', 'junction_45S'), ('4.5S/0', 'rRNA_4.5S'), ('16S/0', 'rRNA_16S'),
                 ('5ETS/0', "ETS5"), ('3ETS/0', "ETS3"), ('ITS1/0', 'ITS_1'), ('ITS2/0', 'ITS_2'),
                 ('5S/0', 'rRNA_5S'), ('28S/16', 'rRNA_28S'), ('18S/16', 'rRNA_18S'), ('5.8S/16', 'rRNA_5_8S'),
                 ('45S_across_region/16', 'junction_45S'), ('4.5S/16', 'rRNA_4.5S'), ('16S/16', 'rRNA_16S'),
                 ('5ETS/16', "ETS5"), ('3ETS/16', "ETS3"), ('ITS1/16', 'ITS_1'), ('ITS2/16', 'ITS_2'),
                 ('5S/16', 'rRNA_5S')])
            flag_refs = [
                '28SrRNA+',
                '18SrRNA+',
                '5.8SrRNA+',
                '45S_rRNA_other+',
                '4.5SrRNA+',
                '',
                '',
                '',
                '',
                '',
                '5SrRNA+',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                'negative chain']
            positive_refs = [
                '28SrRNA+',
                '18SrRNA+',
                '5.8SrRNA+',
                '45S_rRNA_other+',
                '4.5SrRNA+',
                '',
                '',
                '',
                '',
                '',
                '5SrRNA+']
            negative_refs = [
                '28SrRNA-',
                '18SrRNA-',
                '5.8SrRNA-',
                '45S_rRNA_other-',
                '4.5SrRNA-',
                '',
                '',
                '',
                '',
                '',
                '5SrRNA-']
            refs = [
                '28SrRNA',
                '18SrRNA',
                '5.8SrRNA',
                '45S_rRNA_other',
                '4.5SrRNA',
                '',
                '',
                '',
                '',
                '',
                '5SrRNA']
            gene_length = [
                4730,
                1870,
                157,
                540,
                174,
                110,
                4007,
                551,
                1000,
                1088,
                121,
                4730,
                1870,
                157,
                540,
                174,
                110,
                4007,
                551,
                1000,
                1088,
                121]
        positive_referencename = referencename.copy()
        negative_referencename = referencename.copy()

        if chain == 'all':
            for flag_ref in flag_references.keys():
                flag_references[flag_ref] = self.count_different_map(
                    flag=flag_ref.split('/')[1], rname=flag_ref.split('/')[0], threshold=threshold)
            unique_data = list(zip(*flag_references.values()))[0]
            redundant_data = list(zip(*flag_references.values()))[1]
            normalize_unique_data = [
                i / k for i,
                k in zip(
                    unique_data,
                    gene_length)]
            normalize_redundant_data = [
                i / k for i,
                k in zip(
                    redundant_data,
                    gene_length)]
            fig = plt.figure(figsize=(24, 12), tight_layout=True)
            ax1 = fig.add_subplot(121)
            cmapi = plt.get_cmap('tab20')
            colnumbers = [i for i in range(18)]
            last_four_colors = [17, 17, 17, 18]
            colnumbers = colnumbers + last_four_colors
            colors = cmapi(colnumbers)
            plot_pie(ax1, unique_data, redundant_data, flag_refs)
            ax1.set_title('All Mapped Sequence Distribution', fontsize=20)
            ax2 = fig.add_subplot(122)
            plot_pie(
                ax2,
                normalize_unique_data,
                normalize_redundant_data,
                flag_refs)
            ax2.set_title(
                'All Mapped Sequence Normalization Distribution',
                fontsize=20)
        elif chain == 'positive':
            for refname in positive_referencename.keys():
                positive_referencename[refname] = self.count_different_map(
                    flag='0', rname=refname, threshold=threshold)
            unique_data = list(zip(*flag_references.values()))[0]
            redundant_data = list(zip(*flag_references.values()))[1]
            normalize_unique_data = [
                i / k for i,
                k in zip(
                    unique_data,
                    gene_length)]
            normalize_redundant_data = [
                i / k for i,
                k in zip(
                    redundant_data,
                    gene_length)]
            fig = plt.figure(figsize=(24, 12), tight_layout=True)
            ax1 = fig.add_subplot(121)
            cmapi = plt.get_cmap('Set3')
            colnumbers = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1]
            colors = cmapi(colnumbers)
            plot_pie(ax1, unique_data, redundant_data, positive_refs)
            ax1.set_title('Positive Mapped Sequence Distribution', fontsize=20)
            ax2 = fig.add_subplot(122)
            plot_pie(
                ax2,
                normalize_unique_data,
                normalize_redundant_data,
                positive_refs)
            ax2.set_title(
                'Positive Mapped Sequence Normalization Distribution',
                fontsize=20)
        elif chain == 'negative':
            for refname in negative_referencename.keys():
                negative_referencename[refname] = self.count_different_map(
                    flag='16', rname=refname, threshold=threshold)
            unique_data = list(zip(*flag_references.values()))[0]
            redundant_data = list(zip(*flag_references.values()))[1]
            normalize_unique_data = [
                i / k for i,
                k in zip(
                    unique_data,
                    gene_length)]
            normalize_redundant_data = [
                i / k for i,
                k in zip(
                    redundant_data,
                    gene_length)]
            fig = plt.figure(figsize=(24, 12), tight_layout=True)
            ax1 = fig.add_subplot(121)
            cmapi = plt.get_cmap('Set3')
            colnumbers = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1]
            colors = cmapi(colnumbers)
            plot_pie(ax1, unique_data, redundant_data, negative_refs)
            ax1.set_title('Negative Mapped Sequence Distribution', fontsize=20)
            ax2 = fig.add_subplot(122)
            plot_pie(
                ax2,
                normalize_unique_data,
                normalize_redundant_data,
                negative_refs)
            ax2.set_title(
                'Negative Mapped Sequence Normalization Distribution',
                fontsize=20)
        else:
            for refname in referencename.keys():
                referencename[refname] = self.count_different_map(
                    rname=refname, threshold=threshold)
            unique_data = list(zip(*flag_references.values()))[0]
            redundant_data = list(zip(*flag_references.values()))[1]
            normalize_unique_data = [
                i / k for i,
                k in zip(
                    unique_data,
                    gene_length)]
            normalize_redundant_data = [
                i / k for i,
                k in zip(
                    redundant_data,
                    gene_length)]
            fig = plt.figure(figsize=(24, 12), tight_layout=True)
            ax1 = fig.add_subplot(121)
            cmapi = plt.get_cmap('Set3')
            colnumbers = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1]
            colors = cmapi(colnumbers)
            plot_pie(ax1, unique_data, redundant_data, refs)
            ax1.set_title('Mapped Sequence Distribution', fontsize=20)
            ax2 = fig.add_subplot(122)
            plot_pie(
                ax2,
                normalize_unique_data,
                normalize_redundant_data,
                refs)
            ax2.set_title(
                'Mapped Sequence Normalization Distribution',
                fontsize=20)
        if save_directory:
            plt.savefig(save_directory, format="pdf", dpi=200)
        else:
            ask_save_fig = input('Do you want to save the figure?(y/n)')
            if ask_save_fig == 'y':
                self.recursion_savefig()

    def plot_simple_length_distribution(
            self,
            flag='',
            rname='',
            normalize=False,
            threshold=0,
            save_directory=False):
        """
        plot simple length distribution about single flag or rname or both flag and rname
        :param flag: positive chain or negative chain
        :param rname: reference gene
        :param normalize: if True,use gene length to normaliza
        :param threshold: filter sequence use threshold
        :return:length distribution bar picture，you can use plt.savefig() to store the result
        """
        def plot_simple_bar(
                ax,
                unique_array,
                redundant_array,
                ylabel='Read counts'):
            """
            simple bar workflow
            :param ax:subplot name
            :param unique_array
            :param redundant_array
            :return:simple bar picture:unique vs redundant
            """
            x = np.array([i for i in range(16, 47)])
            ax.bar(
                x - 0.225,
                unique_array,
                width=0.45,
                color='peru',
                label='unique')
            ax.bar(
                x + 0.225,
                redundant_array,
                width=0.45,
                color='SkyBlue',
                label='redundant')
            ax.set_xlabel('Length')
            ax.set_ylabel(ylabel)
            ax.legend(loc='best')
            ax.set_xticks([15, 20, 25, 30, 35, 40, 45, 46])
            ax.set_xticklabels(
                ['15', '20', '25', '30', '35', '40', '45', '>45'])
        if self.species == 'human':
            gene_length = collections.OrderedDict([('28S', 5066), ('18S', 1869), ('5.8S', 157), ('5S', 121),
                                                   ('12S', 954), ('16S', 1559), ('45S_across_region', 540),
                                                   ('5ETS', 3654), ('3ETS', 361), ('ITS1', 1077), ('ITS2', 1167),
                                                   ('45S', 13351), ('', 1)])
        elif self.species == 'mouse':
            gene_length = collections.OrderedDict([('28S', 4730), ('18S', 1870), ('5.8S', 157), ('5S', 121),
                                                   ('4.5S', 174), ('16S', 110), ('45S_across_region', 540),
                                                   ('5ETS', 4007), ('3ETS', 551), ('ITS1', 1000), ('ITS2', 1088),
                                                   ('45S', 13400), ('', 1)])
        flag_chain = {
            '0': '(Positive Chain)',
            '16': '(Negative Chain)',
            '': ''}
        unique_referencename, redundant_referencename = self.get_length_distribution(
            flag=flag, rname=rname, threshold=threshold)
        unique_referencename = np.array(list(self.length_distribution(
            unique_referencename).values()))
        redundant_referencename = np.array(list(self.length_distribution(
            redundant_referencename).values()))
        if normalize:
            normalization_unique_referencename = unique_referencename / \
                gene_length[rname]
            normalization_redundant_referencename = redundant_referencename / \
                gene_length[rname]
            fig, axes = plt.subplots(figsize=(16, 18), nrows=2, ncols=1)
            ax0, ax1 = axes.flatten()
            plot_simple_bar(ax0, unique_referencename, redundant_referencename)
            ax0.set_title(
                rname +
                'rsRNA Length distribution:Unique vs Redundant' +
                flag_chain[flag],
                fontsize=15)
            plot_simple_bar(
                ax1,
                normalization_unique_referencename,
                normalization_redundant_referencename)
            ax1.set_title(
                rname +
                'rsRNA Normalization Length distribution:Unique vs Redundant' +
                flag_chain[flag],
                fontsize=15)
        else:
            fig, ax = plt.subplots(figsize=(12, 5))
            plot_simple_bar(ax, unique_referencename, redundant_referencename)
            ax.set_title(
                rname +
                'rsRNA Length distribution:Unique vs Redundant' +
                flag_chain[flag],
                fontsize=15)
        if save_directory:
            plt.savefig(save_directory, format="pdf", dpi=200)
        else:
            ask_save_fig = input('Do you want to save the figure?(y/n)')
            if ask_save_fig == 'y':
                self.recursion_savefig()

    def plot_stack_length_distribution(
            self,
            chain='',
            normalize=False,
            threshold=0,
            aspect=True,
            contrast=False,
            subplot='',
            type='',
            save_directory=False):
        """
        plot stack bar of length distribution
        three type pictures you can choose to draw:
        one is standard stack bar plot:above is unique,below is redundant
        second is stack bar plot and add thumbnails below the pic
        third is contrast unique and redundant in the same axes
        chain == 'all' can't draw normalize figure!
        :param chain:'All':distinct positive and all negative,
                     'positive':distinct positive,
                     'negative':distinct negative,
                     '':all distinct positive+negative(default)
        :param contrast:if True, plot unique and redundant length distribution stack bar in the same axes
        :param normalize:if True, nomalization with reference gene length
        :param subplot:add thumbnails under the stack bar pic,two_type subplot,use
        :param type:when plot second picture,you should choose unique or other specific to paint,
        if 'unique',draw unique,
        otherwise draw redundant(default)
        :return: show picture,you can use plt.savefig to store the pic
        """
        def transfer_array(unique_referencename, redundant_referencename):
            """
            replication code
            """
            for reference in unique_referencename.keys():
                unique_referencename[reference] = np.array(list(self.length_distribution(
                    unique_referencename[reference]).values()))
                redundant_referencename[reference] = np.array(list(self.length_distribution(
                    redundant_referencename[reference]).values()))
            return unique_referencename, redundant_referencename

        def plot(unique_referencename, redundant_referencename, subplots=12):
            """
            replication code
            """
            if subplot:
                stack_bar_subplot(
                    unique_referencename,
                    redundant_referencename,
                    subplots=subplots)
            else:
                stack_bar_plot(unique_referencename, redundant_referencename)

        def stack_bar_plot(unique_referencename, redundant_referencename):
            """
            use two static function and a simple local function
            plot stack bar picture workflow
            :param unique_referencename:unique sequence,a dictionary and its values also are dictionary
            :param redundant_referencename:redundant sequence,a dictionary and its values also are dictionary
            :return: show picture,you can use plt.savefig to store the pic
            """
            unique_referencename, redundant_referencename = transfer_array(
                unique_referencename, redundant_referencename)
            del unique_referencename['45S']
            del redundant_referencename['45S']
            if normalize:
                normalization_unique_referencename = unique_referencename.copy()
                normalization_redundant_referencename = unique_referencename.copy()
                for rname, counts_array in unique_referencename.items():
                    normalization_unique_referencename[rname] = counts_array / \
                        gene_length[rname]
                for rname, counts_array in redundant_referencename.items():
                    normalization_redundant_referencename[rname] = counts_array / \
                        gene_length[rname]
                fig, axes = plt.subplots(figsize=(16, 18), nrows=2, ncols=1)
                ax0, ax1 = axes.flatten()
                self.plot_stack_bar(ax0, normalization_unique_referencename)
                ax0.set_title(
                    'Normalization Unique Length Distribution ' +
                    chain,
                    fontsize=18,
                    pad=35)
                self.plot_stack_bar(ax1, normalization_redundant_referencename)
                ax1.set_title(
                    'Normalization Redundant Length Distribution ' +
                    chain,
                    fontsize=18,
                    pad=35)
            else:
                fig, axes = plt.subplots(figsize=(16, 18), nrows=2, ncols=1)
                ax0, ax1 = axes.flatten()
                self.plot_stack_bar(ax0, unique_referencename)
                ax0.set_title(
                    'Unique Length Distribution ' +
                    chain,
                    fontsize=18,
                    pad=35)
                self.plot_stack_bar(ax1, redundant_referencename)
                ax1.set_title(
                    'Redundant Length Distribution ' +
                    chain,
                    fontsize=18,
                    pad=35)

        def stack_bar_subplot(
                unique_referencename,
                redundant_referencename,
                subplots=12):
            """
            draw stack bar pic and add thumbnails under the figure
            :param unique_referencename:unique sequence,a dictionary and its values also are dictionary
            :param redundant_referencename:redundant sequence,a dictionary and its values also are dictionary
            :return: show picture,you can use plt.savefig to store the pic
            """
            def subplot_stack_bar(
                    stack_length_dic,
                    title='',
                    ylabel='Read counts'):
                """
                subplot workflow
                :param stack_length_dic: a dictionary and its values are arrays
                :param title:eg: Unique
                :return:show pic
                """
                fig = plt.figure(figsize=(15, 9), dpi=200)
                gs = mg.GridSpec(5, 1)
                ax = plt.subplot(gs[:3, :])
                ax.set_title(
                    title +
                    ' Length Distribution',
                    fontsize=18,
                    pad=35)
                stack_length_dic_copy = stack_length_dic.copy()
                stack_length_dic_copy.pop('45S')
                self.plot_stack_bar(ax, stack_length_dic_copy)
                ind = np.arange(16, 47)
                width = 0.645
                colornumber = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 3]
                cmap = plt.get_cmap('Set3')
                counter = [i for i in range(len(stack_length_dic.keys()))]
                if aspect:
                    grid = ImageGrid(
                        fig, 313, (1, subplots), axes_pad=0.05, aspect=True, share_all=False)
                    for i, length_distribution in zip(
                            counter, stack_length_dic.values()):
                        grid[i].set_title(
                            references[i], fontstyle='italic', fontsize=11)
                        grid[i].bar(ind, length_distribution /
                                    (max(length_distribution) +
                                     0.000000000000000000000000000000000000000000001) *
                                    70, width=width, color=cmap(colornumber[i]))
                    grid[0].set_yticklabels([''])
                else:
                    grid = ImageGrid(
                        fig, 313, (1, subplots), axes_pad=0.05, aspect=False, share_all=True)
                    for i, length_distribution in zip(
                            counter, stack_length_dic.values()):
                        grid[i].set_title(
                            references[i], fontstyle='italic', fontsize=11)
                        grid[i].bar(
                            ind,
                            length_distribution,
                            width=width,
                            color=cmap(
                                colornumber[i]))
                grid[5].set_xlabel('Length', fontsize=15)
                grid[0].set_ylabel(ylabel, fontsize=15)
            unique_referencename, redundant_referencename = transfer_array(
                unique_referencename, redundant_referencename)
            if normalize:
                normalization_unique_referencename = unique_referencename.copy()
                normalization_redundant_referencename = unique_referencename.copy()
                for rname, counts_array in unique_referencename.items():
                    normalization_unique_referencename[rname] = counts_array / \
                        gene_length[rname]
                for rname, counts_array in redundant_referencename.items():
                    normalization_redundant_referencename[rname] = counts_array / \
                        gene_length[rname]
                if type == 'unique':
                    subplot_stack_bar(
                        normalization_unique_referencename,
                        title='Normalization Unique')
                else:
                    subplot_stack_bar(
                        normalization_redundant_referencename,
                        title='Normalization Redundant')
            else:
                if type == 'unique':
                    subplot_stack_bar(unique_referencename, title='Unique')
                else:
                    subplot_stack_bar(
                        redundant_referencename, title='Redundant')

        def stack_bar_contrast(unique_referencename, redundant_referencename):
            """
            contrast stack bar plot
            :param unique_referencename:unique sequence,a dictionary and its values also are dictionary
            :param redundant_referencename:redundant sequence,a dictionary and its values also are dictionary
            :return: show picture,you can use plt.savefig to store the pic
            """
            unique_referencename, redundant_referencename = transfer_array(
                unique_referencename, redundant_referencename)
            del unique_referencename['45S']
            del redundant_referencename['45S']
            if normalize:
                normalization_unique_referencename = unique_referencename.copy()
                normalization_redundant_referencename = unique_referencename.copy()
                for rname, counts_array in unique_referencename.items():
                    normalization_unique_referencename[rname] = counts_array / \
                        gene_length[rname]
                for rname, counts_array in redundant_referencename.items():
                    normalization_redundant_referencename[rname] = counts_array / \
                        gene_length[rname]
                fig, axes = plt.subplots(figsize=(16, 18), nrows=2, ncols=1)
                ax0, ax1 = axes.flatten()
                self.plot_stack_bar(ax0,
                                    unique_referencename,
                                    bar_location=(-0.225),
                                    width=0.45)
                self.plot_stack_bar(
                    ax0,
                    redundant_referencename,
                    bar_location=0.225,
                    width=0.45,
                    legend=False)
                ax0.set_title(
                    'Length Distribution:unique vs redundant',
                    fontsize=15,
                    pad=35)
                self.plot_stack_bar(ax1,
                                    normalization_unique_referencename,
                                    bar_location=(-0.225),
                                    width=0.45)
                self.plot_stack_bar(
                    ax1,
                    normalization_redundant_referencename,
                    bar_location=0.225,
                    width=0.45,
                    legend=False)
                ax1.set_title(
                    'Normalization Length Distribution:unique vs redundant',
                    fontsize=15,
                    pad=35)
                plt.tight_layout()
            else:
                fig, ax = plt.subplots(figsize=(12, 5))
                self.plot_stack_bar(ax, unique_referencename,
                                    bar_location=(-0.225), width=0.45)
                self.plot_stack_bar(
                    ax,
                    redundant_referencename,
                    bar_location=0.225,
                    width=0.45,
                    legend=False)
                ax.set_title(
                    'Length Distribution:unique vs redundant',
                    fontsize=15,
                    pad=35)
                plt.tight_layout()
        if self.species == 'human':
            gene_length = collections.OrderedDict([('28S', 5066), ('18S', 1869), ('5.8S', 157), ('5S', 121),
                                                   ('12S', 954), ('16S', 1559), ('45S_across_region', 540),
                                                   ('5ETS', 3654), ('3ETS', 361), ('ITS1', 1077), ('ITS2', 1167),
                                                   ('45S', 13351)])
            if chain == 'all':
                references = [
                    '28SrRNA',
                    '18SrRNA',
                    '5.8SrRNA',
                    '5SrRNA',
                    '12SrRNA',
                    '16SrRNA',
                    '45S_across',
                    "5'ETS",
                    "3'ETS",
                    'ITS1',
                    'ITS2',
                    'negative',
                    '45S']
                unique_referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'),
                                                                ('5.8S', 'rRNA_5_8S'), ('5S', 'rRNA_5S'),
                                                                ('12S', 'rRNA_12S'), ('16S', 'rRNA_16S'),
                                                                ('45S_across_region', 'junction_45S'),
                                                                ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                                ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'),
                                                                ('negative', 'negative_chain'), ('45S', 'rRNA_45S')])
            else:
                references = [
                    '28SrRNA',
                    '18SrRNA',
                    '5.8SrRNA',
                    '5SrRNA',
                    '12SrRNA',
                    '16SrRNA',
                    '45S_across',
                    "5'ETS",
                    "3'ETS",
                    'ITS1',
                    'ITS2',
                    '45S']
                unique_referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'),
                                                                ('5.8S', 'rRNA_5_8S'), ('5S', 'rRNA_5S'),
                                                                ('12S', 'rRNA_12S'), ('16S', 'rRNA_16S'),
                                                                ('45S_across_region', 'junction_45S'),
                                                                ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                                ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'),
                                                                ('45S', 'rRNA_45S')])
        elif self.species == 'mouse':
            gene_length = collections.OrderedDict([('28S', 4730), ('18S', 1870), ('5.8S', 157), ('5S', 121),
                                                   ('4.5S', 174), ('16S', 110), ('45S_across_region', 540),
                                                   ('5ETS', 4007), ('3ETS', 551), ('ITS1', 1000), ('ITS2', 1088),
                                                   ('45S', 13400)])
            if chain == 'all':
                references = [
                    '28SrRNA',
                    '18SrRNA',
                    '5.8SrRNA',
                    '5SrRNA',
                    '4.5SrRNA',
                    '16SrRNA',
                    '45S_across',
                    "5'ETS",
                    "3'ETS",
                    'ITS1',
                    'ITS2',
                    'negative',
                    '45S']
                unique_referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'),
                                                                ('5.8S', 'rRNA_5_8S'), ('5S', 'rRNA_5S'),
                                                                ('4.5S', 'rRNA_4.5S'), ('16S', 'rRNA_16S'),
                                                                ('45S_across_region', 'junction_45S'),
                                                                ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                                ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'),
                                                                ('negative', 'negative_chain'), ('45S', 'rRNA_45S')])
            else:
                references = [
                    '28SrRNA',
                    '18SrRNA',
                    '5.8SrRNA',
                    '5SrRNA',
                    '4.5SrRNA',
                    '16SrRNA',
                    '45S_across',
                    "5'ETS",
                    "3'ETS",
                    'ITS1',
                    'ITS2',
                    '45S']
                unique_referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'),
                                                                ('5.8S', 'rRNA_5_8S'), ('5S', 'rRNA_5S'),
                                                                ('4.5S', 'rRNA_4.5S'), ('16S', 'rRNA_16S'),
                                                                ('45S_across_region', 'junction_45S'),
                                                                ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                                ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'),
                                                                ('45S', 'rRNA_45S')])
        redundant_referencename = unique_referencename.copy()
        if chain == 'all':
            for reference in unique_referencename.keys():
                if reference == 'negative':
                    unique_referencename[reference], redundant_referencename[reference] = self.get_length_distribution(
                        flag='16', threshold=threshold)
                else:
                    unique_referencename[reference], redundant_referencename[reference] = self.get_length_distribution(
                        flag='0', rname=reference, threshold=threshold)
            if contrast:
                stack_bar_contrast(
                    unique_referencename,
                    redundant_referencename)
            else:
                plot(
                    unique_referencename,
                    redundant_referencename,
                    subplots=13)
        elif chain == 'positive':
            for reference in unique_referencename.keys():
                unique_referencename[reference], redundant_referencename[reference] = self.get_length_distribution(
                    flag='0', rname=reference, threshold=threshold)
            if contrast:
                stack_bar_contrast(
                    unique_referencename,
                    redundant_referencename)
            else:
                plot(unique_referencename, redundant_referencename)
        elif chain == 'negative':
            for reference in unique_referencename.keys():
                unique_referencename[reference], redundant_referencename[reference] = self.get_length_distribution(
                    flag='16', rname=reference, threshold=threshold)
            if contrast:
                stack_bar_contrast(
                    unique_referencename,
                    redundant_referencename)
            else:
                plot(unique_referencename, redundant_referencename)
        else:
            for reference in unique_referencename.keys():
                unique_referencename[reference], redundant_referencename[reference] = self.get_length_distribution(
                    rname=reference, threshold=threshold)
            if contrast:
                stack_bar_contrast(
                    unique_referencename,
                    redundant_referencename)
            else:
                plot(unique_referencename, redundant_referencename)
        if save_directory:
            plt.savefig(save_directory, format="pdf", dpi=200)
        else:
            ask_save_fig = input('Do you want to save the figure?(y/n)')
            if ask_save_fig == 'y':
                self.recursion_savefig()

    def location_distribution_pool(
            self,
            flag='',
            type='',
            threshold=0,
            count_format='COUNT',
            ylabel='Read counts',
            save_directory=''):
        """
        gather all location distribution(12 mapping results in one pdf)
        9 choices to plot
        :param flag: positive('0') or negative('16') or positive+negative(else)
        :param type:unique('unique') or redundant('redundant') or unique+redundant(else)
        :param threshold: RPM or Read counts threshold
        :param count_format: RPM or COUNT
        :param ylabel: RPM or Read counts
        :param save_directory:save the result as a pdf
        :return: show figure and ask for whether saving or not
        """
        def transfer_numpy(location_distribution_dic):
            """
            repication code:transfer dic to numpy ndarray
            :param location_distribution_dic
            :return: 2d array 1:location 2:read counts
            """
            location_distribution_dic = location_distribution_dic.items()
            location_distribution_dic = list(zip(*location_distribution_dic))
            location_distribution_dic = np.array(location_distribution_dic)
            return location_distribution_dic
        if self.species == 'human':
            gene_length = collections.OrderedDict([('28S', 5066), ('18S', 1869), ('5.8S', 157), ('5S', 121),
                                                   ('12S', 954), ('16S', 1559), ('45S', 13351),
                                                   ('45S_across_region', 13351), ('5ETS', 3654),
                                                   ('3ETS', 361), ('ITS1', 1077), ('ITS2', 1167)])
        elif self.species == 'mouse':
            gene_length = collections.OrderedDict([('28S', 4730), ('18S', 1870), ('5.8S', 157), ('5S', 121),
                                                   ('4.5S', 174), ('16S', 110), ('45S', 13400),
                                                   ('45S_across_region', 13400),
                                                   ('5ETS', 4007), ('3ETS', 551), ('ITS1', 1000), ('ITS2', 1088)])
        fig, axes = plt.subplots(figsize=(16, 48), nrows=12, ncols=1, dpi=150)
        fig.suptitle('rRNA Mapping Result', fontsize=16, y=0.895)
        i = 0
        if type == 'unique' or type == 'redundant':
            for gene, length in gene_length.items():
                location_dic = self.quick_location_dic(
                    flag=flag,
                    rname=gene,
                    type=type,
                    threshold=threshold,
                    count_format=count_format)
                x = np.array([i for i in range(1, length + 1)])
                patches = axes[i].bar(
                    x, location_dic.values(), width=1, linewidth=0)
                N = np.array(list(location_dic.values()))
                fracs = N / \
                    (N.max() + 0.00000000000000000000000000000000000000000000000000000000001)
                norm = colors.Normalize(fracs.min(), fracs.max())
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.viridis(norm(thisfrac))
                    thispatch.set_facecolor(color)
                axes[i].set_xlabel('Location', fontsize=13)
                axes[i].set_ylabel(ylabel, fontsize=13)
                axes[i].set_title(
                    gene + ' rRNA Mapping Result',
                    fontsize=14,
                    loc='left')
                i += 1
        else:
            for gene, length in gene_length.items():
                unique_location_dic, redundant_location_dic = self.quick_location_dic(
                    flag=flag, rname=gene, type=type, threshold=threshold, count_format=count_format)
                unique_location_distribution_array = transfer_numpy(
                    unique_location_dic)
                redundant_location_distribution_array = transfer_numpy(
                    redundant_location_dic)
                ax_1 = axes[i]
                ax_2 = ax_1.twinx()
                axes[i].set_title(
                    gene + ' rRNA Mapping Result',
                    fontsize=14,
                    loc='left')
                ax_1.plot(
                    unique_location_distribution_array[0],
                    unique_location_distribution_array[1],
                    color='purple',
                    label='unique')
                ax_1.fill_between(
                    unique_location_distribution_array[0],
                    unique_location_distribution_array[1],
                    where=unique_location_distribution_array[1] > 0,
                    facecolor='purple',
                    alpha=0.6)
                ax_2.plot(
                    redundant_location_distribution_array[0],
                    redundant_location_distribution_array[1],
                    color='blue',
                    label='redundant')
                ax_2.fill_between(
                    redundant_location_distribution_array[0],
                    redundant_location_distribution_array[1],
                    where=redundant_location_distribution_array[1] > 0,
                    facecolor='blue',
                    alpha=0.2)
                ax_1.legend(bbox_to_anchor=(0., 1.02, 0.88, .11), fontsize=11)
                ax_2.legend(bbox_to_anchor=(0., 1.02, 1., .11), fontsize=11)
                ax_1.set_xlabel('Location', fontsize=13)
                ax_1.set_ylabel('Unique Read counts', fontsize=13)
                ax_2.set_ylabel('Redundant ' + ylabel, fontsize=13)
                i += 1
        if save_directory:
            plt.savefig(
                save_directory,
                format="pdf",
                dpi=200,
                bbox_inches='tight')
        else:
            ask_save_fig = input('Do you want to save the figure?(y/n)')
            if ask_save_fig == 'y':
                self.recursion_savefig()


def get_location_matrix(
        sam_list='',
        flag='',
        rname='',
        type='redundant',
        threshold=0,
        count_format='RPM',
        save_directory='',
        genome_index_path='',
        species='human',
        dataframe=True):
    """
    get base expression matrix
    if dataframe get dataframe
    else get ndarray and sample information
    :param sam_list: srr_list include path and srr ID
    :param flag: positive or negative or positive+negative and as a part of filename
    :param rname: reference gene
    :param type: unique or redundant not both
    :param threshold: int, rpm or read counts threshold
    :param count_format: RPM or COUNT
    :param save_directory: path+filename
    :param genome_index_path: genome bowtie index path,if 'RPM',you should align the sequence to genome
    and get mapping read counts for RPM normalization
    :param species: human or mouse
    :param dataframe: if dataframe get dataframe else get ndarray and sample information
    :return: if dataframe get dataframe else get ndarray and sample information,finally save the data structure
    dataframe:df1:three indexs ['Tissue','SRP','SRR'] df2:two indexs ['Tissue','Number']
    ndarray:1.ndarray 2.tissue_number ['Tissue','Number'] :a tuple list
    3.sample_info ['Tissue','SRP','SRR'] :a tuple list
    """
    with open(sam_list) as sl:
        sam_lists = sl.readlines()
    sample_info = []
    tissue_number = []
    i = 1
    for sam in sam_lists:
        sam = sam.strip()
        if sam[1:3] == 'RR' or sam[1:3] == 'RX':
            rRNA_sam_path = path + sam + '.sam'
            samclass = SamFile(rRNA_sam_path, species)
            if count_format == 'RPM':
                fastq_path = path + sam + '.fastq'
                samclass.rpm_normalization(fastq_path, genome_index_path)
            if i == 1:
                location_numpy = samclass.quick_location_dic(
                    flag=flag,
                    rname=rname,
                    type=type,
                    threshold=threshold,
                    count_format=count_format)
                location_numpy = np.array(list(location_numpy.values()))
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                i += 1
                j += 1
            else:
                location_series = samclass.quick_location_dic(
                    flag=flag, rname=rname, type=type, threshold=threshold, count_format=count_format)
                location_series = np.array(list(location_series.values()))
                location_numpy = np.vstack((location_numpy, location_series))
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                j += 1
        else:
            path = sam
            j = 1
    if dataframe:
        heatmap_df1 = pd.DataFrame(location_numpy, index=tuple(sample_info), columns=[
                                   i for i in range(1, len(location_numpy[0]) + 1)])
        heatmap_df1.index.names = ['Tissue', 'SRP', 'SRR']
        heatmap_df1.columns.name = 'Location'
        heatmap_df2 = pd.DataFrame(location_numpy, index=tuple(tissue_number), columns=[
                                   i for i in range(1, len(location_numpy[0]) + 1)])
        heatmap_df2.index.names = ['Tissue', 'Number']
        heatmap_df2.columns.name = 'Location'
        heatmap_df1.to_csv(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '_1.csv',
            header=True,
            index=True)
        heatmap_df2.to_csv(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '_2.csv',
            header=True,
            index=True)
        return heatmap_df1, heatmap_df2
    else:
        np.save(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '_data.npy',
            location_numpy)
        np.save(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '_sample_info1.npy',
            np.array(sample_info))
        np.save(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '_sample_info2.npy',
            np.array(tissue_number))
        return location_numpy, tissue_number, sample_info


def get_location_clustermap(
        location_data,
        sample_info='',
        log=False,
        z_score=None,
        save_directory='',
        row_cluster=False,
        col_cluster=True,
        cmap='viridis',
        yticklabels=True,
        flag='',
        rname='',
        type='',
        count_format='RPM',
        standard_scale=None,
        method='average',
        species='human'):
    """
    clustermap
    :param location_data: if ndarray you need deliver sample_info,elif datafrmae needn't deliver
    :param sample_info: ndarray index information Multiindex get from funcation get_matrix
    :param log: whether use np.log (e) to solve data or not True or False
    :param z_score: Either 0 (rows) or 1 (columns).
    Whether or not to calculate z-scores for the rows or the columns.
    Z scores are: z = (x - mean)/std, so values in each row (column)
    will get the mean of the row (column) subtracted, then divided by
    the standard deviation of the row (column).
    This ensures that each row (column) has mean of 0 and variance of 1.
    :param save_directory:path+filename
    :param row_cluster:True or False
    :param col_cluster:True or False
    :param cmap: YlGnBu RdYlBu YlOrRd hot twilight_shifted coolwarm bwr viridis
    :param yticklabels:True or False
    :param flag:as filename part for saving file
    :param rname:reference gene for title
    :param type:unique or redundant for title
    :param count_format:RPM or COUNT for colorbar x-label
    :param standard_scale:Either 0 (rows) or 1 (columns). Whether or not to standardize that dimension,
    meaning for each row or column, subtract the minimum and divide each by its maximum.
    :param method:Linkage method to use for calculating clusters,default:'average'
    :return:None,just show figure and store figure and dataframe
    """
    if log:
        location_data = np.log(location_data + 1)
    if isinstance(location_data, np.ndarray):
        if len(sample_info[0]) == 3:
            tissue, srr, srr_number = zip(*sample_info)
            index_name = ['Tissue', 'SRP', 'SRR']
        elif len(sample_info[0]) == 2:
            tissue, number = zip(*sample_info)
            index_name = ['Tissue', 'Number']
        heatmap_df = pd.DataFrame(location_data, index=tuple(sample_info), columns=[
                                  i for i in range(1, len(location_data[0]) + 1)])
        heatmap_df.index.names = index_name
        heatmap_df.columns.name = 'Location'
    elif isinstance(location_data, pd.core.frame.DataFrame):
        tissue = list(location_data.index.get_level_values('Tissue'))
        heatmap_df = location_data
    if len(set(tissue)) < 13:
        color_palette = 'Paired'
    else:
        color_palette = 'tab20'
    color_pal = sns.color_palette(color_palette, len(set(tissue)))
    single_tissue = list(set(tissue))
    single_tissue.sort(key=tissue.index)
    tissue_color = dict(zip(single_tissue, color_pal))
    side_color_bar = pd.Series(
        heatmap_df.index.get_level_values('Tissue'),
        index=heatmap_df.index).map(tissue_color)
    hm = sns.clustermap(
        heatmap_df,
        row_cluster=row_cluster,
        cmap=cmap,
        row_colors=side_color_bar,
        figsize=(
            12,
            12),
        col_cluster=col_cluster,
        yticklabels=yticklabels,
        z_score=z_score,
        standard_scale=standard_scale,
        method=method)
    hm.ax_heatmap.set_title(
        rname +
        ' rsRNA Location Distribution ' +
        type,
        fontsize=17,
        pad=90)
    hm.ax_heatmap.set_xlabel('Location', fontsize=15)
    hm.ax_heatmap.set_ylabel('Sample', fontsize=15)
    plt.xlabel(count_format)
    for tissue, color in tissue_color.items():
        hm.ax_col_dendrogram.bar(0, 0, color=color,
                                 label=tissue, linewidth=0)
    hm.ax_col_dendrogram.legend(
        loc='best', ncol=4, mode='expand', bbox_to_anchor=(-.08, .4, 1.15, 0.3))
    if save_directory:
        plt.savefig(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '.png',
            dpi=200,
            bbox_inches='tight')
        heatmap_df.to_csv(
            save_directory +
            species +
            '_' +
            rname +
            flag +
            type +
            '.csv',
            header=True,
            index=True)


def location_pool_heatmap(
        sam_list,
        flag='',
        type='redundant',
        threshold=0,
        count_format='RPM',
        save_directory='',
        genome_index_path='',
        species='human',
        log=False,
        z_score=None,
        row_cluster=True,
        col_cluster=False,
        cmap='viridis',
        yticklabels=True,
        standard_scale=None,
        method='average',
        three_index=True):
    """
    batch plot heatmap
    :param sam_list:srr_list include path and srr id
    :param flag: 0 or 16 or None corresponding to positive,negative,positive+negative
    :param type:unique or redundant
    :param threshold:rpm or read count threshold
    :param count_format:RPM or COUNT
    :param save_directory:path+filename
    :param genome_index_path:if RPM ,need map to genome ，so need give a bowtie genome index path
    :param species:human or mouse
    :param log:log e (all data + 1)
    :param z_score:z_score normalization
    :param row_cluster:True or False
    :param col_cluster:True or False
    :param cmap:YlGnBu RdYlBu YlOrRd hot twilight_shifted coolwarm bwr viridis
    :param yticklabels:True or False or a list
    :param standard_scale:1 normalization
    :param method:cluster method
    :param three_index:if True，index is ['Tissue', 'SRP', 'SRR'],otherwise index is ['Tissue','Number']
    :return:save one dataframe and one figure every SRR ID
    """
    if species == 'mouse':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('4.5S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('5S', 'rRNA_5S'),
                                                 ('45S', 'rRNA_45S')])
    elif species == 'human':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('12S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('5S', 'rRNA_5S'),
                                                 ('45S', 'rRNA_45S')])
    with open(sam_list) as sl:
        sam_lists = sl.readlines()
    sample_info = []
    tissue_number = []
    i = 0
    for sam in sam_lists:
        sam = sam.strip()
        if sam[1:3] == 'RR' or sam[1:3] == 'RX':
            rRNA_sam_path = path + sam + '.sam'
            samclass = SamFile(rRNA_sam_path, species)
            if count_format == 'RPM':
                fastq_path = path + sam + '.fastq'
                samclass.rpm_normalization(fastq_path, genome_index_path)
            tissue = rRNA_sam_path.split('/')[-2]
            srr = rRNA_sam_path.split('/')[-1][:6]
            srr_number = rRNA_sam_path.split('/')[-1][:-4]
            sample_name = tissue, srr, srr_number
            sample_info.append(sample_name)
            tissue_number.append((tissue, j))
            i += 1
            j += 1
            for rname in referencename.keys():
                if i == 1:
                    location_numpy = samclass.quick_location_dic(
                        flag=flag, rname=rname, type=type, threshold=threshold, count_format=count_format)
                    referencename[rname] = np.array(
                        list(location_numpy.values()))

                else:
                    location_series = samclass.quick_location_dic(
                        flag=flag, rname=rname, type=type, threshold=threshold, count_format=count_format)
                    location_series = np.array(list(location_series.values()))
                    referencename[rname] = np.vstack(
                        (referencename[rname], location_series))
        else:
            path = sam
            j = 1
    if three_index:
        index = sample_info
        indexname = ['Tissue', 'SRP', 'SRR']
    else:
        index = tissue_number
        indexname = ['Tissue', 'Number']
    for rname, df in referencename.items():
        heatmap_df1 = pd.DataFrame(
            df, index=tuple(index), columns=[
                i for i in range(
                    1, len(
                        df[0]) + 1)])
        heatmap_df1.index.names = indexname
        heatmap_df1.columns.name = 'Location'
        get_location_clustermap(
            heatmap_df1,
            log=log,
            z_score=z_score,
            save_directory=save_directory,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            cmap=cmap,
            yticklabels=yticklabels,
            rname=rname,
            type=type,
            count_format=count_format,
            standard_scale=standard_scale,
            method=method,
            flag=flag,
            species=species)


def image_splice(
        list_or_path,
        col,
        row,
        width,
        height,
        save_directory='',
        same_size_figure=True,
        pic_format='png'):
    """
    splice image
    :param list_or_path: path list or path txt or path
    :param col set col number
    :param row set row number
    :param width every figure width
    :param height every figure height
    :param save_directory: path to store
    :param same_size_figure: the spliced image whether same size or not True or False
    :param pic_format: default 'png'
    :return: show spliced image and if save_directory exist,save it
    """
    import os
    from PIL import Image
    import matplotlib.pyplot as plt

    def walk_directory(path):
        """
        get directory file path
        :param path: directory path
        :return: list about file path
        """
        all_path = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if "tif" or 'png' or 'jpg' or 'jpeg' or 'bmp' or 'svg' or 'gif' in file:
                    all_path.append(os.path.join(root, file))
        all_path.sort()
        return all_path

    try:
        with open(list_or_path) as listpath:
            image_lists = listpath.readlines()
    except BaseException:
        if isinstance(list_or_path, list):
            image_lists = list_or_path
        else:
            image_lists = walk_directory(list_or_path)
    if same_size_figure:
        figure = Image.open(image_lists[0])
        width, height = figure.size
    else:
        width = width
        height = height
    bg_image = Image.new(
        'RGBA', (col * width, row * height), (255, 255, 255, 255))
    num = 0
    for i in range(0, row):
        for j in range(0, col):
            pic = Image.open((image_lists[num]).strip())
            if not same_size_figure:
                pic = pic.resize((width, height))
            loc = (int(j % col * width), int(i % row * height))
            bg_image.paste(pic, loc)
            num = num + 1
            print("第" + str(num) + "存放位置" + str(loc))
            if num >= len(image_lists):
                print("break")
                break
    print(bg_image.size)
    plt.imshow(bg_image)
    if save_directory:
        bg_image.save(save_directory, format=pic_format)


def get_length_matrix(
        sam_list,
        species='human',
        flag='',
        rname='',
        type='',
        threshold=0,
        save_directory=''):
    """
    get length distribution matrix
    :param sam_list: path and srr ID list txt file path
    :param species: 'human' or 'mouse'
    :param flag: positive or negative or positive + negative
    :param rname: reference gene
    :param type: unique and redundant
    :param threshold: threshold of read counts or rpm
    :param save_directory: path to store four dataframes
    :return: four dataframes
    one is three indexs rname length distribution
    two is two indexs rname length distribution
    three is three indexs total reads length distribution
    four is two indexs total reads length distribution
    """
    def save_dataframe(array, index, indexname, suffix=''):
        """
        replication code :save dataframe
        :param array: ndarray for creating dataframe
        :param index: index of dataframe : a tuple list
        :param indexname:a list about index name
        :param suffix:for distinguishing different dataframe
        :return:dataframe about input array
        """
        columns = [str(i) for i in range(16, 46)]
        columns.append('>45')
        length_df = pd.DataFrame(array, index=tuple(index),
                                 columns=columns)
        length_df.index.names = indexname
        length_df.columns.name = 'Length'
        if save_directory:
            length_df.to_csv(
                save_directory +
                species +
                '_' +
                rname +
                flag +
                '_' +
                type +
                suffix,
                header=True,
                index=True)
        return length_df
    with open(sam_list) as sams:
        srr_lists = sams.readlines()
    i = 1
    sample_info = []
    tissue_number = []
    for sam in srr_lists:
        sam = sam.strip()
        if sam[1:3] == 'RR' or sam[1:3] == 'RX':
            rRNA_sam_path = path + sam + '.sam'
            samclass = SamFile(rRNA_sam_path, species)
            if i == 1:
                if type == 'redundant':
                    total_length_array = samclass.length_distribution(
                        samclass.total_redundant_length_distribution)
                    total_length_array = np.array(
                        list(total_length_array.values()))
                elif type == 'unique':
                    total_length_array = samclass.length_distribution(
                        samclass.total_unique_length_distribution)
                    total_length_array = np.array(
                        list(total_length_array.values()))
                length_array = samclass.length_distribution(
                    samclass.get_length_distribution(
                        flag=flag, rname=rname, type=type, threshold=threshold))
                length_array = np.array(list(length_array.values()))
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                i += 1
                j += 1
            else:
                if type == 'redundant':
                    total_length_series = samclass.length_distribution(
                        samclass.total_redundant_length_distribution)
                    total_length_series = np.array(
                        list(total_length_series.values()))
                elif type == 'unique':
                    total_length_series = samclass.length_distribution(
                        samclass.total_unique_length_distribution)
                    total_length_series = np.array(
                        list(total_length_series.values()))
                total_length_array = np.vstack(
                    (total_length_array, total_length_series))
                length_series = samclass.length_distribution(
                    samclass.get_length_distribution(
                        flag=flag, rname=rname, type=type, threshold=threshold))
                length_series = np.array(list(length_series.values()))
                length_array = np.vstack((length_array, length_series))
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                j += 1
        else:
            path = sam
            j = 1
    three_index = ['Tissue', 'SRP', 'SRR']
    two_index = ['Tissue', 'Number']
    length_df1 = save_dataframe(
        length_array,
        sample_info,
        three_index,
        '_length1.csv')
    length_df2 = save_dataframe(
        length_array,
        tissue_number,
        two_index,
        '_length2.csv')
    total_length_df1 = save_dataframe(
        total_length_array,
        sample_info,
        three_index,
        '_total_length1.csv')
    total_length_df2 = save_dataframe(
        total_length_array,
        tissue_number,
        two_index,
        '_total_length2.csv')
    return length_df1, length_df2, total_length_df1, total_length_df2


def get_length_heatmap(
        heatmap_df1,
        heatmap_df2='',
        normalization='',
        cmap='viridis',
        count_format='',
        species='human',
        flag='',
        rname='',
        type='',
        save_directory='',
        del_col='',
        save_fig=''):
    """
    get length distribution heatmap
    :param heatmap_df1: rname heatmap dataframe
    :param heatmap_df2: total read counts heatmap dataframe
    :param normalization: normalization with reads-min/(max-min)
    :param cmap: color map
    :param count_format: color bar label
    :param species: human or mouse for label
    :param flag: positive or negative or positive + negative for label
    :param rname: reference gene for label
    :param type: redundant or unique for label
    :param save_directory:path to save the dataframe csv include '/'
    :param del_col: -1 -2 or other   the last one(>45) the last two other(>45 45)
    :param save_fig: path to save figure png include '/'
    :return: save csv or png and show picture
    """
    if isinstance(heatmap_df2, pd.core.frame.DataFrame):
        if del_col:
            heatmap_df1 = heatmap_df1.iloc[:, :del_col]
            heatmap_df2 = heatmap_df2.iloc[:, :del_col]
        if normalization:
            heatmap_df1 = (
                heatmap_df1.subtract(
                    heatmap_df1.min(
                        axis=1),
                    axis=0)).div(
                heatmap_df1.max(
                    axis=1) -
                heatmap_df1.min(
                    axis=1) +
                0.00000000000000000000000000000000000000000001,
                axis=0)
            heatmap_df2 = (
                heatmap_df2.subtract(
                    heatmap_df2.min(
                        axis=1),
                    axis=0)).div(
                heatmap_df2.max(
                    axis=1) -
                heatmap_df2.min(
                    axis=1) +
                0.00000000000000000000000000000000000000000001,
                axis=0)
        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 24))
        sns.heatmap(
            heatmap_df1,
            ax=ax1,
            cmap=cmap,
            cbar_kws={
                'label': count_format})
        sns.heatmap(
            heatmap_df2,
            ax=ax2,
            cmap=cmap,
            cbar_kws={
                'label': count_format})
        ax1.set_xlabel('Length', fontsize=15)
        ax1.set_ylabel('Sample', fontsize=15)
        ax1.set_title(
            species.capitalize() +
            ' ' +
            rname +
            ' rsRNA Length distribution ' +
            type,
            fontsize=17)
        ax2.set_xlabel('Length', fontsize=15)
        ax2.set_ylabel('Sample', fontsize=15)
        ax2.set_title(
            species.capitalize() +
            ' Total Length distribution ' +
            type,
            fontsize=17)
        if save_directory:
            heatmap_df1.to_csv(
                save_directory +
                species +
                rname +
                '_rsRNA_length_' +
                flag +
                type +
                '1.csv')
            heatmap_df2.to_csv(
                save_directory +
                species +
                rname +
                '_rsRNA_length_' +
                flag +
                type +
                '2.csv')
    else:
        if del_col:
            heatmap_df1 = heatmap_df1.iloc[:, :del_col]
        if normalization:
            heatmap_df1 = (
                heatmap_df1.subtract(
                    heatmap_df1.min(
                        axis=1),
                    axis=0)).div(
                heatmap_df1.max(
                    axis=1) -
                heatmap_df1.min(
                    axis=1) +
                0.000000000000000000000000000000000000000001,
                axis=0)
        fig, ax = plt.subplots(figsize=(12, 12))
        sns.heatmap(
            heatmap_df1,
            cmap=cmap,
            ax=ax,
            cbar_kws={
                'label': count_format})
        ax.set_xlabel('Length', fontsize=15)
        ax.set_ylabel('Sample', fontsize=15)
        ax.set_title(
            species.capitalize() +
            ' ' +
            rname +
            ' rsRNA Length distribution ' +
            type,
            fontsize=17)
        if save_directory:
            heatmap_df1.to_csv(
                save_directory +
                species +
                rname +
                '_rsRNA_length_' +
                flag +
                type +
                '.csv')
    if save_fig:
        plt.savefig(
            save_fig +
            species +
            rname +
            '_rsRNA_length_' +
            flag +
            type +
            '.png',
            dpi=200,
            bbox_inches='tight')


def get_length_clustermap(
        heatmap_df,
        z_score=0,
        save_directory='',
        species='human',
        normalization='',
        row_cluster=False,
        col_cluster=False,
        cmap='viridis',
        yticklabels=True,
        flag='',
        rname='',
        type='',
        count_format='z_score',
        standard_scale=None,
        method='average',
        del_col=0):
    """
    get length distribution clustermap
    :param heatmap_df:dataframe to draw
    :param z_score: normalization method
    :param save_directory: path to save fig and csv
    :param species: human or mouse for label
    :param normalization: normalization with reads-min/(max-min)
    :param row_cluster: True or False
    :param col_cluster: True or False
    :param cmap: color map
    :param yticklabels: True or False
    :param flag: positive or negative for label
    :param rname:reference gene for label
    :param type:unique or redundant for label
    :param count_format:cbar label
    :param standard_scale:standardation method
    :param method:clustering method
    :param del_col:1 or 2 or other to del last columns
    :return:clustermap fig and save fig and csv
    """
    if del_col > 0:
        for i in range(del_col):
            heatmap_df = heatmap_df.drop(heatmap_df.columns[len(
                heatmap_df.columns) - 1], axis=1, inplace=False)
    if normalization:
        heatmap_df = (
            heatmap_df.subtract(
                heatmap_df.min(
                    axis=1),
                axis=0)).div(
            heatmap_df.max(
                axis=1) -
            heatmap_df.min(
                axis=1) +
            0.000000000000000000000000000000000000000000000001,
            axis=0)
    tissue = list(heatmap_df.index.get_level_values('Tissue'))
    if len(set(tissue)) < 13:
        color_palette = 'Paired'
    else:
        color_palette = 'tab20'
    color_pal = sns.color_palette(color_palette, len(set(tissue)))
    single_tissue = list(set(tissue))
    single_tissue.sort(key=tissue.index)
    tissue_color = dict(zip(single_tissue, color_pal))
    side_color_bar = pd.Series(
        heatmap_df.index.get_level_values('Tissue'),
        index=heatmap_df.index).map(tissue_color)
    hm = sns.clustermap(
        heatmap_df,
        row_cluster=row_cluster,
        cmap=cmap,
        row_colors=side_color_bar,
        figsize=(
            12,
            12),
        col_cluster=col_cluster,
        yticklabels=yticklabels,
        z_score=z_score,
        standard_scale=standard_scale,
        method=method)
    hm.ax_heatmap.set_title(
        species.capitalize() +
        ' ' +
        rname +
        ' rsRNA Length Distribution ' +
        type,
        fontsize=17,
        pad=90)
    hm.ax_heatmap.set_xlabel('Length', fontsize=15)
    hm.ax_heatmap.set_ylabel('Sample', fontsize=15)
    plt.xlabel(count_format)
    for tissue, color in tissue_color.items():
        hm.ax_col_dendrogram.bar(0, 0, color=color,
                                 label=tissue, linewidth=0)
    hm.ax_col_dendrogram.legend(
        loc='best', ncol=4, mode='expand', bbox_to_anchor=(-.08, .4, 1.15, 0.3))
    if save_directory:
        plt.savefig(
            save_directory +
            species +
            rname +
            flag +
            type +
            'lengthcluster.png',
            dpi=200,
            bbox_inches='tight')
        heatmap_df.to_csv(
            save_directory +
            species +
            rname +
            flag +
            type +
            'lengthcluster.csv',
            header=True,
            index=True)


def length_pool_heatmap(
        sam_list,
        flag='',
        type='',
        threshold=0,
        count_format='z_score',
        subplot='',
        single=True,
        save_directory='',
        species='human',
        z_score=0,
        row_cluster=False,
        clustermap='',
        del_col1=0,
        col_cluster=False,
        cmap='viridis',
        yticklabels=True,
        standard_scale=None,
        method='average',
        normalization='',
        del_col2='',
        save_fig='',
        three_index=True):
    """
    pool length distribution heatmap:five types
    first:cluster,second:heatmap single figure pool,third:heatmap single figure subplots,
    fourth:double heatmaps pool,fifth:double heatmaps subplots
    :param sam_list:path and SRR ID
    :param flag:0 16 or ''
    :param type:unique or redundant
    :param threshold:int,RPM or read counts value
    :param count_format:color bar label
    :param subplot:True or False
    :param single:True or False
    :param save_directory:path to store
    :param species:human or mouse
    :param z_score:0 or 1
    :param row_cluster:True or False
    :param clustermap:True or False
    :param del_col1:del dataframe columns 0 1 2
    :param col_cluster:True or False
    :param cmap:color map
    :param yticklabels:True or False
    :param standard_scale:0 or 1
    :param method:clustering method
    :param normalization:True or False
    :param del_col2:del dataframe columns -1 -2
    :param save_fig:function get_length_heatmap() get two parameters to store csv file and png,
    save directory is store csv and save fig is store png
    :param three_index:use three indexs or two indexs format,True or False
    :return:store csv and png and show figure
    """
    def length_heatmap_subplot_double(
            heatmap_df1, heatmap_df2, ax1, ax2, rname):
        """
        draw subplots two rows and n columns,first row is reference gene,second is all clean reads
        :param heatmap_df1:first row dataframe
        :param heatmap_df2:second row dataframe
        :param ax1:first dataframe subplot ax
        :param ax2:second dataframe subplot ax
        :param rname:reference gene name for title
        :return:draw two heatmaps subplots in assigned axes position
        """
        if del_col2:
            heatmap_df1 = heatmap_df1.iloc[:, :del_col2]
            heatmap_df2 = heatmap_df2.iloc[:, :del_col2]
        if normalization:
            heatmap_df1 = (
                heatmap_df1.subtract(
                    heatmap_df1.min(
                        axis=1),
                    axis=0)).div(
                heatmap_df1.max(
                    axis=1) -
                heatmap_df1.min(
                    axis=1) +
                0.000000000000000000000000000000000000000000001,
                axis=0)
            heatmap_df2 = (
                heatmap_df2.subtract(
                    heatmap_df2.min(
                        axis=1),
                    axis=0)).div(
                heatmap_df2.max(
                    axis=1) -
                heatmap_df2.min(
                    axis=1) +
                0.000000000000000000000000000000000000000000001,
                axis=0)
        sns.heatmap(
            heatmap_df1,
            ax=ax1,
            cmap=cmap,
            cbar_kws={
                'label': count_format})
        sns.heatmap(
            heatmap_df2,
            ax=ax2,
            cmap=cmap,
            cbar_kws={
                'label': count_format})
        ax1.set_xlabel('Length', fontsize=15)
        ax1.set_ylabel('Sample', fontsize=15)
        ax1.set_title(
            species.capitalize() +
            ' ' +
            rname +
            ' rsRNA Length distribution ' +
            type,
            fontsize=17)
        ax2.set_xlabel('Length', fontsize=15)
        ax2.set_ylabel('Sample', fontsize=15)
        ax2.set_title(
            species.capitalize() +
            ' rsRNA total reads Length distribution ' +
            type,
            fontsize=17)

    def length_heatmap_subplot_single(heatmap_df, ax, rname):
        """
        get subplots,single mean single column to draw subplots
        :param heatmap_df:subplot dataframe
        :param ax:axes position to draw a subplot
        :param rname:reference gene name for title
        :return:draw heatmap subplot in assigned ax position
        """
        if del_col2:
            heatmap_df = heatmap_df.iloc[:, :del_col2]
        if normalization:
            heatmap_df = (
                heatmap_df.subtract(
                    heatmap_df.min(
                        axis=1),
                    axis=0)).div(
                heatmap_df.max(
                    axis=1) -
                heatmap_df.min(
                    axis=1) +
                0.00000000000000000000000000000000000000000001,
                axis=0)
        sns.heatmap(
            heatmap_df,
            cmap=cmap,
            ax=ax,
            cbar_kws={
                'label': count_format})
        ax.set_xlabel('Length', fontsize=15)
        ax.set_ylabel('Sample', fontsize=15)
        ax.set_title(
            species.capitalize() +
            ' ' +
            rname +
            ' rsRNA Length distribution ' +
            type,
            fontsize=17)

    def length_heatmap():
        """
        replication code to get length distribution dic
        :return:a dic key is rname, value is length distribution array
        """
        nonlocal i
        for sam in sam_lists:
            sam = sam.strip()
            if sam[1:3] == 'RR' or sam[1:3] == 'RX':
                rRNA_sam_path = path + sam + '.sam'
                samclass = SamFile(rRNA_sam_path, species)
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                i += 1
                j += 1
                for rname in referencename.keys():
                    if i == 1:
                        if rname == 'all':
                            if type == 'redundant':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_redundant_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = total_length_array
                            elif type == 'unique':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_unique_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = total_length_array
                        else:
                            length_array = samclass.length_distribution(
                                samclass.get_length_distribution(
                                    flag=flag, rname=rname, type=type, threshold=threshold))
                            length_array = np.array(
                                list(length_array.values()))
                            referencename[rname] = length_array
                    else:
                        if rname == 'all':
                            if type == 'redundant':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_redundant_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = np.vstack(
                                    (referencename[rname], total_length_array))
                            elif type == 'unique':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_unique_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = np.vstack(
                                    (referencename[rname], total_length_array))
                        else:
                            length_array = samclass.length_distribution(
                                samclass.get_length_distribution(
                                    flag=flag, rname=rname, type=type, threshold=threshold))
                            length_array = np.array(
                                list(length_array.values()))
                            referencename[rname] = np.vstack(
                                (referencename[rname], length_array))
            else:
                path = sam
                j = 1
        return referencename
    if species == 'mouse':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('4.5S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('45S', 'rRNA_45S'), ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'),
                                                 ('', ''), ('all', 'all_clean'), ('5S', 'rRNA_5S')])
    elif species == 'human':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('12S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('5S', 'rRNA_5S'),
                                                 ('', ''), ('all', 'all_clean'), ('45S', 'rRNA_45S')])
    with open(sam_list) as sl:
        sam_lists = sl.readlines()
    sample_info = []
    tissue_number = []
    i = 0
    columns = [str(i) for i in range(16, 46)]
    columns.append('>45')
    if single and clustermap:
        referencename = length_heatmap()
        if three_index:
            index = sample_info
            indexname = ['Tissue', 'SRP', 'SRR']
        else:
            index = tissue_number
            indexname = ['Tissue', 'Number']
        for rname, df in referencename.items():
            heatmap_df = pd.DataFrame(df, index=tuple(index),
                                      columns=columns)
            heatmap_df.index.names = indexname
            heatmap_df.columns.name = 'Length'
            get_length_clustermap(
                heatmap_df=heatmap_df,
                z_score=z_score,
                save_directory=save_directory,
                species=species,
                row_cluster=row_cluster,
                col_cluster=col_cluster,
                cmap=cmap,
                yticklabels=yticklabels,
                flag=flag,
                rname=rname,
                type=type,
                count_format=count_format,
                standard_scale=standard_scale,
                method=method,
                del_col=del_col1,
                normalization=normalization)
    elif single and not clustermap:
        if not subplot:
            referencename = length_heatmap()
            if three_index:
                index = sample_info
                indexname = ['Tissue', 'SRP', 'SRR']
            else:
                index = tissue_number
                indexname = ['Tissue', 'Number']
            for rname, df in referencename.items():
                heatmap_df = pd.DataFrame(df, index=tuple(index),
                                          columns=columns)
                heatmap_df.index.names = indexname
                heatmap_df.columns.name = 'Length'
                get_length_heatmap(
                    heatmap_df,
                    normalization=normalization,
                    cmap=cmap,
                    species=species,
                    save_fig=save_fig,
                    flag=flag,
                    rname=rname,
                    type=type,
                    save_directory=save_directory,
                    del_col=del_col2,
                    count_format=count_format)
        elif subplot:
            fig, axs = plt.subplots(nrows=14, figsize=(12, 156))
            referencename = length_heatmap()
            if three_index:
                index = sample_info
                indexname = ['Tissue', 'SRP', 'SRR']
            else:
                index = tissue_number
                indexname = ['Tissue', 'Number']
            for rname_df, ax in zip(referencename.items(), axs.flatten()):
                heatmap_df = pd.DataFrame(rname_df[1], index=tuple(index),
                                          columns=columns)
                heatmap_df.index.names = indexname
                heatmap_df.columns.name = 'Length'
                length_heatmap_subplot_single(heatmap_df, ax, rname_df[0])
            if save_directory:
                plt.savefig(
                    save_directory +
                    species +
                    flag +
                    type +
                    '_length.pdf',
                    dpi=150,
                    bbox_inches='tight')
    elif not single:
        if not subplot:
            referencename = length_heatmap()
            if three_index:
                index = sample_info
                indexname = ['Tissue', 'SRP', 'SRR']
            else:
                index = tissue_number
                indexname = ['Tissue', 'Number']
            heatmap_df2 = pd.DataFrame(
                referencename['all'],
                index=tuple(index),
                columns=columns)
            heatmap_df2.index.names = indexname
            heatmap_df2.columns.name = 'Length'
            for rname, df in referencename.items():
                if rname == 'all':
                    continue
                heatmap_df = pd.DataFrame(df, index=tuple(index),
                                          columns=columns)
                heatmap_df.index.names = indexname
                heatmap_df.columns.name = 'Length'
                get_length_heatmap(
                    heatmap_df,
                    heatmap_df2,
                    normalization=normalization,
                    cmap=cmap,
                    species=species,
                    flag=flag,
                    rname=rname,
                    type=type,
                    save_fig=save_fig,
                    save_directory=save_directory,
                    del_col=del_col2,
                    count_format=count_format)
        elif subplot:
            fig, axs = plt.subplots(nrows=2, ncols=13, figsize=(260, 24))
            referencename = length_heatmap()
            if three_index:
                index = sample_info
                indexname = ['Tissue', 'SRP', 'SRR']
            else:
                index = tissue_number
                indexname = ['Tissue', 'Number']
            heatmap_df2 = pd.DataFrame(
                referencename['all'],
                index=tuple(index),
                columns=columns)
            heatmap_df2.index.names = indexname
            heatmap_df2.columns.name = 'Length'
            k = 0
            for rname, df in referencename.items():
                if rname == 'all':
                    continue
                heatmap_df = pd.DataFrame(df, index=tuple(index),
                                          columns=columns)
                heatmap_df.index.names = indexname
                heatmap_df.columns.name = 'Length'
                length_heatmap_subplot_double(
                    heatmap_df, heatmap_df2, axs[0][k], axs[1][k], rname)
                k += 1
            if save_directory:
                plt.savefig(
                    save_directory +
                    species +
                    flag +
                    type +
                    '_length.pdf',
                    dpi=150,
                    bbox_inches='tight')


def tSNE_scatter_plot_location_single(
        tsne_df,
        n_components=2,
        perplexity=30.0,
        early_exaggeration=12.0,
        learning_rate=200.0,
        n_iter=1000,
        n_iter_without_progress=300,
        min_grad_norm=1e-07,
        metric='euclidean',
        init='random',
        verbose=0,
        random_state=None,
        method='barnes_hut',
        angle=0.5,
        save_directory='',
        normalization='',
        flag='',
        species='',
        type='',
        rname=''):
    """
    get location distribution tSNE Reduced dimensional clustering scatter plot
    :param tsne_df:dataframe for tSNE scatter plot
    :param n_components:int, optional (default: 2) Dimension of the embedded space.
    :param perplexity:perplexity : float, optional (default: 30)
    The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms.
    Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50.
    Different values can result in significanlty different results.
    :param early_exaggeration:float, optional (default: 12.0)
    Controls how tight natural clusters in the original space are in the embedded space and
    how much space will be between them. For larger values, the space between natural clusters
    will be larger in the embedded space. Again, the choice of this parameter is not very critical.
    If the cost function increases during initial optimization,
    the early exaggeration factor or the learning rate might be too high.
    :param learning_rate:float, optional (default: 200.0)
    The learning rate for t-SNE is usually in the range [10.0, 1000.0].
    If the learning rate is too high, the data may look like a ‘ball’ with any point approximately equidistant
    from its nearest neighbours. If the learning rate is too low,
    most points may look compressed in a dense cloud with few outliers.
    If the cost function gets stuck in a bad local minimum increasing the learning rate may help.
    :param n_iter:n_iter : int, optional (default: 1000)
    Maximum number of iterations for the optimization. Should be at least 250.
    :param n_iter_without_progress:n_iter_without_progress : int, optional (default: 300)
    Maximum number of iterations without progress before we abort the optimization,
    used after 250 initial iterations with early exaggeration.
    Note that progress is only checked every 50 iterations so this value is rounded to the next multiple of 50.
    :param min_grad_norm:min_grad_norm : float, optional (default: 1e-7)
    If the gradient norm is below this threshold, the optimization will be stopped.
    :param metric:metric : string or callable, optional
    The metric to use when calculating distance between instances in a feature array.
     If metric is a string, it must be one of the options allowed by scipy.
     spatial.distance.pdist for its metric parameter, or a metric listed in pairwise.
     PAIRWISE_DISTANCE_FUNCTIONS. If metric is “precomputed”, X is assumed to be a distance matrix.
     Alternatively, if metric is a callable function, it is called on each pair of instances (rows)
     and the resulting value recorded. The callable should take two arrays from X as input
     and return a value indicating the distance between them. The default is “euclidean”
     which is interpreted as squared euclidean distance.
    :param init:string or numpy array, optional (default: “random”)
    Initialization of embedding. Possible options are ‘random’,
    ‘pca’, and a numpy array of shape (n_samples, n_components).
    PCA initialization cannot be used with precomputed distances and is usually more globally stable
    than random initialization.
    :param verbose: int, optional (default: 0) Verbosity level.
    :param random_state:int, RandomState instance or None, optional (default: None)
    If int, random_state is the seed used by the random number generator;
    If RandomState instance, random_state is the random number generator;
    If None, the random number generator is the RandomState instance used by np.random.
    Note that different initializations might result in different local minima of the cost function.
    :param method:string (default: ‘barnes_hut’)
    By default the gradient calculation algorithm uses Barnes-Hut approximation running in O(NlogN) time.
    method=’exact’ will run on the slower, but exact, algorithm in O(N^2) time.
    The exact algorithm should be used when nearest-neighbor errors need to be better than 3%.
    However, the exact method cannot scale to millions of examples.
    :param angle: float (default: 0.5)
    Only used if method=’barnes_hut’ This is the trade-off between speed and accuracy for Barnes-Hut T-SNE.
    ‘angle’ is the angular size (referred to as theta in [3]) of a distant node as measured from a point.
     If this size is below ‘angle’ then it is used as a summary node of all points contained within it.
     This method is not very sensitive to changes in this parameter in the range of 0.2 - 0.8.
     Angle less than 0.2 has quickly increasing computation time and angle greater 0.8 has quickly increasing error.
    :param save_directory:PATH to store figure
    :param normalization:normalization
    :param flag:for save fig name
    :param species:for title and save fig name
    :param type:for title and save fig name
    :return:show and save tSNE scatter figure
    """
    if normalization:
        tsne_df = (
            tsne_df.subtract(
                tsne_df.min(
                    axis=1),
                axis=0)).div(
            tsne_df.max(
                axis=1) -
            tsne_df.min(
                axis=1) +
            0.00000000000000000000000000000000000000000000000000001,
            axis=0)
    tissues = tsne_df.index.get_level_values('Tissue').unique()
    tissues_series = tsne_df.index.get_level_values(
        'Tissue').value_counts(sort=False)
    sort_tissue_dic = collections.OrderedDict()
    for tissue in tissues:
        sort_tissue_dic[tissue] = tissues_series[tissue]
    tsne_array = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        n_iter=n_iter,
        n_iter_without_progress=n_iter_without_progress,
        min_grad_norm=min_grad_norm,
        metric=metric,
        init=init,
        verbose=verbose,
        random_state=random_state,
        method=method,
        angle=angle).fit_transform(tsne_df)
    n = 0
    k = 0
    palette = np.array(sns.color_palette("tab20"))
    plt.figure(figsize=(9, 9))
    for key, value in sort_tissue_dic.items():
        if n == 0:
            plt.scatter(tsne_array[:value, 0], tsne_array[:value, 1], c=palette[[
                        k]], s=22, alpha=0.8, label=key)
            n += value
            k += 1
        else:
            plt.scatter(tsne_array[n:n +
                                   value, 0], tsne_array[n:n +
                                                         value, 1], c=palette[[k] *
                                                                              value], s=22, alpha=0.8, label=key)
            n += value
            k += 1
    plt.title(
        species.capitalize() +
        ' ' +
        rname +
        ' rsRNA location distribution tSNE cluster ' +
        type,
        fontsize=16)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 11}
    plt.legend(loc=2, bbox_to_anchor=(1.05, 0.85), borderaxespad=0, prop=font)
    plt.tight_layout()
    if save_directory:
        plt.savefig(
            save_directory +
            species +
            rname +
            flag +
            type +
            '_location_tSNE.png',
            dpi=200,
            bbox_inches='tight')
        tsne_df.to_csv(
            save_directory +
            species +
            rname +
            flag +
            type +
            '_location_tSNE.csv')


def tSNE_scatter_plot_location_pool(
    subplot='',
    normalization='',
    sam_list='',
    flag='',
    type='redundant',
    threshold=0,
    count_format='RPM',
    three_index=True,
    genome_index_path='',
    species='human',
    n_components=2,
    perplexity=30.0,
    early_exaggeration=12.0,
    learning_rate=200.0,
    n_iter=1000,
    n_iter_without_progress=300,
    min_grad_norm=1e-07,
    metric='euclidean',
    init='random',
    verbose=0,
    random_state=None,
    method='barnes_hut',
    angle=0.5,
    save_directory='',
):
    """
    location tSNE scatter pool，two types:subplots or single tSNE pool
    :param subplot:draw with subplots
    :param normalization: True or False
    :param sam_list: path and SRR ID
    :param flag: '0' '16' or ''
    :param type: unique or redundant
    :param threshold: RPM or COUNT threshold
    :param count_format: RPM or COUNT
    :param three_index: True or False
    :param genome_index_path: bowtie index path of genome if use RPM to count
    :param species: human or mouse
    :param save_directory: path to store figure and csv
    :return: show figure, save csv and png
    """
    def tSNE_dataframe_dic():
        """
        get tSNE array dictionary
        :return: tSNE array dictionary，key is rname ，value is rname location distribution array
        """
        i = 0
        for sam in sam_lists:
            sam = sam.strip()
            if sam[1:3] == 'RR' or sam[1:3] == 'RX':
                rRNA_sam_path = path + sam + '.sam'
                samclass = SamFile(rRNA_sam_path, species)
                if count_format == 'RPM':
                    fastq_path = path + sam + '.fastq'
                    samclass.rpm_normalization(fastq_path, genome_index_path)
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                i += 1
                j += 1
                for rname in referencename.keys():
                    if i == 1:
                        location_numpy = samclass.quick_location_dic(
                            flag=flag, rname=rname, type=type, threshold=threshold, count_format=count_format)
                        referencename[rname] = np.array(
                            list(location_numpy.values()))

                    else:
                        location_series = samclass.quick_location_dic(
                            flag=flag, rname=rname, type=type, threshold=threshold, count_format=count_format)
                        location_series = np.array(
                            list(location_series.values()))
                        referencename[rname] = np.vstack(
                            (referencename[rname], location_series))
            else:
                path = sam
                j = 1
        return referencename

    def tSNE_subplot(tsne_df, ax, rname):
        """
        draw tSNE subplot
        :param tsne_df: subplot dataframe
        :param ax: subplot ax position
        :param rname: reference gene name for title
        :return: draw dataframe tSNE subplot on assigned ax
        """
        tissues = tsne_df.index.get_level_values('Tissue').unique()
        tissues_series = tsne_df.index.get_level_values(
            'Tissue').value_counts(sort=False)
        sort_tissue_dic = collections.OrderedDict()
        for tissue in tissues:
            sort_tissue_dic[tissue] = tissues_series[tissue]
        tsne_array = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            early_exaggeration=early_exaggeration,
            learning_rate=learning_rate,
            n_iter=n_iter,
            n_iter_without_progress=n_iter_without_progress,
            min_grad_norm=min_grad_norm,
            metric=metric,
            init=init,
            verbose=verbose,
            random_state=random_state,
            method=method,
            angle=angle).fit_transform(tsne_df)
        n = 0
        k = 0
        palette = np.array(sns.color_palette("tab20"))
        for key, value in sort_tissue_dic.items():
            if n == 0:
                ax.scatter(tsne_array[:value, 0], tsne_array[:value, 1], c=palette[[
                           k]], s=22, alpha=0.8, label=key)
                n += value
                k += 1
            else:
                ax.scatter(tsne_array[n:n +
                                      value, 0], tsne_array[n:n +
                                                            value, 1], c=palette[[k] *
                                                                                 value], s=22, alpha=0.8, label=key)
                n += value
                k += 1
        ax.set_title(
            species.capitalize() +
            ' ' +
            rname +
            ' rsRNA location distribution tSNE cluster ' +
            type,
            fontsize=16)
        font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 11}
        ax.legend(
            loc=2,
            bbox_to_anchor=(
                1.05,
                0.85),
            borderaxespad=0,
            prop=font)
    if species == 'mouse':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('4.5S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('5S', 'rRNA_5S'),
                                                 ('45S', 'rRNA_45S')])
        gene_no_replication = ['45S', '5S', '16S', '4.5S']
    elif species == 'human':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('12S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('5S', 'rRNA_5S'),
                                                 ('45S', 'rRNA_45S')])
        gene_no_replication = ['45S', '5S', '16S', '12S']
    with open(sam_list) as sl:
        sam_lists = sl.readlines()
    sample_info = []
    tissue_number = []
    referencename = tSNE_dataframe_dic()
    if three_index:
        index = sample_info
        indexname = ['Tissue', 'SRP', 'SRR']
    else:
        index = tissue_number
        indexname = ['Tissue', 'Number']
    dfs = []
    if subplot:
        fig, axs = plt.subplots(nrows=13, figsize=(11, 143))
        fig.suptitle(
            species +
            ' rRNA Location tSNE ' +
            type,
            fontsize=16,
            y=0.895)
        for rname_df, ax in zip(referencename.items(), axs.flatten()):
            tsne_df = pd.DataFrame(
                rname_df[1], index=tuple(index), columns=[
                    i for i in range(
                        1, len(
                            rname_df[1][0]) + 1)])
            if normalization:
                tsne_df = (
                    tsne_df.subtract(
                        tsne_df.min(
                            axis=1),
                        axis=0)).div(
                    tsne_df.max(
                        axis=1) -
                    tsne_df.min(
                        axis=1) +
                    0.00000000000000000000000000000000000001,
                    axis=0)
            tsne_df.index.names = indexname
            tsne_df.columns.name = 'Location'
            if rname_df[0] in gene_no_replication:
                dfs.append(tsne_df)
            tSNE_subplot(tsne_df, ax, rname_df[0])
            if save_directory:
                tsne_df.to_csv(
                    save_directory +
                    species +
                    rname_df[0] +
                    flag +
                    type +
                    'location_tSNE.csv')
        merge_dataframe = pd.concat(dfs, axis=1, sort=False, join='outer')
        tSNE_subplot(merge_dataframe, axs.flatten()[-1], 'all')
        if save_directory:
            plt.savefig(
                save_directory +
                species +
                flag +
                type +
                'tSNEsubplots.pdf',
                dpi=200,
                bbox_inches='tight')
            merge_dataframe.to_csv(
                save_directory +
                species +
                'allrname' +
                flag +
                type +
                'location_tSNE.csv')
    else:
        for rname, df in referencename.items():
            tsne_df = pd.DataFrame(
                df, index=tuple(index), columns=[
                    i for i in range(
                        1, len(
                            df[0]) + 1)])
            tsne_df.index.names = indexname
            tsne_df.columns.name = 'Location'
            if rname in gene_no_replication:
                dfs.append(tsne_df)
            tSNE_scatter_plot_location_single(
                tsne_df,
                n_components=n_components,
                perplexity=perplexity,
                early_exaggeration=early_exaggeration,
                save_directory=save_directory,
                learning_rate=learning_rate,
                n_iter=n_iter,
                normalization=normalization,
                n_iter_without_progress=n_iter_without_progress,
                type=type,
                min_grad_norm=min_grad_norm,
                metric=metric,
                init=init,
                verbose=verbose,
                random_state=random_state,
                method=method,
                angle=angle,
                flag=flag,
                species=species,
                rname=rname)
        merge_dataframe = pd.concat(dfs, axis=1, sort=False, join='outer')
        tSNE_scatter_plot_location_single(
            merge_dataframe,
            n_components=n_components,
            perplexity=perplexity,
            early_exaggeration=early_exaggeration,
            save_directory=save_directory,
            learning_rate=learning_rate,
            n_iter=n_iter,
            normalization=normalization,
            n_iter_without_progress=n_iter_without_progress,
            type=type,
            min_grad_norm=min_grad_norm,
            metric=metric,
            init=init,
            verbose=verbose,
            random_state=random_state,
            method=method,
            angle=angle,
            flag=flag,
            species=species,
            rname='all')


def tSNE_scatter_plot_length_single(
        tsne_df,
        del_col='',
        normalization='',
        n_components=2,
        perplexity=30.0,
        early_exaggeration=12.0,
        learning_rate=200.0,
        n_iter=1000,
        n_iter_without_progress=300,
        rname='',
        flag='',
        type='',
        species='human',
        min_grad_norm=1e-07,
        metric='euclidean',
        init='random',
        verbose=0,
        random_state=None,
        method='barnes_hut',
        angle=0.5,
        save_directory=''):
    """
    get length distribution tSNE Reduced dimensional clustering scatter plot
    :param tsne_df: dataframe to draw tSNE scatter figure
    :param del_col: del dataframe columns -1 -2 or other
    :param normalization: True or False
    :param rname: for title and saved file name
    :param flag: for saved file name
    :param type: for title and saved file name
    :param species: for title and saved file name
    :param save_directory: path to store csv and png
    :return: show length tSNE scatter figure,store figure and csv
    """
    if del_col:
        tsne_df = tsne_df.iloc[:, :del_col]
    if normalization:
        tsne_df = (
            tsne_df.subtract(
                tsne_df.min(
                    axis=1),
                axis=0)).div(
            tsne_df.max(
                axis=1) -
            tsne_df.min(
                axis=1) +
            0.00000000000000000000000000000000000000000000000000001,
            axis=0)
    tissues = tsne_df.index.get_level_values('Tissue').unique()
    tissues_series = tsne_df.index.get_level_values(
        'Tissue').value_counts(sort=False)
    sort_tissue_dic = collections.OrderedDict()
    for tissue in tissues:
        sort_tissue_dic[tissue] = tissues_series[tissue]
    tsne_array = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        n_iter=n_iter,
        n_iter_without_progress=n_iter_without_progress,
        min_grad_norm=min_grad_norm,
        metric=metric,
        init=init,
        verbose=verbose,
        random_state=random_state,
        method=method,
        angle=angle).fit_transform(tsne_df)
    n = 0
    k = 0
    palette = np.array(sns.color_palette("tab20"))
    plt.figure(figsize=(9, 9))
    for key, value in sort_tissue_dic.items():
        if n == 0:
            plt.scatter(tsne_array[:value, 0], tsne_array[:value, 1], c=palette[[
                        k]], s=22, alpha=0.8, label=key)
            n += value
            k += 1
        else:
            plt.scatter(tsne_array[n:n +
                                   value, 0], tsne_array[n:n +
                                                         value, 1], c=palette[[k] *
                                                                              value], s=22, alpha=0.8, label=key)
            n += value
            k += 1
    plt.title(
        species.capitalize() +
        ' ' +
        rname +
        ' Length distribution tSNE ' +
        type,
        fontsize=16)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 11}
    plt.legend(loc=2, bbox_to_anchor=(1.05, 0.85), borderaxespad=0, prop=font)
    if save_directory:
        plt.savefig(
            save_directory +
            species +
            rname +
            flag +
            type +
            '_length_tSNE.png',
            dpi=200,
            bbox_inches='tight')
        tsne_df.to_csv(
            save_directory +
            species +
            rname +
            flag +
            type +
            '_length_tSNE.csv')


def tSNE_scatter_plot_length_pool(
        sam_list='',
        del_col='',
        normalization='',
        three_index=True,
        subplot='',
        species='human',
        flag='',
        type='',
        threshold=0,
        save_directory='',
        n_components=2,
        perplexity=30.0,
        early_exaggeration=12.0,
        learning_rate=200.0,
        n_iter=1000,
        n_iter_without_progress=300,
        min_grad_norm=1e-07,
        metric='euclidean',
        init='random',
        verbose=0,
        random_state=None,
        method='barnes_hut',
        angle=0.5):
    """
    length distribution tSNE scatter pool，two types:subplots or single tSNE pool
    :param sam_list:path and SRR
    :param del_col:del dataframe columns,-1 -2 or other
    :param normalization:True or False
    :param three_index:if True use three indexs,otherwise use two indexs
    :param subplot:True or False
    :param species:human or mouse
    :param flag:'0' '16'  or ''
    :param type:unique or redundant
    :param threshold:RPM or COUNT threshold
    :param save_directory:path to store png and csv
    :return:show length tSNE scatter figure , save csv and png
    """
    def tsne_length_subplot(tsne_df, ax, rname):
        """
        draw tSNE length distribution subplot
        :param tsne_df: subplot dataframe
        :param ax: subplot ax position
        :param rname: reference gene name for title
        :return: draw subplot on the assigned ax
        """
        if del_col:
            tsne_df = tsne_df.iloc[:, :del_col]
        if normalization:
            tsne_df = (
                tsne_df.subtract(
                    tsne_df.min(
                        axis=1),
                    axis=0)).div(
                tsne_df.max(
                    axis=1) -
                tsne_df.min(
                    axis=1) +
                0.0000000000000000000000000000000000000000000000001,
                axis=0)
        tissues = tsne_df.index.get_level_values('Tissue').unique()
        tissues_series = tsne_df.index.get_level_values(
            'Tissue').value_counts(sort=False)
        sort_tissue_dic = collections.OrderedDict()
        for tissue in tissues:
            sort_tissue_dic[tissue] = tissues_series[tissue]
        tsne_array = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            early_exaggeration=early_exaggeration,
            learning_rate=learning_rate,
            n_iter=n_iter,
            n_iter_without_progress=n_iter_without_progress,
            min_grad_norm=min_grad_norm,
            metric=metric,
            init=init,
            verbose=verbose,
            random_state=random_state,
            method=method,
            angle=angle).fit_transform(tsne_df)
        n = 0
        k = 0
        palette = np.array(sns.color_palette("tab20"))
        for key, value in sort_tissue_dic.items():
            if n == 0:
                ax.scatter(tsne_array[:value, 0], tsne_array[:value, 1], c=palette[[
                           k]], s=22, alpha=0.8, label=key)
                n += value
                k += 1
            else:
                ax.scatter(tsne_array[n:n +
                                      value, 0], tsne_array[n:n +
                                                            value, 1], c=palette[[k] *
                                                                                 value], s=22, alpha=0.8, label=key)
                n += value
                k += 1
        ax.set_title(
            species.capitalize() +
            ' ' +
            rname +
            ' rsRNA length distribution tSNE ' +
            type,
            fontsize=16)
        font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 11}
        ax.legend(
            loc=2,
            bbox_to_anchor=(
                1.05,
                0.85),
            borderaxespad=0,
            prop=font)

    def length_heatmap():
        """
        replication code to get length distribution dictionary
        :return:a dic key is rname, value is length distribution array
        """
        nonlocal i
        for sam in sam_lists:
            sam = sam.strip()
            if sam[1:3] == 'RR' or sam[1:3] == 'RX':
                rRNA_sam_path = path + sam + '.sam'
                samclass = SamFile(rRNA_sam_path, species)
                tissue = rRNA_sam_path.split('/')[-2]
                srr = rRNA_sam_path.split('/')[-1][:6]
                srr_number = rRNA_sam_path.split('/')[-1][:-4]
                sample_name = tissue, srr, srr_number
                sample_info.append(sample_name)
                tissue_number.append((tissue, j))
                i += 1
                j += 1
                for rname in referencename.keys():
                    if i == 1:
                        if rname == 'all':
                            if type == 'redundant':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_redundant_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = total_length_array
                            elif type == 'unique':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_unique_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = total_length_array
                        else:
                            length_array = samclass.length_distribution(
                                samclass.get_length_distribution(
                                    flag=flag, rname=rname, type=type, threshold=threshold))
                            length_array = np.array(
                                list(length_array.values()))
                            referencename[rname] = length_array
                    else:
                        if rname == 'all':
                            if type == 'redundant':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_redundant_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = np.vstack(
                                    (referencename[rname], total_length_array))
                            elif type == 'unique':
                                total_length_array = samclass.length_distribution(
                                    samclass.total_unique_length_distribution)
                                total_length_array = np.array(
                                    list(total_length_array.values()))
                                referencename[rname] = np.vstack(
                                    (referencename[rname], total_length_array))
                        else:
                            length_array = samclass.length_distribution(
                                samclass.get_length_distribution(
                                    flag=flag, rname=rname, type=type, threshold=threshold))
                            length_array = np.array(
                                list(length_array.values()))
                            referencename[rname] = np.vstack(
                                (referencename[rname], length_array))
            else:
                path = sam
                j = 1
        return referencename
    if species == 'mouse':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('4.5S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('', ''), ('all', 'all_clean'),
                                                 ('45S', 'rRNA_45S')])
    elif species == 'human':
        referencename = collections.OrderedDict([('28S', 'rRNA_28S'), ('18S', 'rRNA_18S'), ('5.8S', 'rRNA_5_8S'),
                                                 ('45S_across_region', 'junction_45S'), ('12S', 'rRNA_12S'),
                                                 ('16S', 'rRNA_16S'), ('5ETS', "ETS5"), ('3ETS', "ETS3"),
                                                 ('ITS1', 'ITS_1'), ('ITS2', 'ITS_2'), ('5S', 'rRNA_5S'),
                                                 ('', ''), ('all', 'all_clean'), ('45S', 'rRNA_45S')])
    with open(sam_list) as sl:
        sam_lists = sl.readlines()
    sample_info = []
    tissue_number = []
    i = 0
    columns = [str(i) for i in range(16, 46)]
    columns.append('>45')
    if subplot:
        fig, axs = plt.subplots(nrows=14, figsize=(9, 126))
        fig.suptitle(
            species.capitalize() +
            ' rsRNA Length distribution tSNE ' +
            type,
            fontsize=16,
            y=0.895)
        referencename = length_heatmap()
        if three_index:
            index = sample_info
            indexname = ['Tissue', 'SRP', 'SRR']
        else:
            index = tissue_number
            indexname = ['Tissue', 'Number']
        for rname_df, ax in zip(referencename.items(), axs.flatten()):
            tsne_df = pd.DataFrame(rname_df[1], index=tuple(index),
                                   columns=columns)
            tsne_df.index.names = indexname
            tsne_df.columns.name = 'Length'
            if save_directory:
                tsne_df.to_csv(
                    save_directory +
                    species +
                    rname_df[0] +
                    flag +
                    type +
                    'length_tSNE.csv')
            tsne_length_subplot(tsne_df, ax, rname_df[0])
        if save_directory:
            plt.savefig(
                save_directory +
                species +
                'allrname' +
                flag +
                type +
                'length_tSNE.pdf',
                dpi=200,
                bbox_inches='tight')
    else:
        referencename = length_heatmap()
        if three_index:
            index = sample_info
            indexname = ['Tissue', 'SRP', 'SRR']
        else:
            index = tissue_number
            indexname = ['Tissue', 'Number']
        for rname, df in referencename.items():
            tsne_df = pd.DataFrame(df, index=tuple(index),
                                   columns=columns)
            tsne_df.index.names = indexname
            tsne_df.columns.name = 'Length'
            tSNE_scatter_plot_length_single(
                tsne_df,
                del_col=del_col,
                normalization=normalization,
                n_components=n_components,
                perplexity=perplexity,
                rname=rname,
                flag=flag,
                early_exaggeration=early_exaggeration,
                learning_rate=learning_rate,
                n_iter=n_iter,
                n_iter_without_progress=n_iter_without_progress,
                type=type,
                min_grad_norm=min_grad_norm,
                metric=metric,
                init=init,
                verbose=verbose,
                random_state=random_state,
                method=method,
                angle=angle,
                save_directory=save_directory,
                species=species)
