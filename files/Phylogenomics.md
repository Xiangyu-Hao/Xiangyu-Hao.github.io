# A workflow for phylogenomics

This pipeline was integrated from Xiangyu HAO and Yulu ZOU on Sep 1, 2024.

此指南适用于基于二代高通量测序的低覆盖度基因组、三代测序高质量基因组或二代三代基因组数据相结合的系统发生基因组学分析。
用户可自定义设置分析流程

# NCBI下载基因组

1.

2.

3.

# NCBI下载SRA数据

## (1) sra toolkit下载SRA数据

### Installation

```bash
wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz
tar zxvf sratoolkit.3.1.1-ubuntu64.tar.gz
vim ~/.bashrc
export PATH=$PATH:/path_to_sratoolkit/sratoolkit.3.1.1-ubuntu64/bin
source ~/.bashrc
```

### Usage

```bash
fastq-dump.3.1.0 <accession_number> --split-files
```

aspera下载sra数据https://www.jianshu.com/p/fed19a8821eb
conda安装：conda install -c hcc aspera-cli -y 验证是否安装成功：ascp -h
密钥查找which ascp 输出内容将把bin及bin后面的内容换成etc/asperaweb_id_dsa.openssh
/home/zyl/.aspera/connect/etc/asperaweb_id_dsa.openssh
ENA搜索SRA号，read files勾选fastq_aspera,(click link to copy URL)
ascp -i /home/zyl/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 200M -P 33001 -QT -k 2 era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR266/068/SRR26641668/SRR26641668_2.fastq.gz .[下载目标路径]



# 二代测序原始数据质控

## (1) fastp


**Installation:** `conda install fastp`

| Parameter | Note                                                         |
| :-------: | ------------------------------------------------------------ |
|    -i     | read1输入文件名                                              |
|    -I     | read2输入文件名                                              |
|    -o     | read1输出文件名                                              |
|    -O     | read2输出文件名                                              |
|    -A     | 默认启用接头修剪。如果指定了此选项，则禁用接头修剪。         |
|    -a     | 对于双端数据，通过read1和read2之间的overlap可判断接头，因此通常可以不指定接头序列，但仍然可以人为指定read1和read2的接头序列。当fastp无法找到overlap时(比如低质量碱基)，将根据指定的接头序列进行接头去去除。 |
|    -Q     | 质量过滤默认启用。如果指定了此选项，则禁用质量过滤。         |
|    -q     | 质量值低于此值认为是低质量碱基，默认15，表示质量为Q15。      |
|    -L     | 长度过滤默认启用。如果指定了此选项，则禁用长度过滤。         |
|    -l     | 短于length_required的读取将被丢弃，默认为15。                |
|    -w     | 工作线程数，默认是3。                                        |

```bash
Command:
fastp -i <raw_data_1.fastq> -I <raw_data_2.fastq> -o <clean_data_1.fastq> -O <clean_data_2.fastq>
```

```bash
Example:
fastp -i Npla_1.fq.gz -I Npla_2.fq.gz -o Npla_1_clean.fq.gz -O Npla_2_clean.fq.gz
```



# 基因组组装

## (1) SPAdes

### Installation

```bash
wget -c https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
tar -xzf SPAdes-4.0.0-Linux.tar.gz
cd ./SPAdes-4.0.0-Linux/bin/
vim ~/.bashrc
export PATH=$PATH:/path_to_SPAdes/SPAdes-4.0.0-Linux/bin
source ~/.bashrc
```

### Usage

```bash
Command:
spades.py -1 <trim_1.fq.gz> -2 <trim_2.fq.gz> -o <output_dir> -t <INT>
```

```bash
Example:
spades.py -1 1.trim.fq.gz -2 2.trim.fq.gz -o spades -t 10
```



## (2) GTAB-Minia-Pipeline

**GitHub:** https://github.com/GATB/gatb-minia-pipeline

### Installation

```
python2 -m pip install --user BESST
git clone --recursive https://github.com/GATB/gatb-minia-pipeline
cd gatb-minia-pipeline
make test
```

### Usage

Paired reads assembly: `./gatb -1 read_1.fastq -2 read_2.fastq`

Unpaired reads: `./gatb -s single_reads.fastq`

## (3) RagTag

**GitHub:** https://github.com/malonge/RagTag

### Installation

```bash
conda install -c bioconda ragtag
```

### Dependencies

```bash
Minimap2, Unimap, or Nucmer
Python 3 (with the following auto-installed packages)
	numpy
	intervaltree
	pysam
	networkx
```

```bash
#Minimap2
conda install minimap2

#Unimap
git clone https://github.com/lh3/unimap
cd unimap
make
```

### Usage

#### (1) correct

通过同源基因组进行纠错，该过程不会增加或减少碱基，只是在假定的错误组装点进行打断。

```bash
ragtag.py correct <reference.fa> <query.fa>
```

**参数：**

```bash
<reference.fa> 参考基因组fasta文件 (压缩包也可以)
<query.fa> 所需基因组fasta文件 (压缩包也可以)
correction options:
-f INT 最小唯一比对长度 [默认:1000]
--remove-small 删除比最小唯一比对长度短的唯一比对
-q INT 最小比对质量 (Nucmer比对不适用) [默认:10]
-d INT 最大比对合并距离 [默认:100000]
-b INT 最小距离末端断裂长度[默认:5000]
-e <exclude.txt> 需要剔除的参考序列的列表 [默认:无]
-j <skip.txt> 要保留自己未校正的序列的列表[默认:无]
--inter 只打断参考序列之间的错误组合
--intra 只打破参考序列之内的错误组合
--gff 不在GFF区间内断开序列 [默认:无]
input/output options:
-o PATH 输出目录[./ragtag_output]
-w 覆盖中间文件
-u 为未改变的序列名添加后缀
mapping options:
-t INT minimap2/unimap线程数[默认:1]
–aligner PATH 全基因组比对程序路径(nucmer、unimap或minimap2)[默认:minimap2]
–mm2-params STR 用于minimap2全基因组比对的参数[默认:-x asm5]
–unimap-params STR 用于unimap的参数[默认:-x asm5]
–nucmer-params STR 用于nucmer全基因组比对的参数[默认:–maxmatch -l 100 -c 500]
validation options:
–read-aligner PATH reads比对程序路径（仅允许minimap2）[minimap2]
-R <reads.fasta> 验证用reads（未压缩或gzipped）[默认无]
-F <reads.fofn> 与-R相同，但是是文件列表[默认无]
-T STR 数据类型。sr、ont和corr分别用于Illumina、ont和错误校正后的长reads[默认无]
-v INT 覆盖度验证窗口大小[默认10000]
–max-cov INT 断开高于或等于此覆盖度的区域[默认自动]
–min-cov INT 断开低于或等于此覆盖度的区域[默认自动]
```



# 基因组注释

==基于参考基因组的蛋白编码基因快速同源注释软件：**LiftOn**==

选择近源物种的高质量参考基因组 (带gff格式注释文件) 进行注释

![_images/liftover_illustration.gif](liftover_illustration-3821478.gif)

### Installation

```bash
conda install miniprot
conda install minimap2
pip install lifton
```

### Usage

- **Input**:
  target **Genome** *T* in FASTA format.
  reference **Genome** *R* in FASTA format.
  reference **Annotation** *R~A~* in GFF3 format.
- **Output**:
  LiftOn annotation file, **Annotation** *T~A,~* in GFF3 format.
  Protein sequence identities & mutation types
  Features with extra copies
  Unmapped features

```bash
Command:
lifton -g <reference_genome_annotation.gff> -o <output.gff3> -copies <target_genome.fna> <reference_genome.fna>
```

```bash
Example:
lifton -g d.menogaste_genomic.gff -o lifton.gff3 -copies d.erecta_genomic.fna d.menogaster_genomic.fna
```



# 提取基因组最长CDS

由于基因组中同一个基因可能存在不同的isoform或transcript variant或其他注释问题，所以注释文件中可能存在同一个基因提取出多条冗余序列的情况，提取基因组中所有基因最长的CDS是后续提取单拷贝直系同源基因的前提。

此过程需要**gffread**, **seqkit**软件和一系列脚本来完成，需要使用的脚本及其功能如下：

| Script                       | Function                                                     |
| ---------------------------- | ------------------------------------------------------------ |
| del_seq_length_3.py          | 删除长度小于3的CDS/转录本                                    |
| sort_genelist_from_gff.py    | 从gff注释文件提取gene_id和gene_name，根据gene_name排序，并统计CDS文件中每个gene_id对应的序列的长度，生成list文件gene_list_sort.txt |
| extract_longest_sequences.py | 根据gene_list_sort.txt文件，选择每个基因最长的CDS，过滤冗余的gene_id |
| early_termination_list.py    | 查找无法正常翻译的蛋白序列，即出现翻译提前终止的蛋白序列，创建list文件seq_need_check.txt |
| triplize.py                  | 校正错误序列的读码框：从错误序列的起始位置删除0、1、2个位点，删除后生成对应的一式三份的序列 |
| bestORF.pl                   | 将上一步一式三份的序列进行翻译，并查找校正后可正确翻译的序列，保留校正后可正确翻译的序列，删除校正后仍无法正常翻译的序列，生成.ok文件 |
| del_seq.py                   | 在CDS文件中删除early_termination_list.py脚本产生的list中的序列 |



### Scripts used in this pipeline

You can also download these scripts on FigShare (https://doi.org/10.6084/m9.figshare.26768500.v1)

#### ① **del_seq_length_3.py**

```python
import sys

def read_fasta(file_path):
    sequences = {}
    current_sequence = None
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    return sequences

def filter_sequences(sequences):
    filtered_sequences = {}
    for header, sequence in sequences.items():
        if len(sequence) >= 3:
            filtered_sequences[header] = sequence
    return filtered_sequences

def write_fasta(filtered_sequences, output_file):
    with open(output_file, 'w') as f:
        for header, sequence in filtered_sequences.items():
            f.write(f'>{header}\n{sequence}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    sequences = read_fasta(input_file)
    filtered_sequences = filter_sequences(sequences)
    write_fasta(filtered_sequences, output_file)
```

#### ② sort_genelist_from_gff.py

```python
import sys
from Bio import SeqIO

def extract_info(input_file_gff, input_file_fasta, output_file):
    # Read and extract information from the GFF3 file
    gene_lengths = {}
    gene_info = []  # List to store gene information for sorting
    with open(input_file_gff, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#"):  # Skip comment lines
                continue
            columns = line.strip().split("\t")
            if columns[2] == "mRNA":
                attributes = columns[8].split(";")
                id_value = None
                parent_value = None
                for attribute in attributes:
                    if attribute.startswith("ID="):
                        id_value = attribute.split("=")[1]
                    elif attribute.startswith("Parent="):
                        parent_value = attribute.split("=")[1]
                if id_value and parent_value:
                    # Remove "_" and characters after it from parent_value
                    parent_value = parent_value.split("_")[0]
                    gene_lengths[id_value] = gene_lengths.get(id_value, 0)
                    gene_info.append((parent_value, id_value))

    # Read and calculate gene lengths from the FASTA file
    with open(input_file_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            gene_id = record.id.split()[0]
            if gene_id in gene_lengths:
                gene_lengths[gene_id] += len(record.seq)

    # Sort gene information based on the second column
    gene_info.sort(key=lambda x: x[0])

    # Write sorted results to the output file
    with open(output_file, 'w') as f:
        for parent_value, id_value in gene_info:
            gene_length = gene_lengths.get(id_value, 0)
            f.write(f"{id_value}\t{parent_value}\t{gene_length}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input_file_gff input_file_fasta output_file")
        sys.exit(1)

    input_file_gff = sys.argv[1]
    input_file_fasta = sys.argv[2]
    output_file = sys.argv[3]

    extract_info(input_file_gff, input_file_fasta, output_file)
```

#### ③ extract_longest_sequences.py

```python
import sys

def process_file(input_file):
    gene_sequences = {}

    with open(input_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            gene_id = columns[0]
            gene_name = columns[1]
            sequence_length = int(columns[2])

            if gene_name not in gene_sequences or sequence_length > gene_sequences[gene_name][0]:
                gene_sequences[gene_name] = (sequence_length, gene_id, gene_name)

    for gene_info in gene_sequences.values():
        print('\t'.join(map(str, gene_info)))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    process_file(input_file)
```

#### ④ early_termination_list.py

```python
from Bio import SeqIO
import sys

def find_early_terminations(fasta_file):
    early_termination_sequences = []

    # 使用BioPython读取fasta文件
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        # 检查序列中是否含有提前终止符号
        if '*' in sequence[:-1]:
            early_termination_sequences.append(record.id)

    return early_termination_sequences

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    early_termination_seqs = find_early_terminations(fasta_file)
    if early_termination_seqs:
        print("Sequences with premature termination:")
        for seq_name in early_termination_seqs:
            print(seq_name)
    else:
        print("No sequences with premature termination found.")
```

#### ⑤ triplize.py

```python
import sys

fas = sys.argv[1]

CODE = {
    "standard": {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',  # Alanine
        'TGC': 'C', 'TGT': 'C',  # Cysteine
        'GAC': 'D', 'GAT': 'D',  # Aspartic Acid
        'GAA': 'E', 'GAG': 'E',  # Glutamic Acid
        'TTC': 'F', 'TTT': 'F',  # Phenylalanine
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',  # Glycine
        'CAC': 'H', 'CAT': 'H',  # Histidine
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',  # Isoleucine
        'AAA': 'K', 'AAG': 'K',  # Lysine
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',  # Leucine
        'ATG': 'M',  # Methionine
        'AAC': 'N', 'AAT': 'N',  # Asparagine
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',  # Proline
        'CAA': 'Q', 'CAG': 'Q',  # Glutamine
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R',  # Arginine
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',  # Serine
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',  # Threonine
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',  # Valine
        'TGG': 'W',  # Tryptophan
        'TAC': 'Y', 'TAT': 'Y',  # Tyrosine
        'TAA': 'X', 'TAG': 'X', 'TGA': 'X'  # Stop
    }
    ## more translate table could be added here in future
}


def triplize(fasfile):
    name2seq = fas2hash(fasfile)
    for title in name2seq:
        aex = {}
        seq = name2seq[title].upper()
        for i in range(3):
            j = i
            position = 0
            aex_num = 0
            sequence = []
            while i + 2 < len(seq):
                NNN = seq[i:i + 3]
                sequence.append(NNN)
                if NNN in CODE['standard']:
                    aa = CODE['standard'][NNN]
                    if aa == 'X':
                        aex_num = aex_num + 1
                        position = i
                i = i + 3
            # if aex_num == 0 or aex_num == 1:
            # if aex_num == 0:
            # if aex_num == 1:
            triple_seq = ''.join(sequence)
            print(">%s|%s|%d|%d|%d\n%s" % (title, aex_num, position, len(seq), j, triple_seq))


def fas2hash(fasfile):
    name2seq = {}
    with open(fasfile) as fh:
        for line in fh.readlines():
            line = line.replace('\n', '')
            if line.startswith('>'):
                seq_title = line.replace('>', '')
                name2seq[seq_title] = ''
            else:
                name2seq[seq_title] += line
    return name2seq


triplize(fas)
```

#### ⑥ bestORF.pl

```perl
open IN,$ARGV[0];
while (<IN>){
        chomp;
        if (/^>/){
                s/^>//;
                $title = $_;
        }
        else{$hash{$title} .= $_}
}
close IN;

open OUT, ">$ARGV[0].ok";
foreach (keys %hash){
        @line = split(/\|/,$_);
        if ($line[4] == 0){
                $m++;
                if ($line[3]%3 == 0){
                        $y++;
                        if ($hash{$_} =~ /^ATG/ and $line[1] == 1 and $line[2]+3 == $line[3]){
                                $p++;
                                $temp{$line[0]} = 1;
                                print OUT ">$line[0]\n$hash{$_}\n";
                        }#perfect_start_stop
                        if ($hash{$_} !~ /^ATG/ and $line[1] == 1 and $line[2]+3 == $line[3]){
                                $i++;
                                $temp{$line[0]} = 1;
                                print OUT ">$line[0]\n$hash{$_}\n";
                        }#intact_stop_codon at end
                        if ($line[1] == 0 and $hash{$_} =~ /^ATG/){
                                $j++;
                                $temp{$line[0]} = 1;
                                print OUT ">$line[0]\n$hash{$_}\n";
                        }#intact_with_start_codon
                        if ($line[1] == 0 and $hash{$_} !~ /^ATG/){
                                $n++;
                                $temp{$line[0]} = 1;
                                print OUT ">$line[0]\n$hash{$_}\n";
                        }#no_start_no_stop_codon
                        if ($line[1] == 1 and $line[2]+3 != $line[3]){$q++;}#internal_one_stop
                        if ($line[1] > 1){$x++;}#multiple stop codon

                }
                else {
                        $k++;
                }
        }
}

foreach (keys %hash){
        @line = split(/\|/,$_);
        unless ($temp{$line[0]}){
                if ($line[4] > 0){
                        $length =  length($hash{$_});
                        if ($hash{$_} =~ /^ATG/ and $line[1] == 1 and $line[2]-$line[4]+3 == $length){
                                $z1++;
                                print OUT ">$line[0]\n$hash{$_}\n";
                                $temp{$line[0]} = 1;
                        }
                        if ($hash{$_} =~ /^ATG/ and $line[1] == 0){
                                $z2++;
                                print OUT ">$line[0]\n$hash{$_}\n";
                                $temp{$line[0]} = 1;
                        }
                        if ($hash{$_} !~ /^ATG/ and $line[1] == 1 and $line[2]-$line[4]+3 == $length){
                                $z3++;
                                print OUT ">$line[0]\n$hash{$_}\n";
                                $temp{$line[0]} = 1;
                        }
                        if ($line[1] == 0 and $hash{$_} !~ /^ATG/){
                                $z4++;
                                print OUT ">$line[0]\n$hash{$_}\n";
                                $temp{$line[0]} = 1;
                        }
                        else {
                                print "$_\n";
                        }
                }
        }
}

open IN,$ARGV[1];
while (<IN>){
        chomp;
        @line1 = split;
        $length = length($hash{$line1[0]});
        @info = split /\|/,$line1[0];
        $gene = $info[0];
        $a = $line1[6]-1;
        $b = $a % 3;
        if ($info[4] == 1){
                if (!exists $temp{$gene} && $b == 0 && $info[1] == 0){
                        print OUT ">$gene\n$hash{$line[0]}\n";
                        $temp{$gene} = 1;
                        $zz++;
                        #print "$zz\n";
                }
                if (!exists $temp{$gene} && $b == 0 && $info[1] == 1 && $info[2]-$info[4]+3 == $length){
                        print OUT ">$gene\n$hash{$line[0]}\n";
                        $temp{$gene} = 1;
                        $zzz++;
                }
        }
}

close IN;
print "total_cds_number: $m\n";
print "number of cds triple: $y\n";
print "\tperfect: $p\n";
print "\tintact_with_stop_condon: $i\n";
print "\tintact_with_start_codon: $j\n";
print "\tno_start_no_stop_codon: $n\n";
print "\tinternal_one_stop codon: $q\n";
print "\tmultiple stop codon: $x\n";
print "bad_cds: $k\n";
print "\tperfect in bad cds: z1:$z1\tz2:$z2\tz3:$z3\tz4:$z4\tz5:$zz\tz6:$zzz\n";
```

#### ⑦ del_seq.py

```python
import sys

def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence = ''
        seq_name = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences[seq_name] = sequence
                seq_name = line[1:]
                sequence = ''
            else:
                sequence += line
        if sequence:
            sequences[seq_name] = sequence
    return sequences

def write_fasta(file_path, sequences):
    with open(file_path, 'w') as file:
        for seq_name, sequence in sequences.items():
            file.write('>' + seq_name + '\n')
            file.write(sequence + '\n')

def filter_sequences(fasta_file, list_file):
    sequences = read_fasta(fasta_file)
    with open(list_file, 'r') as file:
        sequences_to_remove = set(line.strip() for line in file)
    sequences = {name: sequence for name, sequence in sequences.items() if name not in sequences_to_remove}
    write_fasta(fasta_file, sequences)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py fasta_file list_file")
        sys.exit(1)
    fasta_file = sys.argv[1]
    list_file = sys.argv[2]
    filter_sequences(fasta_file, list_file)
    print("Filtered sequences saved to", fasta_file)
```



### ==Usage==

1. **gffread**
   使用gffread软件从基因组中提取CDS序列

   ```bash
   gffread genome_annotation.gff -g genome.fna -x cds.fna
   ```

2. **seqkit**
   使用seqkit软件对genome_cds.fna文件中的序列进行多行转一行操作

   ```bash
   seqkit seq cds.fna -w 0 > species_cds.fna
   ```

3. **del_seq_length_3.py**
   使用del_seq_length_3.py脚本删除长度小于3的CDS/转录本

   ```bash
   python3 del_seq_length_3.py species_cds.fna
   ```

4. **sort_genelist_from_gff.py**
   从gff注释文件提取gene_id和gene_name，根据gene_name排序，并统计CDS文件中每个gene_id对应的序列的长度，生成list文件gene_list_sort.txt

   ```bash
   python3 sort_genelist_from_gff.py genome_annotation.gff species_cds.fna species_genelist_sort.txt
   ```

5. **extract_longest_sequence.py**
   根据gene_list_sort.txt文件，选择每个基因最长的CDS，过滤冗余的gene_id

   ```bash
   python3 extract_longest_sequence.py species_genelist_sort.txt > species_genelist_sort_filter.txt
   ```

6. **Shell**
   根据species_genelist_sort_filter.txt结果在CDS文件species_cds.fna提取对应的序列，并生成可执行文件extract_longest_sequence.sh

   ```bash
   gawk '{print"grep -A 1 -w \""$2"\" species_cds.fna"}' species_genelist_sort_filter.txt > extract_longest_sequence.sh
   ```

7. **Shell**
   执行extract_longest_sequence.sh，将结果输入到species_longest_cds.fna文件中，获得每个基因最长CDS序列

   ```bash
   bash extract_longest_sequence.sh > species_longest_cds.fna
   ```

8. **seqkit**
   将最长CDS序列翻译为蛋白序列

   ```bash
   seqkit translate species_longest_cds.fna > species_longest_prot.faa
   ```

9. **early_termination_list.py**
   查找无法正常翻译的蛋白序列，即出现翻译提前终止的蛋白序列，创建list文件seq_need_check.txt

   ```bash
   python3 early-termination.py species_longest_prot.faa > seq_need_check.txt
   ```

10. **Shell**
    根据seq_need_check.txt结果在CDS文件species_cds.fna提取对应的序列，并生成可执行文件seq_need_checked.sh

    ```bash
    gawk '{grep -A 1 -w \""$1"\" species_longest_cds.fna"}' seq_need_check.txt > seq_need_checked.sh
    ```

11. **Shell**
    执行seq_need_checked.sh，将结果输入到seq_need_checked.seq文件中

    ```bash
    bash seq_need_checked.sh > seq_need_checked.seq
    ```

12. **triplize.py**
    校正错误序列的读码框：从错误序列的起始位置删除0、1、2个位点，删除后生成对应的一式三份的序列

    ```bash
    python triplize.py seq_need_checked.seq > seq_need_checked_triplize.seq
    ```

13. **bestORF.pl**
    将上一步一式三份的序列进行翻译，并查找校正后可正确翻译的序列，保留校正后可正确翻译的序列，删除校正后仍无法正常翻译的序列，生成.ok文件

    ```bash
    perl bestORF.pl seq_need_checked_triplize.seq
    ```

14. **del_seq.py**
    在CDS文件中删除early_termination_list.py脚本产生的list中的序列

    ```bash
    python del_seq.py species_longest_cds.fna seq_need_check.txt
    ```

15. **Shell**
    合并删除错误序列后的species_longest_cds.fna文件与校正开放阅读框后的seq_need_checked_triplize.seq.ok文件，获得最终的最长CDS序列

    ```bash
    cat species_longest_cds.fna seq_need_checked_triplize.seq.ok > species_longest_cds_ok.fna
    ```

16. **seqkit**
    通过seqkit翻译最长CDS序列得到最长蛋白序列

    ```bash
    seqkit translate species_longest_cds_ok.fna > species_longest_prot_ok.faa
    ```

    

# OrthoFinder提取直系同源基因

软件主页：https://github.com/davidemms/OrthoFinder

### Installation

```bash
conda install orthofinder
```

### Usage

将所有物种或基因组的最长CDS放入一个文件夹，例如，./orthofinder

```bash
orthofinder -f ./orthofinder
```

### Output

结果文件为Results_Date，例如Results_May06

```bash
└── Results_May06
    ├── Citation.txt
    ├── Comparative_Genomics_Statistics
    ├── Gene_Duplication_Events
    ├── Gene_Trees
    ├── Log.txt
    ├── Orthogroups
    ├── Orthogroup_Sequences
    ├── Orthologues
    ├── Phylogenetically_Misplaced_Genes
    ├── Phylogenetic_Hierarchical_Orthogroups
    ├── Putative_Xenologs
    ├── Resolved_Gene_Trees
    ├── Single_Copy_Orthologue_Sequences
    ├── Species_Tree
    └── WorkingDirectory
```

需要重点关注的文件夹有**Single_Copy_Orthologue_Sequences**, **Orthogroups**和**Orthogroup_Sequences**，其中**Orthogroups**和**Orthogroup_Sequences**文件夹中的内容是重要的。

(1) **(Optional) Single_Copy_Orthologue_Sequences**: 包含了所有物种或基因组的<u>单拷贝直系同源基因</u>，可直接用于后续序列比对和建树。但此文件夹里包含的序列是所有物种共有的，所以无法设置missing data，及matrix的完整性为100%；
例如：

```bash
						species_A		species_B		species_C		species_D		species_E			......
OG0003103				 √					√						√						√						√
OG0003106				 √					√						√						√						√
OG0003107				 √					√						√						√						√
OG0003108				 √					√						√						√						√
OG0003109				 √					√						√						√						√
OG0003110				 √					√						√						√						√
......
```

这一步也可以不做，因为**Orthogroups**文件夹包含了Single_Copy_Orthologue_Sequences的结果。

(2) **Orthogroups**: 主要关注**Orthogroups.tsv**文件 (可以在Excel中打开)。OrthoFinder会在这个文件中记录在所有物种或基因组中搜索出的one to one, one to many, many to one和many to many的gene_id，并且会记录哪些基因在哪些物种或基因组中丢失，因此可以根据这个文件的信息设置missing data，这样可以扩大matrix的数据量；
需要过滤掉one to many, many to one和many to many的结果，只选择one to one的Orthogroup和存在missing gene的Orthogroup；
例如：

```bash
						species_A		species_B		species_C		species_D		species_E			......
OG0003103				 √					√						√						√						√
OG0003106				 √					√						√						√						√
OG0003107				 √					√						√						√						√
OG0003108				 √					√						√						√						√
OG0003109				 √					√						√						√						√
OG0003110				 √					√						√						√						√
OG0003111				 √					×						×						√						√
OG0003112				 √					√						√						√						×
OG0003113				 √					×						√						√						√
OG0003114				 √					√						×						√						×
OG0003115				 √					×						√						×						×
......
```

==这一步推荐使用脚本**orthogroups_filter.py**过滤：==

Usage: `python orthogroups_filter.py <input_file> <output_file>`

**orthogroups_filter.py**

```python
import csv
import sys

def filter_rows(input_file, output_file):
    with open(input_file, 'r', newline='') as f_input, \
         open(output_file, 'w', newline='') as f_output:
        
        reader = csv.reader(f_input, delimiter='\t')
        writer = csv.writer(f_output, delimiter='\t')
        
        # 读取表头
        header = next(reader)
        writer.writerow(header)

        # 找出每列包含2个及以上序列号的行
        for row in reader:
            skip_row = False
            for col in row[1:]:  # 跳过第一列基因名
                if len(col.split(', ')) >= 2:
                    skip_row = True
                    break
            if not skip_row:
                writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    filter_rows(input_file, output_file)
    print("过滤完成，结果保存在", output_file)
```

(3) **Orthogroup_Sequences**: 包含了所有Orthogroup的序列信息，可以根据 (2) 中**Orthogroups**文件夹中需要的Orthogroups提取序列。推荐





# 核苷酸序列 (CDS) 的获取

上一步根据OrthoFinder的输出结果是氨基酸序列，为了设计不同的数据矩阵 (matrix)，我们同样需要获取氨基酸序列对应的CDS核苷酸序列。





# 序列比对

序列比对的软件很多，建议使用MAFFT或prank。

Note: 对于蛋白编码基因，序列比对的策略是先比对氨基酸序列，再根据比对后的氨基酸序列使用脚本back-translate到核苷酸序列；
其中prank可以根据codon模式进行比对 (不需要先比对氨基酸再back-translate到核苷酸序列)，但运行速度慢。

## (1) 氨基酸序列比对

### MAFFT

#### Installation

```bash
conda install mafft
```

#### Usage

##### **氨基酸序列比对**

```bash
mafft --maxiterate 1000 --localpair seq.faa > seq_mafft.faa
```

##### **批量执行序列比对**

###### ① 创建list

例如，orthofinder生成的结果为：

```bash
OG0003103.faa	OG0003106.faa	OG0003107.faa	OG0003108.faa	OG0003109.faa	OG0003110.faa	OG0003111.faa
OG0003112.faa	OG0003113.faa	OG0003114.faa	OG0003115.faa	OG0003116.faa	OG0003117.faa	OG0003118.faa
OG0003119.faa	OG0003120.faa	OG0003121.faa	OG0003122.faa	OG0003123.faa	OG0003125.faa	......
```

对所有需要比对的序列生成list，可删除后缀

```bash
ls *.fa | sed 's/.fa//g' > list 
```

执行此命令后生成的list为：

```bash
OG0003103	OG0003106	OG0003107	OG0003108	OG0003109	OG0003110	OG0003111
OG0003112	OG0003113	OG0003114	OG0003115	OG0003116	OG0003117	OG0003118
OG0003119	OG0003120	OG0003121	OG0003122	OG0003123	OG0003125	......
```



###### ② 使用gawk批量打印命令

**使用gawk批量打印命令，将结果输出到可执行文件mafft.sh中**

```bash
gawk '{print"mafft --maxiterate 1000 --localpair "$1".faa > "$1".mafft.faa"}' list > mafft.sh
```

mafft.sh的结果应为：

```bash
mafft --maxiterate 1000 --localpair OG0003103.faa > OG0003103.mafft.faa
mafft --maxiterate 1000 --localpair OG0003106.faa > OG0003106.mafft.faa
mafft --maxiterate 1000 --localpair OG0003107.faa > OG0003107.mafft.faa
mafft --maxiterate 1000 --localpair OG0003108.faa > OG0003108.mafft.faa
mafft --maxiterate 1000 --localpair OG0003109.faa > OG0003109.mafft.faa
mafft --maxiterate 1000 --localpair OG0003110.faa > OG0003110.mafft.faa
mafft --maxiterate 1000 --localpair OG0003111.faa > OG0003111.mafft.faa
mafft --maxiterate 1000 --localpair OG0003112.faa > OG0003112.mafft.faa
mafft --maxiterate 1000 --localpair OG0003113.faa > OG0003113.mafft.faa
mafft --maxiterate 1000 --localpair OG0003114.faa > OG0003114.mafft.faa
mafft --maxiterate 1000 --localpair OG0003115.faa > OG0003115.mafft.faa
mafft --maxiterate 1000 --localpair OG0003116.faa > OG0003116.mafft.faa
mafft --maxiterate 1000 --localpair OG0003117.faa > OG0003117.mafft.faa
mafft --maxiterate 1000 --localpair OG0003118.faa > OG0003118.mafft.faa
mafft --maxiterate 1000 --localpair OG0003119.faa > OG0003119.mafft.faa
mafft --maxiterate 1000 --localpair OG0003120.faa > OG0003120.mafft.faa
mafft --maxiterate 1000 --localpair OG0003121.faa > OG0003121.mafft.faa
mafft --maxiterate 1000 --localpair OG0003122.faa > OG0003122.mafft.faa
mafft --maxiterate 1000 --localpair OG0003123.faa > OG0003123.mafft.faa
mafft --maxiterate 1000 --localpair OG0003125.faa > OG0003125.mafft.faa
......
```



###### ③ 使用软件ParaFly批量执行序列比对命令

**安装ParaFly：**`conda install -c bioconda parafly`

**运行ParaFly：**`nohup ParaFly -c <command.sh> -CPU <INT> &`

Example: `nohup ParaFly -c mafft.sh -CPU 60 &`
