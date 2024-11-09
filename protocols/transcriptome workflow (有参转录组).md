# transcriptome workflow (有参转录组)

**required software:** fastp, hisat2, samtools, stringtie, gffread, Deseq2 (R包)

**workflow参考:** Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650-1667.



## (1) fastp (对原始数据raw reads进行质控)


**Installation:** `conda install fastp`

| Parameter | Note                                                         |
| :-------: | :----------------------------------------------------------- |
|    -i     | 后接需要过滤的fastq文件(上游)                                |
|    -I     | 后接需要过滤的fastq文件(下游)                                |
|    -o     | 后接过滤完输出的的fastq文件(上游)                            |
|    -O     | 后接过滤完输出的fastq文件(下游)                              |
|    -h     | html报告文件名                                               |
|    -q     | 设置碱基质量阈值，默认阈值为15                               |
|    -u     | 设置允许不合格碱基阈值，超过比例时，整条read将被遗弃，默认为40 |
|    -l     | 设置read的最小长度，默认是15，长度＜15的read将被遗弃         |
|    -w     | 设置软件可用线程数，默认为2                                  |

```bash
Command:
fastp -i <raw_data_1.fastq> -I <raw_data_2.fastq> -o <clean_data_1.fastq> -O <clean_data_2.fastq> -h <html_file> -q <INT> -u <INT> -l <INT> -w <INT>
```

```bash
Example:
fastp -i frog_R1.fastq -I frog_R2.fastq -o frog_R1_clean.fastq -O frog_R2_clean.fastq -l 30 -w 10
```



## (2) hisat2 (比对到参考基因组)


**Installation:** `conda install hisat2`

### ① Index 建立索引

```bash
Command:
mkdir Index
hisat2-build -f genome.fna ./Index/basename
```

```bash
Example:
hisat2-build -f GCA_002284835.2_RCv2.1_genomic.fna /home/Index/frog_index
```



### ② Alignment 使用生成的索引文件对RNA-Seq进行比对拼接

| Parameter | Note                                  |
| :-------: | :------------------------------------ |
|    -p     | 线程数目                              |
|    -x     | 参考基因组索引的basename，即前缀名    |
|    -1     | 双端测序文件的read1                   |
|    -2     | 双端测序文件的read2                   |
|    -S     | SAM写入的文件名，默认写入到标准输出中 |

```bash
Command:
hisat2 -p <INT> -x ./Index/basename -1 <fastq_file_1.fastq> -2 <fastq_file_2.fastq> -S <output_sam_file.sam>
```

```bash
Example:
hisat2 -p 10 -x ./Index/forg_index -1 frog_R1_clean.fastq -2 frog_R2_clean.fastq -S frog.sam
```



## (3) SAMtools (排序及格式转换)

**Installation:** `conda install samtools=1.15.1` (需要安装samtools 1.15.1版本)

```bash
Command:
samtools sort -@ 4 -o <output_bam_file.bam> <input_sam_file.sam>
```

```bash
Example:
samtools sort -@ 4 -o frog.bam frog.sam 
```



## (4) stringtie (提取转录本丰度)

### ① 提取转录本丰度

**Installation:** `conda install stringtie`

| Parameter | Note                                                         |
| :-------: | :----------------------------------------------------------- |
|    -G     | 使用注释好的参考基因组注释文件文件辅助组装，在-e未设置的条件下，输出中包括注释文件中的转录本和预测的新型转录本 |
|    -o     | 输出文件的名字，最好是全路径，默认输出为标准输出             |
|    -l     | 为输出的转录本设置前缀名，默认为STRG                         |
|    -p     | 线程数，默认为1                                              |
|    -A     | 对输出的gff统计基因表达量，并以一个tab分割的文件输出，这里需要提交输出的文件名 |
|    -C     | 对输出的gff中属于-G提交的参考gtf的转录本统一输出到该文件，这里需要提交一个文件名 |
|    -B     | 是否需要输出Ballgown可以识别的文件,在-b设置的情况下，使用-o的路径输出 |
|    -b     | 对Ballgown输出的文件指定一个路径保存                         |
|    -e     | 只统计可以匹配-G提交的参考gff中的转录本，不再对新的转录本做预测，这可以加快程序的运行速度 |
|    -m     | 对预测的转录本设置最小长度，默认为200                        |

```bash
Command:
stringtie -B -p <INT> -G <genome_annotation_file.gff> -o <transcripts.gtf> <bam_file.bam>
#若参考基因组无注释文件: stringtie -p <INT> -o <transcripts.gtf> <bam_file.bam>
```

```bash
Example:
stringtie -B -p 10 -G GCA_002284835.2_RCv2.1_genomic.gff3 -o ./ballgowns/frog.gtf frog.bam
#若参考基因组无注释文件: stringtie -p 10 -o ./ballgowns/frog.gtf frog.bam
```



 **(Optional)** merge合并多组分析转录本 (适用于将多个生物学重复的结果合并为一个转录本)

将多个生物学重复的样本号 (例如frog_1, frog_2, frog_3) 放在一个gtf_list中，每行粘贴一个样本号

```bash
Command:
stringtie --merge -p 10 -G <genome_annotation_file.gff> -o <transcripts.gtf> -l merge <gtf_list>
```

```bash
Example:
stringtie --merge -p 10 -G GCA_002284835.2_RCv2.1_genomic.gff3 -o frog_merged.gtf -l merge gtf_list
```



### ② gffread生成转录本序列

**Installation:** `conda install gffread`

将得到的gtf文件转化为fasta序列文件，即转录组。-g输入参考基因组序列，-w选择输出文件，保存transcripts.fasta文件中。从而得到RNAseq的全部转录本。

```bash
Command:
gffread <transcripts.gtf> -g <genome.fna> -w <transcripts.fasta>
```

```bash
Example:
gffread merged.gtf -g GCA_002284835.2_RCv2.1_genomic.fna -w total_merged.fas
```



### ③ prepDE.py3脚本提取reads count

prepDE.py3下载链接: https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3

**prepDE.py3**

```python
#!/usr/bin/env python3
import re, csv, sys, os, glob, warnings, itertools
from math import ceil
from optparse import OptionParser
from operator import itemgetter

parser=OptionParser(description='Generates two CSV files containing the count matrices for genes and transcripts, using the coverage values found in the output of `stringtie -e`')
parser.add_option('-i', '--input', '--in', default='.', help="a folder containing all sample sub-directories, or a text file with sample ID and path to its GTF file on each line [default: %default/]")
parser.add_option('-g', default='gene_count_matrix.csv', help="where to output the gene count matrix [default: %default")
parser.add_option('-t', default='transcript_count_matrix.csv', help="where to output the transcript count matrix [default: %default]")
parser.add_option('-l', '--length', default=75, type='int', help="the average read length [default: %default]")
parser.add_option('-p', '--pattern', default=".", help="a regular expression that selects the sample subdirectories")
parser.add_option('-c', '--cluster', action="store_true", help="whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)")
parser.add_option('-s', '--string', default="MSTRG", help="if a different prefix is used for geneIDs assigned by StringTie [default: %default]")
parser.add_option('-k', '--key', default="prepG", help="if clustering, what prefix to use for geneIDs assigned by this script [default: %default]")
parser.add_option('-v', action="store_true", help="enable verbose processing")

parser.add_option('--legend', default="legend.csv", help="if clustering, where to output the legend file mapping transcripts to assigned geneIDs [default: %default]")
(opts, args)=parser.parse_args()

samples = [] # List of tuples. If sample list, (first column, path). Else, (subdirectory name, path to gtf file in subdirectory)
if (os.path.isfile(opts.input)):
    # gtfList = True
    try:
        fin = open(opts.input, 'r')
        for line in fin:
            if line[0] != '#':
                lineLst = tuple(line.strip().split(None,2))
                if (len(lineLst) != 2):
                    print("Error: line should have a sample ID and a file path:\n%s" % (line.strip()))
                    exit(1)
                if lineLst[0] in samples:
                    print("Error: non-unique sample ID (%s)" % (lineLst[0]))
                    exit(1)
                if not os.path.isfile(lineLst[1]):
                    print("Error: GTF file not found (%s)" % (lineLst[1]))
                    exit(1)
                samples.append(lineLst)
    except IOError:
        print("Error: List of .gtf files, %s, doesn't exist" % (opts.input))
        exit(1)
else:
    # gtfList = False
    ## Check that opts.input directory exists
    if not os.path.isdir(opts.input):
      parser.print_help()
      print(" ")
      print("Error: sub-directory '%s' not found!" % (opts.input))
      sys.exit(1)
    #####
    ## Collect all samples file paths and if empty print help message and quit
    #####
    samples = []
    for i in next(os.walk(opts.input))[1]:
        if re.search(opts.pattern,i):
         for f in glob.iglob(os.path.join(opts.input,i,"*.gtf")):
            samples.append((i,f)) 

if len(samples) == 0:
  parser.print_help()
  print(" ")
  print("Error: no GTF files found under base directory %s !" % (opts.input))
  sys.exit(1)

RE_GENE_ID=re.compile('gene_id "([^"]+)"')
RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')
RE_COVERAGE=re.compile('cov "([\-\+\d\.]+)"')
RE_STRING=re.compile(re.escape(opts.string))

RE_GFILE=re.compile('\-G\s*(\S+)') #assume filepath without spaces..


#####
## Sort the sample names by the sample ID
#####

samples.sort()

#if opts.v:
#  print "Sample GTFs found:"
#  for s in samples:
#     print s[1]

#####
## Checks whether a given row is a transcript 
## other options: ex. exon, transcript, mRNA, 5'UTR
#####
def is_transcript(x):
  return len(x)>2 and x[2]=="transcript"

def getGeneID(s, ctg, tid):
  r=RE_GENE_ID.search(s)
  #if r: return r.group(1)
  rn=RE_GENE_NAME.search(s)
  #if rn: return ctg+'|'+rn.group(1)
  if r:
    if rn: 
      return r.group(1)+'|'+rn.group(1)
    else:
      return r.group(1)
  return tid

def getCov(s):
  r=RE_COVERAGE.search(s)
  if r:
    v=float(r.group(1))
    if v<0.0: v=0.0
    return v
  return 0.0

def is_overlap(x,y): #NEEDS TO BE INTS!
  return x[0]<=y[1] and y[0]<=x[1]


def t_overlap(t1, t2): #from badGenes: chromosome, strand, cluster, start, end, (e1start, e1end)...
    if t1[0] != t2[0] or t1[1] != t2[1] or t1[5]<t2[4]: return False
    for i in range(6, len(t1)):
        for j in range(6, len(t2)):
            if is_overlap(t1[i], t2[j]): return True
    return False

## Average Readlength
read_len=opts.length

## Variables/Matrices to store t/g_counts
t_count_matrix, g_count_matrix=[],[]

##Get ready for clustering, stuff is once for all samples##
geneIDs={} #key=transcript, value=cluster/gene_id


## For each of the sorted sample paths
for s in samples:
    badGenes=[] #list of bad genes (just ones that aren't MSTRG)
    try:
        ## opts.input = parent directory of sample subdirectories
        ## s = sample currently iterating through
        ## os.path.join(opts.input,s,"*.gtf") path to current sample's GTF
        ## split = list of lists: [[chromosome, ...],...]

        #with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f:
        #    split=[l.split('\t') for l in f.readlines()]
#        if not gtfList:
#            f = open(glob.iglob(os.path.join(opts.input,s[1],"*.gtf")).next())
#        else:
#            f = open(s[1])
        with open(s[1]) as f:
            split=[l.split('\t') for l in f.readlines()]

        ## i = numLine; v = corresponding i-th GTF row
        for i,v in enumerate(split):
            if is_transcript(v):
                t_id=RE_TRANSCRIPT_ID.search(v[8]).group(1)
                try:
                  g_id=getGeneID(v[8], v[0], t_id)
                except:
                  print("Problem parsing file %s at line:\n:%s\n" % (s[1], v))
                  sys.exit(1)
                geneIDs.setdefault(t_id, g_id)
                if not RE_STRING.match(g_id):
                    badGenes.append([v[0],v[6], t_id, g_id, min(int(v[3]),int(v[4])), max(int(v[3]),int(v[4]))]) #chromosome, strand, cluster/transcript id, start, end
                    j=i+1
                    while j<len(split) and split[j][2]=="exon":
                        badGenes[len(badGenes)-1].append((min(int(split[j][3]), int(split[j][4])), max(int(split[j][3]), int(split[j][4]))))
                        j+=1

    except StopIteration:
        warnings.warn("Didn't get a GTF in that directory. Looking in another...")

    else: #we found the "bad" genes!
        break

##THE CLUSTERING BEGINS!##
if opts.cluster and len(badGenes)>0:
    clusters=[] #lists of lists (could be sets) or something of transcripts
    badGenes.sort(key=itemgetter(3)) #sort by start coord...?
    i=0
    while i<len(badGenes): #rather un-pythonic
        temp_cluster=[badGenes[i]]

        k=0
        while k<len(temp_cluster):
            j=i+1
            while j<len(badGenes):
                if t_overlap(temp_cluster[k], badGenes[j]):
                    temp_cluster.append(badGenes[j])
                    del badGenes[j]
                else:
                    j+=1
            k+=1
        if len(temp_cluster)>1:
            clusters.append([t[2] for t in temp_cluster])
        i+=1

    print(len(clusters))

    for c in clusters:
        c.sort()

    clusters.sort(key=itemgetter(0))
    legend=[]
    for u,c in enumerate(clusters):
        my_ID=opts.key+str((u+1))
        legend.append(list(itertools.chain.from_iterable([[my_ID],c]))) #my_ID, clustered transcript IDs
        for t in c:
            geneIDs[t]=my_ID
##            geneIDs[t]="|".join(c) #duct-tape transcript IDs together, disregarding ref_gene_names and things like that

    with open(opts.legend, 'w') as l_file:
        my_writer=csv.writer(l_file)
        my_writer.writerows(legend)

geneDict={} #key=gene/cluster, value=dictionary with key=sample, value=summed counts
t_dict={}
guidesFile='' # file given with -G for the 1st sample
for q, s in enumerate(samples):
    if opts.v:
       print(">processing sample %s from file %s" % s)
    lno=0
    try:
        #with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f: #grabs first .gtf file it finds inside the sample subdirectory
#        if not gtfList:
#            f = open(glob.iglob(os.path.join(opts.input,s[1],"*.gtf")).next())
#        else:
        f = open(s[1])
        transcript_len=0
        
        for l in f:
            lno+=1
            if l.startswith('#'):
                if lno==1:
                    ei=l.find('-e')
                    if ei<0:
                       print("Error: sample file %s was not generated with -e option!" % ( s[1] ))
                       sys.exit(1)
                    gf=RE_GFILE.search(l)
                    if gf:
                       gfile=gf.group(1)
                       if guidesFile:
                          if gfile != guidesFile:
                             print("Warning: sample file %s generated with a different -G file (%s) than the first sample (%s)" % ( s[1], gfile, guidesFile ))
                       else:
                          guidesFile=gfile
                    else:
                       print("Error: sample %s was not processed with -G option!" % ( s[1] ))
                       sys.exit(1)
                continue
            v=l.split('\t')
            if v[2]=="transcript":
                if transcript_len>0:
##                        transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
                    t_dict.setdefault(t_id, {})
                    t_dict[t_id].setdefault(s[0], int(ceil(coverage*transcript_len/read_len)))
                t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
                #g_id=RE_GENE_ID.search(v[len(v)-1]).group(1)
                g_id=getGeneID(v[8], v[0], t_id)
                #coverage=float(RE_COVERAGE.search(v[len(v)-1]).group(1))
                coverage=getCov(v[8])
                transcript_len=0
            if v[2]=="exon":
                transcript_len+=int(v[4])-int(v[3])+1 #because end coordinates are inclusive in GTF

##            transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
        t_dict.setdefault(t_id, {})
        t_dict[t_id].setdefault(s[0], int(ceil(coverage*transcript_len/read_len)))

    except StopIteration:
#        if not gtfList:
#            warnings.warn("No GTF file found in " + os.path.join(opts.input,s[1]))
#        else:
        warnings.warn("No GTF file found in " + s[1])


##        transcriptList.sort(key=lambda bla: bla[1]) #gene_id
    
    for i,v in t_dict.items():
##        print i,v
       try:
          geneDict.setdefault(geneIDs[i],{}) #gene_id
          geneDict[geneIDs[i]].setdefault(s[0],0)
          geneDict[geneIDs[i]][s[0]]+=v[s[0]]
       except KeyError:
          print("Error: could not locate transcript %s entry for sample %s" % ( i, s[0] ))
          raise

if opts.v:
   print("..writing %s " % ( opts.t ))
with open(opts.t, 'w') as csvfile:
   my_writer = csv.DictWriter(csvfile, fieldnames = ["transcript_id"] + [x for x,y in samples])
   my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
   for i in t_dict:
        t_dict[i]["transcript_id"] = i
        my_writer.writerow(t_dict[i])
if opts.v:
   print("..writing %s " % ( opts.g ))
with open(opts.g, 'w') as csvfile:
   my_writer = csv.DictWriter(csvfile, fieldnames = ["gene_id"] + [x for x,y in samples])
##    my_writer.writerow([""]+samples)
##    my_writer.writerows(geneDict)
   my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
   for i in geneDict:
        geneDict[i]["gene_id"] = i #add gene_id to row
        my_writer.writerow(geneDict[i])
if opts.v:
   print("All done.")
```

```bash
Command:
每一个样本创建一个字文件夹，将gtf文件放入子文件夹，再将所有子文件夹放入一个文件夹中，例如ballgowns
python3 prepDE.py3 -i ballgowns

#或者需要准备一个sample_list.txt，该文件为\t制表符分隔的两列，第一列为样本名称，第二列为定量的gtf文件的路径，示例如下：
	forg_1	frog_1.gtf
	frog_2	frog_2.gtf
	frog_3	frog_3.gtf

Command:
python3 prepDE.py3 -i sample_list.txt
```



## (5) 差异表达分析

推荐使用在线分析平台: https://www.omicshare.com/tools/Home/Soft/diffanalysis

## (6) 富集分析

推荐使用在线分析平台

GO富集:
https://omicshare.com/tools/Home/Soft/gogseasenior
https://www.omicshare.com/tools/home/report/goenrich.html

KEGG富集:
https://www.omicshare.com/tools/Home/Soft/pathwaygseasenior
https://www.omicshare.com/tools/home/report/koenrich.html