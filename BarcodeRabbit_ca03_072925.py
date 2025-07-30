#!/usr/bin/env python
## 
## ######################
## BarcodeRabbit-- SSU rRNA-based Reference-agnostic phylogenetic survey of cellular composition from NGS datasets
## Version ca03 072925
## ->Fuction: BarcodeRabbit is designed to make use of a highly conserved "core" motif in defining composition of source material for NGS datasets
##  - The general question is "I have a(some) NGS datasets and want to know how rich the source material was in both known and unknown species.
##  - Two applications are envisioned for BarcodeRabbit
##  - 1.  Provisional assessment of cellular phylogenetic composition that may provide a broader species assessment than standard k-mer tools 
##  - 2.  As an initial tool for discovery of potentially novel species/genus or larger phylogenetic groups
##  - Note that BarcodeRabbit shares considerable conceptual foundation with a number of other phylogenetic tools (e.g. Kraken), with
##  - its unique characteristics centered on the possibility of discovering novel cellular categories in sequencing datasets
##  - A further (major) footnote is that many different processes in molecular evolution and the experimental determination of sequences
##  - can produce sequence reads that would appear novel or otherwise of potential interest.  Thus BarcodeRabbit is intended solely to nominate
##  - candidates rather than to produce output that can be directly interpreted without considerable further work (computational and potentially experimental)
##
## ->Syntax
##  - python BarcodeRabbit<ver>.py  Data=<file,file..>   <Threads=# [e.g.8]>  
##  - BarcodeRabbit uses an ancillary program (VSG_ModuleFV.py).  This should be placed in the same directory that BarcodeRabbit is run from
##  - Two other files are needed (also in the same directory):  
##  -  i) a Reference directory of SSU rRNA 'barcode' sequences.  Currently this is derived from Silva and is 'Rabbit_SSU_Barcodes_aa00_012225.fasta.gz'
##  -    Other datasbases (ncbi, greengenes2, pr2) may have some advantages and formatted versions of these will be available 'soon'
##  -  ii) a list of various illumina linkers 'illuminatetritis1223wMultiN.single.fa'
##  - The 'Threads' variable is optional to speed up analysis of large numbers of files (without providing this, BarcodeRabbit will run with a single core)
##  -    Example:
##  -          pypy3 BarcodeRabbit_ba04_012825.py data=''/Users/rabbit/media/datadrive/**/*.fastq.gz'
##  -      This command will run BarcodeRabbit and analyze all *.fastq.gz files in the indicated directory and any subdirectories
##  -      And you'll need these four files in the directory where BarcodeRabbit is being run from (and where the output will appear)
##  -        'BarcodeRabbit_bb01_040625.py'
##  -        'VSG_ModuleFV.py'
##  -        'illuminatetritis1223wMultiN.single.fa'
##  -        'Rabbit_SSU_Barcodes_bb01_040625.fasta.gz'  
## ->Output
##  - The outputs of RabbitSSUBarcodeCounter are spreadsheets (.tsv files) with a list of data for barcode sequences that may be either
##  -   Barcode_segments from known SSU genes (with some identity information to indicate what species could have provided these
##  -   Barcode_segment-like regions that don't match any known SSU segment (as initially assessed using the Silva database used as input
##  - Many of the "novel" matches are actually not indicative of new species or genera.  Some examples
##  -   SSU rRNA segments that were not in the Silva database but are in other databases such as NCBI nr
##  -   Various chimeras, sequencing artefacts, unannotated variants, and sequencing error/noise
##  - But some matches may be novel SSU genes not previously described
##  - By default, two .tsv files are produced. One includes all SSU rDNA barcodes (novel or not) and one includes strong candidates for
##      novelty based on the provided database. Note again that the majority in this class may still be 'known' rDNAs, as a blast search of nr at NCBI will show.
##  - Output Columns:
##  -    0 Barcode: the sequence of the (novel or known) barcode segment
##  -    1 Count: number of occurences in all datasets
##  -    2 SampleCounts: number of occurences by sample name
##  -    3 DataFileCounts: number of occurences by data file name
##  -    4 CoreSeq: sequence of the core segment used to originally nominate the barcode
##  -    5 CoreScore: mismatches between the core in this barcode and the canonical core (crw0)
##  -    6 UpsLenCount: the number of different upstream lengths observed in all data.  If this number is small (1 or 2, there is limited diversity in the observed reads, which may be indicative of certain artefacts)
##  -    7 maxDwnLen: the longest downstream segment in any individual read.  If this is small, the segment is always near the 3' end of the observed DNA/RNA fragment, somewhat raising concerns of a linker chimera.
##  -    8 Whole_BestMatchScore: Differences between this barcode and the closest match in the reference dataset (Silva or other database).  Single base snps contribute 1 to the score, indels 2
##  -    9 Whole_BestMatchSeq: The sequence of the closest matching reference sequence
##  -   10 Whole_BestMatchTaxa: This is just one example of a reference sequence amongst the 'best match' category (there may be others also)
##  -   11 Ups_BestMatchScore: Best Match score for the sequence of length sups1 upstream of the core.  Small values indicate likely close relationship to a known SSU; a low score (e.g. 0 or 1) here with a higher score in the WholeBestMatch score suggests a chimeric or truncated read, rather than a truly novel SSU segment
##  -   12 Dwn_BestMatchScore: Best Match score for the sequence of length sdwn1 downstream of the core.  Small values indicate likely close relationship to a known SSU; a low score (e.g. 0 or 1) here with a higher score in the WholeBestMatch score suggests a chimeric or truncated read, rather than a truly novel SSU segment
##  -   13 U2_BestMatchScore: Best Match score for a relatively short sequence (sups2) just upstream of the core.  A perfect or near perfect match here with poor matche by the full upstream sequence (col 11) indicates possible chimerism or artefact within the upstream flank.
##  -   14 D2_BestMatchScore: Best Match score for a relatively short sequence (sdwn2) just downstream of the core.  A perfect or near perfect match here with poor matche by the full upstream sequence (col 11) indicates possible chimerism or artefact within the upstream flank.
##  -   15 Conserved_9_Score: A metric indicating the number of bases in the core that match a consensus of 9 residues that are generally conserved between Eukaryotes, Bacteria, Archae, Mitochondria, and Chloroplasts.  A low score here (<3) is strongly indicative of potential SSU rRNA character, while larger values may indicate either a non-rRNA locus that just happened to have the core or a very highly diverged rRNA SSU segment.
##  -   16 Clamp_Affinity_Metric: A metric indicating the quality of any potential duplex just outside the conserved 20-mer that "clamps" this sequence in the ribosome structure.  Higher values (e.g.>20) seem strongly indicative of true ribosomal origin.
##  -   17 WashUSS_Provisional: A diagram of possible secondary structure for a clamp, using the WashUSS conventions (see documentation for Infernal for definitions).
##  -   18 Clamp_Percentile: A percentile of clamp quality (energy) amongst a set of sequences where bases outside the core-20 mer are randomly shuffled (default 1000 different shuffles) to yield a bootstrap comparision (higher number indicates greater significance of clamp structure).
##  -   19 High_Novelty: A boolean flag indicating which sequences match the variety of flags set up to detect novelty amongst barcodes.
##  -   20 Bard_Mnemonic: This is an experimental mnemoic derived from teh sequence and energy properties of each barcode.  In some cases, similar barcodes will have similarities in their mnemonic. Words are from W.Shakespeare with some editing to attempt to avoid profanity.
##  -   21 AssemblyLocal: A very naive local assembly that may be of use in a blast search or other analysis.  Not of quality expected for a program like Spades or Megahit-- this is merely to facilitate interim assessment.
##
## ->Required Parameters
## ->Required Parameter 1: A Reference file with one or many 'reference' sequences in fasta format
##  - Reference=<file[s]> : Datasets to gather kmers for analyzing data
##  -   File lists can be comma separated (no spaces) and/or specified with a wildcard (*) [wildcards require quotes around name]
##  -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.
## ->Required Parameter 2: Data Files and Sample Names
##  - Direct FileName Option: Data=<fasta/q file[s]> : Here you input one or more sequencing data files for analysis
##  -   This can be a single file or a comma-delimited list and can include wildcards (*)
##  -   "Data" can also be a directory (in this case, all fasta/fastq/fastq.gz/fastq.gz in the directory will be analyzed).
##  -   By default the sample name will be the filename and Rabbit "AutoPair" R1/R2 data file pairs as a single ('Rx') entity   
##  - SampleToData MetaFile Option: For more flexibility, provide SampleName/ReadFileNames combinations in your own file
##  -   This is invoked with SampleToData=<your SampleDescriptorFile> instead of Data= in command line
##  -   SampleDescriptor files have simple line-by-line format where each line consists of
##  -   A sample name, followed by tab, followed bya comma-delimited list of relevant the file paths with data from that sample
##
## ->Optional Parameters
##  - Although several different characteristics of the input handling and output can be varied, the program is effectively tuned with
##  - a single set of parameters.  Detailed parameter lists are declared and defined in the first non-commented section of the code.
##  - NoNs=<False> : Setting this to true ignores any read with an N
##  - OutFile=<file name assigned by program if not specified> : Allows user to specify base filenames for data output
##  -   For paired data output, separate R1 and R2 file names with a semicolon, e.g. OutFile='myoutfile_R1.fastq;myoutfile_R2.fastq'
##  - ReportGranularity=<100000> : How often to report progress (default every 100000 reads)
##  - MaxReadsPerFile=<0> : Max number to check for each substrate file (default is zero (all reads)
##  - GzipOut = <True/False>: Gzips the output (default is false, can also be set by providing an output file with extension '.gz')
##  - Threads = <number>: Running with multiple threads may be faster on a system where disk access is not limiting (default=1)
## ->Additional Feature Options
##  - ***Pair/Trimming Function-- Options to pair trim Reads of linker sequences
##  -   Two types of trimming are used
##  -    (a) To avoid reading into linkers, any read that contains the complement of the first k-mer (of len TrimFk1) from its paired read will be trimmed to remove all sequences after then end of that k-mer
##  -    (b) Any read that contains a k-mer of len TrimTk1 matching a linker will be trimmed to remove all sequences from the start of that k-mer
##  -   TrimFk=<16> : sets a value for k-mer length to avoid paired reads that extend into linker
##  -   TrimTk=<13>: sets a value for k-mer to remove tetritis
##  -   TrimTf=<illuminatetritis1223.fa> provides a FastA file name to extract k-mers for tetritis removal. (a file  that can be used for this is available on the GITHub site)
## ->Other Notes
##  - All  features should work on Linux/Mac/Windows.
##  - Requires Python 3.7+.
##  - The pypy interpreter (www.pypy.org) will speed up the program many-fold.
## ###############
## End Help

OriginalReferenceFile1 = ['ncbi_rDNADump.fa.gz','SILVA_138.1_SSURef_tax_silva.fasta.gz','pr2_rDNADump.fa.gz','gbnonbact.fa.gz','genbank_rDNADump.fa.gz','gg2_rDNADump.fa.gz','allgbRDNAI.fasta.gz']
                          ##(originalreference=) Original file to build reference Barcode list from (only needed if pre-compiled barcode file is not available)
ParsedReferenceFile1 = 'Rabbit_SSU_Barcodes_bb01_040625.fasta.gz' ##(reference=)
RecompileReference1 = False ## Will overwrite the existing reference file with a newly compiled one
DataFiles1 = './' ## (data=)
Threads1  = 1  # Multithreading is enabled if Threads>1.  Generally setting this to a value close to the number of cores to be used will be optimal for high-throughput applications
myThread1 = -1 # Internal variable for multithreading
## The parameters below define the core motif searched for and used to generate barcodes
##  Generally this program searches for a defined core motif (e.g. 20bp) that is shared between very divergent species or isoforms
##  Then uses flanking DNA/RNA on each side to generate a barcode that can be used for classification and record keeping
##    The general structure of the barcode is
##        ----upstream--- ----core----  ----downstream---
##        UUUUUUUUUUUUUUU CCCCCCCCCCCC  DDDDDDDDDDDDDDDDD
##    The program is optimized for a single highly-conserved small subunity rRNA motif 'GTGCCAGCAGCCGCGGTAAT' that is present in all known kingdoms of life
##    But can easily be used for other motifs
crw0 = 'GTGCCAGCAGCCGCGGTAAT'  ## Core motif.  Default:'GTGCCAGCAGCCGCGGTAAT', present in the vast majority of SSU rRNAs
maxRefCoreHam1 = 4             ## Maximum hamming distance used to identify motif in reference database /0 implies only reference genes with a perfect match to the input motif will be considered) 
maxSampleCoreHam1 = 2          ## Maximum hamming distance to Core motif in sample datasets (0 implies a perfect match search in the FastQ sample files) 
klen2 = 10                     ## Kmer length used in filtering sample data (default:10).  Affects execution speed but not outcome
sups1 = 22                     ## length of sequence upstream of the core motif to be used in generating the barcode
sdwn1 = 22                     ## length of sequence downstream of the core motif to be used in generating the barcode
sups2 = 12                     ## length of sequence used to identify sequences that match precisely at 5' end of core
sdwn2 = 12                     ## length of sequence used to identify sequences that match precisely at 3' end of core
##   Other input processing features
ReportCadence1 = 1000000       ## How frequently to report progress to user
MaxReadsPerFile1 = 0           ## Set to a positive number to limit the number of reads per data file (sample)
DataIsAlignment1 = False      ## set to true to strip '-' characters in input
## Settings to avoid trim sequences that go beyond R1/R2 overlap or match a preset file of sequencing adaptors
TrimFk1 = 16     ##  sets a value for k-mer length to avoid paired reads that extend into linker
TrimTk1 = 13     ##  sets a value for k-mer to avoid tetritis
TrimTf1 = 'illuminatetritis1223wMultiN.single.fa'     ##  provides a FastA file to extract k-mers for tetritis removal

## Seed sequences to speed up Tax identification: These should generally not need to be changed for a simple search but are listed here since they are
##  technically adjustable by the user in (at present unforseen) circumstances
RawAnchors1 =('AGAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTGTCC',
              'GAATAAGGGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAATACCGGCAGCTCGAGTGGTGGCC',
              'GAAAAAGCCCCGGCTAACTTCGGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTATTC',
              'AGCAATTGGAGGGCAAATTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATT')
## These are a set of anchor sequences used to speed up nearest neighbor searches.  They are a (roughly) arbitrary set of 64b barcodes from
## a survey of Silva138.  Speed might be improved by changing these sequences to ones that are more diverse/representative, but for the moment they seem to work.
inD1 = {'wm':4,'um':2,'dm':2,'um2':1,'dm2':1}  ## allowed indels for fragment aligment
ReportingThresholds1 = {'wm':0,'um':0,'dm':0,'um2':0,'dm2':0} ## Only barcodes with mismatch values at least as large as these will be reported (default is zeros-report everything)
AssemblyThresholds1 = {'wm':0,'um':2,'dm':2,'um2':0,'dm2':0}  ## Only barcodes with mismatch values at least as large as these will be assembled
AssembleAll1 = False ## setting this to true overrides AssemblyThresholds and does the limited local assembly on all barcodes
HNT1 = {'wm':7,'um':3,'dm':3,'um2':1,'dm2':1, 'core':2, 'UpsLenCount':2,'maxDwnLen':6,'conserved-9':7, 'clampE':10}  ## (highnoveltythresholds=)Only barcodes with values that match these (<= for 'core', 'UpsLenCount', 'maxDwnLen', 'conserved-9', >= for others) will be chosen as "HighNovelty"
MinTotalReadCount1 = 3  ## This is the minimum total read count (all samples) below which an individual barcode will not be reported
HighlyNovelOut1 = True ## A separate output of the highly novel sequences
## Five different matches are done with each potential barcode, with mismatch values to the "best match" in the SSU database provided for each
##   'wm' represents the entire sequence -- so a mismatch value of 5 in a 64 base barcode represents a situation with 5 mismatches in the entire sequence
##   'um' and 'dm' represent the upstream and downstream sequences flanking the core respecitively (lengths sups and sdwn).  So if sups and sdwn are each 22nt, these will represent the first and last 22nt of the barcode
##   'um2' and 'dm2' represent a shorter sequence immediately upstream and downstream of the core, allowing some chimeric artefacts to be detected
##        defaults are 12 for each of these.  um2 and/or dm2 values of zero even with larger mismatch values in um or dm raise suspicion that a given barcode is not novel but rather simply fused to an external unrelated sequence
##  As of version bc, these values have been reset based on a random set of barcodes generated and their best matches.  Current code requires a minimum
##  of 3base mismatch in each 22b flanking region and at least one mismatch in each 12b flanking region.  This retained all but
##  one of 100 random barcodes that were constrained (80%) in the consensus 9 base positions.

## Settings for output
OutFileName1 = 'default' ## (out=) File Name for output data
CatchReads1 = False      ## Setting this to true will catch all motif-matching reads into a FastA or FastQ file, with the barcode prepended to the description line to allow later extraction and more complete assembly
Delimiter1 = '\n'              ## line delimiter for tabular output files (default:'\n')

## Settings for clamp folding assessment
CoreTemp1 = 0 ## Equivalent core temperature of ribosome decoding center for deltaG calculation
TryRandom1 = 1000 ## Try a set of sequences shuffled in each core flanking region and get a range of values (bootstrap test)

import sys,os,gzip,pickle
from operator import ne
from itertools import combinations,zip_longest, product
from array import array
from glob import glob
from time import sleep
from random import randint, shuffle

try:
    from VSG_ModuleFV import *
except:
    if myThread1==-1: 
        print()
        print('********** ERROR: Failed to Import Module VSG_ModuleFV.py                          ***********')
        print('********** 1. Get VSG_ModuleFV script (from GitHub.com/FireLabSoftware/CountRabbit)  ***********')
        print('********** 2. Ensure VSG_Module#.py is in the same directory as CountRabbit#.py    ***********')
        print()
vcommand()
if myThread1==-1:
    print()
    print('******************')
    print('******************')
    print('Alpha Version of RabbitSSU_BarcodeCounter-- this version for testing only (no guarantees of anything working or being correct/useful')
    print('If more current versions of CountRabbit are available they will be at: github.com/FireLabSoftware/BarcodeRabbit ')
    print('******************')
    print('******************')
    print()
def pypyreminder1(context):
    if not('pypy' in sys.version.lower()) and myThread1==-1:
        print()
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        if context=='before':
            print('WARNING: BarcodeRabbit will be running in standard Python, not PyPy')
        else:
            print('WARNING: BarcodeRabbit has just run in standard Python, not PyPy')
        print('BarcodeRabbit will be very slow with standard Python (up to 10x slower than PyPy)')
        print('Download current native PyPy from PyPy.org for full-speed function')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print()
pypyreminder1('before')

myvnow1 = vnow.replace('_','').replace('T','_').replace('D','_')

MetaSampleList1 = vMetaFileLists(DataFiles1)
def myDataStr1(s):
    try:
        if type(s)==str:
            myN = ''.join(filter(lambda x:str.isalnum(x) or x in '. _-=',os.path.basename(s.split('*')[0].rstrip('/]').split('[')[-1]).rsplit('.',1)[0].replace('?','x')))
        else:
            myN = ''
            for y in s:
                myN += ''.join(filter(lambda x:str.isalnum(x) or x in '. _-=',os.path.basename(str(DataFiles1).rstrip('/').split('*')[0]).rsplit('.',1)[0].replace('?','x')))
        return myN.split('.fast')[0][:100]
    except:
        return 'Rabbit'

if myThread1==-1 and not(MetaSampleList1):
    print()
    print('************************')
    vLog('   Empty Input Error: ****No input found found for specification',DataFiles1,vFileOut=True)
    raise ValueError('Suggest respecifying file list')
if Threads1>len(MetaSampleList1):
    Threads1 = len(MetaSampleList1)
if Threads1>1:
    MetaSampleList1 = MetaSampleList1[myThread1-1::Threads1]

if OutFileName1 == 'default':
    OutFileName1 = myDataStr1(DataFiles1)+'_RabbitBarcodes'+myvnow1
if myThread1==-1 and not(OutFileName1.endswith('.tsv')):
    OutFileName1 += '.tsv'

klen1 = len(crw0)
mask1 = 4**(klen1-1)-1
ksam1 = 4**(klen1-1)
BaseD1 = Counter({'G':0, 'A':1, 'T':2, 'C':3, 'N':0, 'g':0, 'a':1, 't':2, 'c':3, 'n':0, 'U':2, 'u':2})
aBaseL1 = ['G','A','T','C']
BaseL1 = [0]*256
BaseA1 = [0]*256
for b1 in BaseD1:
    BaseL1[ord(b1)] = BaseD1[b1]
    BaseA1[ord(b1)] = (3-BaseD1[b1])*ksam1
exp4 = [4**z for z in range(klen1)]
mask2 = 4**klen2-1
klen11 = klen1-1
klen0 = klen1+sups1+sdwn1      ## length of full barcode

AnchorD1 = {'wm':[],'um':[],'dm':[],'um2':[],'dm2':[]}
for a1 in RawAnchors1:
    AnchorD1['wm'].append(a1)
    AnchorD1['um'].append(a1[:sups1])
    AnchorD1['dm'].append(a1[-sdwn1:])
    AnchorD1['um2'].append(a1[sups1-sups2:sups1])
    AnchorD1['dm2'].append(a1[-sdwn1:-sdwn1+sdwn2])

def CounterToText1(C):
    myT = ''
    for k in C:
        myT += k+':'+str(C[k])+', '
    return myT[:-2]
def TextToCounter1(T):
    myC = Counter()
    for I in T.split(', '):
        k,n = I.split(':')
        myC[k] = int(n)
    return myC
## Keep ValueToSeq and SeqToValue and dantisense for debugging even if not used in code
def ValueToSeq1(v,l):
    ''' converts numerical representation back to sequence v is value, l is k-mer length'''
    return ''.join([aBaseL1[(v>>(2*i)) & 3] for i in range(l-1,-1,-1)])
def SeqToValue1(s):
    return sum([exp4[i]*BaseL1[ord(c)] for (i,c) in enumerate(s[::-1])])
def dantisense1(v,klen):
    w = 0
    for i in range(klen):
        w += (3-(v>>(i*2))&3)<<(2*(klen-1-i))
    return w
def CleanAlign1(s):
    return s.replace('U','T').replace('-','')
def ham1(a,b):
    return sum(map(ne,a,b))  ## quick hamming distance (with penalties for any mismatch in length

def nw1(s1,s2,f=50,span=4):  ## s1 and s2 are sequence strings or byte arrays, q1 and q2 are quality score Phred Arrays.  No case checking done so all upper or lower to shart
    '''Alignment Quality Check-- input is two sequences plus a distance limit 'span' distance, output is a tuple of score, match number, mismatch number, gap number, total gap length'''
    ## The alignment differs in a few ways from standard Smith=Waterman/Needleman-Wunsch.
    ## First, only base positions within a range of +/- span from each other are considered for alignment
    ## so for span=1, a base is only considered for alignment to the previous or next base.  Thus any gaps can't be longer than length span
    ## and net gap length is similarly limited
    ## PyPy runs this code really really fast (20x faster than standard CPython).
    ## Variables
    ##    s1 is an input reference sequence
    ##    s2 is a proband sequence
    ##    f is a maximum distance; nw1 will break and report back f if the distance gets substantially avove this (avoids wasting time on poorly matching sequences)
    ##    span is how distant matching bases are allowed in the sequences (net indels greater than this will be misaligned)    
    ##    ll1, ll2 are the lengths of l1 and l2 (l11= l1-1)
    ##    internal algorithem variables:
    ##    m is an array with distance metrics for each tested combination of sequence positions
    ##    The first span positions are in a -1 row, corresponding to a boundary value
    ##    i1,i2 are the sequence positions in S1,S2 for the current "focus" of the cellular automaton
    ##    d = i2=i1 and must be in range -span to span (inclusive)
    ##    Y is the width of the resulting array
    ##    x is the index in the m array.  x = Y*(1+i1)+d.  So x=0 corresponds to i1=-1,d=0 and x=Y is the first actual comparison position
    ##    e is the number of elements in the array
    ##    r will hold the distance as eventually reported to the calling program
    MMCost1 = 1
    IndelCost1 = 2
    l1 = len(s1)
    if not(l1): return 0
    l11 = l1-1
    l2 = len(s2)
    Y = 2*span+1
    e = span+1+Y*l1
    m = [0]+[99999999999]*(e-1) ## using a large positive value as a default.  Note that values above 2**32 here enforce use of 64bit data structure for python and speed program due to alignment with word boundaries (with little memoryh cost)
    i1 = -1
    d = 0
    x = 0
    r = 99999999999
    MaxRowScore1 = f+Y*IndelCost1
    while x<e:  
        if x>=Y and l2>i1+d>=0:
            m[x] = m[x-Y]+(s1[i1]!=s2[i1+d])*MMCost1
        if d>-span:
            m[x] = min(m[x-1]+IndelCost1,m[x])
        if d<span:
            m[x] = min(m[x-Y+1]+IndelCost1,m[x])      
        else:
            d = -span-1
            i1 += 1
            if m[x]>=MaxRowScore1:
                return f
        if i1==l11 and m[x]<r:
            r = m[x]
        d += 1
        x += 1
    return r
    
myC2 = set((('G','C'),('C','G'),('G','U'),('U','G'),('A','U'),('U','A')))
def iC(a,b):
    if (a,b) in myC2: return 1
    return 0
def vC(a,b):  ## value of base match in bits
    if (a,b) in (('G','U'),('U','G')): return 1
    elif (a,b) in myC2: return 2
    else: return 0  #maybe this should be negative, but keeping 1 for now
def wMatch1(sa1,sa2):
    for ba1,ba2 in zip(sa1,sa2[::-1]):
        if not(iC(ba1,ba2)):
            return False
    return True
# # Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects
# # Jeffrey Zuber , Susan J Schroeder , Hongying Sun , Douglas H Turner , David H Mathews 
# # NAR 50:9 (2022) 5251–5262, https://doi.org/10.1093/nar/gkac261'''

Zuber22Free37 = Counter({('GC','CG'):-3.46, ('CC','GG'):-3.28, ('GA','CU'):-2.42, ('CG','GC'):-2.33, ('AC','UG'):-2.25,
           ('CA','GU'):-2.07, ('AG','UC'):-2.01, ('UA','AU'):-1.29, ('AU','UA'):-1.09, ('AA','UU'):-0.94,
           ('GC','UG'):-2.23, ('CU','GG'):-1.93, ('GG','CU'):-1.80, ('CG','GU'):-1.05, ('AU','UG'):-0.76,
           ('GA','UU'):-0.60, ('UG','GU'):-0.38, ('UA','GU'):-0.22, ('GG','UU'):-0.20, ('GU','UG'):-0.19,
           ('AG','UU'):+0.02})
Zuber22Entropy = Counter({('GC','CG'):-42.13, ('CC','GG'):-34.41, ('GA','CU'):-36.53, ('CG','GC'):-23.46, ('AC','UG'):-31.37,
                          ('CA','GU'):-27.08, ('AG','UC'):-23.66, ('UA','AU'):-25.40, ('AU','UA'):-25.22, ('AA','UU'):-20.98,
                          ('GC','UG'):-40.32, ('CU','GG'):-23.64, ('GG','CU'):-34.23, ('CG','GU'):-14.83, ('AU','UG'):-27.32,
                          ('GA','UU'):-32.19, ('UG','GU'):-27.04, ('UA','GU'):-8.08, ('GG','UU'):-28.57, ('GU','UG'):-24.11,
                          ('AG','UU'):-16.53})
Zuber22I = Counter()
for x1 in Zuber22Free37:
    Zuber22I[x1] = Zuber22Free37[x1]-(CoreTemp1-37)*Zuber22Entropy[x1]/1000.0
for d1,d2 in list(Zuber22I.keys()):
    Zuber22I[(d2[::-1],d1[::-1])] = Zuber22I[(d1,d2)]
Zuber22E = Counter()
Zuber22S = Counter()
Init22 = 4.1 - 1.78*(CoreTemp1-37)/1000
AUeAU22 = 0.22 - 13.35*(CoreTemp1-37)/1000
AUeCG22 = 0.44 - 8.79*(CoreTemp1-37)/1000
AUeGU22 = -0.71 - 18.96*(CoreTemp1-37)/1000
GUeCG22 = 0.13 - 12.17*(CoreTemp1-37)/1000
GUeAU22 = -0.31 - 12.78*(CoreTemp1-37)/1000
GUeGU22 = -0.74 - 22.47*(CoreTemp1-37)/1000
esym1 = 0.43 - 1.38*(CoreTemp1-37)/1000
for d0f,d0s,d1f,d1s  in product('AGCU',repeat=4):
    if d0f+d1f in ('AU','UA'):
        if d0s+d1s in ('AU','UA'): Zuber22S[(d0f+d0s,d1f+d1s)] = AUeAU22
        if d0s+d1s in ('GC','CG'): Zuber22S[(d0f+d0s,d1f+d1s)] = AUeCG22
        if d0s+d1s in ('GU','UG'): Zuber22S[(d0f+d0s,d1f+d1s)] = AUeGU22
    if d0f+d1f in ('GU','UG'):
        if d0s+d1s in ('AU','UA'): Zuber22S[(d0f+d0s,d1f+d1s)] = GUeAU22
        if d0s+d1s in ('GC','CG'): Zuber22S[(d0f+d0s,d1f+d1s)] = GUeCG22
        if d0s+d1s in ('GU','UG'): Zuber22S[(d0f+d0s,d1f+d1s)] = GUeGU22
    Zuber22E[(d1s+d1f,d0s+d0f)] = Zuber22S[(d0f+d0s,d1f+d1s)]
def mydelG1(w,c,NoInit=False):
    '''determine deltaG provisional for a colinear helix with strands w and c (equal length) w is 5'->3', c is 3'->5'; min length=3'''
    lw = len(w)
    if lw!=len(c) or lw<3: return 0.0
    if NoInit:
        e = 0
    else:
        e = Init22
    for p in range(lw-1):
        mypair1 = w[p:p+2],c[p:p+2]
        e += Zuber22I[mypair1]
        if p==0:
            e += Zuber22S[mypair1]
        if p==lw-2:
            e += Zuber22E[mypair1]
    if w == c[::-1]:
        e += esym1
    return -e
def shuffleMe1(s1):
    x1 = list(s1[:22])
    x2 = list(s1[22:42])
    x3 = list(s1[42:])
    shuffle(x1)
    shuffle(x3)
    return(''.join(x1+x2+x3))
def showMe1(bL,i,n,sc):
    for myB in bL.values():
        if myB.c:
            vtext(text=myB.b,xc=myB.x*18,yc=-myB.y*18-i*136,color=myColor1[myB.c],font='DejaVuMonoBold 18')
            if myB.i:
                vtext(text='-',xc=myB.x*18,yc=-(3-myB.y)*18-i*136,color='black',font='DejaVuMonoBold 18')
        vtext(text='novelC_'+str(i+1)+"  ('"+n.split('$',1)[0][:1000].replace('|','_')+"') ;  score="+str(sc),x2=11*18,yc=-4*18-i*136,color='black',font='TimesNewRoman 14')
def myPercentile1(L,x):
    myL = sorted(L+[x,])
    return 100*myL.index(x)/len(myL)
pnPositions1 = (12,11,13,10,14,9,15,8,16,7,17,6,18,5)

def MakeScore1(s0):
    score1 = []
    s1 = s0.upper().replace('T','U')
    if s0[22:26] == 'UGCC': s1 = 'X'+s1[:25]+s0[26:]
    pn = s1[31:34]
    cAvail = Counter(list(range(5,20))) ## positions available for pseudoknotpairing
    v = 24
    w = 41
    cQ = Counter()
    for dv in range(4):
        for dw in range(4):
            dd = 0
            while v-dv-dd>=0 and w+dw+dd<64 and vC(s1[v-dv-dd],s1[w+dw+dd]): dd += 1
            if dd>=3: cQ[(dv,dw,dd)] = mydelG1(s1[v-dv-dd+1:v-dv+1],s1[w+dw:w+dw+dd][::-1])
    if cQ:            
        mydv,mydw,mydd = cQ.most_common()[0][0]
        score1.append(cQ.most_common()[0][1])
        v -= mydv; w += mydw
        for x1 in range(mydd):
            if not(s1[v]+s1[w] in ('GU','UG')): cAvail[v] = 0
            v-=1; w+=1
    cQ = Counter()
    for dv in range(4,11):
        for dw in range(4):
            dd = 0
            while v-dv-dd>=0 and w+dw+dd<64 and vC(s1[v-dv-dd],s1[w+dw+dd]): dd += 1
            if dd>=3: cQ[(dv,dw,dd)] = mydelG1(s1[v-dv-dd+1:v-dv+1],s1[w+dw:w+dw+dd][::-1])
    if cQ:
        mydv,mydw,mydd = cQ.most_common()[0][0]
        score1.append(cQ.most_common()[0][1])
        v -= mydv; w += mydw
        for x1 in range(mydd):
            if not(s1[v]+s1[w] in ('GU','UG')): cAvail[v] = 0
            v-=1; w+=1
    cQ = Counter()
    for pnv in pnPositions1:
        if wMatch1(pn,s1[pnv:pnv+3]) and cAvail[pnv]+cAvail[pnv+1]+cAvail[pnv+2]==3: ##bList[myv].y*bList[myv-1].y*bList[myv-2].y==0:
            cQ[pnv] = mydelG1(pn[::-1],s1[pnv:pnv+3])
    if cQ:
        score1.append(max(0,cQ.most_common()[0][1]))
    return sum(score1)
def MakeWU1(s0):
    wu = [':',]*64
    score1 = []
    s1 = s0.upper().replace('T','U')
    if s0[22:26] == 'UGCC': s1 = 'X'+s1[:25]+s0[26:]
    pn = s1[31:34]
    cAvail = Counter(list(range(5,20))) ## positions available for pseudoknotpairing
    v = 31
    w = 32
    for dx in range(2):
        wu[v] = '_'; wu[w] = '_'
        v-=1; w+=1
    if vC(s1[v],s1[w]) and vC(s1[v-1],s1[w+1]):
        h2 = True
        for dx in range(2):
            wu[v] = '<'; wu[w] = '>'
            v-=1; w+=1
    else:
        h2 = False
        for dx in range(2):
            wu[v] = '_'; wu[w] = '_'
            v-=1; w+=1
    for dx in range(3):
        if h2:
            wu[v] = '-'
        else:
            wu[v] = '_'
        v-=1
    for dx in range(5):
        if h2:
            wu[w] = '-'
        else:
            wu[w] = '_'
        w+=1
    cQ = Counter()
    for dv in range(4):
        for dw in range(4):
            dd = 0
            while v-dv-dd>=0 and w+dw+dd<64 and vC(s1[v-dv-dd],s1[w+dw+dd]): dd += 1
            if dd>=3: cQ[(dv,dw,dd)] = mydelG1(s1[v-dv-dd+1:v-dv+1],s1[w+dw:w+dw+dd][::-1])
    if cQ:
        h2 = True
        mydv,mydw,mydd = cQ.most_common()[0][0]
        score1.append(cQ.most_common()[0][1])
        for dx in range(mydv):
            if h2:
                wu[v] = '-'
            else:
                wu[v] = '_'
            v-=1
        for dx in range(mydw):
            if h2:
                wu[w] = '-'
            else:
                wu[w] = '_'
            w+=1
        for x1 in range(mydd):
            if not(s1[v]+s1[w] in ('GU','UG')): cAvail[v] = 0
            wu[v] = '<'; wu[w] = '>'
            v-=1; w+=1
    cQ = Counter()
    for dv in range(4,11):
        for dw in range(4):
            dd = 0
            while v-dv-dd>=0 and w+dw+dd<64 and vC(s1[v-dv-dd],s1[w+dw+dd]): dd += 1
            if dd>=3: cQ[(dv,dw,dd)] = mydelG1(s1[v-dv-dd+1:v-dv+1],s1[w+dw:w+dw+dd][::-1])
    if cQ:
        mydv,mydw,mydd = cQ.most_common()[0][0]
        score1.append(cQ.most_common()[0][1])
        for dx in range(mydv):
            if h2:
                wu[v] = '-'
            else:
                wu[v] = '_'
            v-=1
        for dx in range(mydw):
            if h2:
                wu[w] = '-'
            else:
                wu[w] = '_'
            w+=1
        for x1 in range(mydd):
            if not(s1[v]+s1[w] in ('GU','UG')): cAvail[v] = 0
            wu[v] = '<'; wu[w] = '>'
            v-=1; w+=1
    for myv in range(v,-1,-1):
        wu[myv] = '.'
    for myw in range(w,64):
        wu[myw] = '.'
    cQ = Counter()
    for pnv in pnPositions1:
        if wMatch1(pn,s1[pnv:pnv+3]) and cAvail[pnv]+cAvail[pnv+1]+cAvail[pnv+2]==3: ##bList[myv].y*bList[myv-1].y*bList[myv-2].y==0:
            cQ[pnv] = mydelG1(pn[::-1],s1[pnv:pnv+3])
    if cQ:
        pnv,mydd = cQ.most_common()[0]
        for v in range(pnv,pnv+3):
            wu[v] = 'A'
        for w in range(31,34):
            wu[w] = 'a'
        score1.append(max(0,cQ.most_common()[0][1]))
    return sum(score1), ''.join(wu)

pnPositions1 = (12,11,13,10,14,9,15,8,16,7,17,6,18,5)
myColor1 = {1:'red', 2:'blue', 3:'green', 4:'black'}
class base():
    def __init__(self,b,x,y,c,p):
        self.b = b #base [AGCU]
        self.x = x #x position on canvas
        self.y = y #y position on canvas
        self.c = c #color for display (1=red for core, 2=blue for pseudoknot, 3=green for G530 analog, 4=otherwise black)
        self.p = p #position (0-63)
        self.i = 0 #indel?
def MakeBList1(s0):
    score1 = []
    s0 = s0.upper().replace('T','U')
    if s0[22:26] == 'UGCC':
        insC4 = True
        s1 = 'X'+s0[:25]+s0[26:]
    else:
        insC4 = False
        s1 = s0
    bList = {i:base(s1[i],9999,0,4,i) for i in range(64)}
    if insC4:
        bList[24.5] = base('C',3.2,0,1,25.5)
    bList[22].x = 0; bList[22].c = 1; bList[22].y = 0
    bList[23].x = 1; bList[23].c = 1; bList[23].y = 0
    if not(insC4):
        bList[24].x = 2; bList[24].c = 1; bList[24].y = 0
        bList[25].x = 3.5; bList[25].c = 1; bList[25].y = 0
        bList[26].x = 5.0; bList[26].c = 1; bList[26].y = 0
        bList[27].x = 6.5; bList[27].c = 1; bList[27].y = 0
    else:
        bList[24].x = 2; bList[25].c = 1; bList[25].y = 0
        bList[25].x = 4.4; bList[26].c = 1; bList[26].y = 0
        bList[26].x = 5.6; bList[26].c = 1; bList[26].y = 0
        bList[27].x = 6.8; bList[27].c = 1; bList[27].y = 0
    pairedinnerloop = iC(s1[28],s1[35])*iC(s1[29],s1[34])  ## This is 1 if both innerloop residue pairs are complementary
    bList[28].x = 8; bList[28].c = 1; bList[28].y = pairedinnerloop
    bList[29].x = 9; bList[29].c = 1; bList[29].y = pairedinnerloop
    bList[30].x = 10; bList[30].c = 1; bList[30].y = 0
    bList[31].x = 11; bList[31].c = 2; bList[31].y = 1
    bList[32].x = 11; bList[32].c = 2; bList[32].y = 2
    bList[33].x = 10; bList[33].c = 2; bList[33].y = 3
    bList[34].x = 9; bList[34].c = 1; bList[34].y = 3-pairedinnerloop
    bList[35].x = 8; bList[35].c = 1; bList[35].y = 3-pairedinnerloop
    bList[36].x = 7; bList[36].c = 1; bList[36].y = 3
    bList[37].x = 6; bList[37].c = 3; bList[37].y = 3
    bList[38].x = 5; bList[38].c = 1; bList[38].y = 3
    bList[39].x = 4; bList[39].c = 1; bList[39].y = 3
    bList[40].x = 3; bList[40].c = 1; bList[40].y = 3
    bList[41].x = 2; bList[41].c = 1; bList[41].y = 3
    bList[42].x = 1; bList[42].c = 4; bList[42].y = 3-iC(s1[23],s1[42])
    bList[43].x = 0; bList[43].c = 4; bList[43].y = 3-iC(s1[22],s1[43])
    pn = s1[31:34]
    cAvail = Counter(list(range(5,20))) ## positions available for pseudoknotpairing
    v = 24
    w = 41
    myx = 2
    cQ = Counter()
    for dv in range(4):
        for dw in range(4):
            dd = 0
            while v-dv-dd>=0 and w+dw+dd<64 and vC(s1[v-dv-dd],s1[w+dw+dd]):
                dd += 1
            if dd>=3:
                cQ[(dv,dw,dd)] = mydelG1(s1[v-dv-dd+1:v-dv+1],s1[w+dw:w+dw+dd][::-1])
    if cQ:            
        mydv,mydw,mydd = cQ.most_common()[0][0]
        score1.append(cQ.most_common()[0][1])
        myv = v; myw = w
        while v>myv-mydv or w<myw+mydw:
            if v>myv-mydv:
                bList[v].x = myx; bList[v].y = 0
                v-=1
            else:
                bList[w].i = 1
            if w<myw+mydw:
                bList[w].x = myx; bList[w].y = 3
                w+=1
            else:
                bList[v+1].i = 1
            myx-=1
        for x1 in range(mydd):
            bList[v].x = myx; bList[v].y = 1
            bList[w].x = myx; bList[w].y = 2
            if not(s1[v]+s1[w] in ('GU','UG')): cAvail[v] = 0
            v-=1; w+=1; myx-=1
    cQ = Counter()
    for dv in range(4,11):
        for dw in range(4):
            dd = 0
            while v-dv-dd>=0 and w+dw+dd<64 and vC(s1[v-dv-dd],s1[w+dw+dd]):
                dd += 1
            if dd>=3:
                cQ[(dv,dw,dd)] = mydelG1(s1[v-dv-dd+1:v-dv+1],s1[w+dw:w+dw+dd][::-1])
    myv = v; myw = w
    if cQ:            
        mydv,mydw,mydd = cQ.most_common()[0][0]
        score1.append(cQ.most_common()[0][1])        
        while v>myv-mydv or w<myw+mydw:
            if v>myv-mydv:
                bList[v].x = myx; bList[v].y = 0
                v-=1
            else:
                bList[w].i = 1
            if w<myw+mydw:
                bList[w].x = myx; bList[w].y = 3
                w+=1
            else:
                bList[v+1].i = 1
            myx-=1
        for x1 in range(mydd):
            bList[v].x = myx; bList[v].y = 1
            bList[w].x = myx; bList[w].y = 2
            if not(s1[v]+s1[w] in ('GU','UG')): cAvail[v] = 0
            v-=1; w+=1; myx-=1
    while v>=0 or w<64:
        if v>=0:
            bList[v].x = myx; bList[v].y = 0
        if w<64:
            bList[w].x = myx; bList[w].y = 3
        v-=1; w+=1; myx-=1
    cQ = Counter()
    for pnv in pnPositions1:
        if wMatch1(pn,s1[pnv:pnv+3]) and cAvail[pnv]+cAvail[pnv+1]+cAvail[pnv+2]==3: ##bList[myv].y*bList[myv-1].y*bList[myv-2].y==0:
            cQ[pnv] = mydelG1(pn[::-1],s1[pnv:pnv+3])
    if cQ:
        myv = cQ.most_common()[0][0]                  
        bList[myv].c = 2
        bList[myv+1].c = 2
        bList[myv+2].c = 2
        score1.append(max(0,cQ.most_common()[0][1]))
    return bList, sum(score1)

BialyL1 = array('b',[0]*(2**(2*klen2)))
ElementC1 = Counter()

scD1 = []
for i1 in range(len(crw0)):
    for c2 in 'GATC':
        scD1.append((i1,c2))
crw4s = set()
for cL1 in combinations(scD1,maxRefCoreHam1):
    crw1 = list(crw0)
    for i1,c1 in cL1:
        crw1[i1] = c1    
    s1 = ''.join(crw1)
    h1 = ham1(s1,crw0)
    v1 = SeqToValue1(s1)
    a1 = dantisense1(v1,klen1)
    if ElementC1[v1] == 0:
        ElementC1[v1] = h1+1
    elif ElementC1[v1]<0:
        ElementC1[v1] = 128
    if ElementC1[a1] == 0:
        ElementC1[a1] = -h1-1
    elif ElementC1[a1]>0:
        ElementC1[v1] = 128 ##  note that 128 is a special value for a sequences matching both sense and antisense
    BialyL1[v1&mask2] = 1
    BialyL1[a1&mask2] = 1

def PhyloList1(pString,dbtype):
    if dbtype==0:
        FirstA = pString.strip().split(' ',1)[-1].split(';')[:-1]
        if FirstA and FirstA[0] and not(FirstA[0][0].isalpha()):  #Fixes some odd special cases in annotation strings
            FirstA = FirstA[1:]
        return FirstA
def LCA1(myL1,myL2):
    for i in range(min(len(myL1),len(myL2))):
        if myL1[i]!=myL2[i]:
            return myL1[:i]
    return myL1
def my2D(f):
    if f>.995:
        return ''
    else:
        return '('+('%.2f'%f)[1:]+')'
class barFrag0():
    def __init__(self,s,t,mn):
        self.s = s
        self.t = t   ##type: um, dm are full upstream and downstream sequence (lengths sups1, sdwn1 respectively
                        ## um2, dm2 are proximal upstream and downstream matches (len sups2 and sdwn2 respectively
                        ## and ws is whole sequence (len sups1+core+sdwn1
        self.pl = []  ##phylogeny list
        self.c  = 0
        self.mn = mn
    def add(self,p):
        self.pl.append(p)
        self.c += 1
    def distance(self):
        self.d = []
        for a in AnchorD1[self.t]:
            self.d.append(nw1(self.s,a,klen0,inD1[self.t]))        
    def phylo(self):
        self.phy = [] ## l
        weightC = Counter()
        maxlen = max([len(x) for x in self.pl])
        myT = 0.000000001
        for j,r in enumerate(self.pl):
            weightC[j] = 1
            if len(r)<maxlen-1:
                weightC[j] += 10
            elif 'uncultured' in r:
                weightC[j] += 100000
                for ri in r:
                    weightC[j] += abs(ham1(ri.lower(),ri)-1)
            myT += 1/weightC[j]
        for i in range(maxlen):
            myC = Counter()
            for j,r in enumerate(self.pl):
                if not(r) or len(r)<=i: continue
                myTax1 = r[i].split(':')[0]
                myC[myTax1] += 1/weightC[j]
            mcc1 = myC.most_common(1)
            self.phy.append(mcc1[0][0]+my2D(mcc1[0][1]/myT))
        if not(self.phy):
            self.phy=['root',]
class barFrag1():
    def __init__(self,s,t,d,phy):
        self.s = s
        self.t = t   
        self.d = d  
        self.phy = phy
myClosestD1 = {}
def myClosest1(s,t):
    global myClosestD1
    if (s,t) in myClosestD1:
        return myClosestD1[(s,t)]
    v = inD1[t]
    ad = [nw1(s,a,klen0,v) for a in AnchorD1[t]]
    minD = 9999
    startflag = False
    for dif  in range(klen0):
        if dif>=minD: break
        for myd in (ad[0]-dif,ad[0]+dif):
            if dif>=minD: break
            if dif==0 and startflag: continue
            startflag = True
            if myd<0 or  myd>=klen0: continue
            for B in barD1[t][myd]:
                if any([abs(d1-d2)>=minD for d1,d2 in zip(ad[1:],B.d[1:])]): continue
                myD = nw1(s,B.s,minD,v)
                if myD<minD:
                    minD = myD
                    myB = B
    myReturn = [minD,myB]
    myClosestD1[(s,t)] = myReturn
    return myReturn
if myThread1==-1 and (RecompileReference1 or not(os.path.isfile(ParsedReferenceFile1)) or os.path.getsize(ParsedReferenceFile1)<50):
    if ParsedReferenceFile1.endswith('.gz'):
        SF1 = gzip.open(ParsedReferenceFile1, mode='wt')
    else:
        SF1 = open(ParsedReferenceFile1, mode='wt')
    if type(OriginalReferenceFile1)==str:
        OriginalReferenceFile1 = [OriginalReferenceFile1,]             
    v0 = 0
    allBC0 = {} ## keys are seq,bctype pair, values are barFrag objects
    for myRefFile1 in OriginalReferenceFile1:
        barD0 = {} ## keys are seq,bctype pair, values are barFrag objects
        F1 = vOpen(myRefFile1,mode='rt',vFileOut=True)
        mn1 = myRefFile1.split('_')[0]
        myIter1 = vFastAToIter(F1)
        for q1,x1 in enumerate(iter(myIter1)):
            if not(x1): break
            (n1,s1) = x1
            if q1%10000==0: vLog('Processing Reference Sequence '+str(q1)+' from input '+myRefFile1,vFileOut=False)
            if not(s1):continue
            BestHam1 = 999
            for j1,c1 in enumerate(s1):
                v0 = ((v0&mask1)<<2)+BaseL1[ord(c1)]
                if BialyL1[v0&mask2]:
                    h1 = ElementC1[v0]
                    if h1>0:
                        BestHam1 = min(BestHam1,abs(h1)-1)
                        p1 = j1
            if BestHam1<=maxRefCoreHam1:
                if p1-klen1>=sups1 and p1+klen1+sdwn1<len(s1):
                    myBarcode1 = s1[p1-klen1-sups1+1:p1+sdwn1+1].upper().replace('U','T')
                    if myBarcode1.count('G')+myBarcode1.count('A')+myBarcode1.count('T')+myBarcode1.count('C')!=len(myBarcode1): continue
                else:
                    continue
                myL1 = PhyloList1(n1,0)
                myS1 = zip((myBarcode1,myBarcode1[:sups1],myBarcode1[-sdwn1:],myBarcode1[sups1-sups2:sups1],myBarcode1[-sdwn1:-sdwn1+sdwn2]),('wm','um','dm','um2','dm2'))
                for s1,t1 in myS1:
                    if not((s1,t1) in allBC0):
                        if not (s1,t1) in barD0:
                            barD0[(s1,t1)] = barFrag0(s1,t1,mn1)
                        barD0[(s1,t1)].add(myL1)
                        barD0[(s1,t1)].c += 1
        for v1 in barD0.values():
            if not(v1 in allBC0):
                 v1.phylo()
        allBC0.update(barD0)
    for j1,(s1,t1) in enumerate(sorted(list(allBC0.keys()),key=lambda x:(allBC0[x].t,allBC0[x].c), reverse=True)):
        B1 = allBC0[(s1,t1)]
        allBC0[(s1,t1)].distance()
        myDString1 = '.'.join([str(a) for a in B1.d])
        mn2 = B1.mn+':'
        SF1.write('>SSU_Barcode_'+t1+'_'+str(j1)+'_'+myDString1+'_ '+mn2+'|'.join(B1.phy)+Delimiter1)
        SF1.write(s1+Delimiter1)
    SF1.close()
#                  0123456789     ..horizontal index
#                  22      28 31
#                  |       |  |
# a      GGCAAGU    UGCCA--  A        ..0
#  GGGCUG       CCGG       GC G       ..1
#  CCCGGC       GGCC       CG C       ..2  
# a      -------    AUAAUGG  C        ..3
#                    |   |  |
#                    41  37 34
#                  0123456789     ..horizontal index
# Numbering is zero based for a 64b barcode
# WashUSS version
# aGGGCUGGGCAAGUCCGGUGCCAGCAGCCGCGGUAAUACCGGCGGCCCa
# :<<<<<<aaa----<<<<----<<_AAA>>------->>>>>>>>>>>:
bardwords1 = open('bardwords.txt', mode='rt').read().split()  ## indeed a list of non obscene words from the bard
d7 = {'G':0,'A':1,'T':2,'C':3}
for i in range(7): d7[i] = 4**i
def SevenBToWord1(s7):
    return bardwords1[sum([d7[c]*d7[i] for i,c in enumerate(s7)])]
class Barcode1():
    def __init__(self,seq):
        self.s = seq ## Barcode
        self.c = 0   ## Count
        self.u = Counter() ## All upstream linked segments (segments from matched reads that go beyond the barcode)
        self.d = Counter() ## All downstream linked segments
        self.ur = Counter() ## All paired reads upstream (reverse comp)
        self.dr = Counter() ## All paired reads downstream (reverse comp)
        self.DataFiles = Counter() ## keys are data file names, values are counts
        self.Samples = Counter() ## keys are sample names, values are counts
        self.CoreSeq = seq[sups1:sups1+klen1]
        self.CoreScore = ham1(self.CoreSeq,crw0)    
    def add(self, other):
        self.c+=other.c; self.ur+=other.ur; self.dr+=other.dr; self.u+=other.u; self.d+=other.d
        self.DataFiles+=other.DataFiles; self.Samples+=other.Samples
    def register(self,ups,dwn,ori,paired,sample,file):
        self.u[ups]+=1
        self.d[dwn]+=1
        self.c += 1
        if ori<0:
            self.ur[paired] += 1
        else:
            self.dr[vantisense(paired)] += 1
        self.Samples[sample] += 1
        self.DataFiles[file] += 1
    def assemble(self):
        ## this is a veryf simple local assembler-- it is not spades or megahit.  It is designed to give a local assembly for relatively
        ## simple (non redundant and low-ish error) sequence groups.  The kmer length (sck) can be adjusted if needed.  
        sck = 12
        taperlen = 6
        dC = [Counter(),Counter()]
        dseqp = self.s[-sck:]
        dseqpC = Counter()
        for x in self.d:
            dseqpC[dseqp+x] = 2*self.d[x]
        useqp = self.s[:sck]
        useqpC = Counter()
        for x in self.u:
            useqpC[x+useqp] = 2*self.u[x]        
        xC = [self.dr+dseqpC,self.ur+useqpC]
        xS = [self.s[-sck:],self.s[:sck][::-1]]
        for j,myC in enumerate(xC):
            for d,w in myC.items():
                if j==1:
                    d=d[::-1]
                ld = len(d)
                for i in range(ld-sck):
                    myweight = w*(1+0.1*i/ld) ## slight upweighting for earlier in the read
                    if ld-sck-i<=taperlen:
                        myweight *= (ld-sck-i)/taperlen  ## downweights the last few bases proportional to distance-to-end
                    ck = d[i:i+sck]
                    if not ck in dC[j]:
                        dC[j][ck] = Counter()
                    dC[j][ck][d[i+sck]] += myweight
            while True:
                lastseq = xS[j][-sck:]
                if not(dC[j][lastseq]): break
                nextbase = dC[j][lastseq].most_common(1)[0][0]
                xS[j] += nextbase
                del(dC[j][lastseq][nextbase])
        self.assembly = xS[1][::-1][:-sck].lower()+self.s+xS[0][sck:].lower()
    def assignbardid(self):
        myletters = 'bdefhjkmpqrsvwyz'
        mynum = '23456789'
        w1 = SevenBToWord1(self.s[42:48]+self.s[17])
        w2 = SevenBToWord1(self.s[48:55])
        w3 = myletters[d7[self.s[16]]*4+d7[self.s[15]]]
        w4 = myletters[d7[self.s[55]]*4+d7[self.s[4]]]
        w5 = str(int(round(self.clampE))).zfill(2)
        self.bardid = w1+'_'+w2+'_'+w3.upper()+w4.upper()+w5
    def takestock(self):
        lenu = set(map(len,self.u))
        self.upslencount = len(lenu)
        lend = set(map(len,self.d))
        self.maxdwnlen = max(lend)
        self.conserved9 = ham1(self.s[5]+self.s[12]+self.s[13]+self.s[14]+self.s[16]+self.s[43]+self.s[54]+self.s[55]+self.s[57],'AGGCACAGG')
        self.wm = myClosest1(self.s,'wm')
        self.um = myClosest1(self.s[:sups1],'um')
        self.um2 = myClosest1(self.s[sups1-sups2:sups1],'um2')
        self.dm = myClosest1(self.s[-sdwn1:],'dm')
        self.dm2 = myClosest1(self.s[-sdwn1:-sdwn1+sdwn2],'dm2')
        self.clampE,self.WashUSS = MakeWU1(self.s)
        self.clampEPercentile = -1
        self.assignbardid()
        if (self.wm[0]>=HNT1['wm'] and self.um[0]>=HNT1['um'] and self.dm[0]>=HNT1['dm']
            and self.um2[0]>=HNT1['um2'] and self.dm2[0]>=HNT1['dm2']
            and self.CoreScore<=HNT1['core'] and self.upslencount>=HNT1['UpsLenCount']
            and self.maxdwnlen>=HNT1['maxDwnLen'] and self.conserved9<=HNT1['conserved-9']
            and self.clampE>=HNT1['clampE']):
            self.highnovelty = 1
            RandomL1 = []
            if TryRandom1:
                for i in range(TryRandom1):
                    s2 = shuffleMe1(self.s)
                    RandomL1.append(MakeScore1(s2))
                self.clampEPercentile = myPercentile1(RandomL1,self.clampE)
        else:
            self.highnovelty = 0
        self.reportme = True
        self.assembleme = True
        for myAttr in ('wm','um','dm','um2','dm2'):
            myC = getattr(self,myAttr)
            if myC[0]<ReportingThresholds1[myAttr]:
                self.reportme = False
            if myC[0]<AssemblyThresholds1[myAttr]:
                self.assembleme = False
        if self.assembleme:
            self.assemble()
        else:
            self.assembly = ' '
        


BarcodeD1 = {}
firsttest1 = sups1+klen1-1
supsMin1 = min(sups1,sdwn1)
def hasCore1(L):
    if L:
        v = 0
        for j in range(supsMin1,len(L)-1):
            v = ((v&mask1)<<2)+BaseL1[ord(L[j])]
            if (j>=firsttest1) and BialyL1[v&mask2]:
                h = ElementC1[v]
                if h and (abs(h)-1<=maxSampleCoreHam1):
                    return (j,h)
BreakHere1 = set()
for s1 in vFastAToDict(TrimTf1).values():
    a1 = vantisense(s1)
    for i1 in range(len(s1)-TrimTk1+1):
        BreakHere1.add(s1[i1:i1+TrimTk1])
        BreakHere1.add(a1[i1:i1+TrimTk1])
if Threads1<2 or myThread1>=0:
    for h11,meta1 in enumerate(MetaSampleList1):
        if len(meta1.files)==2:
            Fn1,Fn2 = meta1.files
        else:
            Fn1,Fn2 = meta1.files[0],''
        if os.path.basename(Fn1).startswith('Rabbit_SSU_Barcodes'): continue
        mySample1 = meta1.sampleName
        if Fn1.lower().endswith('.fastq') or Fn1.lower().endswith('.fastq.gz'):
            DataCadence1 = 4
        else:
            DataCadence1 = 2        
        ReportLineCadence1 = ReportCadence1*DataCadence1
        F1 = vOpen(Fn1, mode='rt',vFileOut=True)
        if Fn2:
            F2 = vOpen(Fn2, mode='rt',vFileOut=True)
            myFile1 = vFileListMnemonic((Fn1,Fn2))
        else:
            F2 = ''
            myFile1 = Fn1
        myFile11 = os.path.basename(myFile1)
        CaughtBarcodes1 = 0
        try:
            for i1,(L1,L2) in enumerate(zip_longest(F1,F2)):
                if i1%DataCadence1 == 1:
                    if i1>1 and i1%ReportLineCadence1 == 1:
                        vLog('Processing Read '+str(i1//DataCadence1)+' from '+str(myFile11)+' BarcodeInstances:'+str(CaughtBarcodes1),vFileOut=False)
                    if MaxReadsPerFile1 and i1>DataCadence1*MaxReadsPerFile1:
                        break
                    if DataIsAlignment1:
                        if L1: L1 = CleanAlign1(L1)
                        if L2: L2 = CleanAlign1(L2)
                    if not(hasCore1(L1) or hasCore1(L2)):
                        continue
                    if not(L1): L1=''
                    if not(L2): L2=''
                    s1 = L1.strip().upper(); s2 = L2.strip().upper()
                    if len(s1)>=TrimTk1:
                        for iA1 in range(len(s1)-TrimTk1):
                            if s1[iA1:iA1+TrimTk1] in BreakHere1:
                                s1 = s1[:iA1]
                                break
                    if len(s2)>=TrimTk1:
                        for iA1 in range(len(s2)-TrimTk1):
                            if s2[iA1:iA1+TrimTk1] in BreakHere1:
                                s2 = s2[:iA1]
                                break
                    if len(s1)>=TrimFk1 and len(s2)>=TrimFk1:
                        p1 = s1.rfind(vantisense(s2[:TrimFk1]))
                        p2 = s2.rfind(vantisense(s1[:TrimFk1]))
                        if p1>=0 or p2>=0:
                            newLen1 = max(p1+TrimFk1,p2+TrimFk1)
                            s1 = s1[:newLen1]
                            s2 = ''
                    x1 = hasCore1(s1)
                    if not(x1):
                        x1 = hasCore1(s2)
                        s1,s2 = s2,s1
                    if not(x1): continue
                    j1,h1 = x1
                    if h1<0:
                        myori1 = -1
                        if (1+j1-klen1-sdwn1<0) or (1+j1+sups1>len(s1)): continue
                        mybar1 = vantisense(s1[1+j1-klen1-sdwn1:1+j1+sups1])
                        myups1 = vantisense(s1[1+j1+sups1:])
                        mydwn1 = vantisense(s1[:1+j1-klen1-sdwn1])
                    elif h1>0:
                        if (1+j1-klen1-sups1<0) or (1+j1+sdwn1>len(s1)): continue
                        myori1 = 1
                        mybar1 = s1[1+j1-klen1-sups1:1+j1+sdwn1]
                        mydwn1 = s1[1+j1+sdwn1:]
                        myups1 = s1[:1+j1-klen1-sups1]
                    if not mybar1 in BarcodeD1:
                        BarcodeD1[mybar1] = Barcode1(mybar1)
                    BarcodeD1[mybar1].register(myups1,mydwn1,myori1,s2,mySample1,myFile11)
                    CaughtBarcodes1 += 1
            vClose(F1,vFileOut=False)
            if Fn2: vClose(F2,vFileOut=False)
            vLog('completing file ' +str(h11+1)+"/"+str(len(MetaSampleList1))+': '+ str(myFile11)+'; Reads:'+str(i1//DataCadence1)+' BarcodeInstances:'+str(CaughtBarcodes1),vFileOut=True)
        except:
            vLog('***Encountered error with file(s)',Fn1,Fn2,'skipping these files',vFileOut=True)
    if Threads1>1:
        pf = open(OutFileName1, mode='wb')
        pp = pickle.Pickler(pf)
        pp.dump(BarcodeD1)
        pf.close()
if Threads1>1 and myThread1==-1:
    AllThreads1 = []
    AllTempFiles1 = []
    for thrn1 in range(Threads1):
        ThreadUID1 = randint(100000000000,999999999999)
        thisFile1 = 'Temp'+str(ThreadUID1)+vnow+str(10000+thrn1)+'.pck'
        AllThreads1.append(subprocess.Popen([sys.executable,]+sys.argv[0:1]+['myThread='+str(thrn1+1),]+sys.argv[1:]+['out='+thisFile1,]))
        AllTempFiles1.append(thisFile1)
    for i1,(thisFile1,thisThread1) in enumerate(zip(AllTempFiles1,AllThreads1)):
        thisThread1.wait()
        pf = open(thisFile1, mode='rb')
        pp = pickle.Unpickler(pf)
        myBarcodeD1 = pp.load()        
        pf.close()
        os.remove(thisFile1)
        if i1==0:
            BarcodeD1 = myBarcodeD1
        else:
            for s00 in myBarcodeD1:
                if not(s00) in BarcodeD1:
                    BarcodeD1[s00] = myBarcodeD1[s00]
                else:
                    BarcodeD1[s00].add(myBarcodeD1[s00])
if myThread1==-1:
    barD1 = {} ## keys here are barcode type, values are hierarchical list: first value is distance to anchor0
    for t1  in ('wm','um','dm','um2','dm2'):
        barD1[t1] = [[] for i in range(100)]
    rD1 = vFastAToDict(ParsedReferenceFile1)
    for n1,s1 in rD1.items():
        t1 = n1.split('_')[2]
        d1 = list(map(int,n1.split('_')[4].split('.')))
        phy = n1.split(' ',1)[-1]
        B1 = barFrag1(s1,t1,d1,phy)
        barD1[t1][d1[0]].append(B1)
        myClosestD1[(t1,s1)] = [0,B1]

    TaskHeader1 =  '<!--Rabbit_Task_Header: '+OutFileName1+'-->'+Delimiter1
    TaskHeader1 += '<!--Command_Line: '+' '.join(sys.argv)+'-->'+Delimiter1
    TaskHeader1 += '<!--PythonVersion: '+','.join(sys.version.splitlines())+'-->'+Delimiter1
    TaskHeader1 += '<!--Rabbit_Version: '+vFileInfo(sys.argv[0])+'-->'+Delimiter1
    TaskHeader1 += '<!--RunTime: '+vnow+'-->'+Delimiter1
    TaskHeader1 += '<!--RunDirectory: '+os.getcwd()+'-->'+Delimiter1
    AbbrevHeader1 = ''.join(TaskHeader1.splitlines()[:-1])+'<!--RabbitTableHeader-->'  ##ending with ':RabbitTableHeader' identifies a line as a row of table headers
    def HeaderTranspose1(hT2):
        hT0 = '<!--\tOutput_Key\t\t-->'+Delimiter1
        hT0 += '<!--\tColumnNumber\tColumnHeader\t-->'+Delimiter1
        for iT2,nT2 in enumerate(hT2):
            nT2 = nT2.strip()
            if not(nT2.startswith('<!')):
                hT0 += '<!--\t'+str(iT2)+'\t'+nT2+'\t-->'+Delimiter1
        return hT0+'<!--\t\t\t-->'+Delimiter1
    hTX1 = ['Barcode','Count', 'SampleCounts','DataFileCounts',
                                                   'CoreSeq','CoreScore','UpsLenCount','maxDwnLen',
                                                   'Whole_BestMatchScore','Whole_BestMatchSeq','Whole_BestMatchTaxa',
                                                   'Ups_BestMatchScore',
                                                   'Dwn_BestMatchScore',
                                                   'U2_BestMatchScore',
                                                   'D2_BestMatchScore',
                                                   'Conserved_9_Score',
                                                   'Clamp_Affiity_Metric',
                                                   'WashUSS_Provisional',
                                                   'Clamp_Percentile',
                                                   'High_Novelty',
                                                   'Bard_mnemonic',
                                                   'AssemblyLocal']

    OutFile1 = open(OutFileName1, mode='wt')
    OutFile1.write(TaskHeader1)
    OutFile1.write(HeaderTranspose1(hTX1))
    OutFile1.write('\t'.join(hTX1)+'\t'+AbbrevHeader1+'\t '+Delimiter1)        
    BCList1 = [k1 for k1 in BarcodeD1.keys() if BarcodeD1[k1].c>=MinTotalReadCount1]
    for k1 in BCList1:
        BarcodeD1[k1].takestock()
    BCList1 = [k1 for k1 in BCList1 if BarcodeD1[k1].reportme]
    BCList1.sort(key=lambda x: BarcodeD1[x].c, reverse=True)
    for s1 in BCList1:
        b1 = BarcodeD1[s1]
        myReport1 = [s1,b1.c,CounterToText1(b1.Samples),CounterToText1(b1.DataFiles),b1.CoreSeq,b1.CoreScore,b1.upslencount,b1.maxdwnlen]
        for myAttr1 in ('wm','um','dm','um2','dm2'):
            myC1 = getattr(b1,myAttr1)
            if myAttr1=='wm':
                myReport1 += [myC1[0],myC1[1].s,myC1[1].phy.replace(' ','_')]
            else:
                myReport1.append(myC1[0])
        myReport1.extend([b1.conserved9,'%.2f'%b1.clampE,b1.WashUSS,'%.1f'%b1.clampEPercentile,b1.highnovelty,b1.bardid,b1.assembly])
        myLine1 = '\t'.join(map(str,myReport1)).replace('Counter','')+Delimiter1
        OutFile1.write(myLine1)
    OutFile1.close()
    if HighlyNovelOut1:
        BCList1 = [k1 for k1 in BCList1 if BarcodeD1[k1].highnovelty]
        BCList1.sort(key=lambda x: (x[-sdwn1:],-BarcodeD1[x].c))
        if BCList1:
            OutFile1 = open('HighNovelty_'+OutFileName1, mode='wt')
            OutFile2 = 'HighNoveltyStickModel_'+OutFileName1[:-4]+'.svg'
            OutFile1.write(TaskHeader1)
            OutFile1.write(HeaderTranspose1(hTX1))
            OutFile1.write('\t'.join(hTX1)+'\t'+AbbrevHeader1+'\t '+Delimiter1)        
            for nd1,s1 in enumerate(BCList1):
                b1 = BarcodeD1[s1]
                myReport1 = [s1,b1.c,CounterToText1(b1.Samples),CounterToText1(b1.DataFiles),b1.CoreSeq,b1.CoreScore,b1.upslencount,b1.maxdwnlen]
                for myAttr1 in ('wm','um','dm','um2','dm2'):
                    myC1 = getattr(b1,myAttr1)
                    if myAttr1=='wm':
                        myReport1 += [myC1[0],myC1[1].s,myC1[1].phy.replace(' ','_')]
                    else:
                        myReport1.append(myC1[0])
                myReport1.extend([b1.conserved9,'%.2f'%b1.clampE,b1.WashUSS,'%.1f'%b1.clampEPercentile,b1.highnovelty,b1.bardid,b1.assembly])
                myLine1 = '\t'.join(map(str,myReport1)).replace('Counter','')+Delimiter1
                OutFile1.write(myLine1)
                bList,score1 = MakeBList1(s1)
                showMe1(bList,nd1,b1.bardid,'%.2f'%b1.clampE+' (bootstrap: %0.1f'%b1.clampEPercentile+'%'+')')
        OutFile1.close()
    vdisplay(OutFile2)
vLog('finished',vFileOut=True)

            
## This somewhat dissheviled glop of code was written by Andrew Fire, Stanford University Departments of Pathology and Genetics
## Major contributions of ideas for conception and impementation: Daniel Chang (Stanford Genetics & Arc Institute)
## Please send additional comments, ideas, complaints, more ideas to us (can send to afire<rat>stanford<rat>edu , replacing the rodent with a useful symbol)
## With thanks for ideas/prompts/encouragement/healthy-speticism (in no particular order) to K.Artiles, C.Benko, S.U.Enam, D. Galls, E.Greenwald, O.Ilbay, D.E. Jeong, J.Kim, D.Lipman, T.Weissman, M.McCoy, M.Shoura, A.Bhatt, J.Salzman, D. Chang, R.N.Hall, I.Zheludev
## Second release version (still alpha, or perhaps sub-alpha) ca03 07-29-2025
## Version Adjustments starting 01-26-25
## 4-11-25 Changed the cutoff for conserved bases to require <= threshold versus >=.  While arguable, this seems to better define
##   criteria for rRNA membership that are the most critical in identifying candidates for novel-but-still-rDNA barcodes.  This only
##   changes the 'high novelty' categorization, nothing else.
## 4-11-25 Also so minor changes to the logging of files used by the program (to avoid huge log files)
## Note 4-11-25: Linux and PyPy3 don't carefully check memory availability, so running a very large job with a lot of threads can result in
## an out of memory system crash.  The system should just reboot, but interim calculations aren't kept.  Recommend limiting thread number on large jobs for this reason.
## This depends on data complexity, database sizes, etc, but at the very least keeping the number of threads used to number of cores on system minus 1 might be optimal
## for very large jobs.
## 6-27-25 Added a routine that estimates clamp energy surrounding the conserved 20-mer.  Very substantial opportunity here to separate
## bone-fide rRNA loci from random hits.  Also added a mnemonic feature that assigns an arbitrary word combination to each novel sequence.  Similar
## sequences will have similar mnemonics(up to a point), which  may be useful in perusing output.
## Copyright 2024-2025 Stanford University


           
                


    


            
    
    