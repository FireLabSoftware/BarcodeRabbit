#!/usr/bin/env python
## 
## ######################
## BarcodeRabbit-- SSU rRNA-based Reference-agnostic phylogenetic survey of cellular composition from NGS datasets
## Version bb01 040625
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
##  -   16 AssemblyLocal: A very naive local assembly that may be of use in a blast search or other analysis.  Not of quality expected for a program like Spades or Megahit-- this is merely to facilitate interim assessment.
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
MaxReadsPerFile1 = 0 ## 0           ## Set to a positive number to limit the number of reads per data file (sample)

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
HNT1 = {'wm':7,'um':3,'dm':3,'um2':1,'dm2':1, 'core':2, 'UpsLenCount':2,'maxDwnLen':6,'conserved-9':4}  ## (highnoveltythresholds=)Only barcodes with values that match these (<= for 'core', 'UpsLenCount', 'maxDwnLen', 'conserved-9', >= for others) will be chosen as "HighNovelty"
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

import sys,os,gzip,pickle
from operator import ne
from itertools import combinations,zip_longest
from array import array
from glob import glob
from time import sleep
from random import randint


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
        if (self.wm[0]>=HNT1['wm'] and self.um[0]>=HNT1['um'] and self.dm[0]>=HNT1['dm']
            and self.um2[0]>=HNT1['um2'] and self.dm2[0]>=HNT1['dm2']
            and self.CoreScore<=HNT1['core'] and self.upslencount>=HNT1['UpsLenCount']
            and self.maxdwnlen>=HNT1['maxDwnLen'] and self.conserved9<=HNT1['conserved-9']):
            self.highnovelty = 1
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
                    if not(hasCore1(L1) or hasCore1(L2)): continue
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
                                                   'High_Novelty',
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
        myReport1.extend([b1.conserved9,b1.highnovelty,b1.assembly])
        myLine1 = '\t'.join(map(str,myReport1)).replace('Counter','')+Delimiter1
        OutFile1.write(myLine1)
    OutFile1.close()
    if HighlyNovelOut1:
        BCList1 = [k1 for k1 in BCList1 if BarcodeD1[k1].highnovelty]
        BCList1.sort(key=lambda x: (x[-sdwn1:],-BarcodeD1[x].c))
        if BCList1:
            OutFile1 = open('HighNovelty_'+OutFileName1, mode='wt')
            OutFile1.write(TaskHeader1)
            OutFile1.write(HeaderTranspose1(hTX1))
            OutFile1.write('\t'.join(hTX1)+'\t'+AbbrevHeader1+'\t '+Delimiter1)        
            for s1 in BCList1:
                b1 = BarcodeD1[s1]
                myReport1 = [s1,b1.c,CounterToText1(b1.Samples),CounterToText1(b1.DataFiles),b1.CoreSeq,b1.CoreScore,b1.upslencount,b1.maxdwnlen]
                for myAttr1 in ('wm','um','dm','um2','dm2'):
                    myC1 = getattr(b1,myAttr1)
                    if myAttr1=='wm':
                        myReport1 += [myC1[0],myC1[1].s,myC1[1].phy.replace(' ','_')]
                    else:
                        myReport1.append(myC1[0])
                myReport1.extend([b1.conserved9,b1.highnovelty,b1.assembly])
                myLine1 = '\t'.join(map(str,myReport1)).replace('Counter','')+Delimiter1
                OutFile1.write(myLine1)
            OutFile1.close()
vLog('finished',vFileOut=True)

            
## Copyright 2024-2025 Andrew Fire, Stanford University, All Rights Reserved
## Please send comments, ideas, complaints, more ideas to me at afire<rat>stanford<rat>edu (replacing the rodent with a useful symbol)
## With thanks for ideas/prompts/encouragement/healthy-speticism (in no particular order) to K.Artiles, C.Benko, S.U.Enam, D. Galls, E.Greenwald, O.Ilbay, D.E. Jeong, J.Kim, D.Lipman, T.Weissman, M.McCoy, M.Shoura, A.Bhatt, J.Salzman, D. Chang, R.N.Hall, I.Zheludev
## First release version (alpha, or perhaps sub-alpha) ba00 01-26-2025
## Version Adjustments starting 01-26-25
## 4-11-25 Changed the cutoff for conserved bases to require <= threshold versus >=.  While arguable, this seems to better define
##   criteria for rRNA membership that are the most critical in identifying candidates for novel-but-still-rDNA barcodes.  This only
##   changes the 'high novelty' categorization, nothing else.
## 4-11-25 Also so minor changes to the logging of files used by the program (to avoid huge log files)
## Note 4-11-25: Linux and PyPy3 don't carefully check memory availability, so running a very large job with a lot of threads can result in
## an out of memory system crash.  The system should just reboot, but interim calculations aren't kept.  Recommend limiting thread number on large jobs for this reason.
## This depends on data complexity, database sizes, etc, but at the very least keeping the number of threads used to number of cores on system minus 1 might be optimal
## for very large jobs.


           
                


    


            
    
    
