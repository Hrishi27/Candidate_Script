#!/tools/bin/python
"""
The pedigree_analysis/candidate script is an annotation script that annotates individual variants at a position in a pedigree.
The current annotation sources include Kaviar,CADD,CMS,UCSC_known_genes,Clinvar and SNPeff(Whenever it is loaded), these annotations 
are pulled from the impala system by querying appropriate database tables.
The candidate script also annotates a variant with a stress score, which is obtained by processing the target inheritance pattern 
and the observed pattern w.r.t to a allele. (Segregation class)
The script currently uses mergedVCF as an input, but in future certain changes will be made to get VCF information from impala.
-----------------------------------------------------------------------------------------------------------------------------------

Usage: python pedigree_analysis.py -i <INPUTVCF> -o <OUTFILE> -c <CONFIG>

Note: The input.config file allows user to control the target inheritance pattern and various annotation cutoff's.

"""
from __future__ import print_function
import re
import sys
import getopt
import vcf
from collections import Counter
import ConfigParser
from annotater import annotater
from segregation import segregation
from impala.dbapi import connect
import time
import logging

LOG_FILENAME = 'Pedigree.log'
logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.ERROR,
                    )

## Connection to impala.
## Since we are using mergedVCF as input it doesn't make sense to have the connection object
## here in the main method. I will either move it to annotater class or keep it here depending
## upon how the genotype representation is handled in the variant table on Impala

conn=connect(host='glados19', port=21050,database="users_hrishi")
DB = conn.cursor()
DB.execute("SELECT VERSION()")
results=DB.fetchone()
if results:
    logging.info ("Connection successful" )
else:
    logging.error ("Failed to make connection")
    sys.exit()

##Vividict allows creation of perl like hashes of hashes or multilevel hash
## One would always have to point their empty hash(dict) to this class before
## storing multiple levels of keys and values 

class Vividict(dict):
    def __missing__(self, key):
        value= self[key] = type(self)()
        return value

##This allows reading and parsing of the configuration file

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


## This function takes all the genotype scores(gqx) from all samples at a query position
## calculates the average and returns a list of max,min and average scores.

def compute_sequence_quality_score_statistics(score_array):
    max_score=max(score_array)
    min_score=min(score_array)
    num_list=[float(x) for x in score_array]
    avg_score=sum(num_list)/len(num_list)
    return (max_score,min_score,avg_score)


##This is currently not used but will be used when a decision is made on what 
## parameter would go into the overall_candidate_score calculation

def quality_score_adjustment(max_score):
    score = max_score
    return_value=0
    if score >=50:
        return 0
    elif score > 35:
        return_value=50-score/15
        return return_value
        #minimum score is 20
    return_value = 1+((35-score)/3.75);  #a score of 20 gets a modification of 4, linearly decreasing to a modification of 1 at a score of 35
    return return_value


def snphash(file):
    hash=Vividict()
    fo = open(file, "r")
    for line in fo:
        line = line.rstrip('\n')
        line = line.strip()
        if line.startswith('#'):
            continue
        array=[]
        array=line.split("\t")
	chromosome=array[0].replace("chr","")
        position=int(array[1])
        hash[chromosome][position][array[3]][array[4]]=array[7]
    return hash
        #print (array[0],array[1],array[3],array[4],array[7],sep="\t",end="\n")
    

## Work in progress


def overall_candidate_score():
    return 1


def usage():
    print ("Usage: python pedigree_analysis.py -i <INPUTVCF> -o <OUTFILE> -c <CONFIG>")

    
## process_genotype is an important function, every heterozygous line [ALT != .] from the 
## mergedVCF is processed here to extract sample specific genotypes(GT) and genotype score(GQX).
## This also translates no-calls (GT=./.) into homozygous reference alleles,skips position and
## reports it in error logs if any of the given sample have only one genoytpe information present
## eg: (GT:GQX:DP:DPF:GQ:AD  0:17:19:11:.:.) The sample output from this function is provided below
## {'LP6005636-DNA_G01': {48: [G, G]}, 'LP6005636-DNA_H01': {137: ['A', G]}, 'LP6005636-DNA_A02': {78: [G, G]}}


def process_genotype(records):
    checker=0
    data_hash=Vividict()
    reference=records.REF
    position=records.POS
    for key in records.samples:
        allele=[]
        if key.data.GT is None:   ## if GT 0/0, I re-create genotypes for this as homozygous reference
            no_call=reference + ', ' +reference
            allele.append(no_call)
            if key.data.GQX is None:
                data_hash[key.sample][1]=allele  ## if GQX is none, I use quality as 1
            else:
                data_hash[key.sample][key.data.GQX]=allele
        else:
            genotype_array=re.split('/',key.data.GT)
            if len(genotype_array) == 1:     ## In this case GT just provides information about single allele. eg GT:GQX:DP:DPF:GQ:AD     0:17:19:11:.:. such positions are skipped. 
                logging.error("skipping" + ' ' + chromosome + ' ' + str(position) +  '' + "due to inconsistent GT representation")
                checker=1
            for genotype in genotype_array:
                genotype=int(genotype)
                if genotype==0:
                    allele.append(reference)
                    if key.data.GQX is None:
                        data_hash[key.sample][0]=allele
                    else:
                        data_hash[key.sample][key.data.GQX]=allele
                else:
                    index=1-genotype
                    allele.append(alt[index])
                    if key.data.GQX is None:
                        data_hash[key.sample][0]=allele
                    else:
                        data_hash[key.sample][key.data.GQX]=allele
    if checker==1:
        return 0
    else:
        return data_hash


configfile= ''
inputfile = ''
outputfile = ''
snpeff_file=''

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:c:h:",["input=","output=","config=","help"])
except getopt.GetoptError:
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print ('Usage:pedigree_analysis.py -i <inputfile> -o <outputfile> -c <configfile>')
        sys.exit()
    elif opt in ("-i", "--input"):
        inputfile = arg
    elif opt in ("-o", "--output"):
        outputfile = arg
    elif opt in ("-c", "--config"):
        configfile = arg
    else:
        print ("Type python pedigree_analysis.py -h to know more")
        sys.exit(2)


if not (inputfile and configfile and outputfile):
    print ("Please provide required input,output and configfile")
    usage()
    sys.exit()

##Reading of the configuration file to
##retrieve users input to the script 


Config = ConfigParser.ConfigParser()
Config.read(configfile)

members=ConfigSectionMap("Inputs")["sample_id"]
target_pattern=ConfigSectionMap("Inputs")["target"]
stress_cutoff=ConfigSectionMap("Inputs")["stress_cutoff"]
freq_cutoff=ConfigSectionMap("Inputs")["qaf_cutoff"]
candidate_cutoff=ConfigSectionMap("Inputs")["score_cutoff"]
snpeff_file=ConfigSectionMap("Inputs")["snpeff_file"]

snpeff_file=snpeff_file.replace('"','')

### Creating objects for annotater and
## segregation classes
segregation_object = segregation()
annotation_object = annotater()

## Reading input mergedVCF file
vcf_Reader=vcf.Reader(open(inputfile),'r')

target = open(outputfile, 'w')

target.write("#Chromosome\tPosition\tReference\tAllele\tcandidate_stress\tmax_score\tmin_score\tavg_score\tCount\tKaviar\tUCSC_gene\tclinvar\tcms_annotation\tcadd_annotation\tsnpeff_annotation\tcandidate_score\n")

## The following reads the mergedVCF using python's VCF tools
## All positions that are homozygous reference [ALT=.] are skipped here.
## Every line is then passed onto the process_genotype function to get information 
## about sample specific genotypes and quality scores.
## Next the script calls the compute_sequence_quality_score_statistics function to get 
## maximum,minimum and average scores. In the further steps unique alleles at every position
## are annotated with kaviar,UCSC_known_gene,CMS,Clinvar,CADD and segregation score/stress scores
## Individual variants with their annotations are printed into tab separated file

snpeff=Vividict()
snpeff=snphash(snpeff_file)
chr_check={}
for records in vcf_Reader:
    if records.ALT[0] is None:
        continue
    reference=records.REF
    position=records.POS
    chromosome=records.CHROM
    chromosome=chromosome.replace("chr","")
    alt=records.ALT
    data_hash=process_genotype(records)
    
    if not chromosome in chr_check:
        print ("Working on chromosome" + chromosome)
        chr_check[chromosome]=1
    
    if data_hash==0:
        continue
    else:
        allele_list=[]
        quality_list=[]
        for line in data_hash:
            for quality in data_hash[line]:
                quality_list.append(quality)
		allele_string=""
                allele_string=str(data_hash[line][quality])
                allele_list.append(allele_string)
        max_score,min_score,avg_score=compute_sequence_quality_score_statistics(quality_list) 
        avg_score="%.2f" % (avg_score)
        allele=""
        counter={}
        for string in allele_list:
            string=string.replace(']',' ').replace('[',' ').replace("'",' ').replace(',',' ') ##This is needs to be be changed.
            allele += string
        allele=re.sub("\s\s+"," ",allele)
        allele=re.sub("^\s","",allele)
        allele=allele.rstrip()
        allele_array=[]
        allele_array=allele.split(" ")
        allele_count={}
        allele_count=Counter(allele_array)
        Unique_array=set(allele_array)
        ### Position Annotations###          
 
        kaviar_annotation=annotation_object.query_kaviar(chromosome,position,reference,DB)
        gene_annotation=annotation_object.get_UCSC_genes(chromosome,position,DB)
        cms_annotation=annotation_object.query_cms(chromosome,position,DB)
	clinvar_annotation=annotation_object.get_clinvar(chromosome,position,DB)
	cadd_phred=annotation_object.get_cadd(chromosome,position,DB)
        
        ### Allele Annotations###

        for candidate in Unique_array:
            genotype_vectors=segregation_object.standardized_genotype_vector_with_reference_to_a_particular_allele(allele_array,candidate)
            candidate_stress=segregation_object.target_test(genotype_vectors,target_pattern)
            if snpeff[chromosome][position][reference][candidate]:
	        snpeff_annotation=snpeff[chromosome][position][reference][candidate]
            else:
	        snpeff_annotation='NA'
	    if candidate_stress=='NA':
		logging.error("length of target_pattern and inheritance pattern are not equal at" + chr + ' ' +pos)
	        continue
            if candidate in kaviar_annotation:
                kaviar=kaviar_annotation[candidate]
	        kaviar= "%.3f" % (kaviar) 
            else:
                kaviar=0

            if candidate in clinvar_annotation:
                clinvar=clinvar_annotation[candidate]
            else:
                clinvar="none"

           ## For cadd annotations if the reference or queried allele is a non SNV
           ## maximum score from the dict is reported in the final output
           
            if len(reference) > 1 or len(candidate)>1:
                maximum=max(cadd_phred.iterkeys(), key=lambda k: cadd_phred[k])
                cadd_annotation=str(cadd_phred[maximum])
            elif candidate==reference: 
                maximum=max(cadd_phred.iterkeys(), key=lambda k: cadd_phred[k])
                cadd_annotation=str(cadd_phred[maximum])
            else:
                cadd_annotation=cadd_phred[candidate]
            target.write (str(chromosome) + "\t" + str(position) + "\t" + str(reference) + "\t" + str(candidate) + "\t" + str(candidate_stress) + "\t" + str(max_score) + "\t" + str(min_score) + "\t" + str(avg_score) +"\t"+ str(allele_count[candidate]) + "\t" + str(kaviar) + "\t" + str(gene_annotation) + "\t" + str(clinvar) + "\t" + str(cms_annotation) + "\t" +str(cadd_annotation) + "\t" + str(snpeff_annotation) + "\t" + '1' +"\n")



