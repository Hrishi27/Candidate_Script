"""
The class annotater includes all the annotation and data processing calls to impala.
Make sure that this class is imported in the pedigree_analysis.py script as 'from annotater import annotater'
The annotater currently uses database handle (cursor object) supplied from the main method. 
Most annotations below are queried with chromosome and position combination.(Kaviar with chr,pos and ref)
A hash/dict with position annotation is returned to the main method where allele/variant
annotation information is retrieved by querying the hash.

"""

from __future__ import print_function
from impala.dbapi import connect
import time
import os

class annotater():

    ##Constructor method##
    def __init__(self):
        self.annotater=""

   ## query_kaviar queries kaviar table on impala with chromosome,position and reference information
   ## Allele frequency(AF) from the each retrieved row is subtracted from 1 and the min of the subtracted value
   ## is considered as reference allele frequency. A dict/hash with alleles as keys and AF's as values is returned
   ## to the main method.  

    def query_kaviar(self,chrom,position,ref,DB):
        impala_query=("select ref,alt,allele_freq from p7_ref_grch37.kaviar where chrom='%s' and pos=%s and ref='%s'" % (chrom,position,ref))
        DB.execute(impala_query)
        kaviar=DB.fetchall()
        QAF_hash={}
        if kaviar:
            for freq in kaviar:
                ref_allele=freq[0]
                alt=freq[1]
                allele_freq=freq[2]
                if alt in QAF_hash:
                    QAF_hash[alt]=QAF_hash[alt] + allele_freq ## This is done if the same allele is repeated with different AF, all the AF's are added together.
                else:
                    QAF_hash[alt]=allele_freq
            reference_array=[]
            for key, value in QAF_hash.iteritems():
                data_ref=1-value
                reference_array.append(data_ref)
            ref_min=min(reference_array)
            QAF_hash[ref_allele]=ref_min
            return(QAF_hash)
        else:
            QAF_hash["none"]=0
            return(QAF_hash)

    ## query_cms queries the cms table from impala based upon chromosome, position combination.
    ## These are range searches if the queried position lies anywhere between the start and stop
    ## position in the database we report 1. Currently there is no score available for this table. 
    
    def query_cms(self,chrom,position,DB):
        chrom='chr'+ str(chr)
        position=int(position)
        cms=""
        DB.execute("select * from p7_itmi.cms_gt1 where %s >= start and %s <= stop and chrom='%s'" %(position,position,chrom))
        cms=DB.fetchall()
        if cms:
            return "1"
        else:
            return "0"

    ## get_UCSC_genes queries the ucsc_genes table from impala to report genes spanning a particular variant position
    ## There are situation where multiple gene names are reported, we report only unique gene names in our final outputs.  

    def get_UCSC_genes(self,chr,pos,DB):
        DB.execute("select gene_name from p7_ref_grch37.ucsc_genes where chrom='%s' and %s>= txstart and %s <= txend" % (chr,pos,pos))
        results=DB.fetchall()
        uniq={}
        gene_annotation=""
        for gene_name in results:
	    gene_name=' '.join(gene_name)
            uniq[gene_name]=1
        if len(uniq) == 0:
            return ("none")
        else:
            for gene in uniq:
                gene_annotation += gene + ' '
        return (gene_annotation)
   
    ## get_cadd queries the cadd table from impala to get CADD scores for the queried chromosome and position combination
    ## A dict/hash is returned for every reported SNV along with their corresponding phred_scores.
    ## Please read the comments in the main method to know how non SNV alleles and reference alleles are handled. 

    def get_cadd(self,chr,pos,DB):
        DB.execute("select phred_a,phred_c,phred_g,phred_t from p7_ref_grch37.cadd where chrom='%s' and pos=%s" %(chr,pos))
        results=DB.fetchall()
        phred={}
        for row in results:
            phred['A']=row[0]
            phred['C']=row[1]
            phred['G']=row[2]
            phred['T']=row[3]
        return phred

    ## get_clinvar queries the clinvar table from impala to report allele specific clinvar significance score.
    ## A hash/dict is returned with allele as keys and significance score as the values.  

    def get_clinvar(self,chr,pos,DB):
        clinvar={}
        DB.execute("select chrom,pos,ref,alt,clin_sig from p7_ref_grch37.clinvar where chrom='%s' and pos=%s" %(chr,pos))
        results=DB.fetchall()
        if len(results) > 0:
            for row in results:
                clinvar[row[3]]=row[4]
            return clinvar
        else:
            clinvar["none"]=1
            return(clinvar)

