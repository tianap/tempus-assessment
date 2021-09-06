import requests, sys, json

"""
Each variant must be annotated with the following pieces of information:
1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent,
intergenic, etc.). If there are multiple effects, annotate with the most deleterious
possibility.

    - Type of variation:
        - substitution (SNP) --> included
        - insertion ( len(ALT) > len(REF) ) --> included
        - deletion ( len(ALT) < len(REF) ) --> included
        - complex --> must be broken down
            - CNV
            - 
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from ExAC API (API documentation is available here:
http://exac.hms.harvard.edu/).
6. Any additional annotations that you feel might be relevant.

"""
def requestVEP(server, request, contentType):
    """
    Uses the Ensembl VEP REST API to find the variant effect. Based on the code from this webinar:
    https://drive.google.com/file/d/1jm8DJOqHxiZCLL0gkMnuyhU3ulLzOUuD/view
    
    """
    myRequest = requests.get(server+request, headers={"Content-Type": contentType})

    if not myRequest.ok:
        myRequest.raise_for_status()
        sys.exit()
    
    if contentType == 'application/json':
        return myRequest.json()
    else:
        return myRequest.text

def fetchVariantEffect(chr, regionStart, allele):
    """
    calls on requestVEP 
    GET vep/:species/region/:region/:allele/
    ignore strand
    """
    species = "human"
    regionEnd = str(int(regionStart) + len(allele) - 1)
    request = "/vep/"+species+"/region/"+chr+":"+regionStart+"-"+regionEnd+":/"+allele+"?"
    server = "http://grch37.rest.ensembl.org/"
    contentType = "content-type=application/json"
    print(requestVEP(server, request, contentType))
    #print(species, regionStart, regionEnd, allele)


def main():
    f = open('Challenge_data_(1).vcf', 'r')
    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "TYPE", "DEPTH"]))
    for line in f.readlines():
        if line.startswith("#"):
            continue
        else:
            chr, pos, id, ref, alt, qual, filter, info, format, normal, va5 = line.split('\t')
            
            # Parse INFO field
            infoDict = dict()
            for x in info.split(';'):
                key,value = x.split('=')
                infoDict[key] = value

            # Find type of variant
            type = infoDict['TYPE']

            # Find sequencing depth
            depth = infoDict['DP']

            # Find effect of variant
            #forwardCount = infoDict['SAF']
            #reverseCount = infoDict['SAR']
            fetchVariantEffect(chr, pos, alt)
            break
            # Write variant and associated information
            #print(chr, pos, id, ref, alt, type, depth, forwardCount, reverseCount)

    f.close()

if __name__ == "__main__":
    main()