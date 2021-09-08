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

def parseVariant():
    """
    """
    return

def requestVEP(server, request):
    """
    Uses the Ensembl VEP REST API to find the variant effect. Based on the code from this webinar:
    https://drive.google.com/file/d/1jm8DJOqHxiZCLL0gkMnuyhU3ulLzOUuD/view
    
    """
    myRequest = requests.get(server+request, headers={"Content-Type" : "application/json"})

    if not myRequest.ok:
        myRequest.raise_for_status()
        sys.exit()
    
    return myRequest.json()

def fetchVariantEffect(chr, startPos, endPos, ref, allele, type):
    """
    calls on requestVEP 
    GET vep/:species/region/:region/:allele/
    ignore strand
    """
    # 1. Get chromosomal positions

    # NOTE: chromosomal positions required for VEP input depend on the type of variant.
    # In particular, the insertion coordinates from the VEP grch37 documentation:
    # The following conditional statements determine whether an end position has been included in
    # the original file. If not, end positions are calculated according to the type of variant.

    if endPos == -1:  # no end position given in the origianl VCF file
        # From the VCF 4.1 documentation:
        # "For precise variants, END is POS + length of REF allele - 1, and for imprecise 
        # variants the corresponding best estimate.
        # NOTE: In the event of a imprecise variant, the END position calculation is the best 
        # estimate.

        if type == "snp":
            endPos = startPos
        else:  # either deletion, insertion, or complex
            endPos = int(startPos) + len(ref) - 1
            
            # Check if insertion or deletion
            if type == "ins":
                # From VEP documentation:
                # "An insertion (of any size) is indicated by start coordinate = end coordinate + 1"

                startPos = endPos + 1  # reset start position, cast as string

    # 2. Parameters for VEP API request
    species = "human"
    request = "/vep/"+species+"/region/"+chr+":"+str(startPos)+"-"+str(endPos)+":/"+allele+"?"
    server = "http://grch37.rest.ensembl.org/"
    #contentType = "content-type=application/json"
    print(server+request)
    # Save the returned JSON

    # 3. Make request to VEP API
    requestJSON = requestVEP(server, request)

    # 4. Grab most sever consequence from JSON and return
    mostSevereEffect = requestJSON[0]['most_severe_consequence']

    return mostSevereEffect
    
def getStartPos(pos, allele, mut):
    """
    """
    if mut == 'ins':
        startCoord = str(int(pos) + len(allele))
    elif mut == 'del':
        startCoord = pos

    return startCoord
    
def main():
    f = open('Challenge_data_(1).vcf', 'r')
    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "TYPE", "DEPTH", "EFFECT"]))
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

            if type == "snp":
                continue
            
            # Find sequencing depth
            depth = infoDict['DP']

            # Check if there are multiple variant alleles
            if "," in alt:
                altAlleles = alt.split(',')
                altAlleleTypes = type.split(',')
                alleleEffects = []
                
                # grab consequences for all variant alelles
                for i in range(len(altAlleles)):
                    allele = altAlleles[i]
                    
                    '''
                    # Check if allele is ins or del
                    if altAlleleTypes[i] == 'ins':
                        # If insertion, startCoord = endCoord + 1
                        # endCoord =
                        # FROM VCF4.1 documentation:
                        # "For precise variants, END is POS + length of REF allele - 1, and the for imprecise variants the corresponding best
                        # estimate."
                        # Check if variant is precise
                        startCoord = str(int(pos) + len(allele))

                    #elif altAlleleTypes[i] == 'del':
                    else:
                        startCoord = pos
                    '''

                    #print(line, startCoord)
                    #lenDiff = len(allele) - len(ref)
                    effect = fetchVariantEffect(chr, pos, -1, ref, allele, type)
                    alleleEffects.append(effect)

                print(alleleEffects)
                sys.exit()
            
            # Find effect of variant
            
            effect = fetchVariantEffect(chr, pos, -1, ref, alt, type)

            # Write variant and associated information
            print(chr, pos, id, ref, alt, type, depth, effect)
            
    f.close()

if __name__ == "__main__":
    main()