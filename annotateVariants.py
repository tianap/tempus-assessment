import requests, sys, json, argparse

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

class CommandLine():
    """
    Handles the command line, usage, and help requests.
    """
    def __init__(self):
        """
        CommandLine constructor. Uses parser to interpret command line using argparse.

        Required arguments:
            --inputFile : ()
            --outputFile : ()
        """

        self.parser = argparse.ArgumentParser(
            description = "Annotates variants from given VCF file.",
            add_help=True,
            prefix_chars="-",
            usage="py script.py -i input.vcf -o output.tsv"
        )

        self.parser.add_argument('-i', metavar="--inputFile", help="input VCF file", required=True)
        self.parser.add_argument('-o', metavar="--outputFile", help="output TSV file", required=True)
    
        self.args = self.parser.parse_args()

class Variant:
    """
    Object that stores associated variant information directly from the input VCF file, 
    (i.e. chromosome, position, ref/alt allele) and information retrieved from external 
    sources (i.e. variant effect, allele frequency).

    Attributes:
        - chr : (str) chromosome in which variant is found
        - pos : (str) chromosomal position of variant
        - id : (str) variant ID from VCF
        - ref : (str) reference allele 
        - alt : (str) variant allele
        - qual : (str) quality of observed variant form VCF
        - type : (str) type of variant
        - depth : (str) read depth of variant
        - varTotalReads : (str) number of reads supporting the variant
        - refTotalReads : (str) number of reads supporting the reference allele
        - percentRatio : (str) percentage of reads supporting variant vs readss upporting ref
        - effect : (str) most severe consequence of variant, retrieved from Ensembl VEP
        - aleleFreq : (str) allele frequency of variant, retrieved from ExAC

    Methods:
        - printVariant(output) : writes formatted variant info to provided output file
        - getPercentageVarVsRef() : calculates ratio of percentages of reads supporting variant 
            vs ref
        - vepFormat() : returns variant info formatted for VEP POST request
        - updateEffect(effect) : updates 'effect' attribute
        - exacFormat() : returns variant info formatted for ExAC POST request
        - updateAlleleFreq(af) : updates 'alleleFreq' attribute

    """
    def __init__(self, vcfLine):
        """
        Class Variant constructor. 

        Parameter(s):
            - vcfLine : (str) line from VCF file
        """
        # Retrieve fields from VCF line
        chr, pos, id, ref, alt, qual, filter, info, format, normal, va5 = vcfLine.split('\t')
        
        # Assign basic variant information
        self.chr = chr
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        
        # Parse INFO field to get variant type, depth, number of reads supporting variant and ref
        infoDict = dict()
        for x in info.split(';'):
            key,value = x.split('=')
            infoDict[key] = value

        self.type = infoDict['TYPE']  # type of variant
        self.depth = infoDict['DP']  # sequencing depth
        self.varTotalReads = infoDict['AO']  # number of reads supporting the variant
        self.refTotalReads = int(infoDict['RO'])  # number of reads supporting the reference

        # Initializes the percentage ratio 
        self.percentRatio = self.getPercentageVarVsRef()
        
        # Updated with variant effect after vep POST request is made
        self.effect = ""

        # Updated with allele frequency after ExAC POST request is made
        self.alleleFreq = ""

    def printVariant(self, output):
        """
        Writes variant information to output file.

        Parameter(s):
            - output : (str) output file

        """
        output.write("\t".join([self.chr, self.pos, self.id, self.ref, self.alt, self.type, 
                                self.effect, self.depth, self.varTotalReads, self.percentRatio, 
                                self.alleleFreq]))
        output.write("\n")
    
    def getPercentageVarVsRef(self):
        """
        Calculates the percentage of reads supporting variant versus reads supporting ref allele.
        If more than one allele, calculates percent ratios for each pair of variant versus ref.

        Returns:
            - percentRatio : (str) string format (var %:ref %) of calculated ratio
        """
        if ',' in self.varTotalReads:  # more than one alternate allele
            ratios = []
            for val in self.varTotalReads.split(','):
                varPercent = int(val) / (int(val) + self.refTotalReads)
                refPercent = self.refTotalReads / (int(val) + self.refTotalReads)
                ratios.append("{:.3f}:{:.3f}".format(varPercent, refPercent))
                percentRatio =",".join(ratios)
        else:
            varPercent = int(self.varTotalReads) / (int(self.varTotalReads) + self.refTotalReads)
            refPercent = self.refTotalReads / (int(self.varTotalReads) + self.refTotalReads)
            percentRatio = "{:.3f}:{:.3f}".format(varPercent, refPercent)
        
        return percentRatio
    
    def vepFormat(self):
        """
        Returns the format of variant information required for VEP POST request.  Called on
        by global method vepPOSTRequest().
        """
        return "\t".join([self.chr, self.pos, self.id, self.ref, self.alt])

    def updateEffect(self, effect):
        """
        Updates the attribute 'effect' with the given consequence. Called on by global method 
        vepPOSTRequest().

        Parameter(s):
            - effect : (str) variant consequence retrieved from VEP
        """
        self.effect = effect 
    
    def exacFormat(self):
        """
        Returns the format of variant information required for ExAC POST request.  Called on
        by global method exacPOSTRequest().
        """
        return "-".join([self.chr, self.pos, self.ref, self.alt])

    def updateAlleleFreq(self, af):
        """
        Updates the attribute 'effect' with the given allele frequency. Called on by global method 
        exacPOSTRequest().

        Parameter(s):
            - af : (str) allele frequency retrieved from ExAC
        """
        self.alleleFreq = af

def vepPOSTRequest(variantInput):
    """
    Makes POST request to Ensembl Variant Effect Predictor (VEP) tool 200 variants at a time.
    Retrieves each variant's most severe consequence and updates the 'effect' attribute of the
    corresponding Variant object. Calls on Variant methods vepFormat() and updateEffect(). Utilizes
    requests and json libraries.

    NOTE: The default human assembly is GRCh37/hg19 and uses this version of VEP.

    Parameter(s):
        - variantInput : (list) list of 200 Variant objects
    
    """

    # Reformat variants in list to acceptable vep format
    vepInput = [var.vepFormat() for var in variantInput]

    # Convert list into json format
    variantJsonInput = json.dumps({"variants" : vepInput})

    # Make POST request
    server = "http://grch37.rest.ensembl.org/"
    request = "vep/homo_sapiens/region/"
    headers = {"Content-Type" : "application/json", "Accept" : "application/json"}
    myRequest = requests.post(server+request, headers=headers, data = variantJsonInput)

    if not myRequest.ok:
        myRequest.raise_for_status()
        sys.exit()
    
    jsonResult = myRequest.json()

    for i in range(len(variantInput)):
        effect = jsonResult[i]['most_severe_consequence']
        variantInput[i].updateEffect(effect)
     
def exacPOSTRequest(variantInput):
    """
    Makes POST request to ExAC search API 200 variants at a time. Retrieves each variant's allele 
    frequency and updates the 'alleleFreq' attribute of the corresponding Variant object. Calls on
    Variant methods vepFormat() and updateEffect(). Utilizes requests and json libraries.
    
    Parameter(s):
        - variantInput : (list) list of 200 Variant objects

    Returns:
        - variantInput : (list) updated list of Variant objects with new 'alleleFreq' attributes
    
    """

    # Reformat variants according to ExAC format requirements
    formattedInput = [var.exacFormat() for var in variantInput]

    # Specify url and headers
    reqUrl = "http://exac.hms.harvard.edu/rest/bulk/variant"
    head = {"Content-Type": "application/json", "Accept":"application/json"}

    # Make request, use json.dumps() on formatted input
    myReq = requests.post(url=reqUrl, data= json.dumps(formattedInput), headers=head)

    if not myReq.ok: 
        print(myReq.raise_for_status)
        sys.exit()
    else:  # request successful
        result = myReq.json()
        
        # Iterate through Variant objects
        for i in range(len(variantInput)):  
            variant = variantInput[i]
            varExacFormat = formattedInput[i]
            try:  # allele frequency is available through ExAC
                varAF = "{:.4}".format(float(result[varExacFormat]["variant"]["allele_freq"]))
            except:  # if no allele frequency provided, use '.'
                varAF = "."
            variant.updateAlleleFreq(varAF)

def main():
    """
    The main() function reads through the provided .vcf file, calls on a variety of functions to 
    annotate variants and write the results to the provided output file. Each variant is made into a 
    Variant object (refer to documentation above) and stores all information pertaining to that 
    variant.
    
    """
    # Create CommandLine object (used for parsing arguments)
    newCommandLine = CommandLine()

    f = open(newCommandLine.args.i, 'r')
    output = open(newCommandLine.args.o, 'w')

    # Write output header
    output.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "TYPE", "EFFECT", "DEPTH", 
                            "SUPPORT", "VAR:REF READS", "ALLELE FREQ"]))
    output.write("\n")

    # Keep track of variants for each post request
    variantList = []

    # Get number of lines
    numLines = len([line.strip("\n") for line in f if line != "\n"])
    
    # Reset readlines to first line
    f.seek(0)

    # Keep track of line count and post requests
    lineCount = 0
    
    # Iterate through lines of file
    for line in f.readlines():
        if line.startswith("#"):
            lineCount += 1
            continue
        else:
            lineCount += 1

            # Create new Variant object
            newVariant = Variant(line)

            # Add new variant to the current list of variants
            variantList.append(newVariant)

            # Make POST requests for VEP (variant effect) and Exac (allele frequency)
            if len(variantList) == 200 or lineCount == numLines:
                
                # make POST request to VEP with current 200 variants
                vepPOSTRequest(variantList)

                # make POST request to ExAC with current 200 variants
                exacPOSTRequest(variantList)

                for var in variantList:
                    var.printVariant(output)

                variantList = []  # reset list
            
    # Close files
    f.close()
    output.close()

if __name__ == "__main__":
    main()