import argparse
import json
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Filter reads based on Taxonomer read classifications of interest")
parser.add_argument("sequence_file_1", type=str, help="fastq file of forward read sequences to filter")
parser.add_argument("sequence_file_2", type=str, help="fastq file of reverse read sequences to filter")
parser.add_argument("json_file", type=str, help="taxonomer json with organisms of interest")
parser.add_argument("taxonomer_file", type=str, help="taxonomer classifier raw output file")
parser.add_argument("-v", "--verbose", help="display every classified read taxonomy and id",
                    action="store_true")
parser.add_argument("-o", "--output", type=str, nargs=2, help="write filtered reads to paired-end fastq files specified")
args = parser.parse_args()

# json file can have many levels of nesting which need to be traversed, using code from:
# 
# https://stackoverflow.com/questions/21028979/
# recursive-iteration-through-nested-json-for-specific-key-in-python/21029414
def item_generator(json_input, lookup_key):
    if isinstance(json_input, dict):
        for k, v in json_input.iteritems():
            if k == lookup_key:
                yield v
            else:
                for child_val in item_generator(v, lookup_key):
                    yield child_val
    elif isinstance(json_input, list):
        for item in json_input:
            for item_val in item_generator(item, lookup_key):
                yield item_val

target_reads = set()
taxids = list()
names = list()
classifications = dict()

# find a better way to get map from json data in the future
# not familiar enough with it now to create a nice solution
with open(args.json_file) as json_data:
    d = json.load(json_data)
    ids = item_generator(d, "id")
    for i in ids:
        taxids.append(i)
    namez = item_generator(d, "name")
    for n in namez:
        names.append(n)
    for i in range(0, len(taxids)):
        classifications[taxids[i]] = names[i]

# taxonomer raw tab-delimited output format as defined at https://www.taxonomer.com/faq
#
# 1. Bin category
# 2. Read classified (C) or unclassified (U) following bin assigment
# 3. Read name
# 4. Taxid of read assignment. For virus this corresponds to NCBI taxid. 
#    For bacteria, this corresponds to an artificial ID to the greengenes database, 
#    likewise for fungi except to the UNITE database. For Human this is an artificial ID to a gene.
# 5. Number of taxonomic levels at which the read assignment was made (more levels is a deeper assignment). 
#    Note this value cannot be compared between bin categories because difference taxonomies are used 
#    for example, one cannot use this to compare the depth of a bacterial assignment to that of a viral 
#    assignment.
# 6. Read length

for line in open(args.taxonomer_file):

    # parse tab delimited raw taxonomer file
    parsed = line.rstrip("\n").split("\t")
    bin_category = parsed[0]
    read_name = parsed[2]
    taxid = parsed[3]
    tid = taxid
 
    # convert taxid to match json format, also prevents collision of ids
    if bin_category == "viral":
        tid = "v"+taxid
    elif bin_category == "bacterial":
        tid = "b"+taxid
    elif bin_category == "fungal":
        tid = "f"+taxid

    if tid in taxids:
        if args.verbose:
            print("%s\t%s\t%s" % (tid, classifications[tid], read_name))
        target_reads.add(read_name)
    
print "Number of target sequences found: %i" % (len(target_reads))

if args.output:
    records_1 = (r for r in SeqIO.parse(args.sequence_file_1, "fastq") if r.id in target_reads)
    count = SeqIO.write(records_1, args.output[0], "fastq")
    print "Saved %i records from %s to %s" % (count, args.sequence_file_1, args.output[0])

    if count < len(target_reads):
        print "Warning %i IDs not found in %s" % (len(target_reads)-count, args.sequence_file_1)

    records_2 = (r for r in SeqIO.parse(args.sequence_file_2, "fastq") if r.id in target_reads)
    count = SeqIO.write(records_2, args.output[1], "fastq")
    print "Saved %i records from %s to %s" % (count, args.sequence_file_2, args.output[1])

    if count < len(target_reads):
        print "Warning %i IDs not found in %s" % (len(target_reads)-count, args.sequence_file_2)
else:
    print "output files not specified so not writing filtered sequences"

