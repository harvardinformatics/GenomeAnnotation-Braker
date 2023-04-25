import argparse
fields = ["seqid", "source", "type", "start",
          "end", "score", "strand", "phase", "attributes"]
          
def ParseBrakerGtfAttributes(attributes):
    attribute_dict = {}
    attribute_list = attributes[:-1].replace('"','').split(';')
    for attribute in attribute_list:
        key,value = attribute.split()
        if key in ['transcript_id','gene_id'] and 'file' in value:
            value = value.split('_')[-1]
        attribute_dict[key] =  value
    return attribute_dict
    

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="Generate CDS transcript and gene interval bed without UTRs")
    parser.add_argument('-supported','--hint-supported-gtf',dest='supported',type=str,help='braker support eval script generated gtf with any support')
    parser.add_argument('-gtf','--braker-gtf',dest='braker',type=str,help='braker support eval script generated gtf with any support')
    parser.add_argument('-o','--output-gtf',dest='output',type=str,help='merge of GeneMark and supported AUGUSTUS annotations')
    
    opts = parser.parse_args()
   
    supported_genes = set()
    supported_tscripts = set()
    supported = open(opts.supported,'r')
    for line in supported:
        if line[0] !='#':
            linedict = dict(zip(fields,line.strip().split('\t')))
            if linedict['type'] in ['intron','start_codon','stop_codon']:
                attribute_dict = ParseBrakerGtfAttributes(linedict['attributes'])
                if attribute_dict['supported'] == 'True':
                    supported_genes.add(attribute_dict['gene_id'])
                    supported_tscripts.add(attribute_dict['transcript_id']) 
    #print(supported_genes)
 
    brakerin = open(opts.braker,'r')
    fout = open(opts.output,'w')
    for line in brakerin:
        if line[0] == '#':
            fout.write(line)
        elif line =='\n':
            pass
        elif "GeneMark" in line:
            fout.write(line)
        else:
            linedict = dict(zip(fields,line.strip().split('\t')))
            if linedict['type'] == 'gene':
                if linedict['attributes'].split('_')[-1] in supported_genes:
                    fout.write(line)
            elif linedict['type'] == 'transcript':
                if linedict['attributes'].split('_')[-1] in supported_tscripts:
                    fout.write(line)
            elif linedict['type'] in ['exon','CDS','stop_codon','start_codon','intron']:
                attribute_dict = ParseBrakerGtfAttributes(linedict['attributes'])
                if attribute_dict['gene_id'] in supported_genes and attribute_dict['transcript_id'] in supported_tscripts:
                    fout.write(line)
            else:
                raise ValueError('%s not in list of valid feature types' % line_dict['type'])

    fout.close() 
