
## A script to extract the major allele, total read count, and major allele read count
## from a pileup file with results from multiple bam files.
## It will determine the major allele and summarize the total reads across all included samples/cells

# import packages
import argparse
import collections
import gzip

# a function to write to the new file
def write_line(out, line, reads, major_allele, major_allele_counts, format = '%s\t%s\t%s\t%i\t%i\n'):
    out.write(format % (
        line[0], line[1], major_allele.upper(), reads, major_allele_counts))

# a function to write a new file with info on minor allele reads across cells for variable sites
def write_cell_line(out, line, cell_string):
    out.write('\t'.join([line[0], line[1]] + cell_string + ['\n']))

def pileup_read_summary(file_path, out_file, cell_file):
    out_header = ['chromosome', 'position', 'pooled_major_allele',
                'read_count', 'pooled_major_allele_read_count']
    counter = 1
    #with open(file_path, 'rt') as pileup_file, open(out_file, 'wt') as out, open(cell_file, 'wt') as cell_out:
    with gzip.open(file_path, 'rt') as pileup_file, open(out_file, 'wt') as out, open(cell_file, 'wt') as cell_out:
        # initiate output file
        out.write('%s\n' % '\t'.join(out_header))
        
        # get the number of lines in the file, and select relevent columns
        firstline = pileup_file.readline().rstrip().split('\t')
        l = [i for i in  range(len(firstline))]
        read_columns = [i for i in l if i > 2 and i % 3 == 0]
        allele_columns = [j for j in l if j > 2 and j % 3 == 1]

        # initiate file that will keep track of minor allele counts in each cell
        cell_header = ['chromosome', 'position'] + [str(i) for i in range(len(read_columns))]
        cell_out.write('%s\n' % '\t'.join(cell_header))

    #### note THIS IS CURRENTLY SKIPPING THE FIRST LINE IN THE OUTPUT FILE

        # iterate through lines and extract the read data
        for line in pileup_file:
            # to write out the first n lines:
            # counter += 1
            # if counter == 1000000:
            #     break
            line = line.strip().split('\t')
            if line[2].upper() == "N":
                print(line)
                continue
            reads = sum([int(line[i]) for i in read_columns])
            string = "".join([line[i] for i in allele_columns]).upper()
            # remove lines with indels:
            if ('+' in string or '-' in string or 'N' in string):
                write_line(out, line, reads, line[2], -1)
            # count up how many reads match the major allele
            else:
                counts = collections.Counter(string)
                matches = counts["."] + counts[","]
                # only include characters that have read meaning
                counts = {key:value for key, value in counts.items() if key in ['A', 'C', 'T', 'G', ',', '.']} # Mon 16 May 2022 11:32:34 AM PDT TODO should probably replace all instances of [,.] with the major allele, as if there are more than 2 alleles this outputs [,.] for the major allele
                # remove lines with have extra or fewer reads due to indels
                if sum(counts.values()) != reads:
                    write_line(out, line, reads, line[2], -1)
                # if major allele is reference allele
                elif matches >= reads / 2:
                    write_line(out, line, reads, line[2], matches)
                    # count up the number of reads to minor allele in each cell
                    if matches != reads:
                        cell_string = ["" if line[i] == '*' else str(line[i].upper().count("A") + line[i].upper().count("T") + line[i].upper().count("G") + line[i].upper().count("C")) for i in allele_columns]
                        write_cell_line(cell_out, line, cell_string)
                else:
                    # otherwise find the major allele and corresponding number of reads
                    # first check if there are multiple nucleotides with the same max value
                    try:
                        major_allele = max(counts, key=counts.get)
                        write_line(out, line, reads, major_allele, counts[major_allele])
                        # again count up the number of reads to minor allele in each cell
                        if counts[major_allele] != reads:
                            cell_string = ["" if line[i] == '*' else str(line[i].upper().count("A") + line[i].upper().count("T") + line[i].upper().count("G") + line[i].upper().count("C") + line[i].upper().count(",") + line[i].upper().count(".") - line[i].upper().count(major_allele)) for i in allele_columns]
                            write_cell_line(cell_out, line, cell_string)

                    # if anything else weird happens, skip that line (write out a -1 for major allele reads)
                    except:
                        write_line(out, line, reads, line[2], -1)

                    
                    

parser = argparse.ArgumentParser(description='Get file path and other info')
parser.add_argument('-f', '--file_path', type=str, 
                    help='path to pileup file', required = True)
parser.add_argument('-i', '--id', type=str, 
                    help='Chromosome', required = False)
parser.add_argument('-t', '--type', type=str, 
                    help='Type (Diploid or Tumor)', required = True)
args = vars(parser.parse_args())

#out_file = ''.join(['/space/s2/sarahj32/rotation/fixed_read_counts/', args['type'], '_chr_', args['id'], '_read_allele_counts.tsv'])
#cell_file = ''.join(['/space/s2/sarahj32/rotation/fixed_read_counts/cell_summary', args['type'], '_chr_', args['id'], '.tsv'])
out_file = ''.join(['/space/s2/sarahj32/rotation/pooledMajorAlleleReadCounts/', args['type'], '_chr_', args['id'], '_read_allele_counts.tsv'])
cell_file = ''.join(['/space/s2/sarahj32/rotation/pooledMajorAlleleReadCounts/cell_summary', args['type'], '_chr_', args['id'], '.tsv'])

pileup_read_summary(args['file_path'], out_file, cell_file)

# testing:
# file_path = '/space/s2/sarahj32/rotation/chromosome_pileup_files/tumor_chr_2_file'
# type = 'tumor'
# id = '2'

# out_file = ''.join(['/space/s2/sarahj32/rotation/fixed_read_counts/', type, '_chr_', id, '_read_allele_counts.tsv'])
# cell_file = ''.join(['/space/s2/sarahj32/rotation/fixed_read_counts/cell_summary', type, '_chr_', id, '.tsv'])

# pileup_read_summary(file_path, out_file, cell_file)
