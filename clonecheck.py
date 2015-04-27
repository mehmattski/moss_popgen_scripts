#	This script checks a file of genotypes for sets of matching clones.
#	!!! IT REQUIRES HAPLOID DATA!!!
#	Missing data is handled by determining which other clones are compatible, and listing them.
#
# Input file must be a comma-delimited file (.csv) with each genotype on a line
#	The first column must be the ID of that genotype, followed by the alleles
#
# One output file (_assign.csv) has each individual in a row, followed by the genotype
#	The genotype is followed by the ID number of each compatible clone
# Another output file (_clones.csv) lists only the unique genotypes and the ID assigned by the script.

# The second argument on the command line is T, the number of mismatched loci allowed.


# 1. As data is read in, put into dictionaries based on the number of blanks
# 	Dictionary Key is the sample ID, genotype is 0 index, empty 1 index
# 2. Data fed into dictionary with zero blanks is checked for uniqueness
# 	Unique genotypes added to a unique dictionary, which is keyed sequentially.
# 3. If a genotype matches one in the unique dictionary, the unique genotype key
# 	is added to a list in the 1 index of the original genotype
# 4. After datafile is completely read, move on to dictionary with one blank, then two, etc.
# 	Unique genotypes are still added to the unique dictionary
# 	May be multiple matches for genotypes, these are listed in the 1 index of that genotype.
# 5. Output is the dictionary keys followed by the 1 index (matching genotypes)
# 	followed by a list of unique genotypes.

import sys
unique_genotypes = {}
all_genotypes= {}

def parse_file(filename):
	class_genotypes = [{}]
	file = open(filename,"U")
	for line in file.xreadlines():
		line = line.rstrip()
		line=line.split(",")
		missing_data = line.count('0')
		while len(class_genotypes) < missing_data+1:
			class_genotypes.append({})
		correctDict = class_genotypes[missing_data]
		correctDict[line[0]] = [0,0]
		correctDict[line[0]][0] =line[1:]
	return class_genotypes
		
		
def uniquecheck(genotype,T):
	match_ID = []
	unique_count = len(unique_genotypes)
	if unique_count == 0:
		unique_genotypes[unique_count]=genotype
	for each in range(unique_count):
		match = 0
		for i in range(len(genotype)):
			if unique_genotypes[each][i] != '0':
				if genotype[i] != '0':
					locus_dist = abs(int(unique_genotypes[each][i]) - int(genotype[i]))
					match += locus_dist
			if i+1 == len(unique_genotypes[each]):
				if match <= T:
					match_ID.append(str(each))
# 			if genotype[i] == unique_genotypes[each][i] or genotype[i] == '0' or unique_genotypes[each][i] == '0':
# 				match += 1
# 			if i+1 == len(unique_genotypes[each]):
# 				if match >= i+1-T:
# 					match_ID.append(str(each))
	if len(match_ID) == 0:
		unique_genotypes[unique_count] = genotype
		match_ID.append(str(unique_count))
	return match_ID
		
def feeder(genotype_list,T):
	missingDataClasses = len(genotype_list)
	for j in range(missingDataClasses):
		if len(genotype_list[j]) > 0:
			for genotype in genotype_list[j]:
				genotype_list[j][genotype][1]= uniquecheck(genotype_list[j][genotype][0],T)
				all_genotypes[genotype] = genotype_list[j][genotype]

def main():
	if len(sys.argv) < 3:
		print "Usage: python clonecheck.py inputfile mismatchLoci"
		sys.exit(1)
	infile = sys.argv[1]
	T = int(sys.argv[2])
	infile_stem = infile.split('.')[0]
	assignfile_stem = infile_stem + '_assign.csv'
	clonefile_stem = infile_stem + '_clones.csv'
	assignfile = open(assignfile_stem, 'w')
	source_genotypes = parse_file(infile)
	feeder(source_genotypes,T)
	for item in sorted(all_genotypes.iterkeys()):
		output_list = [item]
		for each in all_genotypes[item][0]:
			output_list.append(each)
		for thing in all_genotypes[item][1]:
			output_list.append(thing)
		output_string = ','.join(output_list)+'\n'
		assignfile.write(output_string)
	assignfile.close()
	clonefile = open(clonefile_stem,'w')
	for item in sorted(unique_genotypes.iterkeys()):
		outputstring = str(item) + ','+ ','.join(unique_genotypes[item]) + '\n'
		clonefile.write(outputstring)

if __name__ == '__main__': main()		
