from __future__ import division
import sys
import random

'''This script determines the haploid paternity of diploid offspring using the subtraction method.
For heterozygotes, the paternal allele is the one not present in the maternal individual.
For homozygotes, the paternal alelle is assumed to be the same as teh maternal allele.
This method should work for any codominant genetic marker (such as microsatellites or SSRs).

You must provide the multilocus genotypes of the haploid maternal parent and the diploid offspring.

The file is a comma-separated-value (CSV) format
The first line of the document should be comma-separated list of loci names, beginning with "Loci":

	Loci,17a,65a,78a,1a,7a,68a,10a,18a,19a,9a,20a,14a,56a

The multilocus genotypes begin on the next line. The first field is the individual ID. 
The second field is the ID of the maternal parent.
Then there should be two columns per locus. For maternal genotypes, simply list each allele twice.

Example Maternal Genotype:
	MJ872,,158,158,193,193,191,191,254,254,190,190,219,219,249,249,129,129,273,273,180,180,290,290,195,195,198,198

Note that the second field is blank. This identifies this individual as a maternal haploid.

Example Offspring Genotype:
	MJ873,MJ872,158,158,193,193,191,191,254,254,182,190,219,219,249,249,129,129,273,276,180,180,290,293,195,195,198,198

The second field identifies that this individual (MJ873) has individual MJ872 as a maternal haploid parent.

On the command line, enter the genotype CSV file name and the number of loci.

The filestem of the input file (before the .csv) is used to name the output files.

Two output files are generated for each input:

inputfilestem_het.csv 
	outputs the percent heterozygosity for each diploid offspring.
inputfilestem_paternity.csv  
	contains the haploid multilocus genotype of all maternal and inferred paternal individuals.

The script will print an error if a maternal genotype contains heterozygous loci.
The script will also list all offspring for which the maternal allele cannot be found.
'''


class load_matrix():
	def __init__(self,name):
		self.genotype={}
		self.name=name
	def mom(self,DNA_ID):
		self.mom=DNA_ID
	def add_alleles(self,locus,alleles):
		
		self.genotype[locus]=alleles
	def het_loci(self):
		self.num_het_loci = 0
		self.num_loci = 0
		for x in self.genotype:
			if self.genotype[x][0] != '0':		
				if self.genotype[x][0] != self.genotype[x][1]:
					self.num_het_loci+=1
				if self.genotype[x][0] != '0':
					self.num_loci+=1
		if self.num_loci >0:
			self.percent_het = self.num_het_loci/self.num_loci
		else:
			self.percent_het = 0
		self.het_loci_stats = (self.num_loci, self.num_het_loci, self.percent_het)
	
def add_loci(locus_list):
	loci = {}
	for x in range(len(locus_list)):
		loci[x]=locus_list[x]
	return loci
	
def parse_file(filename,nLoci):
	males = False
	ID_List={}
	new_matrix = []
	call_ID_List = []
	line_number=0
	file = open(filename,"U")
	for line in file.xreadlines():
		if "male" in line:
			males = True
		line = line.rstrip()
		line=line.split(",")
		if line[0] == "Loci":
			line.pop(0)
			line = line[0:nLoci]
			loci_IDs=add_loci(line)
		else:
			ID_List[line_number]=line[0]
			callID = ID_List[line_number]
			callID = load_matrix(ID_List[line_number])
			call_ID_List.append(callID)
			if line[1]:
				mom_ID=line[1]
				callID.mom(mom_ID)
			else:
				callID.mom(None)
			alleles = line[2:nLoci*2+2]
			for x in range(0,len(alleles),2):
				callID.add_alleles(loci_IDs[int(x/2)],[alleles[x],alleles[x+1]])
			callID.het_loci()
			new_matrix.append(callID.__dict__)
			line_number +=1
	output_stem = filename.split(".")[0]
	
	if males == True:
		nPops = 3
	else:
		nPops = 2
	#genodive_format(nLoci,len(call_ID_List),nPops,output_stem)
	
	return call_ID_List
	file.close()

def genodive_format(nLoci,nIndiv,nPops,filestem):
	filestem = filestem + "_paternity.csv"
	file=open(filestem,'a')
	file.write(filestem+"\n")
	outputs = [str(nIndiv),str(nPops),str(nLoci),"1\t3","\n"]
	outputstring="\t".join(outputs)
	file.write(outputstring)
	if nPops == 3:
		file.write("dads\nmoms\nmales\n")
	else:
		file.write("dads\nmoms\n")
	file.close()


def find_mom(mom_ID,objects):
	for item in objects:
		if item.name == mom_ID:
			return item

def find_het_loci(matrix):
	for item in matrix:
		individual_genotype = item.genotype
		if item.mom:
			mother = find_mom(item.mom,matrix)
			mother_genotype=mother.genotype
			for locus in mother_genotype:
				if mother_genotype[locus] != individual_genotype[locus]:
					if mother_genotype[locus][0] != '0':
						if individual_genotype[locus][0] != '0':
							print item.name, "is heterozygous at ", locus, " !"
							print mother_genotype[locus],individual_genotype[locus]
		if not item.mom:
			print "female"

def paternity(matrix):
	paternal_genotypes ={}
	mismatch_list = {}
	for item in matrix:
		mismatch_loci = []
		father_genotype = {}
		individual_genotype = item.genotype
		if item.mom:
			mother = find_mom(item.mom,matrix)
			if item.mom != "male":
				mother_genotype=mother.genotype
				
				for locus in mother_genotype:
					if mother_genotype[locus] != individual_genotype[locus]:
						if mother_genotype[locus][0] != '0':
							if individual_genotype[locus][0] != '0':
								for allele in mother_genotype[locus]:
									if allele not in individual_genotype[locus]:
										mismatch_loci.append(locus)
										mismatch_list[item.name]=mismatch_loci
										father_genotype[locus]="0"
										pass
								for allele in individual_genotype[locus]:
									if allele not in mother_genotype[locus]:
										father_genotype[locus]=allele
										#father_genotype[1]="1"
										match = False
								
							else:
							
								father_genotype[locus]="0"
								#father_genotype[1]="1"
						else:
							
							father_genotype[locus]="0"
							#father_genotype[1]="1"
					else: 
						father_genotype[locus]=mother_genotype[locus][0]
						#father_genotype[1]="1"
			else:
				father_genotype = individual_genotype.copy()
				for locus in father_genotype:
					father_genotype[locus]=father_genotype[locus][0]
				#father_genotype[1]="3"
		else:
			father_genotype = individual_genotype.copy()
			for locus in father_genotype:
				father_genotype[locus]=father_genotype[locus][0]
			#father_genotype[1]="2"
		DNAID = item.name
		paternal_genotypes[DNAID]=father_genotype
	for item in sorted(mismatch_list.iterkeys()):
		print item, mismatch_list[item]
	return paternal_genotypes


	
def main():
	if len(sys.argv) < 3:
		print "Usage: python paternity.py genotypes_file numLoci"
		sys.exit(1)
		
	infile = sys.argv[1]
	infile_stem = infile.split('.')[0]
	nLoci = int(sys.argv[2])


	matrix = parse_file(infile,nLoci)	
	hetfilestem = infile_stem + "_het.csv"
	hetfile = open(hetfilestem,'w')
	for item in matrix:
		stats = item.het_loci_stats
		outputstring = item.name + "," + str(stats[0]) + "," + str(stats[1]) + "," + str(stats[2])+"\n"
		hetfile.write(outputstring)
	hetfile.close()
	
	paternal_genotypes =paternity(matrix)
	
	patfilestem = infile_stem + "_paternity.csv"
	patfile = open(patfilestem, 'w')

	loci = open(infile,"U").readline().rstrip().split(",")[1:nLoci+1]
	print "loci to write: ",loci
	#patfile.write(","+",".join(loci)+"\n")

	
	for item in sorted(paternal_genotypes.iterkeys()):
		output_list = [item]
		genotype =  paternal_genotypes[item]
		for locus in loci:
				output_list.append(genotype[locus])
		DNAID = output_list.pop(0)
			
		output_list.insert(0,DNAID)
		output_string = ",".join(output_list)+"\n"
		patfile.write(output_string)
	patfile.close()

if __name__ == '__main__': main()		