import numpy as np 

#create a list of genotypes from the file provided
def construct_genotype(file) :
    with open(file, 'r') as file_text:
        Population_Genotype = []
        for line in file_text:
            line_list = line.strip('\n').split(' ')
            int_line_list = map(int, line_list)
            for each_person in range(0, len(line_list)):
                if len(Population_Genotype) < len(line_list):
                    Population_Genotype.append([])
                Population_Genotype[each_person].append(int_line_list[each_person])
    file_text.close()
    return Population_Genotype


if __name__ == '__main__':
    
    file = "./test.txt"
    print "---------- Loading genotypes from file into memory ----------"
    genome = construct_genotype(file)
    for i in range(0, 50):
        print genome[i]
        
        

