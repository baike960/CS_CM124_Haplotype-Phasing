from Generate_Phases import list_of_phases_for_genotype
from Generate_Phases import finding_haplotype_complement
import numpy as np
import time

def pair_probability(haplotype, haplotype_complement, haplotype_dict):
    probability_of_h = haplotype_dict[' '.join(str(x) for x in haplotype)]
    probability_of_complement = haplotype_dict[' '.join(str(x) for x in haplotype_complement)]
    haplotype_pair_probability = probability_of_h * probability_of_complement
    return haplotype_pair_probability

def construct_all_individual_phase_table(population_genotype):
    geno_list_with_phases = []
    for individuals_genotypes in population_genotype:
        individual_phases = list_of_phases_for_genotype(individuals_genotypes)
        geno_list_with_phases.append(individual_phases)
    return geno_list_with_phases

def initial_haplotype_probability_dict(geno_list_with_phases):
    haplotype_dict = dict()
    for individuals in geno_list_with_phases:
        for phases in individuals:
            for haplotypes in phases:
                haplo_key = ' '.join(str(x) for x in haplotypes)
                if  haplo_key not in haplotype_dict:
                # add this haplotype and its complement in to phase_dict
                    haplotype_dict[haplo_key] = 1.0
    return haplotype_dict

def EM(population_genotype, number_of_iterations):
    start = time.time()
    geno_list_with_phases = construct_all_individual_phase_table(population_genotype)

    #initialize the haplotype_probability_dictionary, phase_probability_table, haplotype 
    haplotype_dict = dict()
    haplotype_location_dict = dict()
    number_of_haplotypes = 0
    phase_probability_table = []
    for individuals in range(0, len(geno_list_with_phases)):
        phase_prob_list = []
        for phases in range(0, len(geno_list_with_phases[individuals])):
            for haplotypes in geno_list_with_phases[individuals][phases]:
                haplo_key = ' '.join(str(x) for x in haplotypes)
                if  haplo_key not in haplotype_dict:
                    haplotype_dict[haplo_key] = 1.0
                    number_of_haplotypes += 1
                if haplo_key not in haplotype_location_dict:
                    location = []
                    location.append(individuals)
                    location.append(phases)
                    haplotype_location_dict[haplo_key] = []
                    haplotype_location_dict[haplo_key].append(location)
                elif haplo_key in haplotype_location_dict:
                    location = []
                    location.append(individuals)
                    location.append(phases)
                    haplotype_location_dict[haplo_key].append(location)
            phase_prob_list.append([1.0])
        phase_probability_table.append(phase_prob_list)
   

    #initialize the probability of haplotypes in haplotype probability table
    for haplotype in haplotype_dict.iterkeys():
        haplotype_dict[haplotype] = 1.0/number_of_haplotypes
        

     
    #test
    #print "haplotype_probability_dict:"
    #for k, v in haplotype_dict.iteritems():
    #    print k, v

    #print "phase_probability table:"
    #print phase_probability_table

    #print "haplotype_location_dict"
    #for k, v in haplotype_location_dict.iteritems():
    #    print k, v

    for i in range(0, number_of_iterations):
        number_of_individual = len(geno_list_with_phases)
        for individuals in range(0, number_of_individual):
            sum_of_phases_probability = 0
            list_of_haplotype_pair_prob = []
            for phases in range(0, len(geno_list_with_phases[individuals])):
                list_phases = geno_list_with_phases[individuals][phases]
                phase_haplotype_probability = pair_probability(list_phases[0], list_phases[1], haplotype_dict)
                list_of_haplotype_pair_prob.append(phase_haplotype_probability)
                sum_of_phases_probability += phase_haplotype_probability
            for phases in range(0, len(geno_list_with_phases[individuals])):
                phase_probability = list_of_haplotype_pair_prob[phases] / sum_of_phases_probability
                phase_probability_table[individuals][phases] = phase_probability
        
        
        for haplotypes in haplotype_location_dict.iterkeys():
            haplo_location = haplotype_location_dict[haplotypes]
            expectation = 0
            for coordinate in range(0, len(haplo_location)):
                individual_pos = haplo_location[coordinate][0]
                phase_pos = haplo_location[coordinate][1]
                expectation += phase_probability_table[individual_pos][phase_pos]
            haplotype_dict[haplotypes] = expectation/(2*number_of_individual)

        #print "haplotype_probability_dict:"
        #for k, v in haplotype_dict.iteritems():
        #    print k, v

        #print "phase_probability table:"
        #print phase_probability_table

        #print "haplotype_location_dict"
        #for k, v in haplotype_location_dict.iteritems():
        #    print k, v


    phase_with_max_probability = []
    for individuals in range (0, number_of_individual):
        phase_highest = geno_list_with_phases[individuals][0]
        highest_probability = phase_probability_table[individuals][0]
        for phase in range(0, len(phase_probability_table[individuals])):
            current_probability = phase_probability_table[individuals][phase]
            if current_probability > highest_probability:
                phase_highest = geno_list_with_phases[individuals][phase]
                highest_probability = current_probability
        phase_with_max_probability.append(phase_highest)# new
    
    
    
        #phase_with_max_probability.append(phase_highest[0])
        #phase_with_max_probability.append(phase_highest[1])
    #phase_with_max_probability_transpose = np.array(phase_with_max_probability).transpose()

    
    end = time.time()
    print "EM took", end - start, "seconds."
    return phase_with_max_probability






if __name__ == "__main__":
    genotypesToPhase = [[2, 1, 1, 1, 0], [1, 0, 0, 0, 1], [2, 2, 2, 2, 1]]
    result =construct_all_individual_phase_table(genotypesToPhase)
    result = EM(genotypesToPhase, 16)
    print result


    

    
    




