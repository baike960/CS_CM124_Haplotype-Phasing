from play import construct_genotype
from Generate_Phases import list_of_phases_for_genotype
from Generate_Phases import finding_haplotype_complement
from Stephanie_EM import initial_haplotype_probability_dict
from Stephanie_EM import EM
from Stephanie_EM import construct_all_individual_phase_table
from Stephanie_EM import pair_probability
import numpy as np
import time


def overlap(chunk_phases, lastChunk,  windowSize):
    knownHaplo=[]
    #print "size of the chunk is: ", len(chunk_phases)
    #print "size of the last chunk is: ", len(lastChunk)
    phase_with_max_probability=[]
    for i in range(len(chunk_phases)):
        # the first part of current chuck of each individual
        individual=chunk_phases[i]
        haplo=individual[0]
        haplo = haplo[0:windowSize / 2]
        haplo_comp=individual[1]
        haplo_comp=haplo_comp[0:windowSize / 2]

        # the last part of last chunk of each individual

        lastIndividual=lastChunk[i]
        lastHaplo = lastIndividual[0]
        lastHaplo = lastHaplo[len(lastHaplo)-(windowSize / 2):len(lastHaplo)]
        lastHaplo_comp = lastIndividual[1]
        lastHaplo_comp = lastHaplo_comp[len(lastHaplo_comp)-(windowSize / 2):len(lastHaplo_comp)]

        #overlap region1
        overlapRegion1=lastHaplo+haplo
        overlapRegion_comp1=lastHaplo_comp+haplo_comp

        # overlap region2
        overlapRegion2 = lastHaplo + haplo_comp
        overlapRegion_comp2 = lastHaplo_comp + haplo

        if i==0:
            knownHaplo.append(overlapRegion1)
            knownHaplo.append(overlapRegion_comp1)
        else:
            maxCount1=0
            maxCount2 = 0
            for haplo in knownHaplo:
                count1=0
                count_comp1=0
                count2 = 0
                count_comp2 = 0
                for j in range(len(haplo)):
                    if overlapRegion1[j]==haplo[j]:
                        count1+=1
                    if overlapRegion_comp1[j]==haplo[j]:
                        count_comp1+=1
                    if overlapRegion2[j]==haplo[j]:
                        count2+=1
                    if overlapRegion_comp2[j]==haplo[j]:
                        count_comp2+=1
                '''
                print "count1 is: ", count1
                print "count_comp1 is: ", count_comp1
                print "count2 is: ", count2
                print "count2_comp is: ", count_comp2
                '''
                if count1> maxCount1:
                    maxCount1=count1
                if count_comp1>maxCount1:
                    maxCount1=count_comp1
                if count2> maxCount2:
                    maxCount2=count2
                if count_comp2>maxCount2:
                    maxCount2=count_comp2

            # if exchanging the haplotypes gives more overlapping region with previous known haplotypes,
            # exchanging them and add them to known haplotypes list for latter individuals
            if maxCount2> maxCount1:
                #print "exchange!"

                temp=individual[0]
                individual[0]=individual[1]
                individual[1]=temp
                knownHaplo.append(overlapRegion2)
                knownHaplo.append(overlapRegion_comp2)
            else:
                knownHaplo.append(overlapRegion1)
                knownHaplo.append(overlapRegion_comp1)
        #print "size of known overlapp region: ", len(knownHaplo)

        phase_with_max_probability.append(individual[0])
        phase_with_max_probability.append(individual[1])
    #print "number of exchange is: ", numCross
    #print "size of phase with max probability: ", len(phase_with_max_probability)
    return phase_with_max_probability







def cut_em_ligation(genotype, cut_size, num_iteration, overlapSize, filename):
    output_file = open(filename, "w")
    start_time = time.time()
    ligation_unfinished = True
    start_pos = 0
    lastChunk=[]#Last window's chunk answer
    
    while ligation_unfinished:
        chunk_end_pos = start_pos + cut_size
        print "end index is: ", chunk_end_pos
        if chunk_end_pos >= len(genotype[0]):
            chunk_end_pos = len(genotype[0])
            ligation_unfinished = False

        # cut the genotypes into chunks for em
        chunk_geno = []
        for ind_geno in genotype:
            chunk_geno.append(ind_geno[start_pos:chunk_end_pos])

        chunk_phases = EM(chunk_geno, num_iteration)
        #print "size of original chunk is: ", len(chunk_phases)
        # Chooes the way of arranging the 2 haplotypes by using Clark algorithm on overlapping region of size k

        phase_with_max_probability=[]
        if start_pos!=0:
            phase_with_max_probability=overlap(chunk_phases, lastChunk, overlapSize)
        else:
            for individual in chunk_phases:
                phase_with_max_probability.append(individual[0])
                phase_with_max_probability.append(individual[1])

        lastChunk=chunk_phases# Update the last chunk
        phase_with_max_probability_transpose=np.array(phase_with_max_probability).transpose()
        
        # each iteration will ligate the next chunk to the ligated chunks in the output file
        for individual in phase_with_max_probability_transpose:
            output_line = ' '.join(str(x) for x in individual)
            output_file.write(output_line + '\n')

        start_pos = chunk_end_pos

    output_file.close()

    end_time = time.time()
    print "This algorithm took", end_time - start_time, "seconds."

        


if __name__ == "__main__":
    transformed_genotypes = construct_genotype("./test_data_2.txt")
    cut_em_ligation(transformed_genotypes, 16, 3, 32, "test_data_2_sol.txt")#window size, iterations, overlap size
