import numpy as np 

def generating_half_of_haplotypes (given_genotype):
    num_of_1 = 0
    half_of_haplotypes = []
    if given_genotype[0] == 2:
        half_of_haplotypes.append([1])
    elif given_genotype[0] == 0:
        half_of_haplotypes.append([0])
    else: 
        half_of_haplotypes.append([1])
        num_of_1 += 1

    for SNP_pos_i in range(1, len(given_genotype)):
        if given_genotype[SNP_pos_i] == 2:
            for h in half_of_haplotypes:
                h.append(1)
        elif given_genotype[SNP_pos_i] == 0:
            for h in half_of_haplotypes:
                h.append(0)
        else:
            num_of_1 += 1
            if num_of_1 == 1:
                for h in half_of_haplotypes:
                    h.append(1)
            else:
                copy_haplotypes = []
                for h in half_of_haplotypes:
                    copy = list(h)
                    copy.append(0)
                    h.append(1)
                    copy_haplotypes.append(copy)
                for copies in copy_haplotypes:
                    half_of_haplotypes.append(copies)
    return half_of_haplotypes

def finding_haplotype_complement(haplotype, given_genotype):
    h_c = []
    for SNPs_i in range(0, len(given_genotype)):
        h_c.append(given_genotype[SNPs_i]-haplotype[SNPs_i])
    return h_c

def list_of_phases_for_genotype (genotype):
    phase_list = []
    haplotypes = generating_half_of_haplotypes(genotype)
    for h in haplotypes:
        phase_pair = []
        phase_pair.append(h)
        phase_pair.append(finding_haplotype_complement(h, genotype))
        phase_list.append(phase_pair)
    return phase_list

        

if __name__ == '__main__':
    haplotypes = generating_half_of_haplotypes([2,0,1,1,0,1,2,0,2,0,0,1])
    for h in haplotypes:
        print h
    print "\n"

    print list_of_phases_for_genotype([2,0,1,1,0,1,2,0,2,0,0,1])
    

