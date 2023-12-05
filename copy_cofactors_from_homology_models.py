"""This script takes homology models and AI-predicted modles with similar sequences and copies fitting ligands from the homology model into the AI-predicted model. Run with >>> pymol -c copy_cofactors_from_homology_models.py"""

import os
import pymol
from Bio.PDB import PDBParser
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

### set maximum RMSD for alignment of ligand-aligning residues
max_rmsd = 3  # max rmsd for high-quality alignment
rmsd_cutoff = 7  # max rmsd for ligand-copying (if RMSD higher than that, the ligand is not incorporated)

input_dir = './'
output_dir = './'

cofactors_to_neglect = []
missing_models = []
failed_to_open = []

logfile = open('./log.txt', 'w')

for items in os.listdir(input_dir): 

    if items.startswith('HM_'):

        identifier = items.split('HM_')[-1].split('.pdb')[0]
        print(identifier)
        logfile.write('Working on {} ...\n'.format(identifier))
        hm = input_dir + 'HM_{}.pdb'.format(identifier)
        try:
            pred = input_dir + '{}_AF2.pdb'.format(identifier)

            ### parse sequences from PDB files

            d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

            parser = PDBParser(QUIET=True)
            try:
                homologymodel = parser.get_structure('struct', hm)
            except ValueError: 
                failed_to_open.append(identifier)
            
            number_of_hetatms = []  # if the it is a homo-multimer and some chains include more ligands than others, we want to align with the chain containing the most ligands
            sequences_hm = []
            chains_hm = []  # list of chains within the homology model
            for model in homologymodel:
                for chain in model:
                    
                    hetatm_count = sum(1 for residue in chain if residue.id[0].startswith('H_') is True)
                    number_of_hetatms.append(hetatm_count)
                    chains_hm.append(chain.id)
                    seq = []
                    for residue in chain:
                        if residue.resname in d3to1:
                            seq.append(d3to1[residue.resname])
                    sequences_hm.append(''.join(seq))

            
            logfile.write('Chains in HM: {}; lengths of chains: {}; number of cofactors: {}\n'.format(','.join(chains_hm),','.join(map(str, [len(i) for i in sequences_hm])), ','.join(map(str, number_of_hetatms))))
            if max(number_of_hetatms) == 0: 
                logfile.write('No cofactors to add here.\n')
            else: 
            
                prediction = parser.get_structure('struct', pred)
                for model in prediction:
                    for chain in model:
                        seq = []
                        for residue in chain:
                            if residue.resname in d3to1:
                                seq.append(d3to1[residue.resname])
                
                                
                pred_seq = ''.join(seq)  # sequence of the de novo structure prediction (should be only one chain)

                
                ### alignment of each sequence in the homology model with the de novo model sequence

                # set alignment parameters
                gap_open = -10
                gap_extend = -8
                match_score = 5
                mismatch_score = -5

                scores = []  # the scores are used for ranking the sequence alignments later
                chains_hm2 = []
                number_of_hetatms2 = []
                
                
                for seqs, ch, nr in zip(sequences_hm, chains_hm, number_of_hetatms):
                    if seqs != '':
                        all_alignments = []
                        
                        seq2 = Seq('{}'.format(seqs))
                        alignments = pairwise2.align.globalms(pred_seq, seq2, match_score, mismatch_score, gap_open, gap_extend, one_alignment_only=True)
                        

                        all_alignments.append(alignments[0])
                        logfile.write('Aligned chain {} with a score of {}.\n'.format(ch, alignments[0].score))
                        scores.append(alignments[0].score)
                        number_of_hetatms2.append(nr)
                        chains_hm2.append(ch)

                # put alignment scores and the hetatms to dataframe for sorting
                df_alignmentscores = pd.DataFrame(list(zip(chains_hm2, number_of_hetatms2, scores)), columns=['Chain', 'NrHetatm', 'Score'])
                df_alignmentscores.sort_values(by=['Score', 'NrHetatm'], ascending=False, inplace=True)  # first sort for largest nr of hetatm and second sort for best alignment, in order to get, if more than one sequence matches with the best score, the one with the most hetatms

                

                chain_to_align = df_alignmentscores['Chain'].tolist()[0]  # this is the chain which best aligns to the de novo sequence
                logfile.write('Chain of choice is {}.\n'.format(chain_to_align))

                # here comes the pymol part
                pymol.finish_launching()

                output_list = []  # list of hetatms to consider
                pymol.cmd.load(hm, 'reference')  # load the homology model
                pymol.cmd.select('check_if_nucleic_acid', '(byres polymer and name P) or (resn A and (not hetatm)) or (resn C and (not hetatm)) or (resn G and (not hetatm)) or (resn U and (not hetatm)) or (resn DA and (not hetatm)) or (resn DC and (not hetatm)) or (resn DG and (not hetatm)) or (resn DT and (not hetatm))') 
                pymol.cmd.remove('reference and (not chain {}) and not (byres chain {} around 5 and hetatm)'.format(chain_to_align, chain_to_align))  # remove all chains from the homology model which are not part of the de novo model sequence in the de novo model but leave hetatms which align to the homology model atoms
                pymol.cmd.load(pred, 'current')  # load the de novo model


                pymol.cmd.indicate('hetatm and {}'.format('reference')) # indicate the hetatms in the homology model
                pymol.cmd.iterate("hetatm and {}".format('reference'), "output_list.append((resi, resn, chain))", space={"output_list": output_list})  # loop through the hetatms in the homology model
                 
                output_list2 = list(set(output_list))
                
                # check if there are amino acids in the homology model which are defined as hetatm, although they belong to the sequence. Such "hetatms" should not be copied and thus removed from the output_list
                print(output_list2)
                
                pred_seq_dict = dict(zip(list(range(1, len([*pred_seq])+1)), [*pred_seq]))
                
                for item in output_list2:
                    
                    if item[1] in d3to1:
                        key,value = int(item[0]),d3to1[item[1]]
                        if key in pred_seq_dict and value == pred_seq_dict[key]:
                            output_list2.remove(item)
                print(output_list2)
                    

                output_list3 = []
                for entry in output_list2:
                    if entry[1] not in cofactors_to_neglect:
                        output_list3.append(entry)
                
                # remove bases if there where HETATM from list of hetatms, unless they are base analogs (neither A, G, T or C)
                if pymol.cmd.count_atoms('check_if_nucleic_acid') != 0: 
                    out3 = []
                    for i in output_list3: 
                        if i[1] != 'A' and i[1] != 'G' and i[1] != 'C' and i[1] != 'T' and i[1] != 'DA' and i[1] != 'DG' and i[1] != 'DC' and i[1] != 'DT':
                            out3.append(i)
                    output_list3 = out3
                
                pymol.cmd.remove('(byres polymer and name P) or (resn A and (not hetatm)) or (resn C and (not hetatm)) or (resn G and (not hetatm)) or (resn U and (not hetatm)) or (resn DA and (not hetatm)) or (resn DC and (not hetatm)) or (resn DG and (not hetatm)) or (resn DT and (not hetatm))')  # remove the DNA, if it is within the chain_to_align (because it would interfere with the alignment) and do not add it as hetatm, because DNA is converted to "ATOM" if YASARA saves it as pdb and upon upload to the platform -  therefore, we would lose some cavities 
                
                

                if output_list3:  # check if list is empty or not 
                    cofactor_added = False
                    problematic_cofactors = False
               
                    for ind_i, i in enumerate(output_list3):
                    
                        # initial check, if the hetatm is within 5 A of reference
                        pymol.cmd.select('check_distance_hetatms', '{} and (not hetatm) within 5 of (resn {} and resi {} and chain {})'.format('reference', i[1], i[0], i[2]))  # here not byres, because any interaction with the protein is counted
                        if pymol.cmd.count_atoms('check_distance_hetatms') == 0:  # hetatm is too far away
                            pass
                        else: # now for the alignment
                            logfile.write('Trying to add ligand {}-{}-{}.\n'.format(i[0], i[1], i[2]))
                            rmsds = []
                            cutoffs = [12, 12, 10, 10, 5, 5, 15, 15]
                            align_method = ['align', 'super', 'align', 'super', 'align', 'super', 'align', 'super']
                            max_attempts = len(cutoffs)
                            attempt = 0
                            # check the RMSD and go through several options, if RMSD too high
                            while attempt < max_attempts: 
                                cutoff = cutoffs[attempt]
                                method = align_method[attempt]
                                pymol.cmd.select('near_current_hetatms', 'byres {} and (not hetatm) within {} of (resn {} and resi {} and chain {})'.format('reference', cutoff, i[1], i[0], i[2]))
                                        
                                if method == 'align':
                                    pymol.cmd.align('reference and near_current_hetatms', 'current')
                                elif method == 'super':
                                    pymol.cmd.super('reference and near_current_hetatms', 'current')
                                
                                rmsd = (pymol.cmd.rms_cur('reference and near_current_hetatms', 'current'))
                                rmsds.append(rmsd)
                                if rmsd <= max_rmsd:
                                    logfile.write('Final RMSD = {}, with method {} and {} A cutoff.\n'.format(rmsd, method, cutoff))
                                    break
                                
                                attempt += 1
                                    
                            if min(rmsds) > max_rmsd and min(rmsds) > rmsd_cutoff:   
                                logfile.write('When copying this ligand, the max rmsd cutoff of {} A was exceeded: RMSD = {} A.\n'.format(rmsd_cutoff, min(rmsds))) 
                            if min(rmsds) > max_rmsd and min(rmsds) <= rmsd_cutoff:
                                problematic_cofactors = True
                                min_index = rmsds.index(min(rmsds))
                                pymol.cmd.select('near_current_hetatms', 'byres {} and (not hetatm) within {} of (resn {} and resi {} and chain {})'.format('reference', cutoffs[min_index], i[1], i[0], i[2]))
                                
                                if align_method[min_index] == 'align':
                                    pymol.cmd.align('reference and near_current_hetatms', 'current')
                                elif align_method[min_index] == 'super':
                                    pymol.cmd.super('reference and near_current_hetatms', 'current')
                                rmsd = (pymol.cmd.rms_cur('reference and near_current_hetatms', 'current'))
                                logfile.write('Final problematic RMSD = {} with method {} and cutoff {} A\n'.format(rmsd, align_method[min_index], cutoffs[min_index]))
                                
                            if min(rmsds) <= rmsd_cutoff:      
                                pymol.cmd.create('current', '(resn {} and resi {} and chain {}) or current'.format(i[1], i[0], i[2]))
                                cofactor_added = True
                                
                    if cofactor_added is True: 
                        if problematic_cofactors is True: 
                            pymol.cmd.save(output_dir + 'cofactor_problematic_' + pred.split('/')[-1], 'current')  # save the model with cofactors
                        else:
                            pymol.cmd.save(output_dir + 'cofactor_' + pred.split('/')[-1], 'current')  # save the model with cofactors
                    
                        logfile.write('Cofactors sucessfully added to model.\n')
                    else: 
                        logfile.write('Cofactor was present but too far away from de novo model after alignment.\n')

                else:
                    logfile.write('No cofactors to add here.\n') 
                    pass
                    
                pymol.cmd.delete('all')

        except FileNotFoundError: 
            logfile.write('The respective model was not found. You can find its ID in missing_models.txt.\n')
            missing_models.append(identifier)

logfile.close()

    
out = open(output_dir + 'missing_models.txt', 'w')

for item in missing_models: 
    out.write(item + '\n')
out.close()

out = open(output_dir + 'HMs_failed_to_open.txt', 'w')

for item in failed_to_open: 
    out.write(item + '\n')
out.close()
