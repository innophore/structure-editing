"""This script removes low confidence regions (pLDDT < 70) from AI-predicted models."""

import shutil
import os
from Bio.PDB import *
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
import traceback

parser = PDBParser()

directory = './'
bfactor_cutoff = 70
average_over = 3
tempfolder= directory + 'temp' + str(bfactor_cutoff) +'/'
os.mkdir(tempfolder)

outdir=directory + 'cut' + str(bfactor_cutoff) +'/'
os.mkdir(outdir)

if os.path.exists(directory + 'logfile.txt'):
    os.remove(directory + 'logfile.txt')


# Loop through all files in the directory
for filename in os.listdir(directory):
    # Check if the file has a .pdb extension
    if filename.endswith('.pdb'):
        print(filename)

        input_file_path=directory + filename
        try:
            lines_without_hetatm = []
            with open(input_file_path, "r") as f:
                for line in f:
                    if not (line.startswith("HETATM") or line.startswith("END")):
                        lines_without_hetatm.append(line)

            lines_with_hetatm = []
            with open(input_file_path, "r") as f:
                for line in f:
                    if line.startswith("HETATM"):
                        lines_with_hetatm.append(line)

            with open(tempfolder + filename, "w") as f:
                f.writelines(lines_without_hetatm)

             # read in the PDB file
            structure = parser.get_structure('PDB', tempfolder + filename)


            #residues below threshold
            resnbelow=[]

            #set cutoff and number of aa to average over
            delete_residues=[]

            for i in range(1, len(structure[0]['A'])):
                if 0 > i - average_over:
                    bfactorlist = [structure[0]['A'][x]['CA'].get_bfactor() for x in range(1, i + average_over)]
                    extendedlist = [res for res in range(1, i + average_over)]

                    if np.average(bfactorlist) < bfactor_cutoff:
                        delete_residues.extend(extendedlist)

                elif 0 < i - average_over and i + average_over < len(structure[0]['A']):

                    bfactorlist = [structure[0]['A'][x]['CA'].get_bfactor() for x in range(i - average_over, i + average_over)]
                    extendedlist = [res for res in range(i - average_over, i + average_over)]

                    if np.average(bfactorlist) < bfactor_cutoff:
                        delete_residues.extend(extendedlist)#[0:5])

                elif i + average_over > len(structure[0]['A']):
                    bfactorlist = [structure[0]['A'][y]['CA'].get_bfactor() for y in range(i - average_over, len(structure[0]['A']) + 1)]
                    extendedlist = [res for res in range(i - average_over, len(structure[0]['A']) + 1)]


                    if np.average(bfactorlist) < bfactor_cutoff:
                        delete_residues.extend(extendedlist)

            delete_residues_u = []
            for i in delete_residues:
                if i not in delete_residues_u:
                    delete_residues_u.append(i)


            ns = NeighborSearch(list(structure.get_atoms()))
            #delete the residues from the chain
            seqlength_old = len(structure[0]['A'])
            for model in structure:
               for chain in model:
                   for id in delete_residues_u:
                       chain.detach_child((' ', id, ' '))


            previous_residue_number = None
            current_residue_list = []
            residue_lists = []

            for model in structure:
                for chain in model:
                    for residue in chain:
                        residue_number = residue.get_id()[1]
                        if previous_residue_number is not None and residue_number - previous_residue_number > 1:
                            # There is a gap, so start a new list
                            residue_lists.append(current_residue_list)
                            current_residue_list = []
                        current_residue_list.append(residue.get_id()[1])
                        previous_residue_number = residue_number

            # Append the last list to the residue lists
            residue_lists.append(current_residue_list)

            longest_list = max(residue_lists, key=len)
            shorter_lists = [lst for lst in residue_lists if lst != longest_list]


            cutoff = 10
            cutoff_short= 5

            nearby_residues_list = []
            for lst in shorter_lists:
                nearby_residues = []
                for residue_id in lst:
                    residue = structure[0]["A"][residue_id]
                    atom = residue['CA']
                    nearby_residues += [r.get_id()[1] for r in ns.search(atom.get_coord(), cutoff, 'R')]
                nearby_residues_list.append(nearby_residues)

            nearby_short=[]
            for lst in shorter_lists:
                nearby_residues = []
                for residue_id in lst:
                    residue = structure[0]["A"][residue_id]
                    atom = residue['CA']
                    nearby_residues += [r.get_id()[1] for r in ns.search(atom.get_coord(), cutoff, 'R')]
                nearby_short.append(nearby_residues)

            residues_to_keep = []
            for i, lst in enumerate(nearby_short):
                exclude_lst = shorter_lists[i]  # exclude the current lst
                for item in lst:
                    for other_lst in shorter_lists:
                        if item in other_lst and other_lst != exclude_lst:
                            for i in other_lst:
                                residues_to_keep.append(i)

            residues_to_keep_u = []
            for i in residues_to_keep:
                if i not in residues_to_keep_u:
                    residues_to_keep_u.append(i)

            residues_todelete = []

            for i, lst in enumerate(nearby_residues_list):
                if not any(item in longest_list for item in lst):
                    residues_todelete += shorter_lists[i]

            residues_todelete_u = []
            for i in residues_todelete:
                if i not in residues_todelete_u:
                    residues_todelete_u.append(i)
            residues_todelete_u = [r for r in residues_todelete_u if r not in residues_to_keep]

            for model in structure:
               for chain in model:
                   for id in residues_todelete_u:
                       chain.detach_child((' ', id, ' '))

            seq_lengthnew = len(structure[0]['A'])
            with open(directory + 'logfile.txt' , "a" ) as file:
                 file.write(f"{filename.strip('.pdb')},  {seqlength_old}, {seq_lengthnew}, {seq_lengthnew/seqlength_old}\n")

            # #save the modified PDB file
            io = PDBIO()
            io.set_structure(structure)
            io.save(outdir + filename.strip('.pdb') + '_pLDDT_{}.pdb'.format(bfactor_cutoff))

            os.system(
                f"sed -i '/{'END'}/d' {outdir + filename.strip('.pdb')  + '_pLDDT_{}.pdb'.format(bfactor_cutoff)}")

            f = open(outdir + filename.strip('.pdb') + '_pLDDT_{}.pdb'.format(bfactor_cutoff),'a')
            f.writelines(lines_with_hetatm)
            f.writelines('END')
            f.close()

        except Exception:
            print(filename + 'not working')
            print(traceback.print_exc())

shutil.rmtree(tempfolder)
