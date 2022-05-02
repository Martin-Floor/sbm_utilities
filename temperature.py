
import sbmOpenMM
import os
import numpy as np
#import bsc_calculations

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcol
import matplotlib.cm as cm

import mdtraj as md
from Bio import PDB

from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit

from collections import OrderedDict
import numpy as np

from sbmOpenMM.core.system import system
from sbmOpenMM.core.models import models


class smog2:
    """
    Class for automating the generation of SBM contact maps using the SMOG2 executable.

    The Executable can be obtained at:

    http://smog-server.org/smog2/

    The program needs to be installed and the executables avaiable through the linux PATH variable.

    Methods
    -------
    generateSBMTopology(pdb_file, kwarg**)
        Generates an all-heavy atom (AA) or alpha-carbon (CA) SBM topology and a protein contact file
        using SMOG2.
    """

    def correctPDBformat(pdb_file, output_file):

        parser = PDB.PDBParser()
        structure = parser.get_structure(pdb_file.replace('.pdb',''), pdb_file)
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        with open(pdb_file, 'r') as pdbf1:
            lines = pdbf1.readlines()
            if lines[-2].startswith('TER'):
                del lines[-2]
        pdbf = open(output_file, 'w')
        for line in lines:
            if line.startswith('ATOM') or line.startswith('TER'):
                pdbf.write(line)
        pdbf.write('END')
        pdbf.close()
        return output_file



    def generateSBMTopology(pdb_file, model_type='CA', contact_file=None, adjust_pdb=False, output_path=None, change_terminal_names=False):
        """
        Generates an AA or CA SBM topology and a protein contact file using SMOG2.

        Parameters
        ----------
        pdb_file : string
            Path to the input PDB file.
        model_type : string ('CA')
            Whether to generate an all atom 'AA' or a alpha-carbon 'CA' SBM model.
        contact_file : string
            Path to the output contact file.
        adjust_pdb : boolean (False)
            Wheter to adjust the format of the pdb file to be read by smog2. If given True this will
            overwrite the input file.
        output_path : str
            Path to save the adjusted file
        Returns
        -------
        None
        """

        # Executes the adjust pdb command and replaces the input pdb file
        if adjust_pdb:
            adjust_command = 'smog_adjustPDB -i '+pdb_file
            os.system(adjust_command)
            if output_path == None:
                output_path = pdb_file.split('/')[-1]
            else:
                output_path = output_path+'/'+pdb_file.split('/')[-1]
            shutil.copyfile('adjusted.pdb', output_path)
            os.remove('adjusted.pdb')
            pdb_file = output_path

        # Executes the SBM topology command for smog2
        if model_type == 'AA':
            command = 'smog2 -i '+pdb_file+' -AA -dname '+pdb_file.replace('.pdb','')
        elif model_type == 'CA':
            command = 'smog2 -i '+pdb_file+' -CA -dname '+pdb_file.replace('.pdb','')
        if contact_file != None:
            command.replace('smog2','smog2 -c '+contact_file)
        print(command)
        os.system(command)

        if change_terminal_names:
            with open(pdb_file) as pf:
                lines = pf.readlines()
            with open(pdb_file, 'w') as of:
                for l in lines:
                    if l.startswith('ATOM'):
                        if len(l.split()) < 9:
                            nname = l.split()[3][:3]+' '+l.split()[3][-1]
                            of.write(l.replace(l.split()[3], nname))
                        else:
                            of.write(l)
                    else:
                        of.write(l)


    def renameTerminalResidues(pdb):

        pdbf = open(pdb,"r")
        lines = pdbf.readlines()
        pdbf.close()
        pdbc = open(pdb[:-4]+"_fixed.pdb","wt")
        for line in lines:
            if line.startswith('ATOM') and line[20]=="T":
                line = line[:20] + " " + line[21:]
                pdbc.write(line)
            else:
                pdbc.write(line)
        pdbc.close()


class system_to_simulate:
    def generate_contact_files(PDB_FILE,OUTPUT):

        smog2.correctPDBformat(PDB_FILE,output_file=OUTPUT)
        smog2.generateSBMTopology(OUTPUT,adjust_pdb=True)
        smog2.renameTerminalResidues(OUTPUT)



    def factor_generator(OUTPUT,CONTACTS,NAME,Lid_START,Lid_END,Prot_START,
                     Prot_END,RANGE_START=0,
                     RANGE_END=1.1,RANGE_STEP=0.1):

        ff_file = {}
        for factor in np.arange(0,1.1,0.1):
            sbm_model = sbmOpenMM.models.getAllAtomModel(OUTPUT, CONTACTS,
                                                         masses_per_element=True, radii_per_atom_type=True,
                                                        ff_radii='oplsaa')

            lid_domain = set([*range(Lid_START,Lid_END+1)])
            not_lid_domain = set(range(Prot_START,Prot_END+1))
            not_lid_domain = not_lid_domain - not_lid_domain.intersection(lid_domain)

            topology = sbm_model.structure.getTopology()

            for c in sbm_model.contacts:
                if int(c[0].residue.id) in lid_domain and int(c[1].residue.id) in not_lid_domain:
                    sbm_model.contacts[c] = (sbm_model.contacts[c][0], sbm_model.contacts[c][1]*factor)
                if int(c[0].residue.id) in not_lid_domain and int(c[1].residue.id) in lid_domain:
                    sbm_model.contacts[c] = (sbm_model.contacts[c][0], sbm_model.contacts[c][1]*factor)
            sbm_model.dumpForceFieldData(NAME + '_'+str('%.2f' % factor)+'.ff')
            ff_file['%.2f' % factor] = NAME + '_'+str('%.2f' % factor)+'.ff'
