import os
import shutil
import subprocess

class SequenceProcessor:
    def __init__(self, sequence, ring_len, loop_len, fold_direction, wk_dir=None, main_script_dir=None):
        self.module_dir = os.path.dirname(os.path.realpath(__file__))
        self.ring_len = ring_len
        self.loop_len = loop_len
        self.sequence = sequence
        self.tail_len = len(self.sequence) - self.ring_len - self.loop_len
        self.fold_direction = fold_direction
        self.wk_dir = os.path.abspath(wk_dir) if wk_dir else os.path.abspath(main_script_dir)
        if not os.path.exists(self.wk_dir):
            os.makedirs(self.wk_dir, exist_ok=True)  
        self.raw_output_dir = os.path.join(self.wk_dir, 'raw_output')
        self.final_output_dir = os.path.join(self.wk_dir, 'final_output')
        self.create_directories()

    def create_directories(self):
        os.makedirs(self.raw_output_dir, exist_ok=True)
        os.makedirs(self.final_output_dir, exist_ok=True)

    def check_iso_peptide(self):
        # Check if the sequence length is sufficient
        if len(self.sequence) < self.ring_len + self.loop_len + 2:
            print("Sequence is too short.")
            return
        # Get the residue at the ring_len position
        self.iso_res = self.sequence[self.ring_len - 1]  # Indexing starts from 0, so subtract 1
        # Full name mapping for residues
        aa_full_name = {'E': 'Glu', 'D': 'Asp', 'N': 'Asn'}
        if self.iso_res in aa_full_name:
            print(f"Isopeptide is {aa_full_name[self.iso_res]}")
        else:
            print("No isopeptide found.")
    
    def copy_scaffold_file(self):
        scaffold_file = f'r{self.ring_len}_l{self.loop_len}.pdb'
        # Choose directory based on fold_direction
        if self.fold_direction == 'right':
            directory = 'scaf_gly_right'
        elif self.fold_direction == 'left': 
            directory = 'scaf_gly_left'
        
        source_path = os.path.join(self.module_dir, directory, scaffold_file)
        destination_path = os.path.join(self.raw_output_dir, scaffold_file)
        shutil.copyfile(source_path, destination_path)

    def generate_pymol_script(self):  
        tail_start = self.ring_len + self.loop_len + 1
        # Generate PyMOL script
        pymol_script = f"""
    # Step 1: Load the scaffold
    load {self.wk_dir}/raw_output/r{self.ring_len}_l{self.loop_len}.pdb
    
    # Step 2: Create a new peptide segment
    fab {'A' * self.tail_len}, my_tail, resi={tail_start}, chain=A, ss=4
    
    # Step 3: Align the new peptide segment to the scaffold
    select mobile, resi {tail_start} & name C+N+O+CA & my_tail
    select target, resi {tail_start} & name C+N+O+CA & r{self.ring_len}_l{self.loop_len}
    pair_fit mobile, target
    
    # Step 4: Remove the tail residues from the scaffold
    remove resi {tail_start} & r{self.ring_len}_l{self.loop_len}
    remove resi {tail_start + 1} & r{self.ring_len}_l{self.loop_len}
    
    # Step 5: Merge scaffold and tail and create bond
    create merged_object, r{self.ring_len}_l{self.loop_len} or my_tail
    bond resi {tail_start - 1} and name C and merged_object, resi {tail_start} and name N and merged_object
    """
        # Amino acid conversion dictionary
        aa_dict = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 
            'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 
            'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 
            'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }
        # Adding mutation instructions
        for i, amino_acid in enumerate(self.sequence):
            pymol_script += f"""
    # Step 6: Mutation for residue {i+1}
    select mut_res, resi {i+1} & merged_object
    cmd.wizard('mutagenesis')
    cmd.do('refresh_wizard')
    cmd.get_wizard().set_mode('{aa_dict[amino_acid]}')  # Use three-letter code
    cmd.get_wizard().do_select('mut_res')
    cmd.get_wizard().apply()
    cmd.set_wizard()
    """
        # Save merged_object
        pymol_script += f"\ncmd.remove('hydrogens')\n"
        pymol_script += f"\nsave {self.wk_dir}/raw_output/r{self.ring_len}_l{self.loop_len}_t{self.tail_len}.pdb, merged_object\n"
        # Save the script
        script_path = os.path.join(self.raw_output_dir, 'pymol_script.pml')
        with open(script_path, 'w') as file:
            file.write(pymol_script)
        print(f"PyMOL script created at {script_path}")
    
    def mut_scaf(self):
        try:
            os.chdir(self.raw_output_dir)
            # Run the PyMOL script
            with open('pymol_script.out', 'w') as output_file:
                subprocess.run(['pymol', '-c', './pymol_script.pml'], stdout=output_file, stderr=subprocess.STDOUT)
            # Check if there are errors
            with open('pymol_script.out', 'r') as output_file:
                output_content = output_file.read()
                if "Error" in output_content or "error" in output_content:
                    print("Error detected in PyMOL script execution. See output in:", os.path.join(self.raw_output_dir, 'pymol_script.out'))
                else:
                    print("Scaffold mutated successfully.")
        finally:
            os.chdir(self.wk_dir)
        print("PyMOL script executed.")

    def pre_tleap1(self):
        # Source PDB file path
        source_pdb_path = os.path.join(self.raw_output_dir, f'r{self.ring_len}_l{self.loop_len}_t{self.tail_len}.pdb')
        # Output PDB file path
        output_pdb_path = os.path.join(self.raw_output_dir, f'r{self.ring_len}_l{self.loop_len}_t{self.tail_len}_rename.pdb')
        # Read the original PDB file
        with open(source_pdb_path, 'r') as file:
            lines = file.readlines()
        new_lines = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Remove hydrogen atoms
                if line[76:78].strip() == 'H' or line[12:14].strip().startswith('H'):
                    continue
                # Modify specific residue names and remove specific atoms if necessary
                if int(line[22:26].strip()) == self.ring_len:
                    if self.iso_res == 'D':
                        line = line[:17] + 'ASX' + line[20:]
                        if line[12:16].strip() == 'OD2':  
                            continue
                    elif self.iso_res == 'E':
                        line = line[:17] + 'GLX' + line[20:]
                        if line[12:16].strip() == 'OE2':  
                            continue
                    elif self.iso_res == 'N':
                        line = line[:17] + 'ANX' + line[20:]
                        if line[12:16].strip() == 'ND2':  
                            continue
            new_lines.append(line)

        with open(output_pdb_path, 'w') as file:
            file.writelines(new_lines)
        print(f"Modified PDB file saved to {output_pdb_path}")


    def lasso_leap(self, step):
        if step == 'min1':
            pdb_file = os.path.join(self.wk_dir, 'raw_output', f'r{self.ring_len}_l{self.loop_len}_t{self.tail_len}_rename.pdb')
            tleap_dir = os.path.join(self.raw_output_dir, 'tleap1')
        elif step == 'min2':
            pdb_file = os.path.join(self.wk_dir, 'raw_output', 'min1', 'min1_clean.pdb')
            tleap_dir = os.path.join(self.raw_output_dir, 'tleap2')
        else:
            raise ValueError("Invalid step. Expected 'min1' or 'min2'.")
        # Create tleap directory
        os.makedirs(tleap_dir, exist_ok=True)
        leap_script = f"""
source {self.module_dir}/iso_ff/leaprc.protein.ff14SB_lasso
source leaprc.water.tip3p
"""
        # Set script content based on iso_res
        if self.iso_res == 'D':
            leap_script += f"""
ASX = loadmol2 {self.module_dir}/iso_ff/ASX.mol2
loadamberparams {self.module_dir}/iso_ff/ASX.frcmod
mol = loadpdb {pdb_file}
bond mol.{self.ring_len}.CG mol.1.N
"""
        elif self.iso_res == 'E':
            leap_script += f"""
GLX = loadmol2 {self.module_dir}/iso_ff/GLX.mol2
loadamberparams {self.module_dir}/iso_ff/GLX.frcmod
mol = loadpdb {pdb_file}
bond mol.{self.ring_len}.CD mol.1.N
"""
        elif self.iso_res == 'N':
            leap_script += f"""
ANX = loadmol2 {self.module_dir}/iso_ff/ANX.mol2
loadamberparams {self.module_dir}/iso_ff/ANX.frcmod
mol = loadpdb {pdb_file}
bond mol.{self.ring_len}.CG mol.1.N
"""
        leap_script += f"""
bond mol.{self.ring_len - 1}.C mol.{self.ring_len}.N
bond mol.{self.ring_len}.C mol.{self.ring_len + 1}.N
solvateoct mol TIP3PBOX 10.0
addions mol Na+ 0
addions mol Cl- 0
charge mol
savepdb mol pro_leap.pdb
saveamberparm mol {tleap_dir}/pro_solv.prmtop {tleap_dir}/pro_solv.inpcrd
quit
"""
        # Set script content based on iso_res
        with open(os.path.join(tleap_dir, 'tleap.in'), 'w') as file:
            file.write(leap_script)
        print(f"tleap input file created at {os.path.join(tleap_dir, 'tleap.in')}")
        # Enter the tleap folder
        os.chdir(tleap_dir)

        with open('tleap.out', 'w') as output_file:
            subprocess.run(['tleap', '-s', '-f', 'tleap.in'], stdout=output_file)
        # Read and print the last line of tleap.out
        with open('tleap.out', 'r') as output_file:
            last_line = output_file.readlines()[-1]
            print("tleap result:", last_line)

        os.chdir(self.wk_dir)

    def md_minimization1(self):
        min1_dir = os.path.join(self.raw_output_dir, 'min1')
        os.makedirs(min1_dir, exist_ok=True)
        shutil.copy(os.path.join(self.module_dir, 'MDin', 'min1.in'), min1_dir)
        # Copy the output files from tleap1
        tleap1_dir = os.path.join(self.raw_output_dir, 'tleap1')
        shutil.copy(os.path.join(tleap1_dir, 'pro_solv.prmtop'), min1_dir)
        shutil.copy(os.path.join(tleap1_dir, 'pro_solv.inpcrd'), min1_dir)
        os.chdir(min1_dir)
        # Run pmemd.cuda
        #subprocess.run(['pmemd.cuda', '-O', '-p', 'pro_solv.prmtop', '-c', 'pro_solv.inpcrd', '-i', 'min1.in', '-o', 'min1.out', '-r', 'min1.rst', '-ref', 'pro_solv.inpcrd', '-AllowSmallBox'])
        subprocess.run(['mpirun','-np','16','pmemd.MPI','-O','-p','pro_solv.prmtop','-c','pro_solv.inpcrd','-i','min1.in','-o','min1.out','-r','min1.rst','-ref','pro_solv.inpcrd','-AllowSmallBox'])
        
        # Run cpptraj
        subprocess.run(['cpptraj', '-p', 'pro_solv.prmtop', '-y', 'min1.rst', '-x', 'min1.pdb'])
        # Process min1.pdb, removing lines containing WAT, Na, TER
        with open('min1.pdb', 'r') as file:
            lines = file.readlines()
        with open('min1_clean.pdb', 'w') as file:
            for line in lines:
                if not any(x in line for x in ['WAT', 'Na+', 'Cl-', 'TER']):
                    file.write(line)
        # Copy to the final_output folder
        final_output_path = os.path.join(self.final_output_dir, 'min1_clean.pdb')
        shutil.copy('min1_clean.pdb', final_output_path)

        os.chdir(self.wk_dir)
        print("MD minimization and cleaning completed. Cleaned file copied to final_output.")

    def md_prod(self):
        mdprod_dir = os.path.join(self.raw_output_dir, 'MDprod')
        os.makedirs(mdprod_dir, exist_ok=True)
        #for file in ['min.in', 'heat.in', 'equi.in', 'prod.in']:
        for file in ['min1.in', 'min2.in', 'min3.in', 'heat.in', 'nvt.in', 'npt.in', 'prod.in','cluster.in','sub_template.sh']:
            shutil.copy(os.path.join(self.module_dir, 'MDin', file), mdprod_dir)
        tleap2_dir = os.path.join(self.raw_output_dir, 'tleap2')
        shutil.copy(os.path.join(tleap2_dir, 'pro_solv.prmtop'), mdprod_dir)
        shutil.copy(os.path.join(tleap2_dir, 'pro_solv.inpcrd'), mdprod_dir)
        os.chdir(mdprod_dir)
        subprocess.run(["cpptraj", "-p", "pro_solv.prmtop", "-y", "pro_solv.inpcrd", "-x", "pro_solv.pdb"])
        # MD simulation
        #subprocess.run(['pmemd.cuda', '-O', '-p', 'pro_solv.prmtop', '-c', 'pro_solv.inpcrd',
        #                '-i', 'min.in', '-o', 'min.out', '-r', 'min.rst', '-x', 'min.mdcrd', '-ref', 'pro_solv.inpcrd'])
        #subprocess.run(['pmemd.cuda', '-O', '-p', 'pro_solv.prmtop', '-c', 'min.rst',
        #                '-i', 'heat.in', '-o', 'heat.out', '-r', 'heat.rst', '-x', 'heat.mdcrd', '-ref', 'min.rst'])
        #subprocess.run(['pmemd.cuda', '-O', '-p', 'pro_solv.prmtop', '-c', 'heat.rst',
        #                '-i', 'equi.in', '-o', 'equi.out', '-r', 'equi.rst', '-x', 'equi.mdcrd', '-ref', 'heat.rst'])
        #subprocess.run(['pmemd.cuda', '-O', '-p', 'pro_solv.prmtop', '-c', 'equi.rst',
        #                '-i', 'prod.in', '-o', 'prod.out', '-r', 'prod.rst', '-x', 'prod.nc', '-ref', 'equi.rst', '-AllowSmallBox'])
        
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'min1.in', '-p', 'pro_solv.prmtop', '-c', 'pro_solv.inpcrd', '-ref', 'pro_solv.inpcrd', '-o', 'min1.out', '-r', 'min1.rst', '-AllowSmallBox'])
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'min2.in', '-p', 'pro_solv.prmtop', '-c', 'min1.rst', '-ref', 'min1.rst', '-o', 'min2.out', '-r', 'min2.rst', '-AllowSmallBox'])
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'min3.in', '-p', 'pro_solv.prmtop', '-c', 'min2.rst', '-ref', 'min2.rst', '-o', 'min3.out', '-r', 'min3.rst', '-AllowSmallBox'])
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'heat.in', '-p', 'pro_solv.prmtop', '-c', 'min3.rst', '-ref', 'min3.rst', '-o', 'heat.out', '-r', 'heat.rst', '-AllowSmallBox'])
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'nvt.in', '-p', 'pro_solv.prmtop', '-c', 'heat.rst', '-ref', 'heat.rst', '-o', 'nvt.out', '-r', 'nvt.rst', '-AllowSmallBox'])
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'npt.in', '-p', 'pro_solv.prmtop', '-c', 'nvt.rst', '-ref', 'nvt.rst', '-o', 'npt.out', '-r', 'npt.rst', '-AllowSmallBox'])
        #subprocess.run(['pmemd.cuda', '-O', '-i', 'prod.in', '-p', 'pro_solv.prmtop', '-c', 'npt.rst', '-ref', 'npt.rst', '-o', 'prod.out', '-r', 'prod.rst', '-x', 'prod.nc', '-AllowSmallBox'])
        os.chdir(self.wk_dir)
        print("MD file copied.")

    def cluster(self):
        mdprod_dir = os.path.join(self.raw_output_dir, 'MDprod')
        # Ensure mdprod_dir exists
        if not os.path.exists(mdprod_dir):
            print(f"Directory {mdprod_dir} does not exist.")
            return
        # Create cluster.in file
        cluster_in_path = os.path.join(mdprod_dir, 'cluster.in')
        with open(cluster_in_path, 'w') as cluster_file:
            cluster_file.write(f"""
parm {mdprod_dir}/pro_solv.prmtop
trajin {mdprod_dir}/prod.nc
strip :WAT,Na+,Cl-
autoimage
cluster C1 @CA clusters 5 out {mdprod_dir}/cnumvtime.out summary {mdprod_dir}/summary.out summaryhalf {mdprod_dir}/summaryhalf.out repout {mdprod_dir}/cluster_result repfmt pdb
go
""")
        print(f"Cluster input file created at {cluster_in_path}")
        os.chdir(mdprod_dir)
        subprocess.run(['cpptraj', '-i', 'cluster.in'])
        os.chdir(self.wk_dir)
        # Copy and rename cluster_result.c0.pdb to the final_output folder
        cluster_output_path = os.path.join(mdprod_dir, 'cluster_result.c0.pdb')
        final_output_path = os.path.join(self.final_output_dir, 'relax_cluster_clean.pdb')
        if os.path.exists(cluster_output_path):
            shutil.copy(cluster_output_path, final_output_path)
            print(f"Cluster file copied to {final_output_path}")
        else:
            print("Cluster file not found.")
        print("Clustering completed.")

    def process_sequence(self):
        self.check_iso_peptide()
        self.copy_scaffold_file()
        self.generate_pymol_script()
        self.mut_scaf()
        self.pre_tleap1()
        self.lasso_leap(step='min1')
        self.md_minimization1()
        self.lasso_leap(step='min2')
        self.md_prod()
        #self.cluster()