"""
Molecular Docking Python Script
A beginner-level guide to molecular docking for bachelor students

This script provides tools for:
- Downloading protein structures from the Protein Data Bank (PDB)
- Preparing proteins and ligands for docking
- Analyzing docking results
- Interpreting binding scores
"""

import os
import json
import urllib.request
import urllib.error
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from pathlib import Path


# ============================================================================
# Data Structures for Molecular Docking
# ============================================================================

@dataclass
class ProteinLigandSystem:
    """Represents a protein-ligand system for docking"""
    name: str
    protein_name: str
    ligand_name: str
    pdb_id: str
    description: str


@dataclass
class DockingResult:
    """Stores docking results"""
    mode: int
    binding_affinity: float  # kcal/mol
    rmsd_lb: float  # Lower bound RMSD
    rmsd_ub: float  # Upper bound RMSD
    
    def __str__(self) -> str:
        return (f"Mode {self.mode}: "
                f"Binding Affinity: {self.binding_affinity:.2f} kcal/mol, "
                f"RMSD: {self.rmsd_lb:.2f}-{self.rmsd_ub:.2f} Å")


# ============================================================================
# Protein-Ligand Systems Database
# ============================================================================

SYSTEMS = {
    "experiment_1": ProteinLigandSystem(
        name="Experiment 1",
        protein_name="Carbonic Anhydrase II",
        ligand_name="Acetazolamide",
        pdb_id="3HS4",
        description="Human enzyme involved in CO2 transport"
    ),
    "experiment_2": ProteinLigandSystem(
        name="Experiment 2",
        protein_name="HIV-1 Protease",
        ligand_name="Amprenavir",
        pdb_id="1HPV",
        description="Viral protease targeted by antiretroviral drugs"
    ),
    "experiment_3": ProteinLigandSystem(
        name="Experiment 3",
        protein_name="COX-2",
        ligand_name="Celecoxib",
        pdb_id="3LN1",
        description="Cyclooxygenase-2, target for anti-inflammatory drugs"
    ),
    "experiment_4": ProteinLigandSystem(
        name="Experiment 4",
        protein_name="Estrogen Receptor Alpha",
        ligand_name="4-Hydroxytamoxifen",
        pdb_id="3ERT",
        description="Nuclear receptor implicated in breast cancer"
    ),
}


# ============================================================================
# PDB Download and File Management
# ============================================================================

class PDBDownloader:
    """Handles downloading protein structures from RCSB PDB"""
    
    BASE_URL = "https://files.rcsb.org/download"
    
    @staticmethod
    def download_pdb(pdb_id: str, output_dir: str = "pdb_files") -> Optional[str]:
        """
        Download a protein structure from RCSB PDB
        
        Args:
            pdb_id: 4-character PDB identifier (e.g., '3HS4')
            output_dir: Directory to save the file
            
        Returns:
            Path to downloaded file, or None if download failed
        """
        # Create output directory if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Construct URL and filename
        pdb_id_lower = pdb_id.lower()
        url = f"{PDBDownloader.BASE_URL}/{pdb_id_lower}.pdb"
        filepath = os.path.join(output_dir, f"{pdb_id_lower}.pdb")
        
        try:
            print(f"Downloading {pdb_id} from RCSB PDB...")
            urllib.request.urlretrieve(url, filepath)
            print(f"✓ Successfully downloaded to {filepath}")
            return filepath
        except urllib.error.URLError as e:
            print(f"✗ Failed to download {pdb_id}: {e}")
            return None
    
    @staticmethod
    def parse_pdb(filepath: str) -> Dict:
        """
        Parse a PDB file and extract basic information
        
        Args:
            filepath: Path to PDB file
            
        Returns:
            Dictionary with PDB information
        """
        info = {
            "title": "",
            "atoms": [],
            "hetatms": [],
            "chains": set()
        }
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    # Extract title
                    if line.startswith("TITLE"):
                        info["title"] = line[10:].strip()
                    
                    # Extract atom information
                    elif line.startswith("ATOM"):
                        atom_data = {
                            "index": int(line[6:11].strip()),
                            "name": line[12:16].strip(),
                            "residue": line[17:20].strip(),
                            "chain": line[21].strip(),
                            "residue_num": int(line[22:26].strip()),
                            "x": float(line[30:38].strip()),
                            "y": float(line[38:46].strip()),
                            "z": float(line[46:54].strip()),
                        }
                        info["atoms"].append(atom_data)
                        info["chains"].add(atom_data["chain"])
                    
                    # Extract ligand/heteroatom information
                    elif line.startswith("HETATM"):
                        hetm_data = {
                            "index": int(line[6:11].strip()),
                            "name": line[12:16].strip(),
                            "residue": line[17:20].strip(),
                            "chain": line[21].strip(),
                        }
                        info["hetatms"].append(hetm_data)
            
            return info
        except FileNotFoundError:
            print(f"✗ File not found: {filepath}")
            return info
    
    @staticmethod
    def get_pdb_info_url(pdb_id: str) -> str:
        """Get the URL for viewing PDB structure information"""
        return f"https://www.rcsb.org/structure/{pdb_id}"


# ============================================================================
# Protein and Ligand Preparation
# ============================================================================

class ProteinPreparer:
    """Handles protein structure preparation for docking"""
    
    @staticmethod
    def remove_water_molecules(filepath: str, output_filepath: str) -> bool:
        """
        Remove water molecules (HOH) from PDB file
        
        Args:
            filepath: Input PDB file
            output_filepath: Output PDB file without water
            
        Returns:
            True if successful, False otherwise
        """
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Filter out water molecules and END statements
            filtered_lines = [
                line for line in lines
                if not line.startswith("HETATM") or "HOH" not in line
            ]
            
            # Ensure file ends properly
            if filtered_lines and not filtered_lines[-1].startswith("END"):
                filtered_lines.append("END\n")
            
            with open(output_filepath, 'w') as f:
                f.writelines(filtered_lines)
            
            print(f"✓ Removed water molecules: {output_filepath}")
            return True
        except Exception as e:
            print(f"✗ Error removing water molecules: {e}")
            return False
    
    @staticmethod
    def extract_chain(filepath: str, chain_id: str, output_filepath: str) -> bool:
        """
        Extract a specific chain from PDB file
        
        Args:
            filepath: Input PDB file
            chain_id: Chain identifier to extract
            output_filepath: Output PDB file
            
        Returns:
            True if successful, False otherwise
        """
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Extract specific chain
            extracted_lines = [
                line for line in lines
                if line.startswith("ATOM") and line[21] == chain_id
            ]
            
            extracted_lines.append("END\n")
            
            with open(output_filepath, 'w') as f:
                f.writelines(extracted_lines)
            
            print(f"✓ Extracted chain {chain_id}: {output_filepath}")
            return True
        except Exception as e:
            print(f"✗ Error extracting chain: {e}")
            return False
    
    @staticmethod
    def get_active_site_info(pdb_info: Dict, center_residue: Optional[int] = None) -> Dict:
        """
        Get information about the protein active site
        
        Args:
            pdb_info: Parsed PDB information
            center_residue: Residue number at center of active site (optional)
            
        Returns:
            Dictionary with active site information
        """
        if not pdb_info["atoms"]:
            return {"error": "No atoms found"}
        
        # Calculate centroid of all atoms
        x_coords = [a["x"] for a in pdb_info["atoms"]]
        y_coords = [a["y"] for a in pdb_info["atoms"]]
        z_coords = [a["z"] for a in pdb_info["atoms"]]
        
        centroid = {
            "x": sum(x_coords) / len(x_coords),
            "y": sum(y_coords) / len(y_coords),
            "z": sum(z_coords) / len(z_coords),
        }
        
        # Find extent for grid size estimation
        x_range = max(x_coords) - min(x_coords)
        y_range = max(y_coords) - min(y_coords)
        z_range = max(z_coords) - min(z_coords)
        
        return {
            "centroid": centroid,
            "grid_center": centroid,
            "grid_size": {
                "x": x_range,
                "y": y_range,
                "z": z_range,
            },
            "num_atoms": len(pdb_info["atoms"]),
            "chains": list(pdb_info["chains"]),
        }


class LigandPreparer:
    """Handles ligand structure preparation for docking"""
    
    @staticmethod
    def extract_ligand(pdb_filepath: str, ligand_name: str, output_filepath: str) -> bool:
        """
        Extract ligand from PDB file
        
        Args:
            pdb_filepath: Input PDB file
            ligand_name: 3-letter ligand code (e.g., 'AAZ' for acetazolamide)
            output_filepath: Output PDB file containing only ligand
            
        Returns:
            True if successful, False otherwise
        """
        try:
            with open(pdb_filepath, 'r') as f:
                lines = f.readlines()
            
            # Extract ligand atoms
            ligand_lines = [
                line for line in lines
                if line.startswith("HETATM") and ligand_name in line
            ]
            
            if not ligand_lines:
                print(f"✗ Ligand '{ligand_name}' not found in {pdb_filepath}")
                return False
            
            ligand_lines.append("END\n")
            
            with open(output_filepath, 'w') as f:
                f.writelines(ligand_lines)
            
            print(f"✓ Extracted ligand '{ligand_name}': {output_filepath}")
            return True
        except Exception as e:
            print(f"✗ Error extracting ligand: {e}")
            return False
    
    @staticmethod
    def get_ligand_properties(pdb_filepath: str, ligand_name: str) -> Dict:
        """
        Extract and calculate ligand properties
        
        Args:
            pdb_filepath: PDB file containing the ligand
            ligand_name: 3-letter ligand code
            
        Returns:
            Dictionary with ligand properties
        """
        properties = {
            "name": ligand_name,
            "num_atoms": 0,
            "atom_types": set(),
            "heavy_atoms": 0,
            "centroid": {"x": 0, "y": 0, "z": 0},
        }
        
        try:
            with open(pdb_filepath, 'r') as f:
                lines = f.readlines()
            
            ligand_atoms = [
                line for line in lines
                if line.startswith("HETATM") and ligand_name in line
            ]
            
            if not ligand_atoms:
                return properties
            
            # Extract coordinates and atom types
            coords = {"x": [], "y": [], "z": []}
            for line in ligand_atoms:
                atom_type = line[76:78].strip()
                if atom_type and atom_type != "H":  # Skip hydrogens for heavy atom count
                    properties["heavy_atoms"] += 1
                
                properties["atom_types"].add(atom_type)
                coords["x"].append(float(line[30:38].strip()))
                coords["y"].append(float(line[38:46].strip()))
                coords["z"].append(float(line[46:54].strip()))
            
            properties["num_atoms"] = len(ligand_atoms)
            properties["atom_types"] = list(properties["atom_types"])
            
            # Calculate centroid
            if coords["x"]:
                properties["centroid"] = {
                    "x": sum(coords["x"]) / len(coords["x"]),
                    "y": sum(coords["y"]) / len(coords["y"]),
                    "z": sum(coords["z"]) / len(coords["z"]),
                }
            
            return properties
        except Exception as e:
            print(f"✗ Error extracting ligand properties: {e}")
            return properties


# ============================================================================
# Docking Analysis and Scoring
# ============================================================================

class DockingAnalyzer:
    """Analyzes and interprets docking results"""
    
    # Binding affinity reference values (kcal/mol)
    AFFINITY_SCALE = {
        "very_strong": (-10.0, -100.0),    # < -10
        "strong": (-8.0, -10.0),           # -10 to -8
        "moderate": (-6.0, -8.0),          # -8 to -6
        "weak": (-4.0, -6.0),              # -6 to -4
        "very_weak": (0.0, -4.0),          # > -4
    }
    
    @staticmethod
    def parse_docking_results(results_file: str) -> List[DockingResult]:
        """
        Parse docking results file
        
        Args:
            results_file: Path to docking results file
            
        Returns:
            List of DockingResult objects
        """
        results = []
        try:
            with open(results_file, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
                if line.startswith("RESULT"):
                    parts = line.split()
                    if len(parts) >= 5:
                        result = DockingResult(
                            mode=int(parts[1]),
                            binding_affinity=float(parts[2]),
                            rmsd_lb=float(parts[3]),
                            rmsd_ub=float(parts[4]),
                        )
                        results.append(result)
            
            return results
        except FileNotFoundError:
            print(f"✗ Results file not found: {results_file}")
            return results
    
    @staticmethod
    def classify_binding_affinity(affinity: float) -> str:
        """
        Classify binding affinity strength
        
        Args:
            affinity: Binding affinity in kcal/mol (typically negative)
            
        Returns:
            Classification string
        """
        for classification, (lower, upper) in DockingAnalyzer.AFFINITY_SCALE.items():
            if lower <= affinity <= upper:
                return classification
        return "unknown"
    
    @staticmethod
    def analyze_results(results: List[DockingResult]) -> Dict:
        """
        Perform comprehensive analysis of docking results
        
        Args:
            results: List of DockingResult objects
            
        Returns:
            Dictionary with analysis results
        """
        if not results:
            return {"error": "No results to analyze"}
        
        affinities = [r.binding_affinity for r in results]
        
        analysis = {
            "num_poses": len(results),
            "best_affinity": min(affinities),
            "worst_affinity": max(affinities),
            "mean_affinity": sum(affinities) / len(affinities),
            "best_pose": results[0],
            "binding_classification": DockingAnalyzer.classify_binding_affinity(
                min(affinities)
            ),
        }
        
        return analysis
    
    @staticmethod
    def print_results_summary(results: List[DockingResult]) -> None:
        """
        Print a formatted summary of docking results
        
        Args:
            results: List of DockingResult objects
        """
        if not results:
            print("No results to display")
            return
        
        analysis = DockingAnalyzer.analyze_results(results)
        
        print("\n" + "="*70)
        print("DOCKING RESULTS SUMMARY")
        print("="*70)
        print(f"Number of poses: {analysis['num_poses']}")
        print(f"Best binding affinity: {analysis['best_affinity']:.2f} kcal/mol")
        print(f"Worst binding affinity: {analysis['worst_affinity']:.2f} kcal/mol")
        print(f"Mean binding affinity: {analysis['mean_affinity']:.2f} kcal/mol")
        print(f"Binding classification: {analysis['binding_classification'].upper()}")
        print("\nDetailed Results:")
        print("-"*70)
        for result in results[:5]:  # Show top 5
            print(f"  {result}")
        if len(results) > 5:
            print(f"  ... and {len(results) - 5} more poses")
        print("="*70 + "\n")


# ============================================================================
# Docking Workflow Coordinator
# ============================================================================

class DockingWorkflow:
    """Coordinates the complete docking workflow"""
    
    def __init__(self, working_dir: str = "docking_workspace"):
        """
        Initialize docking workflow
        
        Args:
            working_dir: Base directory for all files
        """
        self.working_dir = working_dir
        self.pdb_dir = os.path.join(working_dir, "pdb_files")
        self.ligand_dir = os.path.join(working_dir, "ligands")
        self.protein_dir = os.path.join(working_dir, "proteins")
        self.results_dir = os.path.join(working_dir, "results")
        
        # Create directories
        for directory in [self.pdb_dir, self.ligand_dir, self.protein_dir, self.results_dir]:
            Path(directory).mkdir(parents=True, exist_ok=True)
    
    def setup_experiment(self, system_key: str) -> bool:
        """
        Set up a complete docking experiment
        
        Args:
            system_key: Key to SYSTEMS dictionary
            
        Returns:
            True if setup successful, False otherwise
        """
        if system_key not in SYSTEMS:
            print(f"✗ Unknown system: {system_key}")
            return False
        
        system = SYSTEMS[system_key]
        print(f"\n{'='*70}")
        print(f"Setting up: {system.name}")
        print(f"{'='*70}")
        print(f"Protein: {system.protein_name}")
        print(f"Ligand: {system.ligand_name}")
        print(f"PDB ID: {system.pdb_id}")
        print(f"{'-'*70}\n")
        
        # Download PDB structure
        pdb_file = PDBDownloader.download_pdb(system.pdb_id, self.pdb_dir)
        if not pdb_file:
            return False
        
        # Parse PDB file
        pdb_info = PDBDownloader.parse_pdb(pdb_file)
        print(f"✓ Parsed PDB: {len(pdb_info['atoms'])} atoms, Chains: {pdb_info['chains']}")
        
        # Prepare protein
        protein_clean = os.path.join(
            self.protein_dir,
            f"{system.pdb_id}_protein_clean.pdb"
        )
        ProteinPreparer.remove_water_molecules(pdb_file, protein_clean)
        
        # Get active site information
        active_site = ProteinPreparer.get_active_site_info(pdb_info)
        print(f"\nActive site centroid: ({active_site['centroid']['x']:.2f}, "
              f"{active_site['centroid']['y']:.2f}, "
              f"{active_site['centroid']['z']:.2f})")
        print(f"Grid size: {active_site['grid_size']['x']:.1f} × "
              f"{active_site['grid_size']['y']:.1f} × "
              f"{active_site['grid_size']['z']:.1f} Å")
        
        # Save configuration
        config = {
            "system": system_key,
            "pdb_id": system.pdb_id,
            "protein_file": protein_clean,
            "active_site": active_site,
            "ligand_code": "AAZ",  # Default, adjust as needed
        }
        
        config_file = os.path.join(self.results_dir, f"{system_key}_config.json")
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"\n✓ Experiment setup complete!")
        print(f"  Configuration saved to: {config_file}")
        
        return True
    
    def print_workflow_instructions(self) -> None:
        """Print instructions for manual docking steps"""
        print("\n" + "="*70)
        print("NEXT STEPS - MANUAL DOCKING")
        print("="*70)
        print("""
1. EXTRACT LIGAND (if not already in complex):
   - Use PyMOL or PyRx to extract the ligand from the protein structure
   
2. PREPARE LIGAND:
   - Add hydrogens
   - Add Gasteiger charges
   - Adjust protonation state if needed
   
3. PREPARE PROTEIN:
   - Merge non-polar hydrogens
   - Add Gasteiger charges
   - Assign AD4 atom types
   
4. DEFINE DOCKING GRID:
   - Center: Active site centroid (shown above)
   - Size: Based on grid dimensions (shown above)
   - Spacing: 0.375 Å (default)
   
5. RUN DOCKING:
   - Use AutoDock Vina with prepared files
   - Exhaustiveness: 8 (default)
   
6. ANALYZE RESULTS:
   - Load docking results file
   - Visualize top poses
   - Analyze interactions
        """)
        print("="*70 + "\n")


# ============================================================================
# Utility Functions
# ============================================================================

def list_available_systems() -> None:
    """Display all available protein-ligand systems"""
    print("\n" + "="*70)
    print("AVAILABLE PROTEIN-LIGAND SYSTEMS")
    print("="*70)
    for key, system in SYSTEMS.items():
        print(f"\n{key.upper()}")
        print(f"  Name: {system.name}")
        print(f"  Protein: {system.protein_name}")
        print(f"  Ligand: {system.ligand_name}")
        print(f"  PDB ID: {system.pdb_id}")
        print(f"  Description: {system.description}")
    print("="*70 + "\n")


def get_pdb_info(pdb_id: str) -> None:
    """Get information about a PDB entry"""
    print(f"\nPDB Information: {pdb_id}")
    print(f"URL: {PDBDownloader.get_pdb_info_url(pdb_id)}\n")


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    print("\n" + "█"*70)
    print("  MOLECULAR DOCKING - Python Implementation")
    print("  Beginner-Level Practical for Bachelor Students")
    print("█"*70 + "\n")
    
    # Display available systems
    list_available_systems()
    
    # Example: Set up first experiment
    workflow = DockingWorkflow(working_dir="docking_workspace")
    
    # Setup Experiment 1: Carbonic Anhydrase II
    if workflow.setup_experiment("experiment_1"):
        workflow.print_workflow_instructions()
    
    print("\nTo analyze docking results later, use:")
    print("  results = DockingAnalyzer.parse_docking_results('results_file.txt')")
    print("  DockingAnalyzer.print_results_summary(results)")
    print()
