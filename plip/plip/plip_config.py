from pydantic import BaseModel, Field
from typing import Optional, List, Dict

class InferenceRequest(BaseModel):
    """Request model for inference endpoint"""
    pdb_id: Optional[str] = None
    file_content: Optional[str] = None
    output_format: List[str] = ["xml", "txt"]

class PLIPConfig(BaseModel):
    """Configuration model for PLIP analysis"""
    output_format: List[str] = ["xml"]  # xml, txt, pymol
    model: int = Field(default=1, gt=0)  # Model number for multi-model structures
    verbose: bool = False
    peptides: List[str] = []
    nohydro: bool = False  # Do not add hydrogen bonds
    maxthreads: int = Field(default=1, gt=0)  # Maximum number of threads
    breakcomposite: bool = False  # Break up composite ligands
    altlocation: bool = False  # Consider alternate locations
    nofix: bool = False  # Turn off fixing of PDB files
    keepmod: bool = False  # Keep modified residues as ligands
    dnareceptor: bool = False  # Consider DNA as receptor
    chains: Optional[List[List[str]]] = None  # Chain definitions for interactions

    # Threshold settings
    distance_thresholds: Dict[str, float] = {
        "bs_dist": 7.5,  # Max distance for binding site residues
        "hydroph_dist_max": 4.0,  # Max distance for hydrophobic contacts
        "hbond_dist_max": 4.1,  # Max hydrogen bond distance
        "pistack_dist_max": 5.5,  # Max distance for pi-stacking
        "pication_dist_max": 6.0,  # Max distance for pi-cation interactions
        "saltbridge_dist_max": 5.5,  # Max distance for salt bridges
        "halogen_dist_max": 4.0,  # Max distance for halogen bonds
        "metal_dist_max": 3.0  # Max distance for metal complexation
    }

    class Config:
        arbitrary_types_allowed = True

    def to_plip_config(self) -> dict:
        """Convert to dictionary format for PLIP configuration"""
        return {
            "verbose": self.verbose,
            "model": self.model,
            "nohydro": self.nohydro,
            "peptides": self.peptides,
            "output_format": self.output_format,
            "maxthreads": self.maxthreads,
            "breakcomposite": self.breakcomposite,
            "altlocation": self.altlocation,
            "nofix": self.nofix,
            "keepmod": self.keepmod,
            "dnareceptor": self.dnareceptor,
            "chains": self.chains
        }
