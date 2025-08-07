from dataclasses import dataclass
from typing import Dict, Optional, Any

@dataclass
class MoleculeResponse:
    """
    Represents a molecule response from the server

    Attributes:
        graph: Molecular graph structure
        smiles: SMILES representation of the molecule
        inchi_key: InChI key identifier
        molecular_formula: Chemical formula
        molecular_weight: Molecular mass in g/mol (maps to "molecular_weight" in JSON)
        xlogp: Calculated LogP value (maps to "logp" in JSON)
        hbd: Number of hydrogen bond donors
        hba: Number of hydrogen bond acceptors
        tpsa: Topological polar surface area
        rb: Number of rotatable bonds
        qed: Quantitative Estimate of Drug-likeness
        sas: Synthetic Accessibility Score (maps to "sa_score" in JSON)
        lipinski_count: Number of Lipinski's Rule of Five violations
        docking_score: Molecular docking score (maps to "docking" in JSON)
        favorite: User favorite flag
        validation: Validation results
        s3_path: Path to molecule data in S3 storage
        plip_interactions: PLIP interactions information
        status: Molecule status (e.g., "modified", "new")
        activity: Activity status of the molecule (e.g., "active", "inactive")
    """
    graph: Dict[str, Any]
    smiles: Optional[str] = None
    inchi_key: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    xlogp: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    tpsa: Optional[float] = None
    rb: Optional[int] = None
    qed: Optional[float] = None
    sas: Optional[float] = None
    lipinski_count: Optional[float] = None
    docking_score: Optional[float] = None
    favorite: Optional[bool] = None
    validation: Optional[Any] = None
    s3_path: Optional[str] = None
    plip_interactions: Optional[str] = None
    status: Optional[str] = None
    activity: Optional[str] = None

def molecule_response_to_json(molecule: MoleculeResponse) -> Dict[str, Any]:
    """Convert the molecule response to a JSON-compatible dictionary"""
    return {
        "graph": molecule.graph,
        "smiles": molecule.smiles,
        "inchiKey": molecule.inchi_key,
        "molecularFormula": molecule.molecular_formula,
        "molecularWeight": molecule.molecular_weight,
        "xlogp": molecule.xlogp,
        "hbd": molecule.hbd,
        "hba": molecule.hba,
        "tpsa": molecule.tpsa,
        "rb": molecule.rb,
        "qed": molecule.qed,
        "sas": molecule.sas,
        "lipinskiCount": molecule.lipinski_count,
        "dockingScore": molecule.docking_score,
        "favorite": molecule.favorite,
        "validation": molecule.validation,
        "s3Path": molecule.s3_path,
        "plipInteractions": molecule.plip_interactions,
        "status": molecule.status,
        "activity": molecule.activity
    }

def molecule_response_from_json(data: Dict[str, Any]) -> MoleculeResponse:
    """Create a MoleculeResponse instance from a JSON-compatible dictionary"""
    return MoleculeResponse(
        graph=data.get("graph", {}),
        smiles=data.get("smiles"),
        inchi_key=data.get("inchiKey"),
        molecular_formula=data.get("molecularFormula"),
        molecular_weight=data.get("molecularWeight"),
        xlogp=data.get("xlogp"),
        hbd=data.get("hbd"),
        hba=data.get("hba"),
        tpsa=data.get("tpsa"),
        rb=data.get("rb"),
        qed=data.get("qed"),
        sas=data.get("sas"),
        lipinski_count=data.get("lipinskiCount"),
        docking_score=data.get("dockingScore"),
        favorite=data.get("favorite"),
        validation=data.get("validation"),
        s3_path=data.get("s3Path"),
        plip_interactions=data.get("plipInteractions"),
        status=data.get("status"),
        activity=data.get("activity")
    )
