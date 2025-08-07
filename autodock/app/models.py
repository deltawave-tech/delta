from pydantic import BaseModel, Field
from typing import Optional, Union, Dict, Any

class AutoDockConfig(BaseModel):
    # Common parameters
    num_runs: int = 10
    seed: int = 42
    
    # Engine preference - using str instead of Literal for compatibility with pydantic v1
    engine: str = "auto"  # Auto will try AutoDock-GPU first, then fall back to Vina
    
    # AutoDock-GPU specific parameters
    local_search_method: str = "ad"  # ADADELTA (ad), Solis-Wets (sw), Steepest-Descent (sd), FIRE
    heuristics: int = 1  # Use heuristics to set evaluation count based on ligand complexity
    max_evaluations: Optional[int] = None  # Maximum number of energy evaluations, None = use heuristics
    autostop: bool = True  # Automatically stop when converged
    
    # Vina specific parameters
    exhaustiveness: int = 8  # Search exhaustiveness
    num_modes: int = 9  # Number of binding modes to generate
    energy_range: float = 3  # Maximum energy difference between best and worst binding mode
    
    # Search box parameters (used by both engines)
    center_x: float = Field(default=0.0, description="X coordinate of search box center")
    center_y: float = Field(default=0.0, description="Y coordinate of search box center")
    center_z: float = Field(default=0.0, description="Z coordinate of search box center")
    size_x: float = Field(default=40.0, description="Width of search box in X dimension")
    size_y: float = Field(default=40.0, description="Height of search box in Y dimension")
    size_z: float = Field(default=40.0, description="Depth of search box in Z dimension")
    
    # Debug parameters
    timeout: int = 1800  # Timeout in seconds for AutoDock-GPU (30 minutes)
    verbose: bool = False  # Extra logging for debugging
    
    class Config:
        # Allow extra fields for backward compatibility
        extra = "allow"
        
    # For compatibility with different pydantic versions
    def dict(self) -> Dict[str, Any]:
        """Return a dictionary representation of the model."""
        try:
            # Pydantic v1
            return super().dict()
        except AttributeError:
            # Pydantic v2
            return super().model_dump()
