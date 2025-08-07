from pydantic import BaseModel
from typing import Optional

# Empty configuration class that will make us use defaults from the inference scripts
class AutoDockConfig(BaseModel):
    # Removing all default values to ensure we use the defaults in the inference scripts
    pass

class TaskStatus(BaseModel):
    task_id: str
    status: str
    error: Optional[str] = None
    output_dir: Optional[str] = None
