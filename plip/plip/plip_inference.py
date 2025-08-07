from pathlib import Path
from typing import List
import contextlib
import os

from plip.structure.preparation import PDBComplex
from plip.exchange.report import StructureReport
from plip.basic import config
from plip.exchange.webservices import fetch_pdb

class PLIPInference:
    def __init__(self, task_id: str):
        self.task_id = task_id
        self.output_dir = Path(f"storage/{task_id}")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    async def run(self, pdb_file: str, output_format: List[str] = ["xml", "txt"]):
        try:
            error_file = Path(self.output_dir) / "debug.log"
            with open(error_file, "w") as f:
                f.write(f"Input pdb_file: {pdb_file}\n")
                f.write(f"Output dir: {self.output_dir}\n")
                f.write(f"Output format: {output_format}\n")

            self._setup_plip_config(output_format)
            with open(error_file, "a") as f:
                f.write(f"Config set - XML: {config.XML}, TXT: {config.TXT}\n")

            complex = PDBComplex()
            complex.output_path = str(self.output_dir)

            if pdb_file.startswith('pdb:'):
                pdb_id = pdb_file[4:]
                pdb_string, _ = fetch_pdb(pdb_id)
                temp_pdb = Path(self.output_dir) / "input.pdb"
                temp_pdb.write_text(pdb_string)
                complex.load_pdb(str(temp_pdb))
            else:
                complex.load_pdb(pdb_file)

            with open(error_file, "a") as f:
                f.write(f"Found {len(complex.ligands)} ligands\n")

            for i, ligand in enumerate(complex.ligands, 1):
                try:
                    complex.characterize_complex(ligand)
                    report = StructureReport(complex)

                    base_name = f"{ligand.hetid}_{ligand.chain}_{ligand.position}"

                    if config.XML:
                        xml_path = self.output_dir / f"{base_name}.xml"
                        with open(error_file, "a") as f:
                            f.write(f"Writing XML to {xml_path}\n")
                        report.write_xml(as_string=False)

                    if config.TXT:
                        txt_path = self.output_dir / f"{base_name}.txt"
                        with open(error_file, "a") as f:
                            f.write(f"Writing TXT to {txt_path}\n")
                        report.write_txt(as_string=False)

                    with open(error_file, "a") as f:
                        f.write(f"Processed ligand {i}\n")

                except Exception as ligand_error:
                    with open(error_file, "a") as f:
                        f.write(f"Error processing ligand {i}: {str(ligand_error)}\n")
                    raise

            with open(error_file, "a") as f:
                f.write(f"Final directory contents: {list(self.output_dir.glob('*'))}\n")

        except Exception as e:
            error_file = Path(self.output_dir) / "debug.log"
            with open(error_file, "a") as f:
                f.write(f"Error in PLIP analysis: {str(e)}\n")
            raise

    def _setup_plip_config(self, output_format: List[str]):
        config.VERBOSE = False
        config.XML = "xml" in output_format
        config.TXT = "txt" in output_format
        config.OUTPATH = str(self.output_dir)
