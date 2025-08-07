import re
import torch
import selfies as sf
from rdkit import Chem
from transformers import T5Tokenizer, T5EncoderModel, GPT2LMHeadModel, AutoTokenizer
import os 
import gc  # Add garbage collection import
from tqdm import tqdm
## Molecule generation model is taken from https://github.com/HUBioDataLab/Prot2Mol

def get_device():
    """Automatically select GPU if available, otherwise CPU"""
    if torch.cuda.is_available():
        return "cuda"
    elif torch.backends.mps.is_available():
        return "mps"
    else:
        return "cpu"

def cleanup_memory():
    """Comprehensive memory cleanup function"""
    # Force garbage collection
    gc.collect()
    
    # Clear PyTorch caches based on device
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()  # Wait for all operations to complete
    elif torch.backends.mps.is_available():
        torch.mps.empty_cache()
        torch.mps.synchronize()  # Wait for all operations to complete
    
    # Force another garbage collection after cache clearing
    gc.collect()

def init_protein_model(device=None):
    """Initialize the ProtT5 model for protein embedding generation"""
    if device is None:
        device = get_device()
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', 
                                          do_lower_case=False, 
                                          legacy=True, 
                                          clean_up_tokenization_spaces=True)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").to(device)
    return tokenizer, model

def produce_prot_t5_embedding(batch, encoder_max_length, tokenizer, model, device=None):
    """Generate embeddings for a batch of protein sequences"""
    if device is None:
        device = get_device()
    sequence_examples = " ".join(list(re.sub(r"[UZOB]", "X", batch)))
    ids = tokenizer.encode_plus(sequence_examples, 
                              add_special_tokens=True, 
                              padding="max_length", 
                              max_length=encoder_max_length)
    input_ids = torch.tensor(ids['input_ids']).to(device).view(1,-1)
    attention_mask = torch.tensor(ids['attention_mask']).to(device).view(1,-1)
    
    with torch.no_grad():
        embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)
    
    return embedding_repr.last_hidden_state

def init_molecule_model(model_path, device=None):
    """Initialize the GPT molecule decoder model"""
    if device is None:
        device = get_device()
    tokenizer = AutoTokenizer.from_pretrained("HUBioDataLab/Prot2Mol", padding_side="left")
    model = GPT2LMHeadModel.from_pretrained("HUBioDataLab/Prot2Mol").to(device)
    return tokenizer, model

def generate_selfies(model, encoder_hidden_states, num_sequences=2, max_length=202, 
                      do_sample=True, attn_output=False, batch_size=8):
    """
    Generate molecule sequences from protein embeddings, with batch processing.
    
    Args:
        model: The GPT model for molecule generation
        encoder_hidden_states: Protein embeddings to condition on
        num_sequences: Total number of sequences to generate
        max_length: Maximum length of generated sequences
        do_sample: Whether to use sampling for generation
        attn_output: Whether to output attention weights
        batch_size: Maximum number of sequences to generate at once (default: 8)
    
    Returns:
        List of generated sequences
    """
    all_generated_tokens = []
    
    # Calculate total number of batches for progress bar
    total_batches = (num_sequences + batch_size - 1) // batch_size
    
    # Process in batches with progress bar
    remaining = num_sequences
    with tqdm(total=total_batches, desc="Generating molecules") as pbar:
        while remaining > 0:
            current_batch_size = min(batch_size, remaining)
            
            batch_tokens = model.generate(
                encoder_hidden_states=encoder_hidden_states,
                num_return_sequences=current_batch_size,
                do_sample=do_sample,
                max_length=max_length,
                pad_token_id=1,
                bos_token_id=1,
                output_attentions=attn_output
            )
            
            all_generated_tokens.append(batch_tokens)
            remaining -= current_batch_size
            pbar.update(1)
    
    # Concatenate all batches
    if len(all_generated_tokens) > 1:
        # Concatenate along the first dimension (batch dimension)
        return torch.cat(all_generated_tokens, dim=0)
    else:
        return all_generated_tokens[0]

def process_generated_molecules(generated_sequences, tokenizer):
    """Convert generated SELFIES to SMILES and calculate validity"""
    valid_smiles = []
    valid_mols = []
    valid_count = 0
    
    for seq in generated_sequences:
        selfies_str = tokenizer.decode(seq, skip_special_tokens=True)
        try:
            # Convert SELFIES to SMILES
            smiles = sf.decoder(selfies_str)
            # Check validity using RDKit
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                valid_smiles.append(smiles) 
                valid_mols.append(mol)
                valid_count += 1
        except:
            continue
    
    validity_rate = valid_count / len(generated_sequences) if len(generated_sequences) > 0 else 0
    return valid_smiles, valid_mols, validity_rate

def prot2mol(protein_sequence, encoder_max_length=1000, num_molecules=1, max_length=202, output_dir=None):
    """
    Main function to generate molecules from a protein sequence
    
    Args:
        protein_sequence (str): Input protein sequence
        model_path (str): Path to the pretrained molecule generation model
        encoder_max_length (int): Maximum length for protein encoder
        num_molecules (int): Number of molecules to generate
        max_length (int): Maximum length of generated molecule sequences
    
    Returns:
        list: Valid generated molecules in SMILES format
        float: Validity rate of generated molecules
    """
    try:
        # Get device
        device = get_device()
        print("getting device")
        # Initialize protein model
        print(f"Device: {device}")
        prot_tokenizer, prot_model = init_protein_model(device)
        print("initializing protein model")
        
        # Generate protein embedding
        protein_embedding = produce_prot_t5_embedding(
            protein_sequence, 
            encoder_max_length, 
            prot_tokenizer, 
            prot_model, 
            device
        )
        
        # Clean up protein model immediately after use
        del prot_model
        del prot_tokenizer
        cleanup_memory()  # Comprehensive cleanup after protein model
    
        print("generating protein embedding")
        # Initialize molecule model
        mol_tokenizer, mol_model = init_molecule_model(device)
        print("initializing molecule model")
        
        # Generate molecule sequences
        generated_sequences = generate_selfies(
            mol_model,
            protein_embedding,
            num_sequences=num_molecules,
            max_length=max_length,
            batch_size=8  # Add the batch_size parameter
        )
        
        # Clean up after molecule generation
        del protein_embedding
        cleanup_memory()
            
        print("generating molecule sequences")
        # Process and validate generated molecules
        valid_smiles, valid_mols, validity_rate = process_generated_molecules(generated_sequences, mol_tokenizer)
        
        # Final cleanup
        del generated_sequences
        del mol_model
        del mol_tokenizer
        cleanup_memory()  # Final comprehensive cleanup
            
        print("processing generated molecules")
        if output_dir:
            # Convert output_dir to string if it's a Path object
            output_dir = str(output_dir)
            os.makedirs(output_dir, exist_ok=True)
            output_path = output_dir + "/generated_molecules.txt"
            print(f"Writing results to: {output_path}")
            try:
                with open(output_path, "w") as f:
                    for smiles in valid_smiles:
                        f.write(smiles + ':' + str(validity_rate) + "\n")
                for idx, mol in enumerate(valid_mols):
                    sdf_path = os.path.join(output_dir, f"molecule_{idx}.sdf")
                    writer = Chem.SDWriter(sdf_path)
                    writer.write(mol)
                    writer.close()
            except Exception as e:
                print(f"Error writing to file: {str(e)}")
                print(f"Output directory type: {type(output_dir)}")
                raise

        return valid_smiles
    except Exception as e:
        print(f"Error in prot2mol: {str(e)}")
        print(f"Error type: {type(e)}")
        print(f"Output directory type: {type(output_dir)}")
        raise
    finally:
        # Ensure cleanup happens even if there's an error
        cleanup_memory()

from src.tools.tool_definitions import tool_registry

@tool_registry.register("generate_molecules")
def generate_molecules(arguments: dict, output_dir_id: str) -> dict:
    print("Prot2Mol is being used")
    sequence = arguments.get("protein_sequence") or tool_registry.previous_sequence
    if not sequence:
        return {"error": "No protein sequence provided"}

    try:
        result = prot2mol(
            protein_sequence=sequence,
            output_dir=tool_registry.output_dir / output_dir_id / "mol_gen",
            num_molecules=arguments.get("num_molecules", 1)
        )
        return result
    finally:
        # Enhanced cleanup specifically for the tool function
        cleanup_memory()
        import time
        time.sleep(60)
        print("Final memory cleanup completed.")