import pandas as pd
import os
import time

## Adapted from https://github.com/aolgac/plipinteractor

# path to the input directory
input_directory = os.getcwd()  # current working directory ./
# path to the output file
output_file = "residue_all_interactions.csv"


# Define functions to extract interaction information
def extract_hydrophobic_int(file_path, source):
    data = []
    in_hydrophobic_int_section = False

    with open(file_path, "r") as file:
        # extract only relevant part (from hyd int to next **)
        for line in file:
            line = line.strip()
            if line.startswith("**Hydrophobic Interactions"):
                in_hydrophobic_int_section = True
                continue
            elif line.startswith("**"):
                in_hydrophobic_int_section = False
                continue
            if in_hydrophobic_int_section and line.startswith("|"):
                # split lines based on delimiter
                columns = [col.strip() for col in line.split("|") if col.strip()]
                if len(columns) < 9:  # Ensure there are enough columns
                    continue
                    
                # extract the relevant columns
                RESNR = columns[0]
                RESTYPE = columns[1]
                RESCHAIN = columns[2]
                RESNR_LIG = columns[3]
                RESTYPE_LIG = columns[4]
                RESCHAIN_LIG = columns[5]
                LIGCARBONIDX = columns[7]
                PROTCARBONIDX = columns[8]
                
                # Extract distance if available
                DIST = columns[6] if len(columns) > 6 else ""
                
                # append data with source and distance info
                data.append([RESNR, RESTYPE, RESCHAIN, RESNR_LIG, RESTYPE_LIG, RESCHAIN_LIG, 
                            LIGCARBONIDX, PROTCARBONIDX, DIST, source])
    return data


def extract_hydrogen_bonds(file_path, source):
    # store data in a list
    data = []
    in_hydrogen_bonds_section = False

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # extract only relevant part (from hydrogen bond to next **)
            if line.startswith("**Hydrogen Bonds**"):
                in_hydrogen_bonds_section = True
                continue
            elif line.startswith("**"):
                in_hydrogen_bonds_section = False
                continue

            if in_hydrogen_bonds_section and line.startswith("|"):
                # split the line into columns based on the delimiter "|"
                columns = [col.strip() for col in line.split("|") if col.strip()]
                if len(columns) < 14:  # Ensure there are enough columns
                    continue
                    
                # extract the relevant columns
                RESNR = columns[0]
                RESTYPE = columns[1]
                RESCHAIN = columns[2]
                RESNR_LIG = columns[3]
                RESTYPE_LIG = columns[4]
                RESCHAIN_LIG = columns[5]
                DONORIDX = columns[11]
                ACCEPTORIDX = columns[13]
                
                # Extract distance information if available - correct indices based on the report format
                # In the report, columns are: SIDECHAIN | DIST_H-A | DIST_D-A | DON_ANGLE | ...
                DIST_H_A = columns[7] if len(columns) > 7 else ""  # DIST_H-A
                DIST_D_A = columns[8] if len(columns) > 8 else ""  # DIST_D-A
                DON_ANGLE = columns[9] if len(columns) > 9 else ""  # DON_ANGLE
                
                # append the extracted data to the list with distance info
                data.append([RESNR, RESTYPE, RESCHAIN, RESNR_LIG, RESTYPE_LIG, RESCHAIN_LIG,
                            DONORIDX, ACCEPTORIDX, DIST_H_A, DIST_D_A, DON_ANGLE, source])
    return data


def extract_salt_bridges(file_path, source):
    data = []
    in_salt_bridges_section = False

    with open(file_path, "r") as file:
        # extract only relevant part (from salt bridges to next **)
        for line in file:
            line = line.strip()
            if line.startswith("**Salt Bridges**"):
                in_salt_bridges_section = True
                continue
            elif line.startswith("**"):
                in_salt_bridges_section = False
                continue
            if in_salt_bridges_section and line.startswith("|"):
                # split lines based on delimiter
                columns = [col.strip() for col in line.split("|") if col.strip()]
                if len(columns) < 11:  # Ensure there are enough columns
                    continue
                    
                # extract the relevant columns
                RESNR = columns[0]
                RESTYPE = columns[1]
                RESCHAIN = columns[2]
                RESNR_LIG = columns[4]
                RESTYPE_LIG = columns[5]
                RESCHAIN_LIG = columns[6]
                PROT_IDX_LIST = columns[3]
                LIG_IDX_LIST = columns[10]
                
                # Extract distance if available
                DIST = columns[8] if len(columns) > 8 else ""
                
                # append data with distance info
                data.append([RESNR, RESTYPE, RESCHAIN, RESNR_LIG, RESTYPE_LIG, RESCHAIN_LIG,
                            PROT_IDX_LIST, LIG_IDX_LIST, DIST, source])
    return data


def extract_pi_cation_int(file_path, source):
    data = []
    in_pi_cation_section = False

    with open(file_path, "r") as file:
        # extract only relevant part (from pi-cation to next **)
        for line in file:
            line = line.strip()
            if line.startswith("**pi-Cation Interactions**"):
                in_pi_cation_section = True
                continue
            elif line.startswith("**"):
                in_pi_cation_section = False
                continue
            if in_pi_cation_section and line.startswith("|"):
                # split lines based on delimiter
                columns = [col.strip() for col in line.split("|") if col.strip()]
                if len(columns) < 12:  # Ensure there are enough columns
                    continue
                    
                # extract the relevant columns
                RESNR = columns[0]
                RESTYPE = columns[1]
                RESCHAIN = columns[2]
                RESNR_LIG = columns[4]
                RESTYPE_LIG = columns[5]
                RESCHAIN_LIG = columns[6]
                PROT_IDX_LIST = columns[3]
                LIG_IDX_LIST = columns[11]
                
                # Extract distance if available
                DIST = columns[8] if len(columns) > 8 else ""
                
                # append data with distance info
                data.append([RESNR, RESTYPE, RESCHAIN, RESNR_LIG, RESTYPE_LIG, RESCHAIN_LIG,
                            PROT_IDX_LIST, LIG_IDX_LIST, DIST, source])
    return data


def extract_pi_stacking(file_path, source):
    data = []
    in_pi_stacking_section = False

    with open(file_path, "r") as file:
        # extract only relevant part (from pi-Stacking to next **)
        for line in file:
            line = line.strip()
            if line.startswith("**pi-Stacking**"):
                in_pi_stacking_section = True
                continue
            elif line.startswith("**"):
                in_pi_stacking_section = False
                continue
            if in_pi_stacking_section and line.startswith("|"):
                # split lines based on delimiter
                columns = [col.strip() for col in line.split("|") if col.strip()]
                if len(columns) < 12:  # Ensure there are enough columns
                    continue
                    
                # extract the relevant columns
                RESNR = columns[0]
                RESTYPE = columns[1]
                RESCHAIN = columns[2]
                RESNR_LIG = columns[3]
                RESTYPE_LIG = columns[4]
                RESCHAIN_LIG = columns[5]
                PROT_IDX_LIST = columns[6]
                LIG_IDX_LIST = columns[11]
                
                # Extract distance if available
                DIST = columns[8] if len(columns) > 8 else ""
                
                # append data with distance info
                data.append([RESNR, RESTYPE, RESCHAIN, RESNR_LIG, RESTYPE_LIG, RESCHAIN_LIG,
                            PROT_IDX_LIST, LIG_IDX_LIST, DIST, source])
    return data


def extract_halogen_bonds(file_path, source):
    data = []
    in_halogen_bonds_section = False

    with open(file_path, "r") as file:
        # extract only relevant part (from halogen bonds to next **)
        for line in file:
            line = line.strip()
            if line.startswith("**Halogen Bonds**"):
                in_halogen_bonds_section = True
                continue
            elif line.startswith("**"):
                in_halogen_bonds_section = False
                continue
            if in_halogen_bonds_section and line.startswith("|"):
                # split lines based on delimiter
                columns = [col.strip() for col in line.split("|") if col.strip()]
                if len(columns) < 13:  # Ensure there are enough columns
                    continue
                    
                # extract the relevant columns
                RESNR = columns[0]
                RESTYPE = columns[1]
                RESCHAIN = columns[2]
                RESNR_LIG = columns[3]
                RESTYPE_LIG = columns[4]
                RESCHAIN_LIG = columns[5]
                DON_IDX = columns[10]
                ACC_IDX = columns[12]
                
                # Extract distance and angle information if available
                DIST = columns[7] if len(columns) > 7 else ""
                DON_ANGLE = columns[8] if len(columns) > 8 else ""
                ACC_ANGLE = columns[9] if len(columns) > 9 else ""
                
                # append data with distance and angle info
                data.append([RESNR, RESTYPE, RESCHAIN, RESNR_LIG, RESTYPE_LIG, RESCHAIN_LIG,
                            DON_IDX, ACC_IDX, DIST, DON_ANGLE, ACC_ANGLE, source])
    return data


# Store dataframes in lists
dfs_hyd = []
dfs_hb = []
dfs_sb = []
dfs_pc = []
dfs_ps = []
dfs_hal = []  # New list for halogen bonds

# Start time
start_time = time.time()

# Iterate through each file in the input directory
for file_name in os.listdir(input_directory):
    if file_name.endswith(".txt"):
        file_path = os.path.join(input_directory, file_name)

        # Hydrophobic interactions
        hydrophobic_int_data = extract_hydrophobic_int(file_path, file_name)
        if hydrophobic_int_data:
            df_hyd = pd.DataFrame(hydrophobic_int_data, 
                                columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                        "RESCHAIN_LIG", "LIGCARBONIDX", "PROTCARBONIDX", "DIST", "SOURCE"])
            df_hyd = df_hyd[~((df_hyd['LIGCARBONIDX'] == 'LIGCARBONIDX') & (df_hyd['PROTCARBONIDX'] == 'PROTCARBONIDX'))]
            df_hyd["INTERACTION TYPE"] = "Hydrophobic Interaction"
            dfs_hyd.append(df_hyd)

        # Hydrogen bonds
        hydrogen_bonds_data = extract_hydrogen_bonds(file_path, file_name)
        if hydrogen_bonds_data:
            df_hb = pd.DataFrame(hydrogen_bonds_data, 
                                columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                        "RESCHAIN_LIG", "DONORIDX", "ACCEPTORIDX", "DIST_H_A", "DIST_D_A", "DON_ANGLE", "SOURCE"])
            df_hb = df_hb[~((df_hb['DONORIDX'] == 'DONORIDX') & (df_hb['ACCEPTORIDX'] == 'ACCEPTORIDX'))]
            df_hb["INTERACTION TYPE"] = "Hydrogen Bond"
            dfs_hb.append(df_hb)

        # Salt bridges
        salt_bridges_data = extract_salt_bridges(file_path, file_name)
        if salt_bridges_data:
            df_sb = pd.DataFrame(salt_bridges_data, 
                                columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                        "RESCHAIN_LIG", "PROT_IDX_LIST", "LIG_IDX_LIST", "DIST", "SOURCE"])
            df_sb = df_sb[~((df_sb['PROT_IDX_LIST'] == 'PROT_IDX_LIST') & (df_sb['LIG_IDX_LIST'] == 'LIG_IDX_LIST'))]
            df_sb["INTERACTION TYPE"] = "Salt Bridge"
            dfs_sb.append(df_sb)

        # Pi-Cation interactions
        pi_cation_data = extract_pi_cation_int(file_path, file_name)
        if pi_cation_data:
            df_pc = pd.DataFrame(pi_cation_data, 
                                columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                        "RESCHAIN_LIG", "PROT_IDX_LIST", "LIG_IDX_LIST", "DIST", "SOURCE"])
            df_pc = df_pc[~((df_pc['PROT_IDX_LIST'] == 'PROT_IDX_LIST') & (df_pc['LIG_IDX_LIST'] == 'LIG_IDX_LIST'))]
            df_pc["INTERACTION TYPE"] = "Pi-Cation Interaction"
            dfs_pc.append(df_pc)

        # Pi-Stacking interactions
        pi_stacking_data = extract_pi_stacking(file_path, file_name)
        if pi_stacking_data:
            df_ps = pd.DataFrame(pi_stacking_data, 
                                columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                        "RESCHAIN_LIG", "PROT_IDX_LIST", "LIG_IDX_LIST", "DIST", "SOURCE"])
            df_ps = df_ps[~((df_ps['PROT_IDX_LIST'] == 'PROT_IDX_LIST') & (df_ps['LIG_IDX_LIST'] == 'LIG_IDX_LIST'))]
            df_ps["INTERACTION TYPE"] = "Pi-Stacking"
            dfs_ps.append(df_ps)

        # Halogen bonds
        halogen_bonds_data = extract_halogen_bonds(file_path, file_name)
        if halogen_bonds_data:
            df_hal = pd.DataFrame(halogen_bonds_data, 
                                columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                        "RESCHAIN_LIG", "DON_IDX", "ACC_IDX", "DIST", "DON_ANGLE", "ACC_ANGLE", "SOURCE"])
            df_hal = df_hal[~((df_hal['DON_IDX'] == 'DON_IDX') & (df_hal['ACC_IDX'] == 'ACC_IDX'))]
            df_hal["INTERACTION TYPE"] = "Halogen Bond"
            dfs_hal.append(df_hal)

# Process results only if we have data
all_dfs = []

# Process hydrophobic interactions if any found
if dfs_hyd:
    hydrophobic_bond = pd.concat(dfs_hyd, ignore_index=True)
    all_dfs.append(hydrophobic_bond)

# Process hydrogen bonds if any found
if dfs_hb:
    hydrogen_bond = pd.concat(dfs_hb, ignore_index=True)
    all_dfs.append(hydrogen_bond)

# Process salt bridges if any found
if dfs_sb:
    salt_bridge = pd.concat(dfs_sb, ignore_index=True)
    all_dfs.append(salt_bridge)

# Process pi-cation interactions if any found
if dfs_pc:
    pi_cation = pd.concat(dfs_pc, ignore_index=True)
    all_dfs.append(pi_cation)

# Process pi-stacking interactions if any found
if dfs_ps:
    pi_stacking = pd.concat(dfs_ps, ignore_index=True)
    all_dfs.append(pi_stacking)

# Process halogen bond interactions if any found
if dfs_hal:
    halogen_bond = pd.concat(dfs_hal, ignore_index=True)
    all_dfs.append(halogen_bond)

# Combine all interactions if we have any data
if all_dfs:
    # Concatenate all interaction dataframes
    all_interactions = pd.concat(all_dfs, ignore_index=True)
    
    # Count occurrences by each source & add as a new column
    all_interactions["COUNT_by_source"] = all_interactions.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE", "SOURCE"]
    ).transform("size")
    
    # Group by residue & reset index by count_by_source
    a_i = all_interactions.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE", "SOURCE"]
    ).size().reset_index(name='COUNT_by_source')
    
    # Count overall occurrences
    a_i["COUNT"] = a_i.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE"]
    ).transform("size")
    
    # Source + number of occurrences = sources => add as a new column
    a_i['SOURCES'] = a_i['SOURCE'] + ': ' + a_i['COUNT_by_source'].astype(str)
    
    # Merge duplicate sources
    final_a_i = a_i.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE", "COUNT"]
    )['SOURCES'].apply(lambda x: '|'.join(x)).reset_index()
    
    # Make unique list of sources and sort
    final_a_i['SOURCES'] = final_a_i['SOURCES'].apply(lambda x: '|'.join(sorted(set(x.split('|')))))
    
    # For debug/processing only, not in final output
    final_a_i["PERCENT %"] = (final_a_i["COUNT"] * 100) / 1002
    
    # Extract distance information from original data frames
    # For each type of interaction, we'll gather the distance data differently
    
    # Create a copy of the final dataframe to store distance info
    final_with_dist = final_a_i.copy()
    
    # Add distance information columns
    final_with_dist["DIST"] = ""
    final_with_dist["DIST_H_A"] = ""
    final_with_dist["DIST_D_A"] = ""
    final_with_dist["DON_ANGLE"] = ""  # Will be dropped later
    final_with_dist["ACC_ANGLE"] = ""  # Will be dropped later
    
    # For each interaction in final_with_dist, find corresponding rows in all_interactions
    for idx, row in final_with_dist.iterrows():
        # Filter matching interactions
        mask = ((all_interactions["RESNR"] == row["RESNR"]) & 
                (all_interactions["RESTYPE"] == row["RESTYPE"]) & 
                (all_interactions["RESCHAIN"] == row["RESCHAIN"]) & 
                (all_interactions["RESNR_LIG"] == row["RESNR_LIG"]) & 
                (all_interactions["RESTYPE_LIG"] == row["RESTYPE_LIG"]) & 
                (all_interactions["RESCHAIN_LIG"] == row["RESCHAIN_LIG"]) & 
                (all_interactions["INTERACTION TYPE"] == row["INTERACTION TYPE"]))
        
        matching_rows = all_interactions[mask]
        
        if not matching_rows.empty:
            # Get average distances if they exist
            if "DIST" in matching_rows.columns:
                dist_values = matching_rows["DIST"].astype(str).str.replace("Å", "").str.strip()
                dist_values = pd.to_numeric(dist_values, errors='coerce')
                avg_dist = dist_values.mean()
                if not pd.isna(avg_dist):
                    final_with_dist.at[idx, "DIST"] = f"{avg_dist:.2f}Å"
            
            # Only for hydrogen bonds
            if row["INTERACTION TYPE"] == "Hydrogen Bond":
                # Handle DIST_H_A - this is the hydrogen-acceptor distance
                if "DIST_H_A" in matching_rows.columns:
                    dist_h_a_values = matching_rows["DIST_H_A"].astype(str).str.replace("Å", "").str.strip()
                    dist_h_a_values = pd.to_numeric(dist_h_a_values, errors='coerce')
                    avg_dist_h_a = dist_h_a_values.mean()
                    if not pd.isna(avg_dist_h_a):
                        final_with_dist.at[idx, "DIST_H_A"] = f"{avg_dist_h_a:.2f}Å"
                
                # Handle DIST_D_A - this is the donor-acceptor distance
                if "DIST_D_A" in matching_rows.columns:
                    dist_d_a_values = matching_rows["DIST_D_A"].astype(str).str.replace("Å", "").str.strip()
                    dist_d_a_values = pd.to_numeric(dist_d_a_values, errors='coerce')
                    avg_dist_d_a = dist_d_a_values.mean()
                    if not pd.isna(avg_dist_d_a):
                        final_with_dist.at[idx, "DIST_D_A"] = f"{avg_dist_d_a:.2f}Å"
                
                # Handle DON_ANGLE - this is an angle, not a distance
                if "DON_ANGLE" in matching_rows.columns:
                    don_angle_values = matching_rows["DON_ANGLE"].astype(str).str.strip()
                    don_angle_values = pd.to_numeric(don_angle_values, errors='coerce')
                    avg_don_angle = don_angle_values.mean()
                    if not pd.isna(avg_don_angle):
                        final_with_dist.at[idx, "DON_ANGLE"] = f"{avg_don_angle:.2f}°"
            
            # Only for halogen bonds
            if row["INTERACTION TYPE"] == "Halogen Bond":
                # Handle DON_ANGLE for halogen bonds
                if "DON_ANGLE" in matching_rows.columns:
                    don_angle_values = matching_rows["DON_ANGLE"].astype(str).str.strip()
                    don_angle_values = pd.to_numeric(don_angle_values, errors='coerce')
                    avg_don_angle = don_angle_values.mean()
                    if not pd.isna(avg_don_angle):
                        final_with_dist.at[idx, "DON_ANGLE"] = f"{avg_don_angle:.2f}°"
                
                # Handle ACC_ANGLE for halogen bonds
                if "ACC_ANGLE" in matching_rows.columns:
                    acc_angle_values = matching_rows["ACC_ANGLE"].astype(str).str.strip()
                    acc_angle_values = pd.to_numeric(acc_angle_values, errors='coerce')
                    avg_acc_angle = acc_angle_values.mean()
                    if not pd.isna(avg_acc_angle):
                        final_with_dist.at[idx, "ACC_ANGLE"] = f"{avg_acc_angle:.2f}°"
    
    # Remove the DON_ANGLE column if it's not needed in the final output for hydrogen bonds
    if "DON_ANGLE" in final_with_dist.columns:
        # Keep DON_ANGLE for halogen bonds, remove it for other interactions
        for idx, row in final_with_dist.iterrows():
            if row["INTERACTION TYPE"] != "Halogen Bond":
                final_with_dist.at[idx, "DON_ANGLE"] = ""
    
    # Remove the ACC_ANGLE column if it's not needed for non-halogen bonds
    if "ACC_ANGLE" in final_with_dist.columns:
        for idx, row in final_with_dist.iterrows():
            if row["INTERACTION TYPE"] != "Halogen Bond":
                final_with_dist.at[idx, "ACC_ANGLE"] = ""
    
    # Remove empty columns from the output
    final_with_dist = final_with_dist.loc[:, (final_with_dist != '').any()]
    
    # Explicitly remove COUNT, PERCENT %, and SOURCES columns from final output
    final_with_dist = final_with_dist.drop(columns=['COUNT', 'PERCENT %', 'SOURCES'], errors='ignore')
    
    # Write output to a csv file
    final_with_dist.to_csv(output_file, index=False)
else:
    print("No interactions found in the input files.")

# Calculate execution time
execution_time = time.time() - start_time

print("Execution time:", execution_time, "seconds")

def process_plip_file(file_path, file_name):
    """
    Process a single PLIP output file and extract all interaction data.
    
    Args:
        file_path (str): Path to the PLIP output file
        file_name (str): Name of the file (used as source identifier)
        
    Returns:
        dict: Dictionary containing DataFrames for each interaction type
    """
    interaction_data = {}
    
    # Hydrophobic interactions
    hydrophobic_int_data = extract_hydrophobic_int(file_path, file_name)
    if hydrophobic_int_data:
        df_hyd = pd.DataFrame(hydrophobic_int_data, 
                            columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                    "RESCHAIN_LIG", "LIGCARBONIDX", "PROTCARBONIDX", "DIST", "SOURCE"])
        df_hyd = df_hyd[~((df_hyd['LIGCARBONIDX'] == 'LIGCARBONIDX') & (df_hyd['PROTCARBONIDX'] == 'PROTCARBONIDX'))]
        df_hyd["INTERACTION TYPE"] = "Hydrophobic Interaction"
        interaction_data['hydrophobic'] = df_hyd

    # Hydrogen bonds
    hydrogen_bonds_data = extract_hydrogen_bonds(file_path, file_name)
    if hydrogen_bonds_data:
        df_hb = pd.DataFrame(hydrogen_bonds_data, 
                            columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                    "RESCHAIN_LIG", "DONORIDX", "ACCEPTORIDX", "DIST_H_A", "DIST_D_A", "DON_ANGLE", "SOURCE"])
        df_hb = df_hb[~((df_hb['DONORIDX'] == 'DONORIDX') & (df_hb['ACCEPTORIDX'] == 'ACCEPTORIDX'))]
        df_hb["INTERACTION TYPE"] = "Hydrogen Bond"
        interaction_data['hydrogen'] = df_hb

    # Salt bridges
    salt_bridges_data = extract_salt_bridges(file_path, file_name)
    if salt_bridges_data:
        df_sb = pd.DataFrame(salt_bridges_data, 
                            columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                    "RESCHAIN_LIG", "PROT_IDX_LIST", "LIG_IDX_LIST", "DIST", "SOURCE"])
        df_sb = df_sb[~((df_sb['PROT_IDX_LIST'] == 'PROT_IDX_LIST') & (df_sb['LIG_IDX_LIST'] == 'LIG_IDX_LIST'))]
        df_sb["INTERACTION TYPE"] = "Salt Bridge"
        interaction_data['salt_bridge'] = df_sb

    # Pi-Cation interactions
    pi_cation_data = extract_pi_cation_int(file_path, file_name)
    if pi_cation_data:
        df_pc = pd.DataFrame(pi_cation_data, 
                            columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                    "RESCHAIN_LIG", "PROT_IDX_LIST", "LIG_IDX_LIST", "DIST", "SOURCE"])
        df_pc = df_pc[~((df_pc['PROT_IDX_LIST'] == 'PROT_IDX_LIST') & (df_pc['LIG_IDX_LIST'] == 'LIG_IDX_LIST'))]
        df_pc["INTERACTION TYPE"] = "Pi-Cation Interaction"
        interaction_data['pi_cation'] = df_pc

    # Pi-Stacking interactions
    pi_stacking_data = extract_pi_stacking(file_path, file_name)
    if pi_stacking_data:
        df_ps = pd.DataFrame(pi_stacking_data, 
                            columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                    "RESCHAIN_LIG", "PROT_IDX_LIST", "LIG_IDX_LIST", "DIST", "SOURCE"])
        df_ps = df_ps[~((df_ps['PROT_IDX_LIST'] == 'PROT_IDX_LIST') & (df_ps['LIG_IDX_LIST'] == 'LIG_IDX_LIST'))]
        df_ps["INTERACTION TYPE"] = "Pi-Stacking"
        interaction_data['pi_stacking'] = df_ps

    # Halogen bonds
    halogen_bonds_data = extract_halogen_bonds(file_path, file_name)
    if halogen_bonds_data:
        df_hal = pd.DataFrame(halogen_bonds_data, 
                            columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
                                    "RESCHAIN_LIG", "DON_IDX", "ACC_IDX", "DIST", "DON_ANGLE", "ACC_ANGLE", "SOURCE"])
        df_hal = df_hal[~((df_hal['DON_IDX'] == 'DON_IDX') & (df_hal['ACC_IDX'] == 'ACC_IDX'))]
        df_hal["INTERACTION TYPE"] = "Halogen Bond"
        interaction_data['halogen'] = df_hal
        
    return interaction_data


def combine_interactions(all_interaction_data):
    """
    Combine all interaction DataFrames and process them.
    
    Args:
        all_interaction_data (list): List of dictionaries containing interaction DataFrames
        
    Returns:
        DataFrame: Processed DataFrame with all interactions
    """
    # Store all dataframes by interaction type
    all_dfs_by_type = {
        'hydrophobic': [],
        'hydrogen': [],
        'salt_bridge': [],
        'pi_cation': [],
        'pi_stacking': [],
        'halogen': []
    }
    
    # Combine all interaction data by type
    for interaction_data in all_interaction_data:
        for key, df in interaction_data.items():
            if df is not None and not df.empty:
                all_dfs_by_type[key].append(df)
    
    # Concatenate by interaction type
    all_dfs = []
    for key, dfs in all_dfs_by_type.items():
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            all_dfs.append(combined_df)
    
    # Return None if no data found
    if not all_dfs:
        return None
    
    # Combine all interactions into one DataFrame
    all_interactions = pd.concat(all_dfs, ignore_index=True)
    
    # Count occurrences by each source & add as a new column
    all_interactions["COUNT_by_source"] = all_interactions.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE", "SOURCE"]
    ).transform("size")
    
    # Group by residue & reset index by count_by_source
    a_i = all_interactions.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE", "SOURCE"]
    ).size().reset_index(name='COUNT_by_source')
    
    # Count overall occurrences
    a_i["COUNT"] = a_i.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE"]
    ).transform("size")
    
    # Source + number of occurrences = sources => add as a new column
    a_i['SOURCES'] = a_i['SOURCE'] + ': ' + a_i['COUNT_by_source'].astype(str)
    
    # Merge duplicate sources
    final_a_i = a_i.groupby(
        ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "INTERACTION TYPE", "COUNT"]
    )['SOURCES'].apply(lambda x: '|'.join(x)).reset_index()
    
    # Make unique list of sources and sort
    final_a_i['SOURCES'] = final_a_i['SOURCES'].apply(lambda x: '|'.join(sorted(set(x.split('|')))))
    
    # For debug/processing only, not in final output
    final_a_i["PERCENT %"] = (final_a_i["COUNT"] * 100) / 1002
    
    # Extract distance information from original data frames
    # For each type of interaction, we'll gather the distance data differently
    
    # Create a copy of the final dataframe to store distance info
    final_with_dist = final_a_i.copy()
    
    # Add distance information columns
    final_with_dist["DIST"] = ""
    final_with_dist["DIST_H_A"] = ""
    final_with_dist["DIST_D_A"] = ""
    
    # For each interaction in final_with_dist, find corresponding rows in all_interactions
    for idx, row in final_with_dist.iterrows():
        # Filter matching interactions
        mask = ((all_interactions["RESNR"] == row["RESNR"]) & 
                (all_interactions["RESTYPE"] == row["RESTYPE"]) & 
                (all_interactions["RESCHAIN"] == row["RESCHAIN"]) & 
                (all_interactions["RESNR_LIG"] == row["RESNR_LIG"]) & 
                (all_interactions["RESTYPE_LIG"] == row["RESTYPE_LIG"]) & 
                (all_interactions["RESCHAIN_LIG"] == row["RESCHAIN_LIG"]) & 
                (all_interactions["INTERACTION TYPE"] == row["INTERACTION TYPE"]))
        
        matching_rows = all_interactions[mask]
        
        if not matching_rows.empty:
            # Get average distances if they exist
            if "DIST" in matching_rows.columns:
                dist_values = matching_rows["DIST"].astype(str).str.replace("Å", "").str.strip()
                dist_values = pd.to_numeric(dist_values, errors='coerce')
                avg_dist = dist_values.mean()
                if not pd.isna(avg_dist):
                    final_with_dist.at[idx, "DIST"] = f"{avg_dist:.2f}Å"
            
            # Only for hydrogen bonds
            if row["INTERACTION TYPE"] == "Hydrogen Bond":
                # Handle DIST_H_A - hydrogen-acceptor distance
                if "DIST_H_A" in matching_rows.columns:
                    dist_h_a_values = matching_rows["DIST_H_A"].astype(str).str.replace("Å", "").str.strip()
                    dist_h_a_values = pd.to_numeric(dist_h_a_values, errors='coerce')
                    avg_dist_h_a = dist_h_a_values.mean()
                    if not pd.isna(avg_dist_h_a):
                        final_with_dist.at[idx, "DIST_H_A"] = f"{avg_dist_h_a:.2f}Å"
                
                # Handle DIST_D_A - donor-acceptor distance
                if "DIST_D_A" in matching_rows.columns:
                    dist_d_a_values = matching_rows["DIST_D_A"].astype(str).str.replace("Å", "").str.strip()
                    dist_d_a_values = pd.to_numeric(dist_d_a_values, errors='coerce')
                    avg_dist_d_a = dist_d_a_values.mean()
                    if not pd.isna(avg_dist_d_a):
                        final_with_dist.at[idx, "DIST_D_A"] = f"{avg_dist_d_a:.2f}Å"
    
    # Remove columns that have no values
    final_with_dist = final_with_dist.loc[:, (final_with_dist != '').any()]
    
    # Explicitly remove COUNT, PERCENT %, and SOURCES columns from final output
    final_with_dist = final_with_dist.drop(columns=['COUNT', 'PERCENT %', 'SOURCES'], errors='ignore')
    
    return final_with_dist


def summarize_plip_output(input_path, output_file=None, sample_ratio=1.0):
    """
    Main function to extract interactions from PLIP output files and save to CSV.
    
    Args:
        input_path (str): Path to a single file or directory containing PLIP output files
        output_file (str, optional): Path to save the output CSV. Defaults to same name as input with ".csv" extension
                                    in the same directory as the input file/directory.
        sample_ratio (float, optional): Ratio of files to process (for large directories). Defaults to 1.0 (all files).
    
    Returns:
        bool: True if processing was successful, False otherwise
    """
    import os
    import time
    import random
    
    # Start time
    start_time = time.time()
    
    # Set default output file if not provided
    if output_file is None:
        if os.path.isdir(input_path):
            output_file = os.path.join(input_path, "residue_all_interactions.csv")
        else:
            # Use the same name as input file but with .csv extension
            input_basename = os.path.basename(input_path)
            input_name_without_ext = os.path.splitext(input_basename)[0]
            output_file = os.path.join(os.path.dirname(input_path), input_name_without_ext + ".csv")
    
    # List to store interaction data from all files
    all_interaction_data = []
    
    # Process a single file
    if os.path.isfile(input_path):
        if input_path.endswith(".txt"):
            file_name = os.path.basename(input_path)
            print(f"Processing file: {file_name}")
            interaction_data = process_plip_file(input_path, file_name)
            if interaction_data:
                all_interaction_data.append(interaction_data)
        else:
            print(f"Error: {input_path} is not a text file.")
            return False
    
    # Process a directory of files
    elif os.path.isdir(input_path):
        txt_files = [f for f in os.listdir(input_path) if f.endswith(".txt")]
        
        # Sample files if needed
        if sample_ratio < 1.0:
            num_files = int(len(txt_files) * sample_ratio)
            txt_files = random.sample(txt_files, num_files)
        
        if not txt_files:
            print(f"No text files found in {input_path}")
            return False
        
        print(f"Processing {len(txt_files)} files...")
        for file_name in txt_files:
            file_path = os.path.join(input_path, file_name)
            interaction_data = process_plip_file(file_path, file_name)
            if interaction_data:
                all_interaction_data.append(interaction_data)
    
    else:
        print(f"Error: {input_path} is not a valid file or directory.")
        return False
    
    # Combine and process all interaction data
    if all_interaction_data:
        final_df = combine_interactions(all_interaction_data)
        
        if final_df is not None:
            # Save to CSV
            final_df.to_csv(output_file, index=False)
            print(f"Results saved to {output_file}")
        else:
            print("No interactions found in the input files.")
            return False
    else:
        print("No interactions found in the input files.")
        return False
    
    # Calculate and print execution time
    execution_time = time.time() - start_time
    print(f"Execution time: {execution_time:.2f} seconds")
    
    return True


if __name__ == "__main__":
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Extract molecular interactions from PLIP output files and save to CSV.")
    parser.add_argument("input", help="Input file or directory containing PLIP output files")
    parser.add_argument("-o", "--output", help="Output CSV file path")
    parser.add_argument("-s", "--sample", type=float, default=1.0, 
                        help="Sample ratio for processing (0.0-1.0, for large directories)")
    
    args = parser.parse_args()
    
    # Run the main function with provided arguments
    summarize_plip_output(args.input, args.output, args.sample) 