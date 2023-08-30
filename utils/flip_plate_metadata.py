import pandas as pd

def flip_number(num: int) -> int:
    """
    Flips the given number according to a pre-defined dictionary.

    Args:
        num: The number to flip.

    Returns:
        The flipped number.
    """
    flip_dict = {
        1: 12, 2: 11, 3: 10, 4: 9, 5: 8, 6: 7,
        7: 6, 8: 5, 9: 4, 10: 3, 11: 2, 12: 1
    }
    return flip_dict.get(num, num)

def flip_str_combination(s: str) -> str:
    """
    Flips the end combination of a string like Chrom_40_1_b_A1.

    Args:
        s: The string to flip.

    Returns:
        The flipped string.
    """
    if pd.isna(s) or not isinstance(s, str):
        return s
    parts = s.split('_')[-1]
    letter = parts[:-2] if len(parts) > 2 else parts[:-1]
    num = int(parts[-2:]) if len(parts) > 2 else int(parts[-1])
    new_num = flip_number(num)
    return s.rsplit('_', 1)[0] + '_' + letter + str(new_num)

# Load the CSV file
df = pd.read_csv('acute_metadata.csv')

# Separate the values based on _h_ and _b_
df['rna_body'] = df['rna_head'].apply(lambda x: x if '_b_' in str(x) else None)
df['rna_head'] = df['rna_head'].apply(lambda x: x if '_h_' in str(x) else None)

# Apply the flip operation on 'rna_body' column
df['rna_head'] = df['rna_head'].apply(flip_str_combination)

# Remove NaNs or None and shift values up in 'rna_body' and 'rna_head' columns
df['rna_body'] = df['rna_body'].dropna().reset_index(drop=True)
df['rna_head'] = df['rna_head'].dropna().reset_index(drop=True)

# Save the modified dataframe
df.to_csv('modified_file.csv', index=False)

