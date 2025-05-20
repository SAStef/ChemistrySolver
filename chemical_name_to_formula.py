import pubchempy as pcp

def get_formula_from_name(compound_name):
    """
    Retrieves chemical formula for a given compound name using PubChemPy
    
    Args:
        compound_name (str): Name of the chemical compound
        
    Returns:
        dict: Dictionary containing formula and other information if successful,
              or error message if not successful
    """
    try:
        compounds = pcp.get_compounds(compound_name, 'name')
        if compounds:
            compound = compounds[0]  # Take the first match
            return {
                'success': True,
                'name': compound_name,
                'formula': compound.molecular_formula,
                'iupac_name': compound.iupac_name,
                'weight': compound.molecular_weight
            }
        else:
            return {
                'success': False,
                'error': f"Could not find compound: {compound_name}"
            }
    except Exception as e:
        return {
            'success': False,
            'error': f"Error: {str(e)}"
        }

def main():
    print("Chemical Name to Formula Converter")
    print("----------------------------------")
    
    while True:
        compound_name = input("\nEnter chemical name (or 'q' to quit): ")
        
        if compound_name.lower() == 'q':
            break
            
        result = get_formula_from_name(compound_name)
        
        if result['success']:
            print(f"\nFormula: {result['formula']}")
            print(f"IUPAC Name: {result['iupac_name']}")
            print(f"Molecular Weight: {result['weight']} g/mol")
        else:
            print(f"\nError: {result['error']}")

if __name__ == "__main__":
    main()