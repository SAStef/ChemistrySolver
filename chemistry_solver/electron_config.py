"""
Electron Configuration Solver
A tool for calculating electron configurations, orbital analysis, and related properties of elements.
"""

# Dictionary of elements with atomic number and symbol
ELEMENTS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109,
    "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
}

# Create reverse lookup dictionary
ATOMIC_NUMBERS = {v: k for k, v in ELEMENTS.items()}

# Standard orbital filling order for electron configurations
ORBITAL_FILLING_ORDER = [
    "1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", 
    "5f", "6d", "7p", "8s", "5g", "6f", "7d", "8p", "9s"
]

# Maximum electrons per orbital type
ORBITAL_CAPACITIES = {
    "s": 2,
    "p": 6,
    "d": 10,
    "f": 14,
    "g": 18
}

# Elements with special/anomalous electron configurations due to half-filled or filled subshell stability
SPECIAL_CONFIGURATIONS = {
    "Cr": {"expected": "[Ar] 4s2 3d4", "actual": "[Ar] 4s1 3d5"},
    "Cu": {"expected": "[Ar] 4s2 3d9", "actual": "[Ar] 4s1 3d10"},
    "Nb": {"expected": "[Kr] 5s2 4d3", "actual": "[Kr] 5s1 4d4"},
    "Mo": {"expected": "[Kr] 5s2 4d4", "actual": "[Kr] 5s1 4d5"},
    "Ru": {"expected": "[Kr] 5s2 4d6", "actual": "[Kr] 5s1 4d7"},
    "Rh": {"expected": "[Kr] 5s2 4d7", "actual": "[Kr] 5s1 4d8"},
    "Pd": {"expected": "[Kr] 5s2 4d8", "actual": "[Kr] 4d10"},
    "Ag": {"expected": "[Kr] 5s2 4d9", "actual": "[Kr] 5s1 4d10"},
    "La": {"expected": "[Xe] 6s2 5d1", "actual": "[Xe] 6s2 5d1"},  # Not anomalous, but included for completeness
    "Ce": {"expected": "[Xe] 6s2 4f1 5d1", "actual": "[Xe] 6s2 4f1 5d1"},
    "Gd": {"expected": "[Xe] 6s2 4f7 5d1", "actual": "[Xe] 6s2 4f7 5d1"},
    "Au": {"expected": "[Xe] 6s2 4f14 5d9", "actual": "[Xe] 6s1 4f14 5d10"},
    "Ac": {"expected": "[Rn] 7s2 6d1", "actual": "[Rn] 7s2 6d1"},
    "Th": {"expected": "[Rn] 7s2 6d2", "actual": "[Rn] 7s2 6d2"},
    "Pa": {"expected": "[Rn] 7s2 5f2 6d1", "actual": "[Rn] 7s2 5f2 6d1"},
    "U": {"expected": "[Rn] 7s2 5f3 6d1", "actual": "[Rn] 7s2 5f3 6d1"},
    "Np": {"expected": "[Rn] 7s2 5f4 6d1", "actual": "[Rn] 7s2 5f4 6d1"},
    "Cm": {"expected": "[Rn] 7s2 5f7 6d1", "actual": "[Rn] 7s2 5f7 6d1"}
}

# Noble gas configurations (for shorthand notation)
NOBLE_GAS_CONFIGS = {
    2: "He",
    10: "Ne",
    18: "Ar",
    36: "Kr",
    54: "Xe",
    86: "Rn",
    118: "Og"
}

def get_element_by_symbol(symbol):
    """Get atomic number from element symbol."""
    symbol = symbol.strip()
    if symbol in ELEMENTS:
        return {
            "symbol": symbol,
            "atomic_number": ELEMENTS[symbol],
            "success": True
        }
    return {
        "success": False,
        "error": f"Element symbol '{symbol}' not found."
    }

def get_element_by_number(atomic_number):
    """Get element symbol from atomic number."""
    try:
        atomic_number = int(atomic_number)
        if atomic_number in ATOMIC_NUMBERS:
            return {
                "symbol": ATOMIC_NUMBERS[atomic_number],
                "atomic_number": atomic_number,
                "success": True
            }
        return {
            "success": False,
            "error": f"Atomic number {atomic_number} not found."
        }
    except ValueError:
        return {
            "success": False,
            "error": f"Invalid atomic number: {atomic_number}"
        }

def get_nearest_noble_gas(atomic_number):
    """Get the nearest noble gas with atomic number less than the given element."""
    noble_gases = sorted(NOBLE_GAS_CONFIGS.keys())
    for gas in noble_gases:
        if gas >= atomic_number:
            if gas == atomic_number:  # If the element itself is a noble gas
                if gas > 2:  # For noble gases after He
                    return noble_gases[noble_gases.index(gas) - 1]
                return None  # For He
            else:
                return noble_gases[noble_gases.index(gas) - 1]
    return noble_gases[-1]  # Return the highest noble gas if none found

def get_full_electron_configuration(atomic_number):
    """Generate the full electron configuration for an element."""
    if atomic_number <= 0:
        return {
            "success": False,
            "error": "Atomic number must be positive."
        }
    
    symbol = ATOMIC_NUMBERS.get(atomic_number, f"Element {atomic_number}")
    
    # Check for special configuration
    if symbol in SPECIAL_CONFIGURATIONS:
        return {
            "success": True,
            "element": symbol,
            "atomic_number": atomic_number,
            "full_configuration": SPECIAL_CONFIGURATIONS[symbol]["actual"],
            "standard_configuration": SPECIAL_CONFIGURATIONS[symbol]["expected"],
            "is_anomalous": True
        }
    
    # Standard electron configuration
    remaining_electrons = atomic_number
    configuration = []
    
    for orbital in ORBITAL_FILLING_ORDER:
        if remaining_electrons <= 0:
            break
        
        orbital_type = orbital[-1]  # Last character is the orbital type (s, p, d, f)
        max_electrons = ORBITAL_CAPACITIES[orbital_type]
        electrons_in_orbital = min(remaining_electrons, max_electrons)
        
        if electrons_in_orbital > 0:
            configuration.append(f"{orbital}{electrons_in_orbital}")
            remaining_electrons -= electrons_in_orbital
    
    full_config = " ".join(configuration)
    
    return {
        "success": True,
        "element": symbol,
        "atomic_number": atomic_number,
        "full_configuration": full_config,
        "is_anomalous": False
    }

def get_shorthand_electron_configuration(atomic_number):
    """Generate shorthand (noble gas) electron configuration for an element."""
    result = get_full_electron_configuration(atomic_number)
    
    if not result["success"]:
        return result
    
    # Check if this is a noble gas
    if atomic_number in NOBLE_GAS_CONFIGS:
        result["shorthand_configuration"] = result["full_configuration"]
        return result
        
    # Find nearest noble gas
    noble_gas_number = get_nearest_noble_gas(atomic_number)
    
    if noble_gas_number is None:  # For H and He
        result["shorthand_configuration"] = result["full_configuration"]
        return result
        
    noble_gas_symbol = NOBLE_GAS_CONFIGS[noble_gas_number]
    
    # Get full config and noble gas config
    noble_gas_result = get_full_electron_configuration(noble_gas_number)
    if not noble_gas_result["success"]:
        return noble_gas_result
    
    # If it's an anomalous configuration, use the actual rather than predicted
    if "is_anomalous" in result and result["is_anomalous"]:
        # Extract valence part from special configuration
        actual_config = result["full_configuration"]
        noble_gas_part = f"[{noble_gas_symbol}]"
        if noble_gas_part in actual_config:
            valence_part = actual_config.split(noble_gas_part)[1].strip()
            result["shorthand_configuration"] = f"[{noble_gas_symbol}] {valence_part}"
        else:
            # If anomalous config isn't in noble gas format, we'll construct it
            noble_gas_config = noble_gas_result["full_configuration"]
            full_config = result["full_configuration"]
            # This is a simplification and might need refinement for accurate anomalous configs
            valence_part = " ".join([orb for orb in full_config.split() if orb not in noble_gas_config.split()])
            result["shorthand_configuration"] = f"[{noble_gas_symbol}] {valence_part}"
    else:
        # Standard approach for non-anomalous
        noble_gas_config = noble_gas_result["full_configuration"]
        full_config = result["full_configuration"]
        valence_part = " ".join([orb for orb in full_config.split() if orb not in noble_gas_config.split()])
        result["shorthand_configuration"] = f"[{noble_gas_symbol}] {valence_part}"
    
    return result

def parse_orbital_notation(orbital_str):
    """Parse an orbital notation like '2p3' into principal quantum number, orbital type, and electrons."""
    if not orbital_str:
        return None
    
    # Extract the principal quantum number and orbital type
    i = 0
    while i < len(orbital_str) and orbital_str[i].isdigit():
        i += 1
    
    if i == 0 or i == len(orbital_str):
        return None
    
    n = int(orbital_str[:i])
    
    # Extract the orbital type
    j = i
    while j < len(orbital_str) and orbital_str[j].isalpha():
        j += 1
    
    if j == i:
        return None
    
    orbital_type = orbital_str[i:j]
    
    # Extract the number of electrons
    electrons = None
    if j < len(orbital_str):
        electrons = int(orbital_str[j:])
    
    return {
        "n": n,
        "orbital_type": orbital_type,
        "electrons": electrons
    }

def get_orbital_diagram(orbital_str):
    """Generate an orbital diagram showing electrons as arrows (↑ and ↓)."""
    parsed = parse_orbital_notation(orbital_str)
    if not parsed or parsed["electrons"] is None:
        return "Invalid orbital notation"
    
    n = parsed["n"]
    orbital_type = parsed["orbital_type"]
    electrons = parsed["electrons"]
    
    if orbital_type not in ORBITAL_CAPACITIES:
        return f"Invalid orbital type: {orbital_type}"
    
    capacity = ORBITAL_CAPACITIES[orbital_type]
    if electrons > capacity:
        return f"Too many electrons for {orbital_type} orbital (max: {capacity})"
    
    # Number of orbitals based on orbital type
    num_orbitals = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}[orbital_type]
    
    # Create orbital diagram
    diagram = []
    remaining = electrons
    
    for _ in range(num_orbitals):
        if remaining >= 2:
            diagram.append("↑↓")
            remaining -= 2
        elif remaining == 1:
            diagram.append("↑ ")
            remaining -= 1
        else:
            diagram.append("  ")
    
    return f"{n}{orbital_type}: |" + "| |".join(diagram) + "|"

def calculate_unpaired_electrons(atomic_number):
    """Calculate the number of unpaired electrons in an atom."""
    result = get_full_electron_configuration(atomic_number)
    if not result["success"]:
        return result
    
    # Use actual configuration if anomalous
    config = result["full_configuration"]
    if "is_anomalous" in result and result["is_anomalous"]:
        config = result["full_configuration"]  # Use actual config for anomalous elements
    
    unpaired_count = 0
    paired_count = 0
    steps = []
    
    steps.append(f"Element: {result['element']} (atomic number {atomic_number})")
    steps.append(f"Electron configuration: {config}")
    
    if "is_anomalous" in result and result["is_anomalous"]:
        steps.append(f"Note: {result['element']} has an anomalous configuration due to stability of half-filled or filled subshells.")
        steps.append(f"Standard configuration would be: {result['standard_configuration']}")
    
    # Process each orbital set in the configuration
    orbitals = config.split()
    if orbitals[0].startswith("["):  # Handle noble gas shorthand
        noble_gas = orbitals[0].strip("[]")
        steps.append(f"Using noble gas core: [{noble_gas}]")
        # We only care about the valence electrons outside the noble gas core
        orbitals = orbitals[1:]
    
    # Diagram for each subshell
    diagrams = []
    
    for orbital_str in orbitals:
        parsed = parse_orbital_notation(orbital_str)
        if not parsed or parsed["electrons"] is None:
            continue
        
        n = parsed["n"]
        orbital_type = parsed["orbital_type"]
        electrons = parsed["electrons"]
        
        # Number of orbitals based on orbital type
        num_orbitals = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}[orbital_type]
        
        # Count paired and unpaired electrons in this subshell
        subshell_unpaired = 0
        subshell_paired = 0
        
        # Calculate distribution
        full_orbitals = electrons // 2  # How many orbitals have 2 electrons
        half_orbitals = electrons % 2   # How many orbitals have 1 electron (if odd)
        
        # First distribute 1 electron to each orbital
        first_pass = min(electrons, num_orbitals)
        subshell_unpaired += first_pass
        remaining = electrons - first_pass
        
        # Then pair up electrons starting from the first orbital
        pairs = remaining // 1  # Each additional electron forms a pair with an existing one
        subshell_unpaired -= pairs
        subshell_paired += pairs * 2
        
        # Adjust for Hund's rule (maximize unpaired electrons)
        if electrons <= num_orbitals:
            # All electrons are unpaired
            subshell_unpaired = electrons
            subshell_paired = 0
        else:
            # Some orbitals will be filled (paired)
            filled_orbitals = (electrons - num_orbitals) if electrons <= 2*num_orbitals else num_orbitals
            subshell_paired = 2 * (electrons - num_orbitals) if electrons <= 2*num_orbitals else 2 * num_orbitals - electrons
            subshell_unpaired = electrons - subshell_paired
        
        diagram = get_orbital_diagram(orbital_str)
        diagrams.append(diagram)
        
        unpaired_count += subshell_unpaired
        paired_count += subshell_paired
        
        steps.append(f"\nSubshell {n}{orbital_type} with {electrons} electrons:")
        steps.append(f"  {diagram}")
        steps.append(f"  Unpaired electrons: {subshell_unpaired}")
        steps.append(f"  Paired electrons: {subshell_paired}")
    
    steps.append(f"\nTotal unpaired electrons: {unpaired_count}")
    steps.append(f"Total paired electrons: {paired_count}")
    
    return {
        "success": True,
        "element": result["element"],
        "atomic_number": atomic_number,
        "configuration": config,
        "unpaired_electrons": unpaired_count,
        "paired_electrons": paired_count,
        "steps": steps,
        "diagrams": diagrams
    }

def count_electrons_in_orbital_type(atomic_number, orbital_type):
    """Count the total number of electrons in a specific orbital type (s, p, d, f, g)."""
    result = get_full_electron_configuration(atomic_number)
    if not result["success"]:
        return result
    
    # Use actual configuration if anomalous
    config = result["full_configuration"]
    
    total_electrons = 0
    steps = []
    
    steps.append(f"Element: {result['element']} (atomic number {atomic_number})")
    steps.append(f"Electron configuration: {config}")
    
    if "is_anomalous" in result and result["is_anomalous"]:
        steps.append(f"Note: {result['element']} has an anomalous configuration due to stability.")
        steps.append(f"Standard configuration would be: {result['standard_configuration']}")
    
    # Process each orbital set in the configuration
    orbitals = config.split()
    if orbitals[0].startswith("["):  # Handle noble gas shorthand
        noble_gas_symbol = orbitals[0].strip("[]")
        # We need to expand the noble gas core to count all electrons
        noble_gas_atomic_number = ELEMENTS[noble_gas_symbol]
        noble_gas_result = get_full_electron_configuration(noble_gas_atomic_number)
        if noble_gas_result["success"]:
            # Replace the noble gas notation with its full configuration
            noble_gas_orbitals = noble_gas_result["full_configuration"].split()
            orbitals = noble_gas_orbitals + orbitals[1:]
            steps.append(f"Expanded noble gas core [{noble_gas_symbol}] to: {noble_gas_result['full_configuration']}")
    
    # Count electrons in the specified orbital type
    found_orbitals = []
    
    for orbital_str in orbitals:
        parsed = parse_orbital_notation(orbital_str)
        if not parsed or parsed["electrons"] is None:
            continue
        
        if parsed["orbital_type"] == orbital_type:
            found_orbitals.append(f"{parsed['n']}{orbital_type}{parsed['electrons']}")
            total_electrons += parsed["electrons"]
    
    if found_orbitals:
        steps.append(f"\nFound {orbital_type}-orbitals: {', '.join(found_orbitals)}")
    else:
        steps.append(f"\nNo {orbital_type}-orbitals found in the configuration.")
    
    steps.append(f"Total electrons in {orbital_type}-orbitals: {total_electrons}")
    
    return {
        "success": True,
        "element": result["element"],
        "atomic_number": atomic_number,
        "configuration": config,
        f"{orbital_type}_electrons": total_electrons,
        "steps": steps
    }

def solve_electron_problem(problem_type, element_input, **kwargs):
    """
    Solve various electron configuration problems.
    
    Args:
        problem_type (str): Type of problem to solve ('unpaired', 'orbital_count', etc.)
        element_input (str): Element symbol or atomic number
        **kwargs: Additional problem-specific parameters
    
    Returns:
        dict: Results including steps, answers, etc.
    """
    # Determine if input is symbol or atomic number
    element = None
    if element_input.isdigit():
        element = get_element_by_number(int(element_input))
    else:
        element = get_element_by_symbol(element_input)
    
    if not element["success"]:
        return element
    
    atomic_number = element["atomic_number"]
    
    if problem_type == "unpaired":
        # Calculate unpaired electrons
        return calculate_unpaired_electrons(atomic_number)
    
    elif problem_type == "orbital_count":
        # Count electrons in a specific orbital type
        orbital_type = kwargs.get("orbital_type", "p")
        if orbital_type not in ORBITAL_CAPACITIES:
            return {
                "success": False,
                "error": f"Invalid orbital type: {orbital_type}"
            }
        
        return count_electrons_in_orbital_type(atomic_number, orbital_type)
    
    elif problem_type == "electron_config":
        # Get electron configuration
        full_result = get_full_electron_configuration(atomic_number)
        if not full_result["success"]:
            return full_result
        
        shorthand_result = get_shorthand_electron_configuration(atomic_number)
        
        steps = [
            f"Element: {full_result['element']} (atomic number {atomic_number})",
            f"Full electron configuration: {full_result['full_configuration']}",
            f"Shorthand electron configuration: {shorthand_result['shorthand_configuration']}"
        ]
        
        if "is_anomalous" in full_result and full_result["is_anomalous"]:
            steps.append(f"Note: {full_result['element']} has an anomalous configuration due to stability.")
            steps.append(f"Standard configuration would be: {full_result['standard_configuration']}")
        
        return {
            "success": True,
            "element": full_result['element'],
            "atomic_number": atomic_number,
            "full_configuration": full_result['full_configuration'],
            "shorthand_configuration": shorthand_result['shorthand_configuration'],
            "is_anomalous": full_result.get('is_anomalous', False),
            "steps": steps
        }
    
    else:
        return {
            "success": False,
            "error": f"Unknown problem type: {problem_type}"
        }

def get_valence_electrons(atomic_number):
    """Calculate the number of valence electrons for an element."""
    result = get_shorthand_electron_configuration(atomic_number)
    if not result["success"]:
        return result
    
    element_symbol = result["element"]
    shorthand_config = result["shorthand_configuration"]
    
    # Special cases for transition metals and other complex elements
    if element_symbol in ["Cu", "Ag", "Au"]:  # Group 11 - 1 valence electron
        return {
            "success": True,
            "element": element_symbol,
            "atomic_number": atomic_number,
            "valence_electrons": 1,
            "explanation": f"{element_symbol} is in group 11 with 1 valence electron (nd10 (n+1)s1)"
        }
    elif element_symbol in ["Zn", "Cd", "Hg"]:  # Group 12 - 2 valence electrons
        return {
            "success": True,
            "element": element_symbol,
            "atomic_number": atomic_number, 
            "valence_electrons": 2,
            "explanation": f"{element_symbol} is in group 12 with 2 valence electrons (nd10 (n+1)s2)"
        }
    
    # Extract the valence part from shorthand configuration
    if "[" in shorthand_config:
        valence_part = shorthand_config.split("]")[1].strip()
    else:
        valence_part = shorthand_config
    
    orbitals = valence_part.split()
    valence_count = 0
    
    # Count electrons in valence orbitals
    for orbital_str in orbitals:
        parsed = parse_orbital_notation(orbital_str)
        if not parsed or parsed["electrons"] is None:
            continue
        
        # For main group elements, only s and p orbitals contribute to valence
        if parsed["orbital_type"] in ["s", "p"]:
            valence_count += parsed["electrons"]
    
    return {
        "success": True,
        "element": element_symbol,
        "atomic_number": atomic_number,
        "valence_electrons": valence_count,
        "explanation": f"Valence electrons from configuration {shorthand_config}"
    }

def get_quantum_numbers(element_input, electron_number=None):
    """
    Get possible quantum numbers for a specific electron in an element.
    
    Args:
        element_input (str): Element symbol or atomic number
        electron_number (int): Which electron to analyze (default: last/valence electron)
    
    Returns:
        dict: Results including possible quantum numbers
    """
    # Determine if input is symbol or atomic number
    element = None
    if isinstance(element_input, str) and element_input.isdigit():
        element = get_element_by_number(int(element_input))
    elif isinstance(element_input, int):
        element = get_element_by_number(element_input)
    else:
        element = get_element_by_symbol(element_input)
    
    if not element["success"]:
        return element
    
    atomic_number = element["atomic_number"]
    element_symbol = element["symbol"]
    
    # Default to the last electron if not specified
    if electron_number is None or electron_number > atomic_number:
        electron_number = atomic_number
    
    result = get_full_electron_configuration(atomic_number)
    if not result["success"]:
        return result
    
    # Use actual configuration if anomalous
    config = result["full_configuration"]
    
    # Find which orbital contains the specified electron
    orbitals = config.split()
    if orbitals[0].startswith("["):  # Handle noble gas shorthand
        noble_gas_symbol = orbitals[0].strip("[]")
        noble_gas_atomic_number = ELEMENTS[noble_gas_symbol]
        noble_gas_result = get_full_electron_configuration(noble_gas_atomic_number)
        if noble_gas_result["success"]:
            # Replace the noble gas notation with its full configuration
            noble_gas_orbitals = noble_gas_result["full_configuration"].split()
            orbitals = noble_gas_orbitals + orbitals[1:]
    
    current_electron = 0
    target_orbital = None
    
    for orbital_str in orbitals:
        parsed = parse_orbital_notation(orbital_str)
        if not parsed or parsed["electrons"] is None:
            continue
        
        n = parsed["n"]
        orbital_type = parsed["orbital_type"]
        electrons = parsed["electrons"]
        
        if current_electron + electrons >= electron_number:
            # The target electron is in this orbital
            target_orbital = parsed
            break
        
        current_electron += electrons
    
    if not target_orbital:
        return {
            "success": False,
            "error": f"Could not determine orbital for electron {electron_number}"
        }
    
    # Determine quantum numbers
    n = target_orbital["n"]  # Principal quantum number
    
    # Azimuthal quantum number (l)
    l_values = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
    l = l_values[target_orbital["orbital_type"]]
    
    # Magnetic quantum number (ml)
    ml_possibilities = list(range(-l, l + 1))
    
    # Spin quantum number (ms)
    ms_possibilities = [-0.5, 0.5]
    
    # Determine electron position in subshell
    electron_in_subshell = electron_number - current_electron
    
    # Calculate most likely distribution based on Hund's rule
    # Number of orbitals in this subshell
    num_orbitals = 2*l + 1
    
    steps = [
        f"Element: {element_symbol} (atomic number {atomic_number})",
        f"Analyzing electron #{electron_number}",
        f"Electron configuration: {config}",
        f"This electron is in the {n}{orbital_type} subshell"
    ]
    
    # Apply Hund's rule to determine likely ml and ms values
    likely_quantum_numbers = []
    
    # For simplicity in this example, we'll provide all possible combinations
    # In a real scenario, you'd need to apply Aufbau principle and Hund's rule more carefully
    for ml in ml_possibilities:
        for ms in ms_possibilities:
            likely_quantum_numbers.append({
                "n": n,
                "l": l,
                "ml": ml,
                "ms": ms
            })
    
    steps.append(f"Possible quantum numbers for this electron:")
    steps.append(f"  Principal quantum number (n): {n}")
    steps.append(f"  Azimuthal quantum number (l): {l}")
    steps.append(f"  Magnetic quantum number (ml): possible values {ml_possibilities}")
    steps.append(f"  Spin quantum number (ms): possible values {ms_possibilities}")
    
    return {
        "success": True,
        "element": element_symbol,
        "atomic_number": atomic_number,
        "electron_number": electron_number,
        "n": n,
        "l": l,
        "ml_possibilities": ml_possibilities,
        "ms_possibilities": ms_possibilities,
        "orbital": f"{n}{orbital_type}",
        "steps": steps
    }

def solve_molecular_orbital_problem(formula, problem_type="bond_order"):
    """
    Solve problems related to molecular orbital theory.
    
    Args:
        formula (str): Molecular formula (e.g., 'O2', 'N2')
        problem_type (str): Type of problem to solve ('bond_order', 'paramagnetic', etc.)
    
    Returns:
        dict: Results including bond order, magnetic properties, etc.
    """
    # Dictionary of common diatomic molecules and their MO electron configurations
    DIATOMIC_MO_DATA = {
        "H2": {"bond_order": 1, "paramagnetic": False, "config": "σ(1s)²"},
        "He2": {"bond_order": 0, "paramagnetic": False, "config": "σ(1s)² σ*(1s)²"},
        "Li2": {"bond_order": 1, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)²"},
        "Be2": {"bond_order": 0, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)²"},
        "B2": {"bond_order": 1, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² π(2p)²"},
        "C2": {"bond_order": 2, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² π(2p)⁴"},
        "N2": {"bond_order": 3, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² π(2p)⁴ σ(2p)²"},
        "O2": {"bond_order": 2, "paramagnetic": True, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² σ(2p)² π(2p)⁴ π*(2p)²"},
        "F2": {"bond_order": 1, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² σ(2p)² π(2p)⁴ π*(2p)⁴"},
        "Ne2": {"bond_order": 0, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² σ(2p)² π(2p)⁴ π*(2p)⁴ σ*(2p)²"},
        "NO": {"bond_order": 2.5, "paramagnetic": True, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² σ(2p)² π(2p)⁴ π*(2p)¹"},
        "CO": {"bond_order": 3, "paramagnetic": False, "config": "σ(1s)² σ*(1s)² σ(2s)² σ*(2s)² π(2p)⁴ σ(2p)²"}
    }
    
    if formula not in DIATOMIC_MO_DATA:
        # If molecule not in database, try to calculate based on valence electrons
        return {
            "success": False,
            "error": f"Molecular orbital data for {formula} not available in database."
        }
    
    molecule_data = DIATOMIC_MO_DATA[formula]
    
    if problem_type == "bond_order":
        steps = [
            f"Analyzing molecular orbital structure of {formula}",
            f"Molecular orbital configuration: {molecule_data['config']}",
            f"Bond order = (number of bonding electrons - number of antibonding electrons) / 2",
            f"Bond order = {molecule_data['bond_order']}"
        ]
        
        return {
            "success": True,
            "molecule": formula,
            "bond_order": molecule_data['bond_order'],
            "mo_configuration": molecule_data['config'],
            "steps": steps
        }
    
    elif problem_type == "paramagnetic":
        steps = [
            f"Analyzing magnetic properties of {formula}",
            f"Molecular orbital configuration: {molecule_data['config']}",
            f"Paramagnetic: {molecule_data['paramagnetic']} ({'has' if molecule_data['paramagnetic'] else 'does not have'} unpaired electrons)"
        ]
        
        return {
            "success": True,
            "molecule": formula,
            "paramagnetic": molecule_data['paramagnetic'],
            "mo_configuration": molecule_data['config'],
            "steps": steps
        }
    
    else:
        return {
            "success": False,
            "error": f"Unknown problem type: {problem_type}"
        }

def main():
    """Entry point for command-line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Electron Configuration Solver")
    parser.add_argument("problem_type", help="Type of problem to solve", 
                        choices=["unpaired", "orbital_count", "electron_config", "quantum_numbers", "molecular_orbital"])
    parser.add_argument("element", help="Element symbol or atomic number")
    parser.add_argument("--orbital-type", help="Orbital type (s, p, d, f) for orbital_count problems", default="p")
    parser.add_argument("--electron", help="Electron number for quantum_numbers problems", type=int)
    parser.add_argument("--molecule", help="Molecular formula for molecular_orbital problems")
    parser.add_argument("--mo-problem", help="Type of molecular orbital problem", 
                        choices=["bond_order", "paramagnetic"], default="bond_order")
    
    args = parser.parse_args()
    
    if args.problem_type == "molecular_orbital":
        if not args.molecule:
            print("Error: --molecule is required for molecular_orbital problems")
            return
        
        result = solve_molecular_orbital_problem(args.molecule, args.mo_problem)
    else:
        kwargs = {}
        if args.problem_type == "orbital_count":
            kwargs["orbital_type"] = args.orbital_type
        elif args.problem_type == "quantum_numbers" and args.electron:
            kwargs["electron_number"] = args.electron
        
        result = solve_electron_problem(args.problem_type, args.element, **kwargs)
    
    if result["success"]:
        if "steps" in result:
            for step in result["steps"]:
                print(step)
        else:
            print(result)
    else:
        print(f"Error: {result['error']}")

if __name__ == "__main__":
    main()