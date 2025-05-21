"""
Common Terminal UI Components for Chemistry Problem Solver
"""

def display_main_menu():
    """Display the main application menu."""
    menu = """
    [1] Electron Configuration Problem Solver
    [2] Chemical Kinetics Problem Solver
    [3] Chemical Equation Solver
    [4] Thermodynamics Problem Solver
    [0] Exit
    """
    print(menu)

def get_user_choice():
    """Get and validate user menu choice."""
    return input("\nEnter choice: ").strip()

def display_title(title_text):
    """Display a formatted title."""
    width = 60
    border = "═" * width
    
    print("\n╔" + border + "╗")
    print("║" + title_text.center(width) + "║")
    print("╚" + border + "╝")

def display_results_header():
    """Display a results section header."""
    print("\n==== RESULTS ====")

def display_steps(steps, skip_first=0):
    """Display solution steps with optional skipping of first n steps."""
    if not steps:
        return
        
    print("\nDetailed analysis:")
    for step in steps[skip_first:]:
        print(step)

def wait_for_user():
    """Wait for user to press Enter before continuing."""
    input("\nPress Enter to return to menu...")