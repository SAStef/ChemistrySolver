#!/usr/bin/env python3
"""
Chemistry Problem Solver - Main Entry Point
"""
import sys
from ui.terminal_ui import display_main_menu, get_user_choice
from ui.electron_ui import ElectronConfigUI
from ui.kinetics_ui import KineticsUI
# Import other UI modules as needed

def main():
    """Main application entry point."""
    # Display welcome message
    print("\n" + "="*50)
    print("CHEMISTRY PROBLEM SOLVER".center(50))
    print("="*50 + "\n")
    
    while True:
        # Display the main menu with all available modules
        display_main_menu()
        
        # Get user choice
        choice = get_user_choice()
        
        if choice == "0":
            print("\nExiting program. Goodbye!")
            sys.exit(0)
        elif choice == "1":
            # Launch Electron Configuration module
            electron_ui = ElectronConfigUI()
            electron_ui.run()
        elif choice == "2":
            # Launch Chemical Kinetics module
            kinetics_ui = KineticsUI()
            kinetics_ui.run()
        # Add more modules as needed
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()