import subprocess
import os

def force_update_repo(repo_path="."):
    try:
        os.chdir(repo_path)

        # Ensure it's a git repository
        subprocess.run(["git", "rev-parse", "--is-inside-work-tree"], check=True, stdout=subprocess.DEVNULL)

        # Fetch latest changes
        subprocess.run(["git", "fetch", "--all"], check=True)

        # Reset to latest version from origin/main (or origin/master)
        subprocess.run(["git", "reset", "--hard", "origin/main"], check=True)

        print("Repository forcibly updated to match origin/main.")

    except subprocess.CalledProcessError as e:
        print("Git command failed:")
        print(e)
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    force_update_repo()