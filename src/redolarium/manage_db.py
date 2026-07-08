import os
import sys
import argparse

def link_snapshot(active_snapshot, db_dir):
    source = os.path.abspath(os.path.join(db_dir, active_snapshot))
    target = os.path.abspath(os.path.join(os.path.dirname(db_dir), "active_db"))
    
    if not os.path.exists(source):
        print(f"Error: Database snapshot directory does not exist: {source}")
        sys.exit(1)
        
    if os.path.exists(target) or os.path.islink(target):
        if os.path.islink(target) or os.path.isfile(target):
            os.remove(target)
        else:
            # directory symlink on Windows might require rmdir
            import shutil
            shutil.rmtree(target)
            
    try:
        if os.name == 'nt':
            # On Windows, creating a symlink requires admin or developer mode.
            # Fallback to junction point or copying if symlink fails.
            import subprocess
            subprocess.run(["cmd", "/c", "mklink", "/J", target, source], check=True)
        else:
            os.symlink(source, target)
        print(f"Successfully linked database snapshot: {target} -> {source}")
    except Exception as e:
        print(f"Failed to symlink database snapshot: {e}")
        # Fallback: copy or log warning
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manage database snapshot symlinks.")
    parser.add_argument("--snapshot", required=True, help="Snapshot folder name (e.g., 2026-05)")
    parser.add_argument("--db-dir", default="resources/db_snapshots", help="Path to database snapshots root")
    args = parser.parse_args()
    link_snapshot(args.snapshot, args.db_dir)
