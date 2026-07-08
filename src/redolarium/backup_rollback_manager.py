# type: ignore
import os
import sys
import json
import tarfile
import shutil
import subprocess
from datetime import datetime
import yaml

_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from redolarium.audit_logger import RedolariumAuditLogger

BACKUP_DIR = os.path.join(_project_root, "backup_rollback")
CONFIG_PATH = os.path.join(_project_root, "config", "config.yaml")
TRACKING_FILE = os.path.join(BACKUP_DIR, "last_tracked_version.txt")

class RedolariumBackupRollbackManager:
    @staticmethod
    def create_snapshot(tag: str) -> str:
        os.makedirs(BACKUP_DIR, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        archive_name = f"redolarium_snapshot_{tag}_{timestamp}.tar.gz"
        archive_path = os.path.join(BACKUP_DIR, archive_name)
        
        root_dir = _project_root
        
        with tarfile.open(archive_path, "w:gz") as tar:
            for item in os.listdir(root_dir):
                item_path = os.path.join(root_dir, item)
                # Exclude backup_rollback, virtual environments, git, and large output/cache folders
                if item in ["backup_rollback", "results", "audit", "comparative_genomics_results", "venv", ".venv", ".git", ".snakemake", "__pycache__"]:
                    continue
                if item == "resources":
                    # Add resources directory but filter out transient caches
                    for sub_item in os.listdir(item_path):
                        if sub_item in ["kegg_cache", "ncbi_query_cache.json", "tool_snapshots"]:
                            continue
                        tar.add(os.path.join(item_path, sub_item), arcname=os.path.join("resources", sub_item))
                else:
                    tar.add(item_path, arcname=item)
                
        RedolariumAuditLogger.log_action(f"Snapshot ({tag})", "Success")
        print(f"State snapshot successfully archived to: {archive_path}")
        return archive_path

    @classmethod
    def check_and_auto_snapshot(cls):
        """Reads active config.yaml version and triggers auto-snapshot if incremented."""
        os.makedirs(BACKUP_DIR, exist_ok=True)
        if not os.path.exists(CONFIG_PATH):
            return
            
        try:
            with open(CONFIG_PATH, "r") as f:
                config_data = yaml.safe_load(f) or {}
            current_version = str(config_data.get("version", "2.0.0"))
        except Exception:
            current_version = "2.0.0"
            
        last_version = ""
        if os.path.exists(TRACKING_FILE):
            with open(TRACKING_FILE, "r") as f:
                last_version = f.read().strip()
                
        if current_version != last_version:
            print(f"New major version detected: {current_version} (previous: {last_version}). Triggering automatic snapshot...")
            cls.create_snapshot(f"auto_v{current_version}")
            with open(TRACKING_FILE, "w") as f:
                f.write(current_version)
            print("Automatic version snapshot complete.")

    @classmethod
    def verify_environment(cls) -> bool:
        """Post-restore verification step: compiles python source files and runs unit tests."""
        print("Starting post-restore environment verification step...")
        root_dir = _project_root
        
        # 1. Compile check
        try:
            py_files = []
            for root, dirs, files in os.walk(os.path.join(root_dir, "redolarium")):
                for file in files:
                    if file.endswith(".py"):
                        py_files.append(os.path.join(root, file))
            
            for file in py_files:
                import py_compile
                py_compile.compile(file, doraise=True)
            print("[PASS] Codebase compile check passed.")
        except Exception as e:
            print(f"[FAIL] Post-restore compilation check failed: {e}")
            return False
            
        # 2. Run unit tests
        try:
            test_script = os.path.join(root_dir, "tests", "test_validation.py")
            res = subprocess.run([sys.executable, test_script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if res.returncode == 0:
                print("[PASS] Post-restore test validation check passed.")
                return True
            else:
                print(f"[FAIL] Post-restore test validations failed: {res.stderr}")
                return False
        except Exception as e:
            print(f"[FAIL] Test validation run failed: {e}")
            return False

    @classmethod
    def rollback(cls, archive_name: str, authorized: bool = True) -> bool:
        """Executes rollback sequence: captures failed state, extracts target backup, and runs verification."""
        archive_path = os.path.join(BACKUP_DIR, archive_name)
        if not os.path.exists(archive_path):
            print(f"Error: Snapshot archive not found at {archive_path}")
            RedolariumAuditLogger.log_action(f"Rollback to {archive_name}", "Failure")
            return False
            
        if not authorized:
            print("Rollback aborted: Missing user authorization.")
            return False
            
        # 1. Capture Pre-Rollback Snapshot
        print("Capturing Pre-Rollback Snapshot of the current state before restore...")
        cls.create_snapshot("pre_rollback")
        
        # 2. Extract targeted backup overwriting files
        print(f"Restoring system files from snapshot: {archive_name}...")
        root_dir = _project_root
        
        try:
            # Clean current working files to prevent stale references (except backup directory)
            for item in os.listdir(root_dir):
                item_path = os.path.join(root_dir, item)
                if item in ["backup_rollback", "resources", "results", "audit", "comparative_genomics_results"]:
                    continue
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                else:
                    os.remove(item_path)
                    
            # Extract tarball
            with tarfile.open(archive_path, "r:gz") as tar:
                tar.extractall(path=root_dir)
                
            print("Extraction finished.")
            
            # 3. Post-restore Verification
            verified = cls.verify_environment()
            if verified:
                print("[SUCCESS] Environment is fully operational post-restore.")
                RedolariumAuditLogger.log_action(f"Rollback to {archive_name}", "Success")
                return True
            else:
                print("[FAIL] Post-restore verification checks failed. Environment is unstable.")
                RedolariumAuditLogger.log_action(f"Rollback to {archive_name}", "Failure")
                return False
        except Exception as e:
            print(f"Error during rollback: {e}")
            RedolariumAuditLogger.log_action(f"Rollback to {archive_name}", "Failure")
            return False

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Redolarium Backup & Rollback Manager CLI")
    parser.add_argument("--snapshot", help="Create a manual backup snapshot with specified tag")
    parser.add_argument("--check-version", action="store_true", help="Auto check config version increment and snapshot")
    parser.add_argument("--rollback", help="Initiate rollback to a specific snapshot archive filename")
    parser.add_argument("--verify", action="store_true", help="Run system verification tests")
    
    args = parser.parse_args()
    
    if args.snapshot:
        RedolariumBackupRollbackManager.create_snapshot(args.snapshot)
    elif args.check_version:
        RedolariumBackupRollbackManager.check_and_auto_snapshot()
    elif args.rollback:
        RedolariumBackupRollbackManager.rollback(args.rollback, authorized=True)
    elif args.verify:
        RedolariumBackupRollbackManager.verify_environment()
