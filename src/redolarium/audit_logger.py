# type: ignore
import json
import os
import platform
import subprocess
import hashlib
from datetime import datetime
from redolarium.config import CONFIG

_project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AUDIT_LOG_PATH = os.path.join(_project_root, "audit_log.json")

class RedolariumAuditLogger:
    @staticmethod
    def calculate_system_hash() -> str:
        """Calculate a combined SHA-256 hash over all active source files to verify integrity."""
        root_dir = _project_root
        hasher = hashlib.sha256()
        
        # Sort files to ensure deterministic hashing order
        all_files = []
        for root, dirs, files in os.walk(root_dir):
            # Ignore backup_rollback folder to prevent recursive hashing of backups
            if "backup_rollback" in root or "resources/tool_snapshots" in root or "results" in root:
                continue
            for file in files:
                if file.endswith((".py", ".yaml", ".json", "Snakefile", "Makefile")):
                    if file != "audit_log.json":
                        all_files.append(os.path.join(root, file))
        
        all_files.sort()
        
        for file_path in all_files:
            try:
                with open(file_path, "rb") as f:
                    # Update hasher with relative path to capture structure changes
                    rel_path = os.path.relpath(file_path, root_dir)
                    hasher.update(rel_path.encode('utf-8'))
                    while chunk := f.read(8192):
                        hasher.update(chunk)
            except Exception:
                pass
                
        return hasher.hexdigest()

    @classmethod
    def log_action(cls, action_type: str, status: str = "Success"):
        """Logs an operation inside audit_log.json."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        sys_hash = cls.calculate_system_hash()
        
        log_entry = {
            "Timestamp": timestamp,
            "Action_Type": action_type,
            "Status": status,
            "System_Hash": sys_hash
        }
        
        logs = []
        if os.path.exists(AUDIT_LOG_PATH):
            try:
                with open(AUDIT_LOG_PATH, "r", encoding="utf-8") as f:
                    logs = json.load(f)
            except Exception:
                pass
                
        logs.append(log_entry)
        
        with open(AUDIT_LOG_PATH, "w", encoding="utf-8") as f:
            json.dump(logs, f, indent=4)
        print(f"[{status}] Logged action '{action_type}' to audit_log.json. System Hash: {sys_hash[:10]}...")
