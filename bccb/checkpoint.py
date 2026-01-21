"""Checkpoint management for create_crossbar.py"""
import json
import os
from pathlib import Path
from datetime import datetime, timezone, timedelta

TZ = timezone(timedelta(hours=8))

class CheckpointManager:
    def __init__(self, checkpoint_file="output/.checkpoint.json"):
        self.checkpoint_file = Path(checkpoint_file)
        self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
        self.data = self.load()

    def load(self):
        """Load checkpoint from file"""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {"completed": [], "failed": {}, "last_run": None}

    def save(self):
        """Save checkpoint to file"""
        self.data["last_run"] = datetime.now(TZ).isoformat()
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.data, f, indent=2)

    def is_completed(self, adapter_name):
        """Check if adapter is completed"""
        return adapter_name in self.data["completed"]

    def mark_completed(self, adapter_name):
        """Mark adapter as completed"""
        if adapter_name not in self.data["completed"]:
            self.data["completed"].append(adapter_name)
        # Remove from failed if it was there
        self.data["failed"].pop(adapter_name, None)
        self.save()

    def mark_failed(self, adapter_name, error):
        """Mark adapter as failed"""
        self.data["failed"][adapter_name] = {
            "error": str(error),
            "timestamp": datetime.now(TZ).isoformat()
        }
        self.save()

    def reset(self):
        """Reset checkpoint"""
        self.data = {"completed": [], "failed": {}, "last_run": None}
        if self.checkpoint_file.exists():
            self.checkpoint_file.unlink()

    def should_run(self, adapter_name, skip_until=None, only=None):
        """Determine if adapter should run based on checkpoint and CLI args"""
        # If 'only' is specified, only run those adapters
        if only:
            return adapter_name in only

        # If 'skip_until' is specified, skip until we reach that adapter
        if skip_until:
            if self.data.get("_skip_mode", True):
                if adapter_name == skip_until:
                    self.data["_skip_mode"] = False
                    return True
                return False
            return True

        # Normal mode: skip if already completed
        return not self.is_completed(adapter_name)
