import json
import os
import shutil
import warnings
from os import path
from pathlib import Path
from typing import Any, Dict, Optional

_PACKAGE_DATA = Path(__file__).resolve().parents[1] / "data"
_LEGACY_CONFIG = _PACKAGE_DATA / "config.json"
_SETTINGS_FILENAME = "settings.json"

_PERSISTABLE_KEYS = frozenset({"data_directory", "config"})


class Options_dic:
    """A class to manage and modify configuration paths for gwaslab.

Attributes
        paths (dict): Current configuration paths
        default (dict): Default paths for reset functionality
"""

    def __init__(self, path_dic):
        self.paths = path_dic
        self.default = path_dic.copy()

    def set_option(self, key, path_value, persist: bool = False):
        """Update a specific configuration path.
"""
        self.paths[key] = path_value
        if persist and key in _PERSISTABLE_KEYS:
            save_settings({key: path_value})

    def reset_option(self):
        """Reset all paths to their default values.
"""
        self.paths = self.default.copy()

    def get_option(self):
        """Return the options as a dict.
"""
        return self.paths


def _default_data_directory() -> str:
    return path.expanduser("~/.gwaslab/")


def _settings_path(data_dir: Optional[str] = None) -> Path:
    base = Path(data_dir or _default_data_directory())
    return base / _SETTINGS_FILENAME


def load_settings(data_dir: Optional[str] = None) -> Dict[str, Any]:
    """Load user overrides from ~/.gwaslab/settings.json.
"""
    settings_file = _settings_path(data_dir)
    if not settings_file.is_file():
        return {}
    try:
        with open(settings_file, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        return payload if isinstance(payload, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def save_settings(updates: Dict[str, Any], data_dir: Optional[str] = None) -> None:
    """Merge updates into settings.json under the data directory.
"""
    if not updates:
        return
    current_data_dir = data_dir or options.paths.get("data_directory") or _default_data_directory()
    settings_file = _settings_path(current_data_dir)
    settings_file.parent.mkdir(parents=True, exist_ok=True)
    payload = load_settings(str(current_data_dir))
    payload.update(updates)
    with open(settings_file, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=4, ensure_ascii=False)


def _resolve_default_paths() -> Dict[str, str]:
    """Build default path dict: env > settings.json > built-in defaults.
"""
    data_dir = os.environ.get("GWASLAB_DATA_DIR")
    config_path = os.environ.get("GWASLAB_CONFIG")

    if not data_dir:
        settings = load_settings()
        data_dir = settings.get("data_directory")
    if not data_dir:
        data_dir = _default_data_directory()

    data_dir = path.expanduser(data_dir)
    if not data_dir.endswith("/"):
        data_dir = data_dir + "/"

    if not config_path:
        settings = load_settings(data_dir)
        config_path = settings.get("config")
    if not config_path:
        config_path = path.join(data_dir, "config.json")

    return {
        "config": config_path,
        "reference": str(_PACKAGE_DATA / "reference.json"),
        "formatbook": str(_PACKAGE_DATA / "formatbook.json"),
        "data_directory": data_dir,
    }


def _migrate_legacy_config(new_config_path: str, log=None) -> bool:
    """Copy legacy package-local config.json to the user config path when needed.

    Returns True if migration was performed.
"""
    legacy = _LEGACY_CONFIG
    new_path = Path(new_config_path)
    if new_path.exists():
        return False
    if not legacy.is_file():
        return False
    try:
        with open(legacy, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return False
    if not payload.get("downloaded"):
        return False

    new_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(legacy, new_path)
    msg = (
        f"Migrated legacy config from {legacy} to {new_config_path}. "
        "The old file was kept in the package directory."
    )
    if log is not None:
        log.write(f" -{msg}")
    else:
        warnings.warn(msg, UserWarning, stacklevel=2)
    return True


def ensure_user_layout(log=None) -> None:
    """Create data directory and empty config registry when missing.
"""
    data_dir = options.paths["data_directory"]
    os.makedirs(data_dir, exist_ok=True)

    config_path = options.paths["config"]
    _migrate_legacy_config(config_path, log=log)

    if not path.exists(config_path):
        config_parent = Path(config_path).parent
        config_parent.mkdir(parents=True, exist_ok=True)
        with open(config_path, "w", encoding="utf-8") as handle:
            json.dump({"downloaded": {}}, handle, indent=4)


default_dic = _resolve_default_paths()
options = Options_dic(default_dic)
ensure_user_layout()
