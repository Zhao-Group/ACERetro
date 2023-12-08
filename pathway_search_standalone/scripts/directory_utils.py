from pathlib import Path

def project_directory() -> Path:
    """
    Returns the path to the root directory of the repository
    """
    return Path(__file__).parent.parent.parent.resolve()