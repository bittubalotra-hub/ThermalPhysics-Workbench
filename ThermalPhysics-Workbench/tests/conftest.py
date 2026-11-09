from pathlib import Path
import sys
# Permit pytest without installation in a fresh clone.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
