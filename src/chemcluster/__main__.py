import sys
import streamlit.web.cli as stcli
from pathlib import Path

def run_app():
    app_path = Path(__file__).parent / "app.py"
    sys.argv = ["streamlit", "run", str(app_path)]
    sys.exit(stcli.main())

if __name__ == "__main__":
    run_app()