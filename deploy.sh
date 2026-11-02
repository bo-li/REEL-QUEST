#!/bin/bash
# Deployment helper for the REEL-QUEST dashboard.
# Usage: ./deploy.sh [--port PORT] [--host HOST] [--skip-checks]
#
# Creates a virtualenv, installs dependencies, verifies the database loads,
# then launches Streamlit. Docker users want `docker compose up -d` instead.
#
# Keep this file LF-only: a CRLF shebang makes the kernel look for
# "/bin/bash\r" and the script dies with "bad interpreter". .gitattributes
# pins it, and tests/test_deploy_script.py checks it.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

PORT=8501
HOST="localhost"
RUN_CHECKS=1

# Run from the repo root regardless of where the caller invoked us, so the
# relative paths below and Streamlit's own CWD resolution both hold.
cd "$(dirname "${BASH_SOURCE[0]}")"

while [[ $# -gt 0 ]]; do
    case $1 in
        --port)
            PORT="$2"
            shift 2
            ;;
        --host)
            HOST="$2"
            shift 2
            ;;
        --skip-checks)
            RUN_CHECKS=0
            shift
            ;;
        --help)
            echo "Usage: ./deploy.sh [--port PORT] [--host HOST] [--skip-checks]"
            echo ""
            echo "Options:"
            echo "  --port PORT    Port to run on (default: 8501)"
            echo "  --host HOST    Host to bind to (default: localhost)"
            echo "                 Use 0.0.0.0 for external access"
            echo "  --skip-checks  Skip the database load check and test suite"
            echo "  --help         Show this help message"
            echo ""
            echo "Example:"
            echo "  ./deploy.sh --port 8080 --host 0.0.0.0"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

if ! [[ "$PORT" =~ ^[0-9]+$ ]] || [ "$PORT" -lt 1 ] || [ "$PORT" -gt 65535 ]; then
    echo -e "${RED}ERROR: --port must be a number between 1 and 65535 (got '$PORT')${NC}"
    exit 1
fi

echo -e "${GREEN}REEL-QUEST Dashboard - Deployment Script${NC}"
echo "=============================================="
echo ""

# ---------------------------------------------------------------- Python ---
echo -e "${YELLOW}Checking Python version...${NC}"
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: Python 3 not found${NC}"
    echo "Please install Python 3.10 or higher"
    exit 1
fi

# pyproject.toml declares requires-python >= 3.10.
if ! python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
    echo -e "${RED}ERROR: Python 3.10+ required, found $(python3 --version)${NC}"
    exit 1
fi
echo -e "${GREEN}✓ $(python3 --version) found${NC}"
echo ""

# ----------------------------------------------------------------- Files ---
echo -e "${YELLOW}Checking required files...${NC}"
REQUIRED_FILES=(
    "app.py"
    "requirements.txt"
    "pyproject.toml"
    "data/sx_ligands_merged_v14_flat.jsonl"
    "rel_quest/__init__.py"
    "rel_quest/core.py"
    "vocab/elements.json"
    "vocab/chemical_families.json"
)

MISSING_FILES=()
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        MISSING_FILES+=("$file")
    fi
done

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo -e "${RED}ERROR: Missing required files:${NC}"
    for file in "${MISSING_FILES[@]}"; do
        echo "  - $file"
    done
    exit 1
fi
echo -e "${GREEN}✓ All required files present${NC}"
echo ""

# ------------------------------------------------------------------ venv ---
echo -e "${YELLOW}Setting up virtual environment...${NC}"
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    echo -e "${GREEN}✓ Virtual environment created${NC}"
else
    echo -e "${GREEN}✓ Virtual environment already exists${NC}"
fi

# Older activate scripts read $PS1 unguarded, which trips `set -u`.
set +u
# shellcheck source=/dev/null
source venv/bin/activate
set -u
echo -e "${GREEN}✓ Virtual environment activated${NC}"
echo ""

# ------------------------------------------------------------------ deps ---
# requirements.txt first (that is what the Dockerfile installs), then the
# package itself editable so `import rel_quest` works from any directory.
# --no-deps because requirements.txt already satisfied them.
echo -e "${YELLOW}Installing dependencies...${NC}"
python -m pip install --upgrade pip -q
python -m pip install -r requirements.txt -q
python -m pip install -e . --no-deps -q
echo -e "${GREEN}✓ Dependencies installed${NC}"
echo ""

# --------------------------------------------------------------- verify ----
# `set -e` already aborts on any failure above, so these blocks report the
# specific problem rather than pretending to branch on $? (which is always 0
# by the time it is read).
if [ "$RUN_CHECKS" -eq 1 ]; then
    echo -e "${YELLOW}Testing database loading...${NC}"
    # LigandDatabase() resolves the canonical file relative to the package,
    # so no path is hardcoded here.
    if ! LIGAND_COUNT=$(python -c "
from rel_quest import LigandDatabase
print(len(LigandDatabase().get_all_ligands()))
" 2>db_check.log | tail -1); then
        echo -e "${RED}ERROR: Failed to load database${NC}"
        cat db_check.log
        rm -f db_check.log
        exit 1
    fi
    rm -f db_check.log
    echo -e "${GREEN}✓ Database loaded successfully ($LIGAND_COUNT ligands)${NC}"
    echo ""

    echo -e "${YELLOW}Validating database...${NC}"
    if python scripts/validate_database.py > validate.log 2>&1; then
        echo -e "${GREEN}✓ All records valid${NC}"
    else
        echo -e "${RED}ERROR: Database validation failed${NC}"
        tail -20 validate.log
        rm -f validate.log
        exit 1
    fi
    rm -f validate.log
    echo ""

    if python -c "import pytest" 2>/dev/null; then
        echo -e "${YELLOW}Running tests...${NC}"
        if python -m pytest -q > pytest.log 2>&1; then
            echo -e "${GREEN}✓ All tests passed${NC}"
        else
            echo -e "${YELLOW}⚠ Some tests failed (non-critical, continuing)${NC}"
            tail -10 pytest.log
        fi
        rm -f pytest.log
        echo ""
    fi
fi

# ---------------------------------------------------------------- launch ---
echo -e "${GREEN}==============================================${NC}"
echo -e "${GREEN}🚀 Launching Dashboard${NC}"
echo -e "${GREEN}==============================================${NC}"
echo ""
echo "Configuration:"
echo "  Host: $HOST"
echo "  Port: $PORT"
echo ""

if [ "$HOST" == "0.0.0.0" ]; then
    echo -e "${YELLOW}⚠ Dashboard will be accessible from other machines${NC}"
    echo -e "${YELLOW}  Make sure firewall allows port $PORT${NC}"
    echo ""
fi

echo "Press Ctrl+C to stop the server"
echo ""
echo "=============================================="
echo ""

if [ "$HOST" == "localhost" ]; then
    exec streamlit run app.py --server.port "$PORT"
else
    exec streamlit run app.py --server.port "$PORT" --server.address "$HOST"
fi
