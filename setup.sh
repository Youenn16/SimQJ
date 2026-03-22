#!/bin/bash
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'  


OS="$(uname -s)"

echo "Detected OS: $OS"

install_dependencies() {
    case "$OS" in
        Linux)
            if command -v apt-get &>/dev/null; then
                sudo apt-get update && sudo apt-get install -y \
                    build-essential cmake libopenblas-dev liblapack-dev liblapacke-dev
            elif command -v yum &>/dev/null; then
                # Red Hat / CentOS clusters
                sudo yum install -y cmake openblas-devel lapack-devel
            elif command -v module &>/dev/null; then
                # HPC with environment modules
                module load cmake openblas
            else
                echo "Unknown package manager — install cmake, openblas, lapack manually"
                exit 1
            fi
            ;;
        Darwin)
            if ! command -v brew &>/dev/null; then
                echo "Homebrew not found — install it from https://brew.sh"
                exit 1
            fi
            brew install cmake openblas lapack
            ;;
        *)
            echo "Windows detected — please use WSL2 and rerun this script"
            exit 1
            ;;
    esac
}

build_cpp() {
    echo "Building C++..."
    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DBLA_VENDOR=OpenBLAS
    cmake --build build --parallel $(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
}
 

install_python() {
    echo "Setting up Python environment..."
    
    # Create virtual environment if it doesn't exist
    if [ ! -d "toolbox-binary/.venv" ]; then
        python3 -m venv toolbox-binary/.venv
    fi 
    
    # Activate it
    source toolbox-binary/.venv/bin/activate
    
    # Upgrade pip
    pip install --upgrade pip
    
    # Install dependencies
    if [ -f "requirements.txt" ]; then
        pip install -r requirements.txt
    fi
    
    # Install the project itself
    pip install -e .
    
    source toolbox-binary/.venv/bin/activate
    
}

install_dependencies
build_cpp
install_python 

echo -e "${BLUE}===================================${NC}"
echo -e "        SimQJ ready to use !          "
echo -e "${BLUE}===================================${NC}"
echo -e "${GREEN}[Info]${NC} Virtual environment ready at ./toolbox-binary/.venv"
echo -e "${GREEN}[Info]${NC} Activate it with: source toolbox-binary/.venv/bin/activate to use simQJ command"
