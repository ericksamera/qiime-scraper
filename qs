#!/usr/bin/env bash
# Run the qs CLI using the src/ layout
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${PYTHONPATH:-$REPO_DIR/src}"
exec python -m qs "$@"

