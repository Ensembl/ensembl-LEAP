#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
# Connection values remain inside the JSON, never in the command arguments.
exec python3 "$(dirname "${BASH_SOURCE[0]}")/get_current_annotation.py" "$@"
