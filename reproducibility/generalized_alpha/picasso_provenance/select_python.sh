#!/usr/bin/env bash
set -euo pipefail
candidates=(
  /mnt/home/soft/python/programs/x86_64_sles156/3.14.3/bin/python3
  /mnt/home/soft/python/programs/x86_64/python_3.11.4/bin/python3
  python3
)
for p in "${candidates[@]}"; do
  if command -v "$p" >/dev/null 2>&1 || [[ -x "$p" ]]; then
    if "$p" - <<'PY' >/dev/null 2>&1
import numpy, scipy
PY
    then
      printf '%s\n' "$p"
      exit 0
    fi
  fi
done
echo 'No usable Python with numpy+scipy found' >&2
exit 2
