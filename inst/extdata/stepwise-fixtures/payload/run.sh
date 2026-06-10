#!/bin/bash
# Run the AVC forward-Wald stepwise SAS capture.
#
# Usage:  cd to this directory, then ./run.sh
#
# Set HAZAPPS and MACROS before calling if the auto-discovery below does
# not find the CCF HAZARD install:
#   HAZAPPS=/path/to/hazard/bin  MACROS=/path/to/hazard/macros  ./run.sh

set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

# ---------------------------------------------------------------------------
# Discover HAZAPPS (directory containing hazard.exe / hazard)
# $HAZARD is the CCF-standard env var name on lri-sas-p-* servers
# Prefer ~/hazard/bin (local build) over the system install
# ---------------------------------------------------------------------------
if [ -z "${HAZAPPS:-}" ]; then
  for d in \
      "$HOME/hazard/bin" \
      "${HAZARD:-}" \
      /programs/apps/sas/hazard \
      /opt/hazard/bin \
      /opt/hazard \
      "$HOME/hazard"; do
    [ -z "$d" ] && continue
    if [ -f "$d/hazard.exe" ] || [ -f "$d/hazard" ]; then
      HAZAPPS="$d"
      break
    fi
  done
fi

if [ -z "${HAZAPPS:-}" ]; then
  echo "ERROR: Cannot find the HAZARD binary."
  echo "  Set HAZAPPS=/path/to/dir/containing/hazard.exe before running."
  exit 1
fi
export HAZAPPS

# ---------------------------------------------------------------------------
# Discover MACROS (directory containing hazard.sas for SASAUTOS)
# Always verify -- the ambient $MACROS env var may point to the general
# SAS macro library, not the HAZARD-specific macros.
# ---------------------------------------------------------------------------
_hz_macros=""
for d in \
    "$HAZAPPS" \
    "$HAZAPPS/macros" \
    "$(cd "$HAZAPPS/.." 2>/dev/null && pwd)/macros" \
    /programs/apps/sas/hazard \
    /opt/hazard/macros \
    "$HOME/hazard/macros"; do
  [ -z "$d" ] && continue
  if [ -f "$d/hazard.sas" ]; then
    _hz_macros="$(cd "$d" && pwd)"
    break
  fi
done

if [ -z "$_hz_macros" ]; then
  echo "ERROR: Cannot find hazard.sas (SASAUTOS macros directory)."
  echo "  Looked in: $HAZAPPS, $HAZAPPS/macros, and common paths."
  echo "  Set HZ_MACROS=/path/to/dir/containing/hazard.sas before running."
  exit 1
fi
export HZ_MACROS="$_hz_macros"   # use HZ_MACROS to avoid clobbering ambient $MACROS

echo "HAZAPPS=$HAZAPPS"
echo "HZ_MACROS=$HZ_MACROS"
echo ""

mkdir -p capture

# SAS Foundation overrides $TMPDIR internally to /saswork/, breaking the
# XPORT handoff to the C binary.  Pass HZ_TMPDIR so the SAS script can call
# OPTIONS SET=TMPDIR before %HAZARD runs.  Use $HOME/tmp (short path) because
# the HAZARD C binary likely has a fixed-size path buffer (~64-80 chars).
mkdir -p "$HOME/tmp"
export HZ_TMPDIR="$HOME/tmp"

echo "Running SAS..."
# Run SAS; capture exit code without letting set -e abort before we report
sas_rc=0
sas avc-forward-wald.sas || sas_rc=$?

if [ $sas_rc -ne 0 ]; then
  echo "WARNING: SAS exited with code $sas_rc -- check avc-forward-wald.log for errors."
fi

# Copy the SAS listing to capture/ -- this is the raw output parsed by R
if [ -f avc-forward-wald.lst ]; then
  cp avc-forward-wald.lst capture/
  echo "Listing captured: capture/avc-forward-wald.lst"
else
  echo "WARNING: avc-forward-wald.lst not found -- SAS may not have run."
fi

# Always show last few lines of the log so errors are visible
if [ -f avc-forward-wald.log ]; then
  echo ""
  echo "--- last 20 lines of avc-forward-wald.log ---"
  tail -20 avc-forward-wald.log
fi

echo ""
echo "Files in capture/:"
ls -lh capture/ 2>/dev/null || echo "(empty)"
