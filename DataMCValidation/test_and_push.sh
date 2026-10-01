#!/bin/bash
# Build the program, run a quick test on a few events, and -- only if the test
# succeeds -- copy the sources into a checkout of OPALValidation and push.
#
#   bash test_and_push.sh                 # test on 20000 events, then push
#   TEST_EVENTS=0 bash test_and_push.sh   # skip the test run
#
# Environment variables:
#   REPO        local checkout of OPALValidation [/raid5/data/yjlee/OPAL/OPALValidation]
#   SUBDIR      directory inside the repository  [DataMCValidation]
#   TEST_EVENTS events per tree for the test run [20000]
set -euo pipefail

SRC=$(cd "$(dirname "$0")" && pwd)
REPO=${REPO:-/raid5/data/yjlee/OPAL/OPALValidation}
SUBDIR=${SUBDIR:-DataMCValidation}
TEST_EVENTS=${TEST_EVENTS:-20000}
MERGED=/raid5/data/yjlee/OPAL/converted/1994/merged

# ---- 1. build ---------------------------------------------------------------
make -C "$SRC"

# ---- 2. quick test run ------------------------------------------------------
if [ "$TEST_EVENTS" -gt 0 ]; then
  TESTDIR=$SRC/test_output
  mkdir -p "$TESTDIR"
  "$SRC/CompareDataMC" \
    --data "$MERGED/OPAL_1994_data.root" \
    --mc "$MERGED/OPAL_1994_mc_jt74mh.root" \
    --output "$TESTDIR/OPAL_1994_DataMC_validation.root" \
    --plotDir "$TESTDIR" \
    --maxEvents "$TEST_EVENTS"
  for f in pt_spectrum_linear.pdf pt_spectrum_log.pdf ngood_tracks.pdf \
           phi_kk_mass.pdf phi_kk_subtracted.pdf summary.txt \
           OPAL_1994_DataMC_validation.root; do
    [ -s "$TESTDIR/$f" ] || { echo "Test failed: $TESTDIR/$f missing or empty"; exit 1; }
  done
  echo "Test run OK: outputs in $TESTDIR"
fi

# ---- 3. copy into the repository and push -----------------------------------
if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/yenjie/OPALValidation.git "$REPO"
fi
mkdir -p "$REPO/$SUBDIR"
cp "$SRC/CompareDataMC.cpp" "$SRC/Makefile" "$SRC/README.md" \
   "$SRC/run_grendel01.sh" "$SRC/test_and_push.sh" "$SRC/.gitignore" "$REPO/$SUBDIR/"

cd "$REPO"
git add "$SUBDIR"
if git diff --cached --quiet; then
  echo "Nothing new to commit."
  exit 0
fi
git commit -m "Add OPAL 1994 data/MC validation: pT spectra and phi->KK peak" \
  -m "Compiled ROOT program comparing the converted 1994 data with the
JETSET 7.4 + GOPAL MC: per-event charged-particle pT spectra with
data/MC ratio and generator-level overlay, and the phi(1020) -> K+K-
peak in opposite-sign track pairs fitted with a Voigtian on a threshold
background. Written and commented for undergraduate readers." \
  -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
git push
echo "Pushed $SUBDIR to $(git remote get-url origin)"
