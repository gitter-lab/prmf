#!/bin/sh
# submit job over ssh to optimization cluster
set -e

USAGE="usage: nmf_m3r.sh <outdir>"
OUTDIR="$1"
if [ -z "$OUTDIR" ]; then
  echo $USAGE >&2
  exit 1
fi

BASENAME="${OUTDIR##*/}"
REMOTE_HOME=/mnt/ws/home/abaker/
REMOTE_OUTDIR="${REMOTE_HOME}/ampl/${BASENAME}"

# upload data files
/usr/bin/cp ${HOME}/repo/cs799-f16-nbs/factor/ampl/ampl.sub ${OUTDIR}
/usr/bin/cp ${HOME}/repo/cs799-f16-nbs/factor/ampl/gnmf.mod ${OUTDIR}
/usr/bin/cp ${HOME}/repo/cs799-f16-nbs/factor/ampl/nmf_m3r.cmd ${OUTDIR}
ssh abaker@opt-submit mkdir -p ${REMOTE_OUTDIR}
rsync ${OUTDIR}/{*.dat,ampl.sub,gnmf.mod,nmf_m3r.cmd} abaker@opt-submit:${REMOTE_OUTDIR}

# submit
ssh abaker@opt-submit "${REMOTE_HOME}/repo/cs799-f16-nbs/factor/script/condor_submit_cd" "${REMOTE_OUTDIR}" ampl.sub
