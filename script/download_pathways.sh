#!/bin/sh
# verbose and errexit
set -ve

usage="download_pathways.sh <outdir>"
if [ -z "$1" ]; then
  echo "$usage" 1>&2
  exit 1
fi
OUTDIR="$1"
if [ ! -d "$OUTDIR" ]; then
  echo "$usage" 1>&2
  exit 1
fi

# download KEGG graphml with ENSG as node identifier
KEGG_OUTDIR="$OUTDIR/kegg_ensg"
mkdir -p "$KEGG_OUTDIR"
kegg.R --outdir "$KEGG_OUTDIR"

# download mapping file
# TODO take release version as an option (GRCh38.94)
MAPPING_BASENAME_GZ="Homo_sapiens.GRCh38.94.refseq.tsv.gz"
MAPPING_BASENAME="Homo_sapiens.GRCh38.94.refseq.tsv"
if [ ! -s "$OUTDIR/$MAPPING_BASENAME_GZ" -a ! -s "$OUTDIR/$MAPPING_BASENAME" ]; then
  wget -P "$OUTDIR" "ftp://ftp.ensembl.org/pub/release-94/tsv/homo_sapiens/$MAPPING_BASENAME_GZ"
fi
if [ ! -s "$OUTDIR/$MAPPING_BASENAME" ]; then
  gunzip -c "$OUTDIR/$MAPPING_BASENAME_GZ" > "$OUTDIR/$MAPPING_BASENAME"
fi

# download network file (to help with mapping)
# TODO take release version as an option
NETWORK_VERSION="v10.5"
NETWORK_BASENAME_GZ="9606.protein.links.detailed.$NETWORK_VERSION.txt.gz"
NETWORK_BASENAME="9606.protein.links.detailed.$NETWORK_VERSION.txt"
if [ ! -s "$OUTDIR/$NETWORK_BASENAME_GZ" -a ! -s "$OUTDIR/$NETWORK_BASENAME" ]; then
  wget -P "$OUTDIR" "https://stringdb-static.org/download/protein.links.detailed.$NETWORK_VERSION/$NETWORK_BASENAME_GZ"
fi
if [ ! -s "$OUTDIR/$NETWORK_BASENAME" ]; then
  gunzip -c "$OUTDIR/$NETWORK_BASENAME_GZ" > "$OUTDIR/$NETWORK_BASENAME"
fi

# apply mapping to graphml, resulting in networks with nodes as ENSP identifiers
KEGG_MAPPED_OUTDIR="$OUTDIR/kegg_ensp"
mkdir -p "$KEGG_MAPPED_OUTDIR"
map_ensg_to_ensp.py -m "$OUTDIR/$MAPPING_BASENAME" -n "$OUTDIR/$NETWORK_BASENAME" --infiles "$KEGG_OUTDIR"/*.graphml --outdir "$KEGG_MAPPED_OUTDIR"
