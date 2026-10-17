#!/bin/bash
DIR=/home/data/knee_xray
MANIFEST=/home/data/knee_xray_manifest.tsv
LOG=/home/data/knee_xray_finish.log
BASE=https://ftp.cdc.gov/pub/NHANES/XRays/Nhanes3
MAXPASS=30

cd "$DIR" || exit 1

remaining() {
  awk -F'\t' -v d="$DIR" '{
    f = d "/" $1
    cmd = "stat -c %s \"" f "\" 2>/dev/null"
    sz = 0; cmd | getline sz; close(cmd)
    if (sz != $2) print $1
  }' "$MANIFEST"
}

for pass in $(seq 1 $MAXPASS); do
  remaining > /tmp/knee_remaining2.txt
  n=$(wc -l < /tmp/knee_remaining2.txt)
  echo "=== pass $pass  $(date '+%F %T')  remaining=$n ==="
  [ "$n" -eq 0 ] && { echo "ALL COMPLETE"; break; }

  sed "s|^|$BASE/|" /tmp/knee_remaining2.txt > /tmp/knee_urls2.txt
  xargs -P24 -n1 wget -c -nv -a "$LOG" --no-check-certificate \
        --dns-timeout=30 --connect-timeout=60 --read-timeout=45 \
        -t 20 --waitretry=5 -P . < /tmp/knee_urls2.txt
done

echo "=== finished $(date '+%F %T') ==="
echo "still incomplete: $(remaining | wc -l) / $(wc -l < "$MANIFEST")"
