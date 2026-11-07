#!/usr/bin/env bash
# Encode the frames written by render.py into web video plus a poster image.
# Usage: ./encode.sh <render out dir> <site assets dir>   e.g. ./encode.sh ./out ../../assets/video
set -euo pipefail
out="${1:?render output folder}"; dest="${2:?destination folder}"
mkdir -p "$dest"
fps=30
# H.264 MP4 plays everywhere; +faststart lets it start before it has fully downloaded.
ffmpeg -y -v error -framerate $fps -i "$out/frames/f_%05d.png" \
  -c:v libx264 -preset slow -crf 26 -pix_fmt yuv420p -tune animation -level:v 4.0 -movflags +faststart -an \
  "$dest/denv2-ii.mp4"
# (A VP9 WebM came out larger than this MP4 for this flat, line-art animation, so only the MP4 is made.)
# Poster: the frame at about 2016, shown before playback and to people who prefer reduced motion.
n=$(ls "$out/frames" | wc -l | tr -d ' ')
pick=$(( (2016 - 1998) * fps ))
ffmpeg -y -v error -i "$out/frames/$(printf 'f_%05d.png' $pick)" -q:v 3 "$dest/denv2-ii-poster.jpg"
ls -lh "$dest"
