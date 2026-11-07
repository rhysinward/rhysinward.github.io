# Home-page animation: DENV-2 genotype II

`render.py` draws the frames from the DTA phylogeography MCC tree and map files in
`Documents/grapevne_modules/Beast_Vietnam/visulisation`; `encode.sh` turns them into
`assets/video/denv2-ii.mp4` and its poster. Quarto ignores this folder because its name
starts with an underscore, so none of it is published on the site.

```bash
python3 render.py --data ~/Documents/grapevne_modules/Beast_Vietnam/visulisation --out ./out
./encode.sh ./out ../../assets/video
```

Options: `--start` (default 1998), `--sec-per-year` (default 1), `--fps` (30), `--size` (1080),
and `--only 2008.3,2016.4` to render single stills while adjusting the look.

The approach follows Gytis Dudas's baltic animation examples (curonia,
https://github.com/evogytis/baltic).
