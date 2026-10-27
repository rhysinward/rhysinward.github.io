# rhysinward.github.io

Personal website built with [Quarto](https://quarto.org) and published with GitHub Pages.

## Preview locally

1. Install [Quarto](https://quarto.org/docs/get-started/) 1.5 or later (the announcement bar needs it), plus R with `dplyr`, `ggplot2`, `knitr` and `rmarkdown`.
2. Run `quarto preview` in this folder. The site opens in your browser and reloads as you edit.

## Publish

1. Create a public GitHub repo called `rhysinward.github.io` and push this folder to `main`.
2. Run `quarto publish gh-pages` once. It renders the site and creates the `gh-pages` branch.
3. In the repo's Settings → Pages, set the source to the `gh-pages` branch, `/ (root)`.
4. After that, every push to `main` republishes the site through `.github/workflows/publish.yml`.

R code only runs on your machine. `freeze: auto` saves the results in `_freeze/`, so commit that folder and the GitHub Action never needs R.

## Where things live

| To change | Edit |
|---|---|
| Colours and fonts | `assets/theme.scss`, the `:root` block near the top of the rules |
| Navigation, footer, announcement bar | `_quarto.yml` |
| Home page | `index.qmd` |
| Timeline | `timeline.yml` |
| Grants and awards | `awards.yml` (shown under the timeline) |
| Talks and posters | `talks.yml`, with slides and posters in `files/` (optional `role` and `note` fields add a speaker role and a prize label) |
| Publications | `publications.yml` |
| PPIE projects | `ppie.qmd` |
| A new tutorial or primer | copy a `.qmd` in `tutorials/`; its first category picks the hub, its second labels it (`Primer` or `Tutorial`) |
| A new pathogen hub | add it to the list at the top of `_templates/hubs.ejs` |
| Your photo | replace `assets/img/rhys-inward.jpg` (square, 480 x 480) |
| Selected publications (home page) | add `featured: <number>` and `summary:` to the entry in `publications.yml`; lowest number first; it needs a Paper or Preprint link |
| Timeline highlights | `key: true` in `timeline.yml` or `awards.yml` |

Text in `[square brackets]` (and `20XX` in the timeline) shows as a dashed placeholder. Delete the `announcement` block in `_quarto.yml` once everything is filled in.
