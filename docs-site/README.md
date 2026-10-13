# Flan Documentation

Source for the Flan documentation site, published at
**https://shawnzamperini.github.io/flan/**.

The site is built with [MkDocs](https://www.mkdocs.org/) and the
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme. It
covers installing and running Flan, the full list of input options, the physics
models, developer notes on the code structure and parallelism, and a list of
Flan publications and presentations.

## Layout

```
docs-site/
├── mkdocs.yml      # Site config: theme, extensions and the navigation menu
├── requirements.txt  # Python packages needed to build the docs
├── docs/           # Markdown source for every page
│   ├── index.md
│   ├── getting-started.md
│   ├── input-options.md
│   ├── publications.md
│   ├── installation/
│   ├── running/
│   ├── physics/
│   ├── developer/
│   ├── assets/     # Images, including the browser-tab icon (favicon)
│   └── styles/     # Custom CSS, selected with extra_css in mkdocs.yml
├── .mkdoxy/        # Generated API Reference Markdown (not committed)
└── site/           # Built HTML (generated, not committed)
```

## Building the docs

Install MkDocs, the Material theme and the MkDoxy plugin, preferably in a
virtual environment. `requirements.txt` keeps MkDocs below 2.0, which drops
support for plugins and themes like the ones this site uses:

```bash
pip install -r requirements.txt
```

[Doxygen](https://www.doxygen.nl/) must also be installed. MkDoxy runs it
during every build to generate the API Reference section (under Developer)
from the Doxygen comments in `../src` and `../include`, including the CUDA
files. The generated Markdown is cached in `.mkdoxy/`, which is not committed.

Then run all of the following commands from inside `docs-site/`.

### Preview locally

```bash
mkdocs serve
```

Open http://127.0.0.1:8000 in a browser. The page reloads automatically as you
edit files. On WSL, if the browser can't reach it, use
`mkdocs serve -a 0.0.0.0:8000`.

### Build

```bash
mkdocs build --strict
```

This writes the HTML to `site/`. `--strict` turns warnings (such as broken
links or pages missing from the nav) into errors, which makes it a quick check
before deploying.

### Deploy

```bash
mkdocs gh-deploy
```

This builds the site, commits it to the `gh-pages` branch and pushes it to
GitHub, which updates https://shawnzamperini.github.io/flan/. Don't edit the
`gh-pages` branch by hand, because each deploy overwrites it.

## Adding a page

1. Create a Markdown file under `docs/`, e.g. `docs/developer/new-page.md`.
2. Add it to the `nav:` section of `mkdocs.yml`:

    ```yaml
    - Developer:
        - New Page: developer/new-page.md
    ```

   Pages left out of `nav` still build, but they don't appear in the menu.
