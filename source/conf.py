import os
import re
import sys
from datetime import date
from urllib.parse import ParseResult, urljoin, urlparse

from pygit2 import GitError, Repository
from sphinx.errors import ConfigError

# needed for Sphinx to load extensions properly
sys.path.append(os.path.abspath("_ext"))

project = "TTS Research Technology Guides"
author = "Tufts University"
email = "tts-research@tufts.edu"

language = "en"

try:
    repo = Repository(".")
    repo_url: str = repo.remotes["origin"].url
    match = re.match(r"(?:https://|git@)([^/:]+)[/:]([^/]+)/(.+?)(?:\.git)?$", repo_url)
    github_user, github_repo = match.group(2), match.group(3)
    if repo.head.shorthand == "main":
        release = "1.0.0"
    else:
        release = "1.0.0-dev"
except (AttributeError, GitError, IndexError):
    raise ConfigError(
        "Unable to automatically derive repository information.\n"
        'Please set "github_user", "github_repo", and "release" manually.'
    ) from None

version = date.today().strftime("%Y%m%d")
copyright = f"{date.today().year}, {author}"

announcement_branch = "announcement"
announcement_file = "announcement.html"
switcher_branch = "switcher"
switcher_file = "switcher.json"

exclude_patterns = ["**/README*"]

extensions = [
    "gallery_directive",
    "myst_nb",
    "notfound.extension",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_external_toc",
    "sphinx.ext.graphviz",
    "tags",
]

external_toc_path = "_toc.yaml"

html_baseurl = os.environ.get(key="BASEURL", default="./")
html_css_files = [
    "style/bugfix.css",
    "style/footer.css",
    "style/gallery.css",
    "style/navbar.css",
    "style/sidebar.css",
    "style/switcher.css",
]
html_favicon = "_static/favicon.ico"
html_last_updated_fmt = ""
html_js_files = ["script/dynamic-nav-dropdown.js"]
html_logo = "_static/jumbo.png"
html_static_path = ["_static"]
html_theme = "pydata_sphinx_theme"
html_title = project

disclaimer = f"Linked external resources not affiliated with or endorsed by {author}."

footer_links = {
    "Accessibility": "https://access.tufts.edu/digital-accessibility-policy",
    "Non-Discrimination": "https://oeo.tufts.edu/policies-procedures/non-discrimination-statement/",
    "Privacy": "https://www.tufts.edu/about/privacy",
}

html_context = {
    "disclaimer": disclaimer,
    "doc_path": "source",
    "footer_links": footer_links,
    "github_repo": github_repo,
    "github_user": github_user,
    "github_version": repo.head.shorthand,
    "short_title": "TTS RT Guides",
}

icon_links = [
    {
        "name": "Tags",
        "url": urljoin(
            html_baseurl,
            f"{'dev/' if 'dev' in release and html_baseurl != './' else ''}"
            "tags/index.html",
        ),
        "icon": "fa-solid fa-tags",
        "attributes": {"target": "_self"},
    },
    {
        "name": "GitHub",
        "url": f"https://github.com/{github_user}/{github_repo}",
        "icon": "fa-brands fa-github",
    },
    {
        "name": "Email",
        "url": f"mailto:{email}",
        "icon": "fa-solid fa-envelope",
    },
]

html_theme_options = {
    "analytics": {
        # cSpell:disable-next-line
        "google_analytics_id": "G-ZJ5LXFRJ2G",
    },
    "announcement": (
        f"https://raw.githubusercontent.com/{github_user}/{github_repo}/"
        f"refs/heads/{announcement_branch}/{announcement_file}"
    ),
    "collapse_navigation": True,
    "footer_center": ["disclaimer", "footer-links"],
    "footer_end": ["switcher-with-label"],
    "footer_start": ["copyright", "last-updated"],
    "icon_links": icon_links,
    "logo": {"text": project},
    "navigation_with_keys": False,
    "search_bar_text": "",
    "secondary_sidebar_items": [
        "page-toc",
        "tags",
        "edit-this-page",
        "sourcelink",
        "report-error",
    ],
    "show_version_warning_banner": True,
    "switcher": {
        "json_url": (
            f"https://raw.githubusercontent.com/{github_user}/{github_repo}/"
            f"refs/heads/{switcher_branch}/{switcher_file}"
        ),
        "version_match": release,
    },
    "use_edit_page_button": True,
}

myst_enable_extensions = [
    "amsmath",
    "attrs_block",
    "attrs_inline",
    "deflist",
    "dollarmath",
    "fieldlist",
    "linkify",
    "replacements",
    "strikethrough",
    "substitution",
    "tasklist",
]

myst_heading_anchors = 2
myst_linkify_fuzzy_links = False

myst_substitutions = {
    "email": f"<{email}>",
    "mailto:": f"mailto:{email}",
}

nb_custom_formats = {
    ".Rmd": "rmd.convert",
}

baseurl_obj: ParseResult = urlparse(
    urljoin(
        html_baseurl, f"{'dev/' if 'dev' in release and html_baseurl != './' else ''}"
    )
)
notfound_urls_prefix = baseurl_obj.path.strip(".")

templates_path = ["_templates"]
