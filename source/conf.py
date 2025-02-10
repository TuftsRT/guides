import os
import sys
from datetime import date
from urllib.parse import ParseResult, urljoin, urlparse

from pygit2 import GitError, Repository
from sphinx.errors import ConfigError

# needed for Sphinx to load extensions properly
sys.path.append(os.path.abspath("_ext"))

project = "TTS Research Computing Guides"
author = "Tufts University"
email = "tts-research@tufts.edu"

try:
    repo = Repository(".")
    repo_url: str = repo.remotes["origin"].url
    github_user, github_repo = repo_url.split(":")[1].split(".")[0].split("/")
    if repo.head.shorthand == "main":
        release = "pub"
    else:
        release = "dev"
except (GitError, IndexError):
    raise ConfigError(
        "Unable to automatically derive repository information.\n"
        'Please set "github_user", "github_repo", and "release" manually.'
    )

version = date.today().strftime("%Y%m%d")
copyright = f"{date.today().year}, {author}"

announcement_branch = "announcement"
announcement_file = "announcement.html"
switcher_branch = "switcher"
switcher_file = "switcher.json"

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_external_toc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.extlinks",
    "sphinx.ext.graphviz",
    "sphinx.ext.intersphinx",
    "tags",
]

autosectionlabel_prefix_document = True

html_baseurl = os.environ.get(key="BASEURL", default="./")
html_css_files = ["custom.css"]
html_favicon = "_static/favicon.ico"
html_last_updated_fmt = ""
html_logo = "_static/jumbo.png"
html_static_path = ["_static"]
html_theme = "pydata_sphinx_theme"
html_title = project

html_context = {
    "github_user": github_user,
    "github_repo": github_repo,
    "github_version": "main",
    "doc_path": "source",
}

if release == "published":
    tags_url = urljoin(html_baseurl, "tags/index.html")
else:
    tags_url = urljoin(html_baseurl, "dev/tags/index.html")

icon_links = [
    {
        "name": "Tags",
        "url": tags_url,
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
    "announcement": (
        f"https://raw.githubusercontent.com/{github_user}/{github_repo}/"
        f"refs/heads/{announcement_branch}/{announcement_file}"
    ),
    "collapse_navigation": True,
    "footer_center": ["last-updated"],
    "footer_end": ["switcher-label", "version-switcher"],
    "footer_start": ["copyright"],
    "header_links_before_dropdown": 4,
    "icon_links": icon_links,
    "logo": {"text": project},
    "navbar_align": "right",
    "navigation_depth": 1,
    "navigation_with_keys": False,
    "search_bar_text": "",
    "secondary_sidebar_items": [
        "page-toc",
        "tags",
        "edit-this-page",
        "sourcelink",
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
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "linkify",
    "replacements",
    "substitution",
]

nb_custom_formats = {
    ".Rmd": "rmd.convert",
}

baseurl_obj: ParseResult = urlparse(html_baseurl)
notfound_urls_prefix = baseurl_obj.path

templates_path = ["_templates"]
