import json
import os
import sys
from datetime import date
from urllib.parse import urljoin, urlparse

from pygit2 import GitError, Repository

sys.path.append(os.path.abspath("_ext"))

project = "Tufts RT Guides"
author = "Tufts University"
email = "tts-research@tufts.edu"

github_user = "tuftsrt"
github_repo = "guides"

dev_label = "dev"
pub_label = "pub"
switcher_json_file = "_static/switcher.json"

version = date.today().strftime("%Y%m%d")
release = dev_label

try:
    if Repository(".").head.shorthand == "main":
        release = pub_label
except GitError:
    pass

copyright = f"{date.today().year}, {author}"

extensions = [
    "gallery_directive",
    "myst_nb",
    "notfound.extension",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_external_toc",
    "sphinx_togglebutton",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.extlinks",
    "sphinx.ext.graphviz",
    "sphinx.ext.intersphinx",
    "tags",
]

autosectionlabel_prefix_document = True

html_baseurl = os.environ.get(key="BASEURL", default="/")
html_favicon = "_static/favicon.ico"
html_last_updated_fmt = ""
html_logo = "_static/jumbo.png"
html_static_path = ['_static']
html_css_files = ["custom.css"]
html_theme = "pydata_sphinx_theme"
html_title = project

html_context = {
    "github_user": github_user,
    "github_repo": github_repo,
    "github_version": "main",
    "doc_path": "source",
}

icon_links = [
    {
        "name": "Tags",
        "url": urljoin(html_baseurl, "tags/index.html"),
        "icon": "fa-solid fa-tags",
        "attributes": {"target": "_self"},
    },
    {
        "name": "GitHub",
        "url": "https://github.com/{}/{}".format(github_user, github_repo),
        "icon": "fa-brands fa-github",
    },
    {
        "name": "Email",
        "url": "mailto:{}".format(email),
        "icon": "fa-solid fa-envelope",
    },
]

html_theme_options = {
    "announcement": (
        "Incomplete and under active development. Subject to change without notice."
    ),
    "collapse_navigation": True,
    "footer_center": ["last-updated"],
    "footer_end": ["version-switcher"],
    "footer_start": ["copyright"],
    "header_links_before_dropdown": 8,
    "icon_links": icon_links,
    "logo": {"text": project},
    "navbar_align": "content",
    "navigation_depth": 1,
    "navigation_with_keys": False,
    "search_bar_text": "",
    "secondary_sidebar_items": [
        "page-toc",
        "tags",
        "edit-this-page",
        "sourcelink",
    ],
    "show_nav_level": 0,
    "switcher": {
        "json_url": urljoin(html_baseurl, switcher_json_file),
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

notfound_urls_prefix = urlparse(html_baseurl).path

templates_path = ["_templates"]

with open(file=switcher_json_file, mode="w") as f:
    json.dump(
        obj=[
            {
                "version": pub_label,
                "url": html_baseurl,
                "preferred": True,
            },
            {
                "version": dev_label,
                "url": urljoin(html_baseurl, "dev/"),
            },
        ],
        fp=f,
        indent=4,
    )
