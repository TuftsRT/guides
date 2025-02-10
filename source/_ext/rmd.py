import re

import frontmatter
from jupytext import reads


def convert(text, **kwargs):
    post = frontmatter.loads(text)
    kernelspec = {"display_name": "R", "language": "r", "name": "ir"}
    if "jupyter" in post:
        if "kernelspec" not in post["jupyter"]:
            post["jupyter"]["kernelspec"] = kernelspec
    else:
        post["jupyter"] = {"kernelspec": kernelspec}
    if "tags" in post:
        post["jupyter"]["tags"] = post.metadata.pop("tags")
    node = reads(text=frontmatter.dumps(post), fmt="Rmd", **kwargs)
    for cell in node["cells"]:
        if cell["cell_type"] == "code":
            cell["source"] = re.sub(pattern="^%%.+\n", repl="", string=cell["source"])
    return node
