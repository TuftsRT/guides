import nbformat as nbf

notebooks = [
    "./deep-learning/decoder-pytorch.ipynb",
    "./deep-learning/encoder-pytorch.ipynb",
    "./deep-learning/ner-spacy.ipynb",
    "./deep-learning/rnn-pytorch.ipynb",
    "./deep-learning/w2v-from-scratch.ipynb",
]

text_search_dict = {
    "# HIDDEN": "remove-cell",  # Remove the whole cell
    "# NO CODE": "remove-input",  # Remove only the input
    "# HIDE CODE": "hide-input",  # Hide the input w/ a button to show
}

for path in notebooks:
    ntbk = nbf.read(path, nbf.NO_CONVERT)
    ntbk.get("metadata", {}).get("mystnb", {})["execution"] = "off"
    nbf.write(ntbk, path)
