import json

notebooks = [
    "./deep-learning/decoder-pytorch.ipynb",
    "./deep-learning/encoder-pytorch.ipynb",
    "./deep-learning/ner-spacy.ipynb",
    "./deep-learning/rnn-pytorch.ipynb",
    "./deep-learning/w2v-from-scratch.ipynb",
    "./deep-learning/tokenizers.ipynb",
    "./deep-learning/w2v-gensim.ipynb",
    "./introduction/gentleintroto-nlp-1.ipynb",
    "./introduction/regex-intro.ipynb",
    "./pretrained-models/introducing-semantic-search.ipynb",
    "./pretrained-models/language-hacking-workshop.ipynb",
    "./pretrained-models/ner-pretrained.ipynb",
    "./text-proc/textual-feature-extraction-using.ipynb",
    "./text-proc/traditional-topic-modeling-using.ipynb",
    "./text-proc/whoosh-search-engine.ipynb",
    "./webscraping/dynamic-webscraping.ipynb",
    "./webscraping/static-webscraping-with.ipynb",
    "./webscraping/webscraping-iii-scraping-reddit.ipynb",
]


def update_notebook_metadata(notebook_path, new_metadata):
    try:
        with open(notebook_path, "r", encoding="utf-8") as f:
            notebook = json.load(f)

        if "metadata" not in notebook:
            notebook["metadata"] = {}

        notebook["metadata"].update(new_metadata)

        with open(notebook_path, "w", encoding="utf-8") as f:
            json.dump(notebook, f, indent=1, ensure_ascii=False)

        print(f"✓ Updated: {notebook_path}")
        return True

    except Exception as e:
        print(f"✗ Error updating {notebook_path}: {e}")
        return False


def update_notebooks(paths):
    new_metadata = {"mystnb": {"execution_mode": "off"}}

    success_count = 0
    fail_count = 0

    for path in paths:
        if update_notebook_metadata(path, new_metadata):
            success_count += 1
        else:
            fail_count += 1

    print(f"✓ Successfully updated {success_count} notebooks")
    print(f"✗ Failed to update {fail_count} notebooks")


if __name__ == "__main__":
    update_notebooks(notebooks)
    print("All notebooks have been processed.")
